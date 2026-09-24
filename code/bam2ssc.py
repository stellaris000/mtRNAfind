#!/usr/bin/env python3
"""
BAM Processor v8
================
功能：
  1. 默认过滤 secondary / supplementary（--keep-secondary/--keep-supplementary 关闭）
  2. 剪切位点 indel 检测：CIGAR 中 N 紧邻 I/D → is_indel=1
  3. _flnc.ssc  第9列：is_indel (0/1)
     _ssc.count 第6-8列：mean_coverage, mean_identity, indel_ratio

加速优化：
  - has_adjacent_indel 内联到 CIGAR 主循环，消除函数调用和二次遍历
  - ec 累加 sum_identity/sum_coverage，合并时再除以 count 得均值
  - 其余 v7 已有优化保持不变
"""
import argparse
import os
import sys
import logging
import multiprocessing as mp
import tempfile
import subprocess
from collections import defaultdict
from functools import partial

import pysam

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger()

RC_MAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def get_rc(seq): return seq.translate(RC_MAP)[::-1]


def parse_args():
    parser = argparse.ArgumentParser(description='BAM processor v8')
    parser.add_argument("--reference", "-r", required=True)
    parser.add_argument("--bam", "-b", nargs='+', required=True)
    parser.add_argument("--output", "-o", default="isatools_output")
    parser.add_argument("--threads", "-t", type=int, default=4)
    parser.add_argument("--tmpdir", default=None, help="临时目录，建议 /dev/shm")
    parser.add_argument("--target-reads", type=int, default=0)
    parser.add_argument("--keep-secondary", action="store_true")
    parser.add_argument("--keep-supplementary", action="store_true")
    return parser.parse_args()


def get_bam_contigs(bam):
    read_counts = {}
    lengths = {}
    try:
        result = subprocess.run(['samtools', 'idxstats', bam],
                                capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split('\n'):
            if not line.strip(): continue
            parts = line.split('\t')
            contig = parts[0]
            if contig == '*': continue
            lengths[contig] = int(parts[1])
            read_counts[contig] = int(parts[2])
    except Exception as e:
        logger.warning(f"idxstats failed: {e}, using header")
        with pysam.AlignmentFile(bam, 'rb') as bf:
            for ref, length in zip(bf.references, bf.lengths):
                lengths[ref] = length
                read_counts[ref] = 0
    return read_counts, lengths


def build_regions(read_counts, lengths, target_reads_per_chunk):
    regions = []
    for contig, total in read_counts.items():
        if total == 0: continue
        chrom_len = lengths.get(contig, 0)
        if total <= target_reads_per_chunk or chrom_len == 0:
            regions.append((contig, 0, chrom_len, total))
        else:
            n_splits = max(2, round(total / target_reads_per_chunk))
            step = (chrom_len + n_splits - 1) // n_splits
            for i in range(n_splits):
                start = i * step
                end = min((i + 1) * step, chrom_len)
                est = int(total * (end - start) / chrom_len)
                if est > 0:
                    regions.append((contig, start, end, est))
    return regions


def bin_pack_regions(regions, n_bins):
    sorted_regions = sorted(regions, key=lambda x: -x[3])
    bins = [[] for _ in range(n_bins)]
    bin_totals = [0] * n_bins
    for contig, start, end, est in sorted_regions:
        idx = bin_totals.index(min(bin_totals))
        bins[idx].append((contig, start, end))
        bin_totals[idx] += est
    non_empty = [b for b in bins if b]
    logger.info(f"  Bin loads: min={min(bin_totals):,} max={max(bin_totals):,} "
                f"ratio={max(bin_totals)/max(min(bin_totals),1):.2f}x")
    return non_empty


def process_region_group(bam, fasta_file, temp_dir, group_idx, region_list,
                         filter_secondary, filter_supplementary):
    bam_basename = os.path.splitext(os.path.basename(bam))[0]
    out1_tmp = os.path.join(temp_dir, f'o1_{bam_basename}_{group_idx}.tmp')
    out2_tmp = os.path.join(temp_dir, f'o2_{bam_basename}_{group_idx}.tmp')

    id_count = defaultdict(int)
    # ec[key] = [count, sum_identity, sum_coverage, indel_count]
    ec = defaultdict(lambda: [0, 0.0, 0.0, 0])
    seq_cache = {}
    WBUF = 8 * 1024 * 1024

    with pysam.AlignmentFile(bam, 'rb') as bf, \
         pysam.FastaFile(fasta_file) as fa, \
         open(out1_tmp, 'w', buffering=WBUF) as o1_fh:

        fa_fetch = fa.fetch
        write = o1_fh.write

        for contig, reg_start, reg_end in region_list:
            is_full_contig = (reg_start == 0 and
                              reg_end == bf.get_reference_length(contig))
            try:
                read_iter = (bf.fetch(contig) if is_full_contig
                             else bf.fetch(contig, reg_start, reg_end))
            except ValueError:
                continue

            for read in read_iter:
                if read.is_unmapped: continue
                if filter_secondary     and read.is_secondary:     continue
                if filter_supplementary and read.is_supplementary: continue

                ref_start = read.reference_start
                if not is_full_contig:
                    if ref_start < reg_start or ref_start >= reg_end:
                        continue

                seqlen = read.query_length
                if not seqlen: continue
                cigar = read.cigartuples
                if not cigar: continue

                # ── CIGAR 解析 + indel 检测（单次遍历，内联）─────────────
                # indel 判据：N 的紧邻操作为 I(1) 或 D(2)
                # 用 prev_op 记录上一个操作，遇到 N 时检查 prev_op；
                # 遇到 I/D 时检查下一个操作是否为 N（lookahead）。
                # 两种方向合并，确保 M-I-N 和 M-N-I 都能检出。
                pos = ref_start + 1
                cov = 0; clip = 0; gap = 0
                positions = []
                pos_append = positions.append
                is_indel = 0

                n_cigar = len(cigar)
                for ci in range(n_cigar):
                    op, length = cigar[ci]
                    if op == 0:                          # MATCH
                        cov += length; gap += length
                    elif op == 3:                        # N (intron)
                        e_p = pos + gap - 1
                        pos_append(pos); pos_append(e_p)
                        pos = e_p + length + 1; gap = 0
                        # 检查左侧（上一个 op）
                        if not is_indel and ci > 0 and cigar[ci-1][0] in (1, 2):
                            is_indel = 1
                        # 检查右侧（下一个 op）
                        if not is_indel and ci < n_cigar-1 and cigar[ci+1][0] in (1, 2):
                            is_indel = 1
                    elif op == 1:                        # INS
                        cov += length
                    elif op == 4 or op == 5:             # CLIP
                        clip += length
                    else:                                # DEL/other
                        gap += length

                e_p = pos + gap - 1
                pos_append(pos); pos_append(e_p)

                n_pos = len(positions)
                if n_pos > 2:
                    s1 = positions[0]
                    e1 = positions[n_pos - 1]
                    str_pos = '-'.join(map(str, positions[1:n_pos - 1]))
                else:
                    s1 = positions[0] if positions else 'NA'
                    e1 = positions[1] if n_pos > 1 else 'NA'
                    str_pos = 'NA'
                    is_indel = 0  # 无 intron，强制为 0

                # ── tags ──────────────────────────────────────────────────
                xs = None
                try:
                    xs = read.get_tag('XS')
                except KeyError:
                    try:
                        ts = read.get_tag('ts')
                        if ts in ('+', '-'):
                            xs = ('+' if ts == '-' else '-') if read.is_reverse else ts
                    except KeyError:
                        pass
                strand = xs if xs else ('-' if read.is_reverse else '+')

                try:
                    nm = read.get_tag('NM')
                except KeyError:
                    nm = 0

                identity = 1.0 - nm / cov if cov else 0.0
                coverage = (seqlen - clip) / seqlen

                qname = read.query_name
                id_count[qname] += 1

                write(f'{qname}.m{id_count[qname]}\t{contig}\t{strand}\t'
                      f'{s1}\t{e1}\t{str_pos}\t{identity:.4f}\t{coverage:.4f}\t'
                      f'{is_indel}\n')

                if str_pos != 'NA':
                    key = (contig, strand, str_pos)
                    slot = ec[key]
                    slot[0] += 1
                    slot[1] += identity
                    slot[2] += coverage
                    slot[3] += is_indel

                    if key not in seq_cache:
                        b = str_pos.split('-')
                        nb = len(b)
                        ds_parts = []
                        if strand == '+':
                            for j in range(0, nb - 1, 2):
                                k_e = int(b[j]); k_s = int(b[j+1])
                                ds_parts.append(fa_fetch(contig, k_e, k_e+2) + '-' +
                                                fa_fetch(contig, k_s-3, k_s-1))
                            if nb % 2 == 1:
                                k_v = int(b[nb-1])
                                ds_parts.append(fa_fetch(contig, k_v, k_v+2))
                        else:
                            b_rev = b[::-1]; nb_rev = len(b_rev)
                            for j in range(0, nb_rev - 1, 2):
                                k_e = int(b_rev[j]); k_s = int(b_rev[j+1])
                                ds_parts.append(
                                    get_rc(fa_fetch(contig, k_e-3, k_e-1)) + '-' +
                                    get_rc(fa_fetch(contig, k_s, k_s+2)))
                            if nb_rev % 2 == 1:
                                k_v = int(b_rev[nb_rev-1])
                                ds_parts.append(get_rc(fa_fetch(contig, k_v-3, k_v-1)))
                        seq_cache[key] = ','.join(ds_parts)

    # o2：count \t chrom \t strand \t str_pos \t ds \t sum_identity \t sum_coverage \t indel_count
    with open(out2_tmp, 'w', buffering=WBUF) as o2_fh:
        for k, (cnt, sum_id, sum_cov, indel_c) in ec.items():
            o2_fh.write(f'{cnt}\t{k[0]}\t{k[1]}\t{k[2]}\t{seq_cache[k]}\t'
                        f'{sum_id:.6f}\t{sum_cov:.6f}\t{indel_c}\n')

    return out1_tmp, out2_tmp, bam


def _sendfile_copy(src_path, dst_fd):
    """os.sendfile 零拷贝追加（Linux）；非 Linux 回退到块读写。"""
    try:
        with open(src_path, 'rb') as src:
            src_fd = src.fileno()
            size = os.fstat(src_fd).st_size
            if size > 0:
                offset = 0
                while offset < size:
                    sent = os.sendfile(dst_fd, src_fd, offset, size - offset)
                    if sent == 0:
                        break
                    offset += sent
    except (AttributeError, OSError):
        CBUF = 8 * 1024 * 1024
        with open(src_path, 'rb') as src:
            while True:
                chunk = src.read(CBUF)
                if not chunk:
                    break
                os.write(dst_fd, chunk)


def _read_o2_tmp(o2_tmp):
    """读取单个 o2 临时文件，返回解析行列表；供线程池并行调用。"""
    rows = []
    with open(o2_tmp, 'r', buffering=8*1024*1024) as f:
        for line in f:
            p = line.rstrip('\n').split('\t', 7)
            rows.append(((p[1], p[2], p[3]),
                         int(p[0]), float(p[5]), float(p[6]), int(p[7]), p[4]))
    os.remove(o2_tmp)
    return rows


def merge_single_bam(bam, files, out_dir):
    from concurrent.futures import ThreadPoolExecutor
    bam_basename = os.path.splitext(os.path.basename(bam))[0]
    out1 = os.path.join(out_dir, f'{bam_basename}_flnc.ssc')
    out2 = os.path.join(out_dir, f'{bam_basename}_ssc.count')
    CBUF = 8 * 1024 * 1024

    # ── .ssc：os.sendfile 零拷贝，跳过用户态缓冲区 ─────────────────────────
    dst_fd = os.open(out1, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        for o1_tmp, _ in files:
            _sendfile_copy(o1_tmp, dst_fd)
            os.remove(o1_tmp)
    finally:
        os.close(dst_fd)

    # ── .count：线程池并行读取所有 o2_tmp（I/O 密集，线程比进程开销小）──────
    o2_tmps = [o2_tmp for _, o2_tmp in files]
    global_ec = defaultdict(lambda: [0, 0.0, 0.0, 0])
    global_seq = {}

    with ThreadPoolExecutor(max_workers=min(8, len(o2_tmps))) as ex:
        for rows in ex.map(_read_o2_tmp, o2_tmps):
            for key, cnt, sum_id, sum_cov, indel_c, ds in rows:
                slot = global_ec[key]
                slot[0] += cnt
                slot[1] += sum_id
                slot[2] += sum_cov
                slot[3] += indel_c
                global_seq[key] = ds

    # ── 排序 + 批量写入（writelines 比逐行 write 系统调用少）─────────────────
    lines = []
    for k in sorted(global_ec):
        cnt, sum_id, sum_cov, indel_c = global_ec[k]
        mean_id     = sum_id  / cnt if cnt else 0.0
        mean_cov    = sum_cov / cnt if cnt else 0.0
        indel_ratio = indel_c / cnt if cnt else 0.0
        lines.append(f'{cnt}\t{k[0]}\t{k[1]}\t{k[2]}\t{global_seq[k]}\t'
                     f'{mean_cov:.4f}\t{mean_id:.4f}\t{indel_ratio:.4f}\n')

    with open(out2, 'w', buffering=CBUF) as o2_fh:
        o2_fh.writelines(lines)


def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    temp_dir = args.tmpdir or tempfile.gettempdir()

    filter_secondary     = not args.keep_secondary
    filter_supplementary = not args.keep_supplementary
    logger.info(f"filter_secondary={filter_secondary}, "
                f"filter_supplementary={filter_supplementary}")

    logger.info("Getting contig info...")
    all_read_counts = {}
    all_lengths = {}
    for bam in args.bam:
        rc, ln = get_bam_contigs(bam)
        for c, n in rc.items():
            all_read_counts[c] = all_read_counts.get(c, 0) + n
        all_lengths.update(ln)

    total_reads = sum(all_read_counts.values())
    logger.info(f"Total mapped reads: {total_reads:,} across {len(all_read_counts)} contigs")
    if total_reads == 0:
        logger.warning("No mapped reads found.")
        sys.exit(0)

    target = args.target_reads if args.target_reads > 0 else total_reads // args.threads
    logger.info(f"Target reads/chunk: {target:,}")

    regions = build_regions(all_read_counts, all_lengths, target)
    logger.info(f"Built {len(regions)} regions from {len(all_read_counts)} contigs")

    n_bins = min(args.threads, len(regions))
    region_groups = bin_pack_regions(regions, n_bins)
    logger.info(f"Packed into {len(region_groups)} worker groups")

    tasks = []
    for bam in args.bam:
        for g_idx, rlist in enumerate(region_groups):
            tasks.append((bam, args.reference, temp_dir, g_idx, rlist,
                          filter_secondary, filter_supplementary))

    import time
    logger.info(f"Processing {len(tasks)} tasks on {args.threads} workers...")
    t0 = time.perf_counter()
    with mp.Pool(processes=min(args.threads, len(tasks))) as pool:
        results = pool.starmap(process_region_group, tasks)
    t1 = time.perf_counter()
    logger.info(f"Processing done in {t1-t0:.1f}s")

    logger.info("Merging...")
    groups = defaultdict(list)
    for o1, o2, bam_path in results:
        groups[bam_path].append((o1, o2))

    with mp.Pool(processes=min(len(groups), args.threads)) as pool:
        pool.starmap(partial(merge_single_bam, out_dir=args.output), groups.items())
    t2 = time.perf_counter()
    logger.info(f"Merging done in {t2-t1:.1f}s")

    logger.info("Done.")


if __name__ == '__main__':
    main()
