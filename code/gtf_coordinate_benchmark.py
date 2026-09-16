#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gtf_coordinate_benchmark.py

Independent post-run coordinate benchmark module.

This module is deliberately separate from annotation-free discovery and
annotation-rescue discovery. It reads an already generated GTF display track and
compares its coordinates against a user-specified coordinate reference file.
The reference is never used to generate candidates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd


def _normalize_chrom(chrom: str) -> str:
    value = str(chrom)
    return {'chrM': 'MT', 'M': 'MT', 'chrMT': 'MT'}.get(value, value)


def _parse_tolerances(value: str | Iterable[int]) -> List[int]:
    if isinstance(value, str):
        toks = [x.strip() for x in value.replace(';', ',').split(',') if x.strip()]
        vals = [int(x) for x in toks]
    else:
        vals = [int(x) for x in value]
    vals = sorted(set(vals))
    if not vals:
        raise ValueError('At least one benchmark tolerance must be provided.')
    if any(v < 0 for v in vals):
        raise ValueError('Benchmark tolerances must be non-negative integers.')
    return vals


def _clean_attr(attrs: str, key: str) -> str | None:
    token = f'{key} '
    for part in str(attrs).split(';'):
        part = part.strip()
        if part.startswith(token):
            val = part[len(token):].strip().strip('"')
            return val
    return None


def _detect_format(path: Path, requested: str = 'auto') -> str:
    requested = requested.lower()
    if requested in {'gtf', 'ssc'}:
        return requested
    suffix = path.suffix.lower()
    if suffix in {'.gtf', '.gff', '.gff3'}:
        return 'gtf'
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 9 and cols[2] in {'exon', 'transcript', 'gene', 'CDS'}:
                return 'gtf'
            return 'ssc'
    raise ValueError(f'Could not detect reference format from empty file: {path}')


def read_gtf_intervals(path: Path, id_prefix: str, feature_policy: str = 'exon_or_transcript') -> pd.DataFrame:
    rows = []
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = cols[:9]
            if feature_policy == 'exon_or_transcript' and feature not in {'exon', 'transcript'}:
                continue
            if strand not in {'+', '-'}:
                continue
            try:
                s = int(start); e = int(end)
            except ValueError:
                continue
            if s < 1 or e < s:
                continue
            rid = (_clean_attr(attrs, 'candidate_id') or _clean_attr(attrs, 'transcript_id') or
                   _clean_attr(attrs, 'gene_id') or f'{id_prefix}_{line_no:08d}')
            rows.append({
                'record_id': str(rid),
                'chrom': _normalize_chrom(chrom),
                'strand': str(strand),
                'start': s,
                'end': e,
                'source_file': path.name,
                'source_format': 'gtf',
                'source_feature': feature,
                'source_line': line_no,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=['record_id','chrom','strand','start','end','source_file','source_format','source_feature','source_line'])
    # Avoid double-counting transcript+exon rows with identical coordinates.
    df = df.drop_duplicates(subset=['chrom', 'strand', 'start', 'end', 'record_id']).reset_index(drop=True)
    return df


def _int_like(x: str) -> bool:
    try:
        int(str(x))
        return True
    except Exception:
        return False


def read_ssc_intervals(path: Path, id_prefix: str) -> pd.DataFrame:
    rows = []
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            rec_id = cols[0] if cols else f'{id_prefix}_{line_no:08d}'
            chrom = strand = None
            start = end = None
            # Common SSC2: read_id, chrom, strand, start, end, ...
            if len(cols) >= 5 and cols[2] in {'+', '-'} and _int_like(cols[3]) and _int_like(cols[4]):
                chrom, strand, start, end = _normalize_chrom(cols[1]), cols[2], int(cols[3]), int(cols[4])
            # Lexogen SIRV reference SSC: id, group, NA, chrom, strand, start, end, NA
            elif len(cols) >= 7 and cols[4] in {'+', '-'} and _int_like(cols[5]) and _int_like(cols[6]):
                chrom, strand, start, end = _normalize_chrom(cols[3]), cols[4], int(cols[5]), int(cols[6])
            else:
                # Generic fallback: find first strand column followed by two integer columns.
                for i, val in enumerate(cols[:-2]):
                    if val in {'+', '-'} and _int_like(cols[i+1]) and _int_like(cols[i+2]):
                        strand = val
                        chrom = _normalize_chrom(cols[i-1]) if i >= 1 else 'NA'
                        start, end = int(cols[i+1]), int(cols[i+2])
                        break
            if chrom is None or strand is None or start is None or end is None:
                continue
            if start < 1 or end < start:
                continue
            rows.append({
                'record_id': str(rec_id),
                'chrom': _normalize_chrom(chrom),
                'strand': str(strand),
                'start': int(start),
                'end': int(end),
                'source_file': path.name,
                'source_format': 'ssc',
                'source_feature': 'endpoint_pair',
                'source_line': line_no,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=['record_id','chrom','strand','start','end','source_file','source_format','source_feature','source_line'])
    df = df.drop_duplicates(subset=['chrom', 'strand', 'start', 'end', 'record_id']).reset_index(drop=True)
    return df


def read_coordinate_file(path: Path, fmt: str, id_prefix: str) -> pd.DataFrame:
    fmt2 = _detect_format(path, fmt)
    if fmt2 == 'gtf':
        return read_gtf_intervals(path, id_prefix=id_prefix)
    if fmt2 == 'ssc':
        return read_ssc_intervals(path, id_prefix=id_prefix)
    raise ValueError(f'Unsupported coordinate format: {fmt}')


def _greedy_match(pred: pd.DataFrame, ref: pd.DataFrame, tolerance: int) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if pred.empty or ref.empty:
        return pd.DataFrame(), pred.copy(), ref.copy()
    candidates = []
    for (chrom, strand), psub in pred.groupby(['chrom', 'strand'], sort=False):
        rsub = ref[(ref['chrom'].eq(chrom)) & (ref['strand'].eq(strand))]
        if rsub.empty:
            continue
        # Pre-compute reference arrays as int once per (chrom, strand) group
        r_start_arr = rsub['start'].to_numpy(dtype=int)
        r_end_arr = rsub['end'].to_numpy(dtype=int)
        r_record_ids = rsub['record_id'].to_numpy(dtype=str)
        r_indices = rsub.index.to_numpy(dtype=int)
        for pi, pr in psub.iterrows():
            pr_start = int(pr['start'])
            pr_end = int(pr['end'])
            ds = np.abs(r_start_arr - pr_start)
            de = np.abs(r_end_arr - pr_end)
            ok = (ds <= tolerance) & (de <= tolerance)
            if not ok.any():
                continue
            ok_idx = np.where(ok)[0]
            for ri in ok_idx:
                d_start = int(ds[ri])
                d_end = int(de[ri])
                max_delta = max(d_start, d_end)
                sum_delta = d_start + d_end
                candidates.append({
                    'pred_index': pi,
                    'ref_index': int(r_indices[ri]),
                    'pred_id': pr['record_id'],
                    'ref_id': r_record_ids[ri],
                    'delta_start': d_start,
                    'delta_end': d_end,
                    'max_abs_endpoint_delta': max_delta,
                    'sum_abs_endpoint_delta': sum_delta,
                })
    if not candidates:
        return pd.DataFrame(), pred.copy(), ref.copy()
    cand = pd.DataFrame(candidates)
    cand = cand.sort_values(['max_abs_endpoint_delta','sum_abs_endpoint_delta','pred_id','ref_id']).reset_index(drop=True)
    used_p = set(); used_r = set(); rows = []
    for _, row in cand.iterrows():
        pi = int(row['pred_index']); ri = int(row['ref_index'])
        if pi in used_p or ri in used_r:
            continue
        used_p.add(pi); used_r.add(ri)
        pr = pred.loc[pi]; rr = ref.loc[ri]
        rows.append({
            'tolerance_bp': tolerance,
            'pred_id': pr['record_id'],
            'ref_id': rr['record_id'],
            'chrom': pr['chrom'],
            'strand': pr['strand'],
            'pred_start': int(pr['start']),
            'pred_end': int(pr['end']),
            'ref_start': int(rr['start']),
            'ref_end': int(rr['end']),
            'delta_start': int(abs(int(pr['start']) - int(rr['start']))),
            'delta_end': int(abs(int(pr['end']) - int(rr['end']))),
            'max_abs_endpoint_delta': int(max(abs(int(pr['start']) - int(rr['start'])), abs(int(pr['end']) - int(rr['end'])))),
            'sum_abs_endpoint_delta': int(abs(int(pr['start']) - int(rr['start'])) + abs(int(pr['end']) - int(rr['end']))),
        })
    matched = pd.DataFrame(rows)
    unmatched_pred = pred.loc[[i for i in pred.index if i not in used_p]].copy()
    unmatched_ref = ref.loc[[i for i in ref.index if i not in used_r]].copy()
    return matched, unmatched_pred, unmatched_ref


def nearest_reference_table(pred: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    rows = []
    ref_chroms = set(ref['chrom'].astype(str)) if not ref.empty else set()
    # Pre-index references by (chrom, strand) with numpy arrays
    ref_by_key = {}
    if not ref.empty:
        for key, rsub in ref.groupby(['chrom', 'strand'], sort=False):
            rsub2 = rsub.copy()
            rsub2['start_int'] = rsub2['start'].to_numpy(dtype=int)
            rsub2['end_int'] = rsub2['end'].to_numpy(dtype=int)
            ref_by_key[key] = rsub2
    for _, pr in pred.iterrows():
        row = {
            'pred_id': pr['record_id'], 'pred_chrom': pr['chrom'], 'pred_strand': pr['strand'],
            'pred_start': int(pr['start']), 'pred_end': int(pr['end']),
            'nearest_ref_id': 'NA', 'nearest_ref_chrom': 'NA', 'nearest_ref_strand': 'NA',
            'nearest_ref_start': np.nan, 'nearest_ref_end': np.nan,
            'delta_start': np.nan, 'delta_end': np.nan, 'max_abs_endpoint_delta': np.nan,
            'nearest_status': 'off_reference_contig' if pr['chrom'] not in ref_chroms else 'no_same_strand_reference'
        }
        rsub = ref_by_key.get((pr['chrom'], pr['strand']))
        if rsub is not None and len(rsub):
            pr_start = int(pr['start'])
            pr_end = int(pr['end'])
            ds = np.abs(rsub['start_int'].to_numpy(dtype=int) - pr_start)
            de = np.abs(rsub['end_int'].to_numpy(dtype=int) - pr_end)
            max_delta = np.maximum(ds, de)
            sum_delta = ds + de
            best_idx = int(np.lexsort((sum_delta.astype(str) if hasattr(sum_delta, 'astype') else sum_delta, max_delta))[0])
            # Use numpy argmin with lexsort
            best_idx = int(np.argmin(max_delta * 1e9 + sum_delta))
            best = rsub.iloc[best_idx]
            row.update({
                'nearest_ref_id': best['record_id'], 'nearest_ref_chrom': best['chrom'], 'nearest_ref_strand': best['strand'],
                'nearest_ref_start': int(best['start']), 'nearest_ref_end': int(best['end']),
                'delta_start': int(ds[best_idx]), 'delta_end': int(de[best_idx]),
                'max_abs_endpoint_delta': int(max_delta[best_idx]),
                'nearest_status': 'same_chrom_same_strand',
            })
        rows.append(row)
    return pd.DataFrame(rows)


def add_benchmark_cli(ap: argparse.ArgumentParser) -> None:
    group = ap.add_argument_group('Optional independent post-run GTF coordinate benchmark')
    group.add_argument('--benchmark-enable', action='store_true',
                       help='Enable independent post-run benchmark after GTF generation. Reference is not used for discovery.')
    group.add_argument('--benchmark-reference', type=Path, default=None,
                       help='Coordinate reference file for benchmark. Supported: GTF or SSC/SSC2-like endpoint file.')
    group.add_argument('--benchmark-reference-format', choices=['auto','gtf','ssc'], default='auto',
                       help='Reference parser. Default auto detects GTF vs SSC-like text.')
    group.add_argument('--benchmark-prediction-gtf', type=Path, default=None,
                       help='Optional prediction GTF to compare. If omitted, the pipeline compares its default final output GTF.')
    group.add_argument('--benchmark-tolerances', default='0,5,10,20,30,50',
                       help='Comma-separated endpoint tolerances in bp. Default: 0,5,10,20,30,50.')
    group.add_argument('--benchmark-primary-tolerance', type=int, default=10,
                       help='Primary tolerance in bp for unmatched prediction/reference outputs. Default: 10.')
    group.add_argument('--benchmark-output-dir', type=Path, default=None,
                       help='Optional benchmark output directory. Default: <outdir>/results/benchmark.')


def run_gtf_coordinate_benchmark(
    prediction_gtf: Path,
    reference_file: Path,
    outdir: Path,
    reference_format: str = 'auto',
    tolerances: str | Iterable[int] = '0,5,10,20,30,50',
    primary_tolerance: int = 10,
) -> dict:
    prediction_gtf = Path(prediction_gtf)
    reference_file = Path(reference_file)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not prediction_gtf.exists():
        raise FileNotFoundError(f'Prediction GTF not found: {prediction_gtf}')
    if not reference_file.exists():
        raise FileNotFoundError(f'Benchmark reference not found: {reference_file}')
    tol_values = _parse_tolerances(tolerances)
    if primary_tolerance not in tol_values:
        tol_values = sorted(set(tol_values + [int(primary_tolerance)]))

    pred = read_gtf_intervals(prediction_gtf, id_prefix='pred')
    ref = read_coordinate_file(reference_file, fmt=reference_format, id_prefix='ref')
    pred.to_csv(outdir / 'prediction_intervals.parsed.tsv', sep='\t', index=False)
    ref.to_csv(outdir / 'reference_intervals.parsed.tsv', sep='\t', index=False)

    all_metrics = []
    all_matches = []
    primary_unmatched_pred = pd.DataFrame()
    primary_unmatched_ref = pd.DataFrame()
    ref_chroms = set(ref['chrom'].astype(str)) if not ref.empty else set()
    off_ref_count = int((~pred['chrom'].astype(str).isin(ref_chroms)).sum()) if not pred.empty else 0
    same_ref_contig_count = int((pred['chrom'].astype(str).isin(ref_chroms)).sum()) if not pred.empty else 0

    for tol in tol_values:
        matched, unmatched_pred, unmatched_ref = _greedy_match(pred, ref, tol)
        tp = int(len(matched)); fp = int(len(unmatched_pred)); fn = int(len(unmatched_ref))
        precision = tp / (tp + fp) if (tp + fp) else np.nan
        recall = tp / (tp + fn) if (tp + fn) else np.nan
        f1 = (2 * precision * recall / (precision + recall)) if precision == precision and recall == recall and (precision + recall) else np.nan
        all_metrics.append({
            'prediction_gtf': prediction_gtf.name,
            'reference_file': reference_file.name,
            'reference_format': _detect_format(reference_file, reference_format),
            'tolerance_bp': int(tol),
            'predicted_gtf_count': int(len(pred)),
            'reference_count': int(len(ref)),
            'TP': tp, 'FP': fp, 'FN': fn,
            'precision': precision,
            'recall': recall,
            'F1': f1,
            'off_reference_candidate_count': off_ref_count,
            'same_reference_contig_candidate_count': same_ref_contig_count,
        })
        if not matched.empty:
            all_matches.append(matched)
        if tol == primary_tolerance:
            primary_unmatched_pred = unmatched_pred.copy()
            primary_unmatched_ref = unmatched_ref.copy()

    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(outdir / 'benchmark_metrics.tsv', sep='\t', index=False)
    matches_df = pd.concat(all_matches, ignore_index=True) if all_matches else pd.DataFrame(columns=['tolerance_bp','pred_id','ref_id','chrom','strand','pred_start','pred_end','ref_start','ref_end','delta_start','delta_end','max_abs_endpoint_delta','sum_abs_endpoint_delta'])
    matches_df.to_csv(outdir / 'matched_pairs_by_tolerance.tsv', sep='\t', index=False)
    primary_unmatched_pred.to_csv(outdir / f'unmatched_predictions_tolerance_{primary_tolerance}bp.tsv', sep='\t', index=False)
    primary_unmatched_ref.to_csv(outdir / f'unmatched_references_tolerance_{primary_tolerance}bp.tsv', sep='\t', index=False)
    nearest_reference_table(pred, ref).to_csv(outdir / 'nearest_reference_for_each_prediction.tsv', sep='\t', index=False)

    params = {
        'module': 'gtf_coordinate_benchmark',
        'prediction_gtf': str(prediction_gtf),
        'reference_file': str(reference_file),
        'reference_format_requested': reference_format,
        'reference_format_detected': _detect_format(reference_file, reference_format),
        'tolerances_bp': ','.join(map(str, tol_values)),
        'primary_tolerance_bp': int(primary_tolerance),
        'matching_rule': 'same_chrom_same_strand_start_and_end_within_tolerance_one_to_one_greedy',
        'greedy_priority': 'max_abs_endpoint_delta;sum_abs_endpoint_delta;pred_id;ref_id',
        'reference_used_for_discovery': 'false',
    }
    pd.DataFrame([{'parameter': k, 'value': v} for k, v in params.items()]).to_csv(outdir / 'benchmark_run_parameters.tsv', sep='\t', index=False)
    (outdir / 'benchmark_run_parameters.json').write_text(json.dumps(params, indent=2, ensure_ascii=False, sort_keys=True), encoding='utf-8')

    primary_row = metrics_df[metrics_df['tolerance_bp'].eq(primary_tolerance)]
    result = primary_row.iloc[0].to_dict() if not primary_row.empty else metrics_df.iloc[0].to_dict()
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Standalone coordinate benchmark for a generated GTF against a reference GTF/SSC file.')
    ap.add_argument('--prediction-gtf', type=Path, required=True)
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--reference-format', choices=['auto','gtf','ssc'], default='auto')
    ap.add_argument('--tolerances', default='0,5,10,20,30,50')
    ap.add_argument('--primary-tolerance', type=int, default=10)
    args = ap.parse_args()
    res = run_gtf_coordinate_benchmark(args.prediction_gtf, args.reference, args.outdir, args.reference_format, args.tolerances, args.primary_tolerance)
    print(pd.DataFrame([res]).to_string(index=False))
