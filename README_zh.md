# mtfind 0.311

mtfind 根据长读长 RNA 测序中同一条 read 的 5′ 和 3′ 端点，发现线粒体 RNA 坐标候选。它以比对后生成的 SSC 文件为输入，候选发现阶段不读取注释。

输出的 GTF 是坐标候选展示轨道，不是经过验证的转录本注释。

English version: [README.md](README.md).

## 安装

需要 Python 3.10 或更高版本。

```bash
python3 -m pip install -r requirements.txt
```

## 一键示例

本仓库提供一个来自公开 SG-NEx A549 direct-RNA 样本的 MT-only SSC 派生输入。运行：

```bash
bash examples/sgnex_a549/run_demo.sh
```

脚本复用 mtfind 0.311 对该 A549 输入已记录的参数，并在 `examples/sgnex_a549/output/` 生成结果；若该目录已存在，脚本会拒绝覆盖。输出包括候选 GTF、审计表、运行参数和输入/代码副本。输入来源和引用要求见 [DATA_PROVENANCE_zh.md](examples/sgnex_a549/DATA_PROVENANCE_zh.md)。

除高置信度 GTF 外，`results/tables/annotation_free_supercluster_table.tsv` 还记录全部发现的坐标候选，包括未通过展示阈值的候选。这些记录会保留供使用者复核；未出现在 GTF 中不表示该坐标候选被删除，或在生物学上已被否定。

## 运行自己的 SSC 输入

指定 SSC 路径、新的输出目录和样本名即可运行；其余参数均使用 mtfind 0.311 默认值（`MT`、16,569 bp 和 `full` 模式）。

```bash
python3 code/annotation_free_adaptive_depth.py \
  --ssc examples/sgnex_a549/SGNex_A549_directRNA_replicate6_run1.doradov300.sup.mt.sorted_flnc.ssc \
  --outdir results/a549_default \
  --sample-id SGNex_A549_directRNA_replicate6_run1
```

如只需调整少数参数，可在同一命令中补充对应选项。例如：

```bash
python3 code/annotation_free_adaptive_depth.py \
  --ssc examples/sgnex_a549/SGNex_A549_directRNA_replicate6_run1.doradov300.sup.mt.sorted_flnc.ssc \
  --outdir results/a549_custom \
  --sample-id SGNex_A549_directRNA_replicate6_run1 \
  --run-mode data \
  --adaptive-depth-mode fixed \
  --density-min-raw-pair-support 25 \
  --display-min-cluster-support 1000 \
  --display-score-threshold 0.55
```

`--display-score-threshold` 默认值为 `0.40`，用于设定候选展示评分的最低阈值。提高该值会使筛选更保守，通常减少输出候选；降低该值会放宽筛选，通常增加输出候选。实际有效阈值还受样本内相对评分下限约束，因此调低后不保证每个样本都会增加候选。

## 可选：从 BAM 制作 SSC

已有 SSC 时可跳过此步骤。`code/bam2ssc.py` 是独立的上游转换脚本；除上述 mtfind 依赖外，还需要在同一 Python 环境安装 `pysam`，并确保 `samtools` 位于 `PATH`。输入应为按坐标排序且已建索引的 BAM，以及带 `.fai` 索引的对应参考 FASTA；BAM 中的 contig 名称须与 FASTA 序列名一致。安装和建索引方法参见 [pysam](https://pysam.readthedocs.io/en/latest/) 与 [samtools](https://www.htslib.org/doc/) 官方文档。

在仓库根目录运行，替换为自己的路径，并使用新的输出目录：

```bash
python3 code/bam2ssc.py \
  --reference /path/to/MT.fa \
  --bam /path/to/sample.mt.sorted.bam \
  --output results/ssc_sample \
  --threads 4

python3 code/annotation_free_adaptive_depth.py \
  --ssc results/ssc_sample/sample.mt.sorted_flnc.ssc \
  --outdir results/mtfind_sample \
  --sample-id sample
```

每个 BAM 会生成 `*_flnc.ssc`（逐 read 的端点记录）和 `*_ssc.count`（汇总计数）；默认排除 secondary 和 supplementary 比对。mtfind 只读取 `*_flnc.ssc`，不需要 BAM、参考 FASTA、索引或计数文件。转换脚本会**覆盖已有目录中的同名输出**，因此 `--output` 应指向新目录。mtfind 默认的 `MT` 与 16,569 bp 适用于人类线粒体；其他参考序列需相应设置 `--chrom` 和 `--sequence-length`。

## 内容

- `code/annotation_free_adaptive_depth.py`：mtfind 0.311 候选发现主程序。
- `code/bam2ssc.py`：可选的 BAM 转 SSC 脚本。
- `code/gtf_coordinate_benchmark.py`：可选的运行后坐标比较模块，按实际执行版本原样保留。
- `tests/test_seed_backfill.py`：seed/backfill 聚类的小型回归测试。
- `examples/sgnex_a549/`：一个真实数据 SSC 示例及其准确运行命令。

## 快速检查

```bash
python3 code/annotation_free_adaptive_depth.py --version
python3 -m unittest discover -s tests -v
```
