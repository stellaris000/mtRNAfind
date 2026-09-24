# mtfind 0.311

mtfind discovers mitochondrial RNA coordinate candidates from same-read 5′ and 3′ endpoints in long-read RNA-seq. It operates on an alignment-derived SSC file and does not use an annotation during candidate discovery.

The generated GTF is a coordinate-candidate display track, not a validated transcript annotation.

中文说明见 [README_zh.md](README_zh.md).

## Install

Python 3.10 or later is required.

```bash
python3 -m pip install -r requirements.txt
```

## One-command example

The repository includes one MT-only SSC derived from the public SG-NEx A549 direct-RNA sample. Run:

```bash
bash examples/sgnex_a549/run_demo.sh
```

The script uses the recorded mtfind 0.311 A549 parameters and creates `examples/sgnex_a549/output/`. It does not overwrite an existing output directory. The output includes a candidate GTF, audit tables, run parameters and a copied input/code record. The example input and citation requirement are documented in [DATA_PROVENANCE.md](examples/sgnex_a549/DATA_PROVENANCE.md).

In addition to the high-confidence GTF, `results/tables/annotation_free_supercluster_table.tsv` records all discovered coordinate candidates, including candidates that did not pass the display threshold. Those records are retained for user review; not appearing in the GTF does not mean that a coordinate candidate was deleted or biologically disproven.

## Run your own SSC input

Provide an SSC path, a new output directory and a sample name. All other parameters use mtfind 0.311 defaults (`MT`, 16,569 bp and `full` mode).

```bash
python3 code/annotation_free_adaptive_depth.py \
  --ssc examples/sgnex_a549/SGNex_A549_directRNA_replicate6_run1.doradov300.sup.mt.sorted_flnc.ssc \
  --outdir results/a549_default \
  --sample-id SGNex_A549_directRNA_replicate6_run1
```

To override only selected parameters, add them to the same command. For example:

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

`--display-score-threshold` defaults to `0.40` and sets the minimum candidate display score. Raising it makes output more conservative and usually reduces the number of candidates; lowering it relaxes the filter and usually increases the number of candidates. The effective threshold also has a sample-relative lower bound, so lowering this value does not guarantee additional candidates for every sample.

## Optional: prepare an SSC from BAM

If you already have an SSC, skip this step. `code/bam2ssc.py` is a separate upstream converter. It needs `pysam` in the same Python environment and `samtools` on `PATH`, in addition to the mtfind requirements above. Use a coordinate-sorted, indexed BAM and a matching reference FASTA with a `.fai` index; the BAM contig name must match the FASTA sequence name. See the official [pysam](https://pysam.readthedocs.io/en/latest/) and [samtools](https://www.htslib.org/doc/) documentation for installation and indexing.

From the repository root, replace these paths with your own and choose a new output directory:

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

For each BAM, the converter writes `*_flnc.ssc` (per-read endpoints) and `*_ssc.count` (summary counts). By default it excludes secondary and supplementary alignments. mtfind reads only `*_flnc.ssc`; it does not require the BAM, reference FASTA, indexes, or count file. The converter **overwrites same-named outputs** in an existing directory, so use a fresh `--output` path. The mtfind defaults `MT` and 16,569 bp are for human mitochondrial input; set `--chrom` and `--sequence-length` appropriately for other references.

## Contents

- `code/annotation_free_adaptive_depth.py`: mtfind 0.311 discovery program.
- `code/bam2ssc.py`: optional BAM-to-SSC converter.
- `code/gtf_coordinate_benchmark.py`: optional post-run coordinate-comparison module retained unchanged from the executed release.
- `tests/test_seed_backfill.py`: small regression test for seed/backfill clustering.
- `examples/sgnex_a549/`: one real-data SSC example and its exact command.

## Quick check

```bash
python3 code/annotation_free_adaptive_depth.py --version
python3 -m unittest discover -s tests -v
```
