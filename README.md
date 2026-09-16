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

## Contents

- `code/annotation_free_adaptive_depth.py`: mtfind 0.311 discovery program.
- `code/gtf_coordinate_benchmark.py`: optional post-run coordinate-comparison module retained unchanged from the executed release.
- `tests/test_seed_backfill.py`: small regression test for seed/backfill clustering.
- `examples/sgnex_a549/`: one real-data SSC example and its exact command.

## Quick check

```bash
python3 code/annotation_free_adaptive_depth.py --version
python3 -m unittest discover -s tests -v
```

No licence file is included here; select the repository licence on GitHub before publication.
