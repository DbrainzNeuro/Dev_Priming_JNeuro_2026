# Historical tradeSeq fit

This directory contains the code used to generate the deposited
`tradeseq_fit_combined.rds` object.

## Provenance guarantee

`005_run_tradeseq_biowulf.R` is a byte-for-byte copy of the executed analysis
script. Its SHA-256 checksum is:

```text
4C64E709A1E6CB7BD685BD9FA69951ED6BDC129A585D3E246CBCAEF957041A11
```

No modeling, filtering, chunking, testing, or output code was rewritten. The
public Slurm launcher is also copied from the executed launcher, with changes
limited to scheduler log locations and repository-relative input, script, and
output paths. Personal filesystem paths and account names were removed.

The fitted model used these settings:

| Setting | Value |
|---|---:|
| PCs represented by the input Slingshot object | 1--5 |
| Cluster labels | fine clusters |
| Genes per chunk | 250 |
| Workers | 8 |
| tradeSeq knots | 6 |
| Minimum expressing cells | 10 |
| Random seed | 12345 |

## Included gene lists

The executed job fitted the union of:

- `gene_lists/ALL_TF.rds`: 2,576 transcription factors.
- `gene_lists/HVG_minus_TF.rds`: 1,767 non-TF highly variable genes.

The lists are disjoint and contain 4,343 genes in total. `ALL_TF.rds` contains
the same ordered character vector as the deposited `ALL_TF.Rda`; the `.rds`
copy is retained because the unmodified execution script reads gene lists with
`readRDS()`.

## Slurm execution

From the repository root, place `slingshot_output.rds` in `data/`, create the
log directory, and submit:

```bash
mkdir -p logs
sbatch tradeseq_fit/tradeseq_pc_1_to_5_fineClusters_k6_april13.slurm
```

The launcher writes chunk models, combined statistical tests, manifests,
session information, and `tradeseq_fit_combined.rds` under
`results/tradeseq_fit/`.

The deposited `slingshot_output.rds` is the public filename assigned to the
Slingshot SCE input role. Exact numerical reproduction also depends on using
the deposited input object and the package versions recorded in `renv.lock`.

## Direct execution without Slurm

On a Unix-like system with forked BiocParallel support, the same R script can
be invoked directly:

```bash
Rscript tradeseq_fit/005_run_tradeseq_biowulf.R \
  --sce_rds data/slingshot_output.rds \
  --gene_list_rds tradeseq_fit/gene_lists/ALL_TF.rds,tradeseq_fit/gene_lists/HVG_minus_TF.rds \
  --outdir results/tradeseq_fit \
  --final_outdir results/tradeseq_fit \
  --chunk_size 250 \
  --workers 8 \
  --nknots 6 \
  --min_cells_expressed 10 \
  --seed 12345
```

The full fit is memory- and compute-intensive. Notebook 08 uses the deposited
completed fit and does not require users to rerun this step.
