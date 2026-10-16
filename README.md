# Developmental olfactory epithelium priming analysis

This repository contains the public R analysis workflow for the manuscript
*Progressive Fate Restriction of Ascl1-Expressing Stem Cells in the Olfactory
Epithelial Niche by Intrinsic and Extrinsic Factors*.

The numbered R Markdown notebooks, their shared helper functions, and the fixed
plotting assets are all contained directly in this repository. Large Seurat,
Slingshot, tradeSeq, and CellChat objects are deposited separately on Zenodo
and are intentionally not stored in GitHub.

## Where to start

**Most users should start with notebook 03, not notebook 01.**

Notebooks 01 and 02 document the original raw-count merge, quality control,
and integration workflow. They are included for methodological provenance, but
they are not the recommended entry point for reproducing the public analysis.
The Zenodo deposit provides `WholeOE.rds`, the completed output of notebook 02,
so a normal rerun begins with notebook 03.

The notebooks do not all have to be run from beginning to end. Zenodo contains
canonical checkpoints that let each major analysis begin from its tested input:

| Goal | Start with | Download from Zenodo |
|---|---|---|
| Reproduce lineage extraction and annotation | `03_extract_olfactory_lineage.Rmd` | `WholeOE.rds` |
| Reproduce the GBC subset analysis | `04_subset_and_cluster_gbcs.Rmd` | `olfactory_lineage.rds` |
| Reproduce the HBC subset analysis | `05_subset_and_cluster_hbcs.Rmd` | `olfactory_lineage.rds` |
| Reproduce lineage cleaning | `06_clean_lineage_for_slingshot.Rmd` | `olfactory_lineage.rds` |
| Rerun Slingshot and the rotated trajectory plot | `07_run_slingshot_and_plot.Rmd` | `cleaned_for_slingshot_60265.rds` |
| Reproduce the Figure 4 TF heatmap and tradeSeq plots | `08_plot_precomputed_tradeseq.Rmd` | `gbc_subset.rds`, `slingshot_output.rds`, `tradeseq_fit_combined.rds`, and `ALL_TF.Rda` |
| Reproduce the p5-versus-p22 CellChat analysis | `09_cellchat_gbc_incoming_p5_p22.Rmd` | `cellchat_input.RData` |

Notebook 08 plots the deposited tradeSeq fit; it does **not** rerun the
computationally expensive `fitGAM()` step. Notebook 09 is an independent
CellChat analysis and does not require notebooks 03--08.

## Zenodo files

Zenodo DOI: **PENDING**

Download the following files from Zenodo and place them in a directory named
`data` at the repository root.

| File | Purpose |
|---|---|
| `WholeOE.rds` | QC-filtered and integrated whole olfactory epithelium; starting point for notebook 03 |
| `olfactory_lineage.rds` | Finished lineage object; skips the slow lineage extraction and marker work |
| `gbc_subset.rds` | Finished GBC subset; required directly by notebook 08 |
| `hbc_subset.rds` | Finished HBC subset for users interested in the HBC analysis |
| `cleaned_for_slingshot_60265.rds` | Ready-to-run Slingshot input |
| `slingshot_output.rds` | Completed Slingshot result; required directly by notebook 08 |
| `tradeseq_fit_combined.rds` | Completed tradeSeq GAM fit |
| `ALL_TF.Rda` | TF reference required by notebook 08 |
| `cellchat_input.RData` | Completed per-age CellChat objects used by notebook 09 |

After downloading, the relevant part of the repository should look like this:

```text
Dev_Priming_JNeuro_2026/
|-- data/
|   |-- WholeOE.rds
|   |-- olfactory_lineage.rds
|   |-- gbc_subset.rds
|   |-- hbc_subset.rds
|   |-- cleaned_for_slingshot_60265.rds
|   |-- slingshot_output.rds
|   |-- tradeseq_fit_combined.rds
|   |-- ALL_TF.Rda
|   `-- cellchat_input.RData
|-- assets/
|-- 01_merge_raw_h5.Rmd
|-- 02_qc_and_integration.Rmd
|-- 03_extract_olfactory_lineage.Rmd
|-- 04_subset_and_cluster_gbcs.Rmd
|-- 05_subset_and_cluster_hbcs.Rmd
|-- 06_clean_lineage_for_slingshot.Rmd
|-- 07_run_slingshot_and_plot.Rmd
|-- 08_plot_precomputed_tradeseq.Rmd
|-- 09_cellchat_gbc_incoming_p5_p22.Rmd
|-- publication_helpers.R
|-- DESCRIPTION
|-- renv.lock
`-- README.md
```

The serialized data objects are ignored by Git and should remain in Zenodo,
not be committed to this repository.

## Installation

The validated environment uses R 4.4.2. From the repository root, restore the
pinned package environment recorded in `renv.lock`:

```r
install.packages("renv")
renv::restore()
```

Direct package dependencies are also summarized in [`DESCRIPTION`](DESCRIPTION).

## Running the notebooks

Open the repository as the working project and render the desired notebook from
the repository root. The default input directory is `data`, and generated
Seurat, Slingshot, and tradeSeq results are written under `results`.

For example, start the main public workflow from the deposited whole-epithelium
object with:

```r
rmarkdown::render(
  "03_extract_olfactory_lineage.Rmd",
  params = list(
    download_dir = "data",
    output_dir = "results",
    save_outputs = TRUE
  )
)
```

Run an individual downstream notebook the same way. For example, the deposited
cleaned lineage object can be passed directly to Slingshot:

```r
rmarkdown::render(
  "07_run_slingshot_and_plot.Rmd",
  params = list(
    download_dir = "data",
    output_dir = "results",
    save_outputs = TRUE
  )
)
```

Reproduce the tradeSeq-based plots without refitting the GAM:

```r
rmarkdown::render(
  "08_plot_precomputed_tradeseq.Rmd",
  params = list(
    download_dir = "data",
    output_dir = "results",
    save_outputs = TRUE
  )
)
```

Run the independent CellChat analysis with:

```r
rmarkdown::render(
  "09_cellchat_gbc_incoming_p5_p22.Rmd",
  params = list(
    data_dir = "data",
    save_outputs = TRUE,
    receiver = "GBC",
    early_age = "p5",
    late_age = "p22",
    top_n_senders = 6,
    top_n_pathways = 8,
    use_exported_tables = FALSE
  )
)
```

The notebooks may also be opened and knitted interactively in RStudio. Keep the
repository root as the project root so `here::here()` resolves paths correctly.

## Notebook guide

| Notebook | Role | Main output |
|---|---|---|
| `01_merge_raw_h5.Rmd` | Documents merging the original developmental 10x H5 folders | `001_Dev_merged.rds` |
| `02_qc_and_integration.Rmd` | Documents QC, integration, clustering, and whole-epithelium preparation | `WholeOE.rds` |
| `03_extract_olfactory_lineage.Rmd` | Extracts and annotates the olfactory lineage | `olfactory_lineage.rds` |
| `04_subset_and_cluster_gbcs.Rmd` | Subsets and analyzes Ascl1-positive GBCs | `gbc_subset.rds` |
| `05_subset_and_cluster_hbcs.Rmd` | Subsets and analyzes HBCs | `hbc_subset.rds` |
| `06_clean_lineage_for_slingshot.Rmd` | Removes off-lineage populations and prepares the trajectory input | `cleaned_for_slingshot_60265.rds` |
| `07_run_slingshot_and_plot.Rmd` | Runs Slingshot and recreates the fixed rotated trajectory view | `slingshot_output.rds` and trajectory plot |
| `08_plot_precomputed_tradeseq.Rmd` | Recreates the TF heatmap and plots the deposited tradeSeq fit | Figure 4 plots and marker table |
| `09_cellchat_gbc_incoming_p5_p22.Rmd` | Compares incoming GBC signaling at p5 and p22 | CellChat tables and figures |

## Included fixed plotting assets

Notebook 07 requires a fixed 3D coordinate set, camera state, and saved curve
embedding to reproduce the published rotated trajectory view. These small files
are versioned with the code under [`assets`](assets):

- `04_slingshot_rgl_coordinates.rds`
- `04_slingshot_rgl_view_state.rds`
- `04_slingshot_umap_curves_pcs1_to5_april13.rds`

These are plotting assets, not substitute biological datasets. Project-defined
helper functions are colocated in
[`publication_helpers.R`](publication_helpers.R).

## Repository contents

| Path | Description |
|---|---|
| `01_*.Rmd` through `09_*.Rmd` | Numbered public analysis notebooks |
| `publication_helpers.R` | Shared project-defined helper functions |
| `assets/` | Fixed coordinates, camera state, and curve embedding for notebook 07 |
| `DATA_MANIFEST.csv` | Deposited file sizes, notebook usage, rerun tier, and Zenodo DOI field |
| `DESCRIPTION` | Direct R package dependencies |
| `renv.lock` | Pinned R package environment |
| `CITATION.cff` | Software citation metadata |
| `authors.json` | Contributor metadata used for the Zenodo record |
| `LICENSE` | Repository license |

## Reproducibility notes

The deposited checkpoints are the canonical starting objects used for the
public reruns. UMAP coordinates, clustering, and other stochastic results can
vary with operating system, package builds, thread counts, and numerical
libraries. The fixed notebook 07 projection assets are included specifically to
preserve the published camera view and curve projection.

Exact Zenodo file sizes and the DOI field are recorded in
[`DATA_MANIFEST.csv`](DATA_MANIFEST.csv). Replace `PENDING` with the final DOI
after the Zenodo record is minted.

## Citation and license

Please cite the associated manuscript and software release. Citation metadata
is provided in [`CITATION.cff`](CITATION.cff). See [`LICENSE`](LICENSE) for reuse
terms.
