suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(BiocParallel)
})

parse_args <- function(args) {
  out <- list(
    sce_rds = NULL,
    gene_list_rds = NULL,
    outdir = NULL,
    final_outdir = NULL,
    chunk_size = 250L,
    workers = 8L,
    nknots = 6L,
    min_cells_expressed = 10L,
    seed = 12345L,
    verbose = TRUE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    if (identical(name, "quiet")) {
      out$verbose <- FALSE
      i <- i + 1L
      next
    }
    if (i == length(args)) stop("Missing value for argument: ", key)
    value <- args[[i + 1L]]
    if (name %in% c("chunk_size", "workers", "nknots", "min_cells_expressed", "seed")) {
      value <- as.integer(value)
    }
    out[[name]] <- value
    i <- i + 2L
  }
  out
}

log_msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), ..., "\n", sep = "")
}

read_one_gene_list <- function(path) {
  x <- readRDS(path)
  if (is.character(x)) return(unique(x))
  if (is.factor(x)) return(unique(as.character(x)))
  if (is.data.frame(x)) {
    hit <- c("gene", "genes", "symbol", "symbols", "Gene", "Symbol")
    hit <- hit[hit %in% colnames(x)][1]
    if (!is.na(hit)) return(unique(as.character(x[[hit]])))
    return(unique(as.character(x[[1]])))
  }
  if (is.list(x)) {
    flat <- unlist(x, use.names = FALSE)
    if (is.character(flat)) return(unique(flat))
  }
  stop("Unsupported gene list object: ", path)
}

read_gene_list <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  
  paths <- trimws(strsplit(path, ",", fixed = TRUE)[[1]])
  paths <- paths[nzchar(paths)]
  if (!length(paths)) return(NULL)
  
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths)) {
    stop("Missing gene list file(s): ", paste(missing_paths, collapse = ", "))
  }
  
  unique(unlist(lapply(paths, read_one_gene_list), use.names = FALSE))
}

write_csv_gz <- function(df, path, row.names = FALSE) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(df, con, row.names = row.names)
}

combine_fitted_chunks <- function(fitted_chunks) {
  fitted_chunks_clean <- lapply(fitted_chunks, function(chunk) {
    chunk_copy <- chunk
    if ("tradeSeq" %in% colnames(colData(chunk_copy))) {
      colData(chunk_copy)$tradeSeq <- NULL
    }
    if ("crv" %in% colnames(colData(chunk_copy))) {
      colData(chunk_copy)$crv <- NULL
    }
    chunk_copy
  })

  sce_combined <- do.call(rbind, fitted_chunks_clean)
  colData(sce_combined)$tradeSeq <- colData(fitted_chunks[[1]])$tradeSeq
  if ("crv" %in% colnames(colData(fitted_chunks[[1]]))) {
    colData(sce_combined)$crv <- colData(fitted_chunks[[1]])$crv
  }
  sce_combined
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot(!is.null(args$sce_rds), file.exists(args$sce_rds))
stopifnot(!is.null(args$outdir), nzchar(args$outdir))

final_outdir <- if (!is.null(args$final_outdir) && nzchar(args$final_outdir)) {
  args$final_outdir
} else {
  args$outdir
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
chunk_dir <- file.path(args$outdir, "chunk_models")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

test_dir <- args$outdir
combined_fit_rds <- file.path(args$outdir, "tradeseq_fit_combined.rds")
status_csv <- file.path(args$outdir, "test_status.csv")

set.seed(args$seed)
log_msg("Loading ", args$sce_rds, verbose = args$verbose)
sce <- readRDS(args$sce_rds)

if (!"counts" %in% assayNames(sce)) {
  stop(
    "tradeSeq requires a raw 'counts' assay. Found assays: ",
    paste(assayNames(sce), collapse = ", ")
  )
}
assay_name <- "counts"
counts_mat <- assay(sce, "counts")

if (any(counts_mat < 0, na.rm = TRUE)) {
  stop("counts assay contains negative values.")
}
if (!all(abs(counts_mat - round(counts_mat)) < 1e-8)) {
  warning("counts assay does not look integer-valued; verify this is raw counts.")
}
if (any(Matrix::colSums(counts_mat) == 0)) {
  warning("Some cells have zero total counts in the counts assay.")
}
pt <- slingPseudotime(sce)
cw <- slingCurveWeights(sce)
if (is.null(pt) || is.null(cw)) stop("Slingshot pseudotime/weights missing.")

pt <- as.matrix(pt)
cw <- as.matrix(cw)
cells <- colnames(sce)
pt <- pt[cells, , drop = FALSE]
cw <- cw[cells, , drop = FALSE]

w_sum <- rowSums(cw, na.rm = TRUE)
keep_cells <- w_sum > 0
log_msg("Cells with weight-sum==0: ", sum(!keep_cells), verbose = args$verbose)
if (any(!keep_cells)) {
  sce <- sce[, keep_cells]
  counts_mat <- counts_mat[, keep_cells, drop = FALSE]
  pt <- pt[keep_cells, , drop = FALSE]
  cw <- cw[keep_cells, , drop = FALSE]
}

na_pt <- is.na(pt)
zero_w <- (cw == 0) | is.na(cw)
pt[na_pt & zero_w] <- 0
bad_na <- is.na(pt) & !zero_w
bad_cells <- which(rowSums(bad_na, na.rm = TRUE) > 0)
log_msg("Cells with NA pseudotime in weighted lineages: ", length(bad_cells), verbose = args$verbose)
if (length(bad_cells) > 0) {
  keep2 <- rep(TRUE, nrow(pt))
  keep2[bad_cells] <- FALSE
  sce <- sce[, keep2]
  counts_mat <- counts_mat[, keep2, drop = FALSE]
  pt <- pt[keep2, , drop = FALSE]
  cw <- cw[keep2, , drop = FALSE]
}

stopifnot(!any(is.na(pt)))
stopifnot(ncol(counts_mat) == nrow(pt))
stopifnot(ncol(counts_mat) == nrow(cw))
stopifnot(identical(colnames(counts_mat), rownames(pt)))
stopifnot(identical(colnames(counts_mat), rownames(cw)))

gene_list <- read_gene_list(args$gene_list_rds)
genes_use <- if (is.null(gene_list)) rownames(counts_mat) else intersect(rownames(counts_mat), gene_list)
keep_genes <- Matrix::rowSums(counts_mat[genes_use, , drop = FALSE] > 0) >= args$min_cells_expressed
genes_use <- genes_use[keep_genes]
if (!length(genes_use)) stop("No genes left after filtering.")

chunks <- split(seq_along(genes_use), ceiling(seq_along(genes_use) / args$chunk_size))
bp <- BiocParallel::MulticoreParam(workers = args$workers, progressbar = args$verbose)

chunk_manifest <- vector("list", length(chunks))
fitted_chunks <- vector("list", length(chunks))

for (chunk_id in seq_along(chunks)) {
  chunk_genes <- genes_use[chunks[[chunk_id]]]
  log_msg("Chunk ", chunk_id, "/", length(chunks), " genes=", length(chunk_genes), verbose = args$verbose)

  sce_chunk <- sce[chunk_genes, , drop = FALSE]
  fit <- tradeSeq::fitGAM(
    counts = assay(sce_chunk, assay_name),
    pseudotime = pt,
    cellWeights = cw,
    nknots = args$nknots,
    verbose = args$verbose,
    parallel = TRUE,
    BPPARAM = bp
  )

  fit_rds <- file.path(chunk_dir, sprintf("fitGAM_chunk_%03d.rds", chunk_id))
  saveRDS(fit, fit_rds)
  fitted_chunks[[chunk_id]] <- fit

  chunk_manifest[[chunk_id]] <- data.frame(
    chunk_id = chunk_id,
    n_genes = length(chunk_genes),
    first_gene = chunk_genes[[1]],
    last_gene = chunk_genes[[length(chunk_genes)]],
    fit_rds = fit_rds,
    final_fit_rds = file.path(final_outdir, "chunk_models", basename(fit_rds)),
    stringsAsFactors = FALSE
  )

  rm(sce_chunk, fit)
  gc()
}

chunk_manifest_df <- do.call(rbind, chunk_manifest)
utils::write.csv(chunk_manifest_df, file.path(args$outdir, "chunk_manifest.csv"), row.names = FALSE)

log_msg("Combining fitted chunks", verbose = args$verbose)
sce_combined <- combine_fitted_chunks(fitted_chunks)
saveRDS(sce_combined, combined_fit_rds)

n_lineages <- ncol(pt)
test_status <- list()
run_test <- function(name, expr, csv_name, rds_name) {
  log_msg("Running ", name, verbose = args$verbose)
  t0 <- Sys.time()
  result <- tryCatch(expr, error = function(e) e)
  elapsed_mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (inherits(result, "error")) {
    test_status[[length(test_status) + 1L]] <<- data.frame(
      test = name,
      status = "error",
      minutes = elapsed_mins,
      message = conditionMessage(result),
      stringsAsFactors = FALSE
    )
    log_msg(name, " failed: ", conditionMessage(result), verbose = TRUE)
    return(NULL)
  }

  saveRDS(result, file.path(test_dir, rds_name))
  write_csv_gz(result, file.path(test_dir, csv_name), row.names = TRUE)
  test_status[[length(test_status) + 1L]] <<- data.frame(
    test = name,
    status = "ok",
    minutes = elapsed_mins,
    message = "",
    stringsAsFactors = FALSE
  )
  result
}

assoc_res <- run_test(
  "associationTest",
  tradeSeq::associationTest(sce_combined),
  "associationTest_combined.csv.gz",
  "associationTest_combined.rds"
)

startvend_res <- run_test(
  "startVsEndTest",
  tradeSeq::startVsEndTest(sce_combined),
  "startVsEndTest_combined.csv.gz",
  "startVsEndTest_combined.rds"
)

if (n_lineages > 1) {
  diffend_res <- run_test(
    "diffEndTest",
    tradeSeq::diffEndTest(sce_combined),
    "diffEndTest_combined.csv.gz",
    "diffEndTest_combined.rds"
  )

  pattern_res <- run_test(
    "patternTest",
    tradeSeq::patternTest(sce_combined),
    "patternTest_combined.csv.gz",
    "patternTest_combined.rds"
  )

  earlyde_res <- run_test(
    "earlyDETest",
    tradeSeq::earlyDETest(
      sce_combined,
      knots = c(1, 2),
      global = TRUE,
      pairwise = FALSE
    ),
    "earlyDETest_combined.csv.gz",
    "earlyDETest_combined.rds"
  )
}

status_df <- do.call(rbind, test_status)
utils::write.csv(status_df, status_csv, row.names = FALSE)

plot_manifest <- list(
  created_at = format(Sys.time(), "%F %T"),
  input_sce_rds = normalizePath(args$sce_rds, winslash = "/", mustWork = TRUE),
  sce_rds = file.path(final_outdir, basename(combined_fit_rds)),
  local_sce_rds = normalizePath(combined_fit_rds, winslash = "/", mustWork = TRUE),
  gene_list_rds = if (!is.null(args$gene_list_rds) && nzchar(args$gene_list_rds)) {
    vapply(
      trimws(strsplit(args$gene_list_rds, ",", fixed = TRUE)[[1]]),
      normalizePath,
      character(1),
      winslash = "/",
      mustWork = TRUE
    )
  } else {
    NA_character_
  },
  assay_name = assay_name,
  lineage_names = colnames(pt),
  pseudotime = pt,
  cellWeights = cw,
  sample = if ("sample" %in% colnames(colData(sce_combined))) as.character(colData(sce_combined)$sample) else NULL,
  cluster_name = if ("cluster_name" %in% colnames(colData(sce_combined))) as.character(colData(sce_combined)$cluster_name) else NULL,
  broad_clusters = if ("broad_clusters" %in% colnames(colData(sce_combined))) as.character(colData(sce_combined)$broad_clusters) else NULL,
  phase = if ("Phase" %in% colnames(colData(sce_combined))) as.character(colData(sce_combined)$Phase) else NULL,
  chunk_manifest = chunk_manifest_df,
  genes = genes_use,
  chunk_size = args$chunk_size,
  nknots = args$nknots,
  n_lineages = n_lineages,
  combined_fit_rds = file.path(final_outdir, basename(combined_fit_rds)),
  test_status = status_df
)
saveRDS(plot_manifest, file.path(args$outdir, "plotSmoothers_manifest.rds"))
utils::write.csv(
  data.frame(gene = genes_use, chunk_id = rep(seq_along(chunks), lengths(chunks)), stringsAsFactors = FALSE),
  file.path(args$outdir, "plotSmoothers_gene_to_chunk.csv"),
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(args$outdir, "sessionInfo.txt"))
log_msg("Done. Output dir: ", args$outdir, verbose = args$verbose)
