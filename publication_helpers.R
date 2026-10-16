process_seurat_object <- function(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    selection.method = "vst",
    nfeatures = 2000,
    npcs = 40,
    resolution.range = seq(0.1, 2, by = 0.1),
    clustering.algorithm = 2,
    plot_umap = TRUE,
    label_clusters = TRUE,
    reduction.name = "umap",
    use_harmony = FALSE,
    harmony_dims = 1:30,
    random_seed = 42
) {
  if (!use_harmony) {
    message("Running full normalization, scaling, PCA pipeline")
    
    # Normalization
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = normalization.method, scale.factor = scale.factor)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = selection.method, nfeatures = nfeatures)
    variable.features <- VariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = variable.features)
    
    # PCA
    seurat_obj <- RunPCA(seurat_obj, features = variable.features, npcs = npcs)
    
    # Algorithmically determine number of principal components to use
    pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
    pcs <- min(co1, co2)
    if (is.na(pcs) || pcs < 2) pcs <- min(10, npcs)
    
    # Neighbors, UMAP, Clustering on PCA
    set.seed(random_seed)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs)
    set.seed(random_seed)
    seurat_obj <- RunUMAP(
      seurat_obj,
      dims = 1:pcs,
      reduction = "pca",
      reduction.name = reduction.name,
      umap.method = "uwot",
      metric = "cosine",
      seed.use = random_seed
    )
    set.seed(random_seed)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, resolution = resolution.range, algorithm = clustering.algorithm)
    
  } else {
    message("Skipping normalization and PCA — using Harmony embeddings")
    
    # Neighbors, UMAP, Clustering on Harmony
    set.seed(random_seed)
    seurat_obj <- FindNeighbors(seurat_obj, dims = harmony_dims, reduction = "harmony")
    set.seed(random_seed)
    seurat_obj <- RunUMAP(
      seurat_obj,
      dims = harmony_dims,
      reduction = "harmony",
      reduction.name = paste0(reduction.name, ".harmony"),
      umap.method = "uwot",
      metric = "cosine",
      seed.use = random_seed
    )
    set.seed(random_seed)
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, resolution = resolution.range, algorithm = clustering.algorithm)
  }
  
  # Find cluster columns in metadata
  cluster_cols <- grep("_snn_res", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(cluster_cols) == 0) {
    stop("No cluster columns found in the metadata.")
  }
  
  # Plot UMAP for each resolution
  if (plot_umap) {
    for (res in cluster_cols) {
      p <- DimPlot(seurat_obj, reduction = ifelse(use_harmony, paste0(reduction.name, ".harmony"), reduction.name),
                   group.by = res, label = label_clusters) +
        ggtitle(paste(ifelse(use_harmony, paste0(reduction.name, ".harmony"), reduction.name), "-", res))
      print(p)
    }
  }
  
  return(seurat_obj)
}


######################################################################

add_module_scores_and_plot <- function(seurat_obj,
                                       marker_list,
                                       clustering_colors = NULL,
                                       cluster_resolution = "RNA_snn_res.0.3",
                                       min_genes_required = 2,
                                       ctrl = 5,
                                       feature_plot_color_low = "white",
                                       feature_plot_color_high = "blue",
                                       reduction.name = "umap") {
  meta <- seurat_obj@meta.data
  feature_cols <- grep("Features", colnames(meta), value = TRUE)
  if (length(feature_cols) > 0) {
    seurat_obj@meta.data <- meta[, !colnames(meta) %in% feature_cols, drop = FALSE]
  }
  
  available_genes <- rownames(seurat_obj)
  
  valid <- Filter(function(g) length(g) >= min_genes_required,
                  lapply(marker_list, intersect, available_genes))
  for (nm in setdiff(names(marker_list), names(valid)))
    message(sprintf("Skipping '%s': too few genes present", nm))
  
  computed_features <- names(valid)
  
  # Single AddModuleScore call
  seurat_obj <- AddModuleScore(seurat_obj, features = valid, ctrl = ctrl, name = "MODSCORE_")
  idx <- match(paste0("MODSCORE_", seq_along(valid)), colnames(seurat_obj@meta.data))
  colnames(seurat_obj@meta.data)[idx] <- paste0(computed_features, "_Features1")
  
  cluster_col_name <- cluster_resolution
  if (!cluster_col_name %in% colnames(seurat_obj@meta.data))
    stop("The specified cluster resolution does not exist in the metadata.")
  num_clusters <- length(unique(seurat_obj@meta.data[[cluster_col_name]]))
  if (is.null(clustering_colors) || length(clustering_colors) < num_clusters)
    clustering_colors <- Seurat:::DiscretePalette(num_clusters, "polychrome")
  
  if (is.null(feature_plot_color_low))  feature_plot_color_low  <- "white"
  if (is.null(feature_plot_color_high)) feature_plot_color_high <- "blue"

  for (name in computed_features) {
    feature_col <- paste0(name, "_Features1")
    vln_plot <- VlnPlot(seurat_obj, group.by = cluster_col_name, features = feature_col, cols = clustering_colors) +
      theme(legend.position = "none")
    feature_plot <- FeaturePlot(seurat_obj, reduction = reduction.name, features = feature_col,
                                cols = c(feature_plot_color_low, feature_plot_color_high), order = TRUE)
    print(vln_plot / (feature_plot))
    rm(vln_plot, feature_plot)
  }
  
  invisible(seurat_obj)
}

find_markers <- function(seurat_obj, resolution) {
  # Set the identity class to the specified resolution
  Idents(object = seurat_obj) <- resolution
  cat("Identities set to resolution:", resolution, "\n")
  
  # Find all markers
  cat("Starting to find all markers...\n")
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  cat("Markers found. Number of markers:", nrow(markers), "\n")
  
  # Filter markers based on adjusted p-value
  cat("Filtering markers by adjusted p-value...\n")
  top_markers <- markers %>% filter(p_val_adj <= 0.01)
  cat("Filtered markers by adjusted p-value. Number of markers after filter:", nrow(top_markers), "\n")
  
  # Calculate pct.diff and filter based on it
  cat("Calculating pct.diff...\n")
  top_markers$pct.diff <- top_markers$pct.1 - top_markers$pct.2
  top_markers <- top_markers %>% filter(pct.diff > 0.1)
  cat("Filtered markers by pct.diff. Number of markers after filter:", nrow(top_markers), "\n")
  
  # Arrange the top markers
  cat("Arranging top markers...\n")
  top_markers <- top_markers %>%
    dplyr::arrange(cluster, desc(pct.diff))
  
  cat("Finished arranging top markers.\n")
  
  # Return results
  return(list(markers = markers, top_markers = top_markers))
}

################################################
score_and_plot_cell_cycle <- function(seurat_obj, cc_genes = cc.genes, reduction = "umap") {
  # Helper function: convert to sentence case (Mouse gene format)
  to_sentence_case <- function(text) {
    text_lower <- tolower(text)
    text_sentence_case <- paste0(
      toupper(substring(text_lower, 1, 1)), 
      substring(text_lower, 2)
    )
    return(text_sentence_case)
  }
  
  # Transform gene lists
  s.genes <- sapply(cc_genes$s.genes, to_sentence_case)
  g2m.genes <- sapply(cc_genes$g2m.genes, to_sentence_case)
  
  # Cell cycle scoring
  seurat_obj <- CellCycleScoring(
    seurat_obj,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )
  
  # Plot S phase score
  s_plot <- FeaturePlot(seurat_obj, features = "S.Score", reduction = reduction, cols = c("lightgrey", "blue"))
  # Plot G2M phase score
  g2m_plot <- FeaturePlot(seurat_obj, features = "G2M.Score", reduction = reduction, cols = c("lightgrey", "red"))
  # Plot Phase
  phase_plot <- DimPlot(seurat_obj, reduction = reduction, group.by = "Phase")
  
  # Print plots
  print(s_plot)
  print(g2m_plot)
  print(phase_plot)
  
  # Return all plots and scored object
  return(list(
    scored_seurat = seurat_obj,
    s_plot = s_plot,
    g2m_plot = g2m_plot,
    phase_plot = phase_plot
  ))
}


qc_plots <- function(seurat_object){
  
  metadata <- seurat_object@meta.data
  
  # Number of cells per sample
  p1 <- ggplot(metadata, aes(x=seq_folder, fill=sample)) +
    geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, size=3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) + 
    ggtitle("NCells")
  
  # UMI counts per cell
  p2 <- ggplot(metadata, aes(x=nUMI, fill= seq_folder)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() + 
    ylab("Cell density") +
    geom_vline(xintercept = 1000)+
    facet_wrap(~seq_folder)

  # Genes per cell
  p3 <- ggplot(metadata, aes(x=nGene, fill= seq_folder)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = 1000)+
    facet_wrap(~seq_folder)
  
  # Genes per UMI
  p4 <- ggplot(metadata, aes(x=log10GenesPerUMI, fill= seq_folder)) + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)+
    facet_wrap(~seq_folder)

  # Mito ratio
  p5 <- ggplot(metadata, aes(x=mitoRatio, fill= seq_folder)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)+
    facet_wrap(~seq_folder)
  
  # Genes vs UMIs
  p6 <- suppressWarnings(
    ggplot(metadata, aes(x = nUMI, y = nGene, color = mitoRatio)) +
      geom_point() +
      stat_smooth(method = lm) +
      scale_x_log10() +
      scale_y_log10() +
      theme_classic() +
      geom_vline(xintercept = 500) + 
      geom_hline(yintercept = 250) +
      facet_wrap(~ seq_folder)
  )
  
  # Arrange plots
  plots <- list(p1, p2, p3, p4, p5, p6)
  
  return(plots)
  
}

lighten <- function(color, factor = 0.2) {
  rgb_col <- col2rgb(color) / 255
  rgb_light <- rgb_col + (1 - rgb_col) * factor
  rgb(rgb_light[1], rgb_light[2], rgb_light[3])
}

generate_shades <- function(base_color, n) {
  colorRampPalette(c(lighten(base_color, 0.4), base_color))(n)
}

create_gene_pair_plot <- function(seurat_obj, genes) {
    
    plots <- lapply(genes, function(gene) {
        
        FeaturePlot(
            seurat_obj,
            reduction = "umap",
            features = gene,
            order = TRUE,
            min.cutoff = 0.5
        ) +
            theme_umap_arrows(
                x_label = "UMAP_1",
                y_label = "UMAP_2"
            ) +
            ggtitle(gene) +
            coord_cartesian(
                xlim = c(umap_limits$xmin, umap_limits$xmax),
                ylim = c(umap_limits$ymin, umap_limits$ymax)
            ) +
            theme(
                plot.title = element_text(
                    size = 12,
                    face = "bold",
                    hjust = 0.5
                ),
                legend.position = "right",
                legend.key.size = unit(0.3, "cm"),
                legend.text = element_text(size = 8),
                legend.title = element_text(size = 9)
            )
    })
    
    cowplot::plot_grid(
        plotlist = plots,
        ncol = 2,
        align = "hv",
        axis = "tblr"
    )
}

process_lineage_object <- function(seurat_obj, resolutions) {
  DefaultAssay(seurat_obj) <- "RNA"

  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000
  )

  variable_features <- VariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(
    seurat_obj,
    features = variable_features
  )
  seurat_obj <- RunPCA(
    seurat_obj,
    features = variable_features,
    npcs = 40,
    seed.use = pca_seed
  )

  # PC-selection 
  pct <- Stdev(seurat_obj, reduction = "pca") /
    sum(Stdev(seurat_obj, reduction = "pca")) * 100
  cumulative_pct <- cumsum(pct)
  cutoff_1 <- which(cumulative_pct > 90 & pct < 5)[1]
  cutoff_2 <- sort(
    which((pct[-length(pct)] - pct[-1]) > 0.1),
    decreasing = TRUE
  )[1] + 1
  pcs <- min(cutoff_1, cutoff_2)
  if (is.na(pcs) || pcs < 2) pcs <- min(10, length(pct))

  message("Using PCs 1-", pcs, " for neighbors and the 2D UMAP.")

  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = "pca",
    dims = seq_len(pcs),
    k.param = 20,
    compute.SNN = TRUE,
    prune.SNN = 1 / 15,
    nn.method = "annoy",
    n.trees = 50,
    annoy.metric = "euclidean",
    nn.eps = 0,
    l2.norm = FALSE
  )
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca",
    dims = seq_len(pcs),
    umap.method = "uwot",
    n.neighbors = 30,
    n.components = 2,
    metric = "cosine",
    learning.rate = 1,
    min.dist = 0.3,
    spread = 1,
    set.op.mix.ratio = 1,
    local.connectivity = 1,
    repulsion.strength = 1,
    negative.sample.rate = 5,
    uwot.sgd = FALSE,
    seed.use = umap_2d_seed,
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )
  seurat_obj <- FindClusters(
    seurat_obj,
    graph.name = "RNA_snn",
    resolution = resolutions,
    algorithm = 2,
    modularity.fxn = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = cluster_seed,
    group.singletons = TRUE,
    verbose = FALSE
  )

  seurat_obj
}

add_umap3d <- function(seurat_obj) {
  RunUMAP(
    seurat_obj,
    reduction = "pca",
    dims = 1:20,
    umap.method = "uwot",
    return.model = TRUE,
    n.neighbors = 30,
    n.components = 3,
    metric = "cosine",
    learning.rate = 1,
    min.dist = 0.3,
    spread = 1,
    set.op.mix.ratio = 1,
    local.connectivity = 1,
    repulsion.strength = 1,
    negative.sample.rate = 5,
    uwot.sgd = FALSE,
    seed.use = umap_3d_seed,
    reduction.name = "umap3d",
    reduction.key = "UMAP3D_"
  )
}

compute_mean_marker_scores <- function(seurat_obj, marker_list, prefix) {
  if (inherits(seurat_obj[["RNA"]], "Assay5") &&
      length(Layers(seurat_obj[["RNA"]])) > 1L) {
    seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
  }

  expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  for (marker_name in names(marker_list)) {
    genes_present <- intersect(marker_list[[marker_name]], rownames(expr))
    score_name <- paste0(prefix, marker_name)
    seurat_obj[[score_name]] <- if (length(genes_present) == 0L) {
      rep(0, ncol(seurat_obj))
    } else {
      Matrix::colMeans(expr[genes_present, , drop = FALSE])
    }
  }
  seurat_obj
}

assign_top_scoring_label <- function(
    seurat_obj,
    score_prefix = "score_",
    out_col = "refined_celltype"
) {
  score_cols <- paste0(score_prefix, names(lineage_marker_list))
  missing_score_cols <- setdiff(score_cols, colnames(seurat_obj@meta.data))
  if (length(missing_score_cols) > 0L) {
    stop("Missing score columns: ", paste(missing_score_cols, collapse = ", "))
  }

  score_matrix <- as.matrix(
    seurat_obj@meta.data[, score_cols, drop = FALSE]
  )
  top_index <- max.col(score_matrix, ties.method = "first")
  sorted_scores <- t(apply(score_matrix, 1, sort, decreasing = TRUE))

  seurat_obj[[out_col]] <- sub(
    paste0("^", score_prefix),
    "",
    score_cols[top_index]
  )
  seurat_obj[[paste0(out_col, "_score")]] <- sorted_scores[, 1]
  seurat_obj[[paste0(out_col, "_gap")]] <-
    sorted_scores[, 1] - sorted_scores[, 2]
  seurat_obj
}

score_and_label_lineage <- function(seurat_obj) {
  seurat_obj <- compute_mean_marker_scores(
    seurat_obj,
    lineage_marker_list,
    prefix = "score_"
  )
  seurat_obj <- compute_mean_marker_scores(
    seurat_obj,
    off_lineage_marker_list,
    prefix = "offscore_"
  )
  seurat_obj <- assign_top_scoring_label(
    seurat_obj,
    score_prefix = "score_",
    out_col = "refined_celltype"
  )
  seurat_obj$broad_clusters <- as.character(seurat_obj$refined_celltype)
  seurat_obj
}

summarise_sus <- function(seurat_obj, stage) {
  seurat_obj@meta.data %>%
    mutate(sample = tolower(as.character(sample))) %>%
    group_by(sample) %>%
    summarise(
      total_cells = n(),
      sus_cells = sum(cluster_name == "Sus"),
      sus_percent = 100 * sus_cells / total_cells,
      .groups = "drop"
    ) %>%
    mutate(stage = stage)
}

get_assay_matrix <- function(seurat_object, assay_name, layer_name) {
  assay_object <- seurat_object[[assay_name]]
  if (inherits(assay_object, "Assay5")) {
    layer_names <- SeuratObject::Layers(assay_object)
    if (!(layer_name %in% layer_names)) {
      seurat_object <- SeuratObject::JoinLayers(seurat_object, assay = assay_name)
    }
  }
  SeuratObject::GetAssayData(
    seurat_object,
    assay = assay_name,
    layer = layer_name
  )
}

lineage_table <- function(sce_object) {
  paths <- slingshot::slingLineages(sce_object)
  data.frame(
    lineage = paste0("Lineage_", seq_along(paths)),
    endpoint = vapply(paths, function(path) tail(as.character(path), 1), character(1)),
    path = vapply(paths, function(path) paste(as.character(path), collapse = " -> "), character(1)),
    stringsAsFactors = FALSE
  )
}

project_rgl_view <- function(coords_3d) {
  projected <- rgl::rgl.user2window(
    coords_3d[, 1], coords_3d[, 2], coords_3d[, 3]
  )
  projected <- as.matrix(projected)
  if (nrow(projected) != nrow(coords_3d) && ncol(projected) == nrow(coords_3d)) {
    projected <- t(projected)
  }
  if (nrow(projected) != nrow(coords_3d) || ncol(projected) < 2) {
    stop("Unexpected rgl.user2window output shape.")
  }
  data.frame(
    Dim_1 = projected[, 1],
    Dim_2 = projected[, 2],
    row.names = rownames(coords_3d)
  )
}

format_sender_label <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

load_cellchat_list <- function(path) {
  env <- new.env(parent = baseenv())
  loaded_names <- load(path, envir = env)
  preferred_names <- c("object.list", "obj.list")
  object_name <- preferred_names[preferred_names %in% loaded_names][1]

  if (is.na(object_name)) {
    stop(
      "The RData file does not contain `object.list` or `obj.list`. Found: ",
      paste(loaded_names, collapse = ", ")
    )
  }

  object_list <- get(object_name, envir = env, inherits = FALSE)
  if (!is.list(object_list)) {
    stop("Loaded object `", object_name, "` is not a list.")
  }
  object_list
}

get_incoming_lr <- function(cellchat, dataset_name, receiver_name) {
  CellChat::subsetCommunication(cellchat, slot.name = "net") %>%
    filter(target == receiver_name) %>%
    transmute(
      dataset = dataset_name,
      sender = as.character(source),
      receiver = as.character(target),
      pathway_name = as.character(pathway_name),
      ligand = as.character(ligand),
      receptor = as.character(receptor),
      prob = as.numeric(prob)
    )
}

get_incoming_pathways <- function(cellchat, dataset_name, receiver_name) {
  CellChat::subsetCommunication(cellchat, slot.name = "netP") %>%
    filter(target == receiver_name) %>%
    transmute(
      dataset = dataset_name,
      sender = as.character(source),
      receiver = as.character(target),
      pathway_name = as.character(pathway_name),
      prob = as.numeric(prob)
    )
}

make_sender_pathway_panel <- function(sender_name, add_panel_tag = FALSE) {
  panel_data <- focused_pathways %>%
    filter(sender == sender_name) %>%
    arrange(desc(abs_delta), pathway_name) %>%
    slice_head(n = params$top_n_pathways) %>%
    mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))

  sender_color <- sender_palette[[sender_name]]
  if (is.null(sender_color) || is.na(sender_color)) sender_color <- "#555555"

  panel <- ggplot(
    panel_data,
    aes(x = delta_p5_minus_p22, y = pathway_name)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.35
    ) +
    geom_point(
      aes(size = abs_delta),
      color = sender_color,
      alpha = 0.95
    ) +
    scale_size(range = c(1.5, 4.3), guide = "none") +
    scale_x_continuous(
      limits = c(-pathway_x_limit, pathway_x_limit),
      breaks = breaks_pretty(n = 3),
      labels = label_number(accuracy = 0.001)
    ) +
    labs(
      title = format_sender_label(sender_name),
      x = NULL,
      y = NULL,
      tag = if (add_panel_tag) "B" else NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 11, face = "bold"),
      plot.tag.position = c(0.01, 1.10),
      axis.text.x = element_text(size = 6.5),
      axis.text.y = element_text(size = 7.2),
      axis.ticks.y = element_blank(),
      plot.margin = margin(4, 7, 4, 7)
    )

  panel
}

