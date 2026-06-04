#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(Seurat)
  library(ggplot2)
  library(ggExtra)
})

# --------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------

#checks if directory ends in a backslash and adds one if not
add_back_slash <- function(directory) {
  return(ifelse(endsWith(directory, "/"), directory, paste0(directory, "/")))
}

#checks if file exists
#option to exit program when the directory is missing
check_file_exists <- function(file, force_stop = FALSE) {
  if (!force_stop) {
    return(file.exists(file))
  } else {
    if (!file.exists(file)) {
      stop(sprintf("File %s is missing! Exiting program!", file))
    }
  }
}

create_bygene_df <- function(seurat_obj,
                             outcol_suffix = "",
                             layer = "counts",
                             assay = "all") {
  # Get the assays to process
  if (assay == "all") {
    count_assays <- names(seurat_obj@assays)
  } else {
    if (!assay %in% names(seurat_obj@assays)) {
      stop(sprintf("Assay %s is not in the seurat object provided!", assay))
    }
    count_assays <- assay
  }

  # Extract expression matrices and get the union of all gene names
  exp_mats <- lapply(count_assays, function(assay_name) {
    Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = layer)
  })

  # Get the union of all genes
  genes <- unique(unlist(lapply(exp_mats, rownames)))

  # Initialize bygene_df with all genes
  bygene_df <- data.frame(matrix(ncol = 0, nrow = length(genes)))
  rownames(bygene_df) <- genes

  # Sum counts and cells per gene across assays
  counts_list <- lapply(exp_mats, function(x) {
    x_genes <- rownames(x)
    common_genes <- intersect(genes, x_genes)
    # Subset x to common genes
    x_common <- x[common_genes, , drop = FALSE]
    # Create zero counts for missing genes
    missing_genes <- setdiff(genes, x_genes)
    if (length(missing_genes) > 0) {
      zero_mat <- Matrix::Matrix(
        0,
        nrow = length(missing_genes),
        ncol = ncol(x),
        dimnames = list(missing_genes, colnames(x))
      )
      x_common <- rbind(x_common, zero_mat)
    }
    # Ensure the genes are in the same order
    x_common <- x_common[genes, , drop = FALSE]
    rowSums(x_common)
  })
  counts_mat <- do.call(cbind, counts_list)
  total_counts <- rowSums(counts_mat)

  cells_list <- lapply(exp_mats, function(x) {
    x_genes <- rownames(x)
    common_genes <- intersect(genes, x_genes)
    x_common <- x[common_genes, , drop = FALSE]
    missing_genes <- setdiff(genes, x_genes)
    if (length(missing_genes) > 0) {
      zero_mat <- Matrix::Matrix(0,
                                 nrow = length(missing_genes), ncol = ncol(x),
                                 dimnames = list(missing_genes, colnames(x)))
      x_common <- rbind(x_common, zero_mat)
    }
    x_common <- x_common[genes, , drop = FALSE]
    rowSums(x_common != 0)
  })
  cells_mat <- do.call(cbind, cells_list)
  total_cells <- rowSums(cells_mat)

  # Assign counts and cells to bygene_df
  col_count <- paste0("nCount", outcol_suffix)
  col_ncell <- paste0("nCell", outcol_suffix)
  bygene_df[[col_count]] <- total_counts[rownames(bygene_df)]
  bygene_df[[paste0(col_count, "_log10")]] <- log10(bygene_df[[col_count]] + 1)
  bygene_df[[col_ncell]] <- total_cells[rownames(bygene_df)]
  bygene_df[[paste0(col_ncell, "_log10")]] <- log10(bygene_df[[col_ncell]] + 1)
  return(bygene_df)
}

# The transfer_labels function from your code:
transfer_labels <- function(seurat_obj,
                            reference_rds_path,
                            bygene_df = NULL,
                            label_col,
                            sample_name,
                            outdir,
                            norm_meth = NULL,  # Auto-detect from reference
                            pca_dims = 1:30,
                            future_mem_limit,
                            reference_assay = NULL,  # Auto-detect from reference
                            reference_layer = "counts",
                            min_probability = .50,
                            seed = 42) {

  outdir <- add_back_slash(outdir)
  options(future.globals.maxSize = future_mem_limit)

  reference_seurat <- readRDS(reference_rds_path)
  
  # Auto-detect reference assay if not provided
  if (is.null(reference_assay)) {
    reference_assay <- DefaultAssay(reference_seurat)
    message(sprintf("Using default assay from reference: %s", reference_assay))
  }
  
  # Auto-detect normalization method based on assay type
  if (is.null(norm_meth)) {
    # Check if the assay is an SCT assay
    assay_class <- class(reference_seurat[[reference_assay]])[1]
    if (assay_class == "SCTAssay" || grepl("SCT", assay_class, ignore.case = TRUE)) {
      norm_meth <- "SCT"
      message("Auto-detected SCTransform normalization from reference")
    } else {
      norm_meth <- "LogNormalize"
      message("Using LogNormalize for standard RNA assay")
    }
  }
  
  # Auto-detect reference layer if the default doesn't exist
  available_layers <- Layers(reference_seurat[[reference_assay]])
  if (!reference_layer %in% available_layers) {
    # Try common alternatives: data, counts, scale.data
    for (alt_layer in c("data", "counts", "scale.data")) {
      if (alt_layer %in% available_layers) {
        message(sprintf("Layer '%s' not found, using '%s' instead", reference_layer, alt_layer))
        reference_layer <- alt_layer
        break
      }
    }
  }
  
  if (!label_col %in% colnames(reference_seurat[[]])) {
    stop(
      sprintf(
        "Column name %s does not exist in the reference RDS %s!",
        label_col,
        reference_rds_path
      )
    )
  } else if (
    !reference_layer %in% Layers(reference_seurat[[reference_assay]]) # nolint
  ) {
    stop(
      sprintf(
        "Layer %s does not exist in the rds %s!",
        reference_layer,
        reference_rds_path
      )
    )
  }

  # Print available assays
  print("Available assays in seurat_obj:")
  print(Seurat::Assays(seurat_obj))
  print("Available assays in reference_seurat:")
  print(Seurat::Assays(reference_seurat))

  # Reduce reference to features in the query
  common_features <- intersect(
    Features(seurat_obj, assay = "Xenium"), # nolint
    rownames(reference_seurat[[reference_assay]])
  )
  if (length(common_features) == 0) {
    stop("No common features found between the reference and query objects.")
  }
  reference_seurat <- subset(reference_seurat, features = common_features)

  # If bygene_df not provided, create a minimal one using create_bygene_df
  if (is.null(bygene_df)) {
    bygene_df <- create_bygene_df(seurat_obj)
  }

  print("bygene_df after creation:")
  print(head(bygene_df))

  bygene_ref <- create_bygene_df(
    reference_seurat,
    outcol_suffix = "_scRef",
    layer = reference_layer,
    assay = reference_assay
  )

  print("bygene_ref after creation:")
  print(head(bygene_ref))

  bygene_df <- merge(bygene_df, bygene_ref, by = "row.names", all.x = TRUE)
  rownames(bygene_df) <- bygene_df$Row.names
  bygene_df$Row.names <- NULL

  print("bygene_df after merging with bygene_ref:")
  print(head(bygene_df))

  # Filter bygene_df to remove rows with non-finite values
  print("bygene_df before filtering for lm_data:")
  print(head(bygene_df))

  lm_data <- bygene_df[
    is.finite(bygene_df$nCount_log10) & is.finite(bygene_df$nCount_scRef_log10),
  ]

  print("lm_data after filtering:")
  print(head(lm_data))

  # Check if there is enough data to fit the model
  if (nrow(lm_data) == 0) {
    stop("No finite data available for linear modeling.")
  }

  # Fit the linear model using the filtered data
  mod <- lm(nCount_log10 ~ nCount_scRef_log10, data = lm_data)

  resids <- data.frame(
    scRef_residual = abs(residuals(mod)),
    row.names = rownames(lm_data)
  )

  # Merge residuals back into bygene_df
  bygene_df <- merge(bygene_df, resids, by = "row.names", all.x = TRUE)
  rownames(bygene_df) <- bygene_df$Row.names
  bygene_df$Row.names <- NULL

  # Handle any NA values in scRef_residual
  bygene_df$scRef_residual[is.na(bygene_df$scRef_residual)] <- 0

  # Detect if reference has a pre-computed PCA reduction
  reference_reduction <- NULL
  available_reductions <- names(reference_seurat@reductions)
  
  # Check for common PCA reduction names: pca, refDR (Azimuth), spca, integrated.dr
  for (pca_name in c("refDR", "pca", "spca", "integrated.dr")) {
    if (pca_name %in% available_reductions) {
      reference_reduction <- pca_name
      message(sprintf("Using pre-computed reference reduction: %s", pca_name))
      break
    }
  }
  
  anchors <- Seurat::FindTransferAnchors(
    reference = reference_seurat,
    query = seurat_obj,
    normalization.method = norm_meth,
    reference.reduction = reference_reduction,  # Use pre-computed PCA if available
    features = Features(seurat_obj, assay = "Xenium"), # nolint
    seed.use = seed
  )

  pred <- Seurat::TransferData(
    anchorset = anchors,
    refdata = reference_seurat[[]][[label_col]],
    prediction.assay = TRUE,
    weight.reduction = seurat_obj[["pca"]],
    dims = pca_dims
  )

  pred_df <- data.frame(t(LayerData(pred))) # nolint

  celltype_columns <- names(pred_df)[!grepl("max", names(pred_df))]
  pred_df$predicted_celltype <- apply(
    pred_df[celltype_columns],
    1,
    function(x) celltype_columns[which.max(x)]
  )
  pred_df$probability_celltype <- pred_df$max
  pred_df$max <- NULL

  # Ensure no column name conflicts during merge
  colnames(pred_df) <- make.names(colnames(pred_df), unique = TRUE)

  seurat_obj@meta.data <- merge(
    seurat_obj@meta.data,
    pred_df[c("predicted_celltype", "probability_celltype")],
    by = "row.names",
    all.x = TRUE
  )
  rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$Row.names
  seurat_obj@meta.data$Row.names <- NULL

  # Ensure the column names are correctly aligned
  if ("predicted_celltype.x" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data$predicted_celltype <- seurat_obj@meta.data$predicted_celltype.x # nolint
    seurat_obj@meta.data$probability_celltype <- seurat_obj@meta.data$probability_celltype.x # nolint
    seurat_obj@meta.data$predicted_celltype.x <- NULL
    seurat_obj@meta.data$probability_celltype.x <- NULL
  }

  seurat_obj$filtered_celltype <- seurat_obj$predicted_celltype
  seurat_obj$filtered_celltype[seurat_obj$probability_celltype < min_probability] <- "Ambiguous" # nolint

  plot_bar(
    seurat_obj@meta.data,
    "predicted_celltype",
    figsize = c(14, 7),
    savefig = paste0(outdir, sample_name, "_predicted_celltype_frequencies.png")
  )

  scatter_with_hist(
    seurat_obj@meta.data,
    "nCount_Xenium",
    "probability_celltype",
    title = sample_name,
    savefig = paste0(
      outdir,
      sample_name,
      "_probability_celltype_by_nCount_Xenium.png"
    )
  )

  scatter_with_hist(
    seurat_obj@meta.data,
    "nFeature_Xenium",
    "probability_celltype",
    title = sample_name,
    savefig = paste0(
      outdir,
      sample_name,
      "_probability_celltype_by_nFeature_Xenium.png"
    )
  )
  #plot the predicted cell clusters on the umap
  predcell_umap <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "predicted_celltype",
    cols = "polychrome"
  ) + ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    predcell_umap,
    width = 14,
    filename = paste0(outdir, sample_name, "_umap_predicted_cellclusters.png")
  )


  #plot data on the tissue itself
  spat <- Seurat::ImageDimPlot(
    seurat_obj,
    fov = gsub("-", "_", paste0("sample_", sample_name, "_fov")),
    group.by = "predicted_celltype",
    cols = "polychrome"
  ) + ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    spat,
    filename = paste0(
      outdir,
      sample_name,
      "_insitu_predicted_cellclusters.png"
    ),
    width = 21,
    bg = "white"
  )

  #plot the filtered celltype frequencies
  #plot the predicted cell clusters on the umap
  predcell_umap <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "filtered_celltype",
    cols = "polychrome"
  ) + ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    predcell_umap,
    width = 21,
    filename = paste0(outdir, sample_name, "_umap_filtered_celltype.png")
  )


  #plot data on the tissue itself
  spat <- Seurat::ImageDimPlot(
    seurat_obj,
    fov = gsub("-", "_", paste0("sample_", sample_name, "_fov")),
    group.by = "filtered_celltype",
    cols = "polychrome"
  ) + ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    spat,
    filename = paste0(outdir, sample_name, "_insitu_filtered_celltype.png"),
    width = 21,
    bg = "white"
  )

  plot_bar(
    seurat_obj@meta.data,
    "filtered_celltype",
    figsize = c(14, 7),
    savefig = paste0(outdir, sample_name, "_filtered_celltype_frequencies.png")
  )


  # Save the results
  return(list(
    "seurat_obj" = seurat_obj,
    "bygene_df" = bygene_df,
    "pred_df" = pred_df
  ))
}

#########Plotting functions##########
plot_bar <- function(data_tbl,
                     categorical_column,
                     figsize = c(7, 7),
                     savefig = NULL) {

  #get counts of each category for plotting
  plot_df <- data.frame(table(data_tbl[categorical_column]))

  bar <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[categorical_column]], y = Freq) # nolint
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = categorical_column, y = "Count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  if (!is.null(savefig)) {
    ggplot2::ggsave(
      filename = savefig,
      plot = bar,
      bg = "white",
      width = figsize[1],
      height = figsize[2]
    )
  }
}


#plots a scatterplot with marginal histograms
scatter_with_hist <- function(plot_df, # nolint
                              x_col,
                              y_col,
                              figsize = c(6, 6),
                              scatter_marker_size = .25,
                              hist_bins = 50,
                              color_column = NULL,
                              marker_size_column = NULL,
                              add_fit_line = FALSE,
                              hline = NULL,
                              vline = NULL,
                              xlabel = NULL,
                              ylabel = NULL,
                              ylim = NULL,
                              xlim = NULL,
                              xscale = "linear",
                              yscale = "linear",
                              title = NULL,
                              legend_title = NULL,
                              savefig = NULL) {

  # Set default labels if not provided
  if (is.null(xlabel)) xlabel <- x_col
  if (is.null(ylabel)) ylabel <- y_col

  # Create scatter plot
  #the .data unquotes the strings so ggplot will be happy
  scatter_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[x_col]], # nolint
      y = .data[[y_col]], # nolint
      color = if (!is.null(color_column)) .data[[color_column]] else NULL,
      size = if (!is.null(marker_size_column)) {
        .data[[marker_size_column]] # nolint
      } else {
        NULL
      }
    )
  ) +
    ggplot2::labs(x = xlabel, y = ylabel) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "left")

  #adjust scatter marker size unless its set to vary
  if (!is.null(marker_size_column)) {
    scatter_plot <- scatter_plot + ggplot2::geom_point(show.legend = TRUE)
  } else {
    scatter_plot <- scatter_plot +
      ggplot2::geom_point(
        size = scatter_marker_size,
        show.legend = !is.null(color_column)
      )
  }

  #add linear fit line if requested
  if (add_fit_line) {
    scatter_plot <- scatter_plot +
      ggplot2::stat_smooth(method = "lm", col = "red")
  }

  # Add vline/hline if requested
  if (!is.null(vline)) {
    scatter_plot <- scatter_plot +
      ggplot2::geom_vline(xintercept = vline, color = "red")
  }
  if (!is.null(hline)) {
    scatter_plot <- scatter_plot +
      ggplot2::geom_hline(yintercept = hline, color = "red")
  }

  # Adjust scales
  if (xscale == "log") scatter_plot <- scatter_plot + ggplot2::scale_x_log10()
  if (yscale == "log") scatter_plot <- scatter_plot + ggplot2::scale_y_log10()
  if (!is.null(xlim)) scatter_plot <- scatter_plot + xlim(xlim)
  if (!is.null(ylim)) scatter_plot <- scatter_plot + ylim(ylim)

  #adjust legend title if asked
  #I think it will complain if we don't need a legend (ie no color column)
  if (!is.null(color_column) && !is.null(legend_title)) {
    scatter_plot <- scatter_plot + ggplot2::labs(color = legend_title)
  }
  if (!is.null(marker_size_column) && !is.null(legend_title)) {
    scatter_plot <- scatter_plot + ggplot2::labs(size = legend_title)
  }

  # Add title if provided
  if (!is.null(title)) scatter_plot <- scatter_plot + ggplot2::ggtitle(title)

  #add marginal histograms
  #color by group if color_column provided
  scatter_plot <- ggExtra::ggMarginal(
    scatter_plot,
    type = "histogram",
    bins = hist_bins,
    groupColour = !is.null(color_column),
    groupFill = !is.null(color_column)
  )

  # Save the figure if a path is provided
  if (!is.null(savefig)) {
    ggplot2::ggsave(
      filename = savefig,
      plot = scatter_plot,
      width = figsize[1],
      height = figsize[2]
    )
  }
}


# --------------------------------------------------------------------------
# Main Script
# --------------------------------------------------------------------------

p <- arg_parser("Run label transfer and generate label-based QC plots")
p <- add_argument(p, "--inrds", help = "Path to processed Seurat RDS", type = "character") # nolint
p <- add_argument(p, "--reference_rds", help = "Path to reference single cell RDS for label transfer", type = "character") # nolint
p <- add_argument(p, "--label_col", help = "Column in reference RDS to use for labeling", type = "character") # nolint
p <- add_argument(p, "--outdir", help = "Output directory for results", type = "character") # nolint
p <- add_argument(p, "--sample_name", help = "Sample name", type = "character")
p <- add_argument(p, "--min_celltype_probability", help = "Minimum probability for cell type assignment", type = "double", default = 0.5) # nolint
p <- add_argument(p, "--future_mem_limit", help = "Max memory for future package", type = "double", default = 8e9) # nolint
p <- add_argument(p, "--seed", help = "Random seed for reproducibility (anchor finding, PCA)", type = "integer", default = 42)

args <- parse_args(p)

# Set RNG seed for reproducibility
set.seed(args$seed)

inrds <- args$inrds
reference_rds <- args$reference_rds
label_col <- args$label_col
outdir <- args$outdir
sample_name <- args$sample_name
min_prob <- args$min_celltype_probability
future_mem_limit <- args$future_mem_limit

check_file_exists(inrds, force_stop = TRUE)
check_file_exists(reference_rds, force_stop = TRUE)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

seurat_obj <- readRDS(inrds)
# If bygene_df is needed and not available, create a minimal one:
bygene_df <- create_bygene_df(seurat_obj)

res <- transfer_labels(
  seurat_obj = seurat_obj,
  reference_rds_path = reference_rds,
  bygene_df = bygene_df,
  label_col = label_col,
  sample_name = sample_name,
  outdir = outdir,
  min_probability = min_prob,
  future_mem_limit = future_mem_limit,
  seeding = args$seed
)



# Save updated Seurat object and predictions
saveRDS(
  res$seurat_obj,
  file = sprintf("%s/%s_labeled.rds", outdir, sample_name)
)
write.csv(
  res$pred_df,
  file = sprintf("%s/%s_celltype_predictions.csv", outdir, sample_name),
  row.names = TRUE
)

message("Label transfer completed successfully.")
