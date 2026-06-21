set.seed(1234)

# Libraries --------------------------------------------------------------------
library(yaml)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(hdf5r)
library(magrittr)
library(viridis)
library(scDblFinder)
library(harmony)
library(Seurat)
library(SoupX)
library(SeuratObject)
library(patchwork)
library(DropletUtils)
library(glmGamPoi)
library(org.Hs.eg.db)
library(reshape2)

# Read config file and create output folder
yml <- read_yaml("configs/preprocessing_workflow.yml")

# Default-fallback helper for optional YAML keys
`%||%` <- function(a, b) if (is.null(a)) b else a

# Resolve configurable manifest column names (default to current names if absent)
participant_col <- yml$participant_col %||% "participant"
visit_col       <- yml$visit_col       %||% "visit"
arm_col         <- yml$arm_col         %||% "arm"
sex_col         <- yml$sex_col         %||% "sex"
age_col         <- yml$age_col         %||% "age"

# Resolve allowed values for the visit and arm columns
baseline_visit <- yml$baseline_visit %||% "bl"
followup_visit <- yml$followup_visit %||% "eot"
reference_arm  <- yml$reference_arm  %||% "placebo"
treatment_arm  <- yml$treatment_arm  %||% "treatment"

# Read manifest file to go through all samples
sample_table <- read.table(yml$ManifestPath, header = TRUE)

# Validate manifest: required columns present, and visit/arm contain exactly the two configured values
required_cols <- c("key", "filepath", participant_col, visit_col, arm_col, sex_col, age_col)
missing_cols  <- setdiff(required_cols, colnames(sample_table))
if (length(missing_cols) > 0) {
  stop("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
}

visit_vals <- unique(sample_table[[visit_col]])
arm_vals   <- unique(sample_table[[arm_col]])
if (!setequal(visit_vals, c(baseline_visit, followup_visit))) {
  stop(sprintf(
    "Column '%s' must contain exactly the two values '%s' and '%s'; found: %s",
    visit_col, baseline_visit, followup_visit, paste(visit_vals, collapse = ", ")
  ))
}
if (!setequal(arm_vals, c(reference_arm, treatment_arm))) {
  stop(sprintf(
    "Column '%s' must contain exactly the two values '%s' and '%s'; found: %s",
    arm_col, reference_arm, treatment_arm, paste(arm_vals, collapse = ", ")
  ))
}

sample_ids <- sample_table$key
sample_paths <- sample_table$filepath

# 1. Read data --------------------------------------------------------------
read_data <- function(sample_path, sample_id) {
  # Create results path
  resultspath <- paste0(yml$QCPath, "/", sample_id)
  dir.create(resultspath, showWarnings = FALSE, recursive = TRUE)

  # Read filtered_bc_ft_matrix or h5
  if (yml$matrix_barcode_feats) {
    counts <- Seurat::ReadMtx(
      cells = paste0(sample_path, "/sample_filtered_feature_bc_matrix/barcodes.tsv"),
      mtx = paste0(sample_path, "/sample_filtered_feature_bc_matrix/matrix.mtx"),
      features = paste0(sample_path, "/sample_filtered_feature_bc_matrix/features.tsv"), feature.column = 1
    )

    # Extract current ENSG IDs from the matrix and map
    ensembl_ids <- rownames(counts)
    ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)
    gene_symbols <- mapIds(org.Hs.eg.db,
      keys = ensembl_ids_clean,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    # Keeping the Ensembl ID for not found genes
    new_gene_names <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
    # Check for duplicates among the new gene names
    duplicate_indices <- which(duplicated(new_gene_names))
    if (length(duplicate_indices) > 0) {
      message(paste("Found", length(duplicate_indices), "duplicate gene symbols after mapping. Making them unique."))
      # Make them unique by appending suffixes like "_1", "_2", etc.
      new_gene_names <- make.unique(new_gene_names)
    }
    # Assign the new gene names to your matrix rownames
    rownames(counts) <- new_gene_names

    # Check how many IDs mapped
    message(paste("Mapped", sum(!is.na(gene_symbols)), "out of", length(ensembl_ids), "Ensembl IDs."))

    if (yml$DoEmptyDrop) {
      cat("Performing Empty Drop ... \n")
      ret <- DropletUtils::emptyDropsCellRanger(counts)
      counts <- counts[, which(ret$FDR <= yml$FDR_cutoff_EmptyDrop), drop = FALSE]
      df <- as.data.frame(ret) |>
        dplyr::filter(!is.na(FDR)) |>
        dplyr::mutate(Group = factor(ifelse(PValue < yml$FDR_cutoff_EmptyDrop, "Cell", "EmptyDrop")))

      pdf(file = paste0(resultspath, "/EmptyDrop_Diag.pdf"), width = 10, height = 8)
      print(ggplot2::ggplot(df, aes(x = Total, y = -LogProb, color = Group)) +
              geom_point(alpha = 0.5) +
              labs(x = "Total UMI Count", y = "-Log Probability"))
      dev.off()
    }
  }

  # Input filtered h5 matrix
  if (yml$h5_filtered) {
    hdf5 <- Seurat::Read10X_h5(
      filename = paste0(sample_path, "/sample_filtered_feature_bc_matrix.h5"),
      use.names = TRUE,
      unique.features = TRUE
    )

    cat("Merging counts into matrix...\n")
    if (yml$modality) {
      counts <- hdf5[[yml$modality_type]]
    } else {
      counts <- h5df5
    }

    if (yml$DoEmptyDrop) {
      cat("Performing Empty Drop ... \n")
      ret <- DropletUtils::emptyDropsCellRanger(counts)
      counts <- counts[, which(ret$FDR <= yml$FDR_cutoff_EmptyDrop), drop = FALSE]

      df <- as.data.frame(ret) |>
        dplyr::filter(!is.na(FDR)) |>
        dplyr::mutate(Group = factor(ifelse(PValue < yml$FDR_cutoff_EmptyDrop, "Cell", "EmptyDrop")))

      pdf(file = paste0(resultspath, "/EmptyDrop_Diag.pdf"), width = 10, height = 8)
      print(ggplot2::ggplot(df, aes(x = Total, y = -LogProb, color = Group)) +
              geom_point(alpha = 0.5) +
              labs(x = "Total UMI Count", y = "-Log Probability"))
      dev.off()
    }
  }

  # Input h5 filtered and raw matrixes
  if (yml$h5_raw_filtered) {
    # Read matrix filtered and raw matrix
    filt.matrix <- Read10X_h5(paste0(sample_path, "/filtered_feature_bc_matrix.h5"))
    raw.matrix <- Read10X_h5(paste0(sample_path, "/raw_feature_bc_matrix.h5"))

    if (yml$modality) {
      filt.matrix <- filt.matrix[[yml$modality_type]]
      raw.matrix <- raw.matrix[[yml$modality_type]]
    }

    # Calculate and plot Barcode Ranks
    barcode_umis <- Matrix::colSums(filt.matrix)
    barcode_umis_sorted <- sort(barcode_umis, decreasing = TRUE)
    brank_df <- data.frame(
      rank = seq_along(barcode_umis_sorted),
      total = barcode_umis_sorted
    )
    pdf(file = paste0(resultspath, "/barcode_rank_before_ed.pdf"), width = 10, height = 8)
    print(
      ggplot(brank_df, aes(x = rank, y = total)) +
        geom_point(size = 0.5) +
        scale_x_log10() +
        scale_y_log10() +
        labs(title = "Barcode Rank Plot",
             x = "Barcode Rank",
             y = "Total UMI Count") +
        theme_bw()
    )
    dev.off()

    ## Remove ambient RNA with SoupX
    cat("Removing ambient RNA with SoupX...\n")
    soup.channel <- SoupChannel(tod = raw.matrix, toc = filt.matrix)

    srat <- CreateSeuratObject(counts = filt.matrix)
    srat <- SCTransform(srat, verbose = FALSE)
    srat <- RunPCA(srat, verbose = FALSE)
    srat <- RunUMAP(srat, dims = 1:30, verbose = FALSE)
    srat <- FindNeighbors(srat, dims = 1:30, verbose = FALSE)
    srat <- FindClusters(srat, verbose = TRUE)

    # Adding clusters to the soup channel
    meta <- srat@meta.data
    umap <- srat@reductions$umap@cell.embeddings
    soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
    soup.channel <- setDR(soup.channel, umap)

    # Begin calculating ambient RNA profile
    cat("Calculating ambient RNA profile \n")

    # Use tryCatch to handle potential errors
    soup.channel <- tryCatch(
      {
        autoEstCont(soup.channel) # Attempt to calculate ambient RNA profile
      },
      error = function(e) {
        message("An error occurred during autoEstCont: ", e$message) # Print error message
        return(NULL) # Return NULL or any fallback option if an error occurs
      }
    )

    # Proceed with plotting only if soup.channel is not NULL
    if (!is.null(soup.channel)) {
      # Plot ambient RNA
      pdf(file = paste0(resultspath, "/SoupX_contamination_densityPlot.pdf"), width = 10, height = 8)
      print(autoEstCont(soup.channel))
      dev.off()

      # Get counts after SoupX and remove unnecessary columns
      top_genes <- head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = TRUE), ], n = 20)
      print(top_genes)

      counts <- adjustCounts(soup.channel, roundToInt = TRUE) # Adjusted counts
    } else {
      # If an error occurred, use the alternative option
      counts <- filt.matrix # Assign counts to filt.matrix if an error occurred
      message("Using filtered matrix as backup due to an error in autoEstCont.")
    }
  }

  cat("Creating Seurat Object...\n")
  s <- CreateSeuratObject(counts = counts, min.cells = yml$nCells_cutoff, min.features = yml$nFeatures_min)
  DefaultAssay(s) <- "RNA"
  
  s[["percent_mito"]]    <- PercentageFeatureSet(s, pattern = "^MT-")
  s[["percent_ribo"]]  <- PercentageFeatureSet(s, pattern = "^RPS|^RPL")
  s[["percent_plat"]] <- PercentageFeatureSet(s, "PECAM1|PF4")
  s[["complexity"]] <- s$nFeature_RNA / s$nCount_RNA

  # Subset seurat object, according to min and max nfeatures and percent of mt rna
  s <- subset(s, subset = nFeature_RNA < yml$nFeatures_max & percent_mito < yml$percentMT_cutoff)

  # Add metadata
  cat("Adding metadata and QC metrics to the seurat object...\n")
  s$patient_id <- sample_table[[participant_col]][sample_table$key == sample_id[1]]
  s$visit      <- sample_table[[visit_col]][sample_table$key == sample_id[1]]
  s$arm        <- sample_table[[arm_col]][sample_table$key == sample_id[1]]
  s$age        <- sample_table[[age_col]][sample_table$key == sample_id[1]]
  s$sex        <- sample_table[[sex_col]][sample_table$key == sample_id[1]]

  # Return the Seurat object and resultspath as a list
  rm(counts)
  return(list(seurat_object = s, resultspath = resultspath))
}

# 2. Run QC --------------------------------------------------------------------
run_QC <- function(s, resultspath) {
  # Add QC metrics


  # Plotting QC metrics
  pdf(file = paste0(resultspath, "/QC_vln_scCustom.pdf"), width = 12, height = 12)
  p0 <- VlnPlot(s, features = "nFeature_RNA", pt.size = 0)
  p1 <- VlnPlot(s, features = "nCount_RNA", pt.size = 0)
  p2 <- VlnPlot(s, features = "percent_mito", pt.size = 0)
  p3 <- VlnPlot(s, features = "percent_ribo", pt.size = 0)
  p4 <- VlnPlot(s, features = "percent_plat", pt.size = 0)
  p5 <- VlnPlot(s, features = "nCount_RNA", pt.size = 0)
  p6 <- VlnPlot(s, features = "complexity", pt.size = 0)
  print(wrap_plots(p0, p1, p2, p3, p4, p5, p6, nrow = 2))
  dev.off()

  # Plot scatter plots
  pdf(file = paste0(resultspath, "/Scatter_plots.pdf"), width = 12, height = 10)
  p1 <- print(FeatureScatter(s, feature1 = "percent_mito", feature2 = "percent_ribo"))
  p2 <- print(FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  print(p1 + p2)
  dev.off()

  # Plot histogram of QC metrics
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_plat")
  df <- s[[]]
  df_long <- melt(df, variable.name = "Feature", value.name = "Value")
  pdf(file = paste0(resultspath, "/Histogram_Seurat.pdf"), width = 16, height = 12)
  p_hist <- ggplot(df_long, aes(x = Value)) +
    geom_histogram(fill = "blue", color = "black", bins = 40) +
    facet_wrap(~Feature, scales = "free") +
    theme_minimal()
  print(p_hist)
  dev.off()

  # Filter seurat object
  s <- s[!grepl("MALAT1", rownames(s)), ]
  s <- s[!grepl("^MT-", rownames(s)), ]
  s <- s[!grepl("^RP[SL]", rownames(s)), ]

  # Gender plot
  gender_genes <- c("XIST", "DDX3Y", "TTTY14")
  pdf(file = paste0(resultspath, "/GenderCheck_VlnPlot.pdf"), width = 6, height = 10)
  p_vln <- VlnPlot(s, features = gender_genes, pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  print(p_vln)
  dev.off()

  return(s)
}

# 3. Normalize -----------------------------------------------------------------
run_normalization <- function(s, resultspath) {
  cat("Normalizing data...\n")
  s <- NormalizeData(s)

  # Cell cycle scoring, cells in S and G2m phase
  g2m <- cc.genes$g2m.genes
  g2m <- g2m[g2m %in% rownames(s)]

  s_features <- cc.genes$s.genes
  s_features <- s_features[s_features %in% rownames(s)]

  s <- CellCycleScoring(
    object = s, g2m.features = g2m,
    s.features = s_features
  )

  # Cell cycle plot
  gender_feats <- c("S.Score", "G2M.Score")
  pdf(file = paste0(resultspath, "/CellCyle_VlnPlot.pdf"), width = 12, height = 8)
  p_vln <- VlnPlot(s, features = gender_feats, pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  print(p_vln)
  dev.off()

  return(s)
}

# 4. Identify highly variable features------------------------------------------
find_variables <- function(s, resultspath) {
  cat("Identifying highly variable features...\n")
  s <- FindVariableFeatures(s, verbose = FALSE)

  # View and plot top10 if wanted (he does not do it)
  top10_features <- head(VariableFeatures(s), 10)
  pdf(file = paste0(resultspath, "/Top10_highVariablesPlot.pdf"), width = 12, height = 8)
  plot <- VariableFeaturePlot(s)
  print(LabelPoints(plot = plot, points = top10_features, repel = TRUE, xnudge = 0, ynudge = 0))
  dev.off()

  return(s)
}
# 5. Scale and PCA -------------------------------------------------------------
run_scaling_clustering <- function(s, resultspath) {
  cat("Scaling...")
  s <- ScaleData(s,
    vars.to.regress = c("percent_mito", "S.Score", "G2M.Score"),
    verbose = FALSE
  )

  cat("Performing linear dimensionality reduction...\n")
  s <- RunPCA(s, verbose = FALSE, npcs = 30)

  # Visualize PCA results using heatmap and elbowplot
  pdf(file = paste0(resultspath, "/Heatmap_PCAPlot.pdf"))
  print(DimHeatmap(s, dims = 1:30, cells = 500, balanced = TRUE))
  dev.off()

  # Choosing right number of PCAs
  # Determine percent of variation associated with each PC
  pct <- s[["pca"]]@stdev / sum(s[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC, last point where change of % of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  # Minimum of the two calculation
  pcs <- min(co1, co2)

  pdf(file = paste0(resultspath, "/PCA_ElbowPlot.pdf"), width = 12, height = 8)
  p <- ElbowPlot(s, ndims = pcs + 5)
  print(p + ggtitle(paste0("PCAs selected:", pcs)) + theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  s <- FindNeighbors(s, dims = 1:pcs)
  s <- FindClusters(s, resolution = yml$cluster_resolution)
  s <- RunUMAP(s, dims = 1:pcs, verbose = FALSE)
  head(Idents(s), 5)

  pdf(file = paste0(resultspath, "/Dimplot_umap_scCustom.pdf"), width = 12, height = 8)
  print(DimPlot(s, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE, ncol = 2))
  dev.off()

  #Plot genes of interest
  genes_of_interest <- yml$genes_of_int
  pdf(file = paste0(resultspath, "/CheckGenesofInterest.pdf"), width = 16, height = 24)
  print(FeaturePlot(s, features = genes_of_interest, ncol = 2))
  dev.off()

  return(s)
}

# 6. Doublet detection ---------------------------------------------------------
run_doubletdetection <- function(s, resultspath) {
  sce <- scDblFinder(GetAssayData(s, layer = "counts"), clusters = Idents(s))

  # Import the scDblFinder class from sce object to the Seurat Object
  s$scDblFinder.class <- sce$scDblFinder.class

  pdf(file = paste0(resultspath, "/Umap_Doublet_vs_Singlet.pdf"), width = 14, height = 8)
  par(mfrow = c(1, 2))
  p1 <- DimPlot(s, reduction = "umap", group.by = "scDblFinder.class", raster = FALSE)
  p2 <- DimPlot(s, reduction = "umap", split.by = "scDblFinder.class", raster = FALSE)
  print(wrap_plots(p1, p2, ncol = 2))
  dev.off()

  Idents(s) <- s@meta.data$scDblFinder.class
  pdf(file = paste0(resultspath, "/VlnPlot_Doublet_vs_singlet.pdf"), width = 12, height = 8)
  p1 <- VlnPlot(s, features = "nFeature_RNA", split.by = "scDblFinder.class", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p2 <- VlnPlot(s, features = "nCount_RNA", split.by = "scDblFinder.class", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1 / p2)
  dev.off()

  # Keep only singlets (re-think plotting)
  s <- subset(s, idents = "singlet")
  # Save and zip seurat object
  cat("Saving seurat object...")
  saveRDS(s, file = paste0(resultspath, "/seurat_object.rds"))

  return(s)
}

# 7. Run functions -------------------------------------------------------------
process_sample <- function(sample_path, sample_id) {
  results <- read_data(sample_path = sample_path, sample_id = sample_id)
  # Save output
  s <- results$s
  resultspath <- results$resultspath

  s <- run_QC(s, resultspath)
  s <- run_normalization(s, resultspath)
  s <- find_variables(s, resultspath)
  s <- run_scaling_clustering(s, resultspath)
  s <- run_doubletdetection(s, resultspath)
  return(list(s = s, resultspath = resultspath))
}

for (i in seq_along(sample_paths)) {
  process_sample(sample_path = sample_paths[i], sample_id = sample_ids[i])
}

# 8.Data Integration -----------------------------------------------------------
if (yml$DoIntegration) {
  # Initialize an empty list so save the seurat objects
  nn <- vector("list", length(sample_paths))
  nn <- lapply(seq_along(sample_paths), function(i) {
    readRDS(paste0(yml$QCPath, "/", sample_ids[i], "/seurat_object.rds"))
  })

  if (!all(sapply(nn, is, "Seurat"))) {
    stop("All elements in nn must be Seurat objects")
  }
  ## RPCA
  if (yml$integration_method == "rpca") {
    cat("Performing RPCA...")
    # normalize and identify variable features for each dataset independently
    NN.list <- lapply(X = nn, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = NN.list)
    NN.list <- lapply(X = NN.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })

    # RPCA Integration
    # With the dataset probably going to be as large as I anticipate, in here, we want to add reference datasets. The formatting for that will be after the commented line after here...
    NNanchors <- FindIntegrationAnchors(object.list = NN.list, reference = yml$rpca_anchors_features, anchor.features = features, reduction = "rpca")

    NNcombined <- IntegrateData(anchorset = NNanchors, dims = 1:30)
    DefaultAssay(NNcombined) <- "integrated"

    # Run the standard workflow for visualization and clustering
    NNcombined <- ScaleData(NNcombined)
    NNcombined <- RunPCA(NNcombined, features = VariableFeatures(object = NNcombined))

    # Determine percent of variation associated with each PC
    pct <- NNcombined[["pca"]]@stdev / sum(NNcombined[["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
    pcs <- min(co1, co2)

    NNcombined <- RunUMAP(NNcombined, dims = 1:pcs)
    NNcombined <- FindNeighbors(NNcombined, dims = 1:pcs)
    # For running louvain algorithm through clustering...
    NNcombined <- FindClusters(NNcombined, resolution = yml$cluster_resolution)

    # Plot umap
    pdf(file = paste0(yml$QCPath, "/rpca_UMAP.pdf"), width = 12, height = 8)
    p1 <- DimPlot(NNcombined, reduction = "umap", group.by = "patient_id", split.by = "visit", pt.size = .1, raster = FALSE)
    p2 <- DimPlot(NNcombined, reduction = "umap", label = TRUE, pt.size = .1, raster = FALSE)
    print(wrap_plots(p1, p2, nrow = 2))
    dev.off()

    saveRDS(NNcombined, yml$SeuratPath)
  }

  ## Harmony -----------------------------------------------------------------
  if (yml$integration_method == "harmony") {
    cat("Performing Harmony...")
    # Merge datasets
    nn1 <- nn[[1]]
    nn_left <- nn[-1]
    AllNN <- merge(nn1, y = nn_left, add.cell.ids = sample_ids, project = "NNIntegrate")
    rm(nn, nn1, nn_left)

    # Normalize and scale the merged object
    AllNN <- NormalizeData(AllNN)
    AllNN <- FindVariableFeatures(AllNN, selection.method = "vst", nfeatures = 2000)
    AllNN <- ScaleData(AllNN)

    # Run PCA
    AllNN <- RunPCA(AllNN, features = VariableFeatures(object = AllNN))
    AllNN$orig.ident <- paste0(AllNN$patient_id, "_", AllNN$visit)

    # Run Harmony
    AllNN_harmony <- AllNN |>
      RunHarmony(group.by.vars = "orig.ident", plot_convergence = FALSE, nclust = 50, max_iter = 10, early_stop = TRUE)
    rm(AllNN)

    pdf(file = paste0(yml$QCPath, "/Harmonization_Quality.pdf"), width = 16, height = 8)
    p1 <- DimPlot(object = AllNN_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident", raster = FALSE)
    p2 <- VlnPlot(AllNN_harmony, features = "harmony_1", group.by = "orig.ident", pt.size = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p1 / p2)   
    dev.off()

    # Plot genes correlated with the harmonized PCs
    pdf(file = paste0(yml$QCPath, "/Harmony_PCsvsGenes.pdf"), width = 12, height = 8)
    print(DimHeatmap(object = AllNN_harmony, reduction = "harmony", cells = 500, dims = 1:5))
    dev.off()

    # Proceed with downstream analysis
    AllNN_harmony <- FindNeighbors(AllNN_harmony, reduction = "harmony", dims = 1:30)
    AllNN_harmony <- FindClusters(AllNN_harmony, resolution = 0.5)
    AllNN_harmony <- RunUMAP(AllNN_harmony, reduction = "harmony", dims = 1:30)
    cat("Save seurat object")
    saveRDS(AllNN_harmony, yml$SeuratPath)

    # Plot UMAP
    pdf(file = paste0(yml$QCPath, "/Harmony_UMAP.pdf"), width = 12, height = 8)
    p1 <- DimPlot(AllNN_harmony, reduction = "umap", group.by = "patient_id", split.by = "visit", pt.size = .1, raster = FALSE)
    p2 <- DimPlot(AllNN_harmony, reduction = "umap", label = TRUE, pt.size = .1, raster = FALSE)
    print(wrap_plots(p1, p2, nrow = 2))
    dev.off()
  }
}
