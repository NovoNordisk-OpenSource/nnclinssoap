# Load in libraries
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(scales)
library(grid)
library(dplyr)
library(parallel)
library(Seurat)
library(Matrix)
library(data.table)
library(edgeR)
library(tibble)


# Initialize -------------------------------------------------------------------
# Read config, create results path, load data and call source functions
yml <- yaml::read_yaml("configs/preprocessing_workflow.yml")
s <- readRDS(yml$SeuratPath)

# Default-fallback helper for optional YAML keys
`%||%` <- function(a, b) if (is.null(a)) b else a

# Resolve allowed values for the visit and arm columns (defaults match prior hardcoded behaviour)
baseline_visit <- yml$baseline_visit %||% "bl"
followup_visit <- yml$followup_visit %||% "eot"
reference_arm  <- yml$reference_arm  %||% "placebo"
treatment_arm  <- yml$treatment_arm  %||% "treatment"

resultspath <-  yml$DEPath
dir.create(resultspath, showWarnings = FALSE, recursive = TRUE)

# Source your libraries.R script to load all functions

# DE analysis ------------------------------------------------------------------
# Subset seurat object for specific cell type
run_edgeR <- function(seurat.obj, cell_type, resultspath) {
  cell_type_label <- gsub(" ", "_", cell_type)
  
  # Subset cell type
  cat(paste0("Running DE analysis for: ", cell_type))
  s_celltype <- subset(s, subset = type == cell_type)
  
  # Keep samples with min 5 cells
  sample_counts <- table(s_celltype$patient_id, s_celltype$visit)
  keep_samples <- names(which(apply(sample_counts, 1, min) >= 5))

  if (length(keep_samples) < 1) {
    print("Not enough cells for this specific celltype")
    return(NULL)  # Skip the rest of the code
  }

  # Subset remaining samples
  s_subset <- subset(s_celltype, subset = patient_id %in% keep_samples)

  if (length(unique(s_subset$patient_id[s_subset$arm == reference_arm])) < 2 ||
      length(unique(s_subset$patient_id[s_subset$arm == treatment_arm])) < 2) {
    print("Not enough paired samples for this specific celltype")
    return(NULL)  # Skip the rest of the code
  }
  
  # Pseudo bulk 
  s_bulk <- AggregateExpression(
    s_subset,
    return.seurat = TRUE,
    assays = "RNA",
    group.by = c("arm", "patient_id", "visit", "sex", "age"))
  
  # Counts matrix
  pseudo_mat <- as.matrix(GetAssayData(s_bulk, layer = "counts"))
  meta <- s_bulk@meta.data
  
  # Split the row names by underscore and make it a df
  split_names <- strsplit(rownames(meta), "_")
  split_df <- do.call(rbind, split_names)
  
  # Create vectors for each component
  meta$arm <- split_df[, 1]   # First element is arm
  meta$patient_id <- split_df[, 2]  # Second element is patient
  meta$visit <- split_df[, 3]    # Third element is visit
  meta$sex <- split_df[, 4]       # Fourth element is sex
  meta$age <- split_df[, 5]       # Fifth element is age

  meta$arm <- factor(meta$arm)
  meta$arm <- relevel(meta$arm, reference_arm)
  meta$visit <- factor(meta$visit)
  meta$visit <- relevel(meta$visit, baseline_visit)
  
  y <- DGEList(
    counts = pseudo_mat,
    samples = meta,  # Your full metadata data.frame
    group = interaction(meta$arm, meta$visit)  # Create grouping factor
  )
  
  # Create the design matrix starting with patient effects (to account for between-patient variation)
  design <- model.matrix(~patient_id, data = meta)
  
  # Define the treatment effects for each arm
  # These will capture the Pre vs Post effect within each arm
  placebo_post <- meta$arm == reference_arm & meta$visit == followup_visit
  trt_post <- meta$arm == treatment_arm & meta$visit == followup_visit
  
  # Add these effects to the design matrix
  design <- cbind(design, placebo_post, trt_post)
  
  # Determine genes expressed in at least 25% of samples at 1 CPM threshold
  keep <- filterByExpr(y, design)
  y <- y[keep, ]
  
  # Normalize
  y <- normLibSizes(y)
  
  # Estimate dispersion factors and run lm
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  
  pdf(file = paste0(resultspath, "/", cell_type_label, '_bcv.pdf'), height = 6, width = 10)
  print(plotBCV(y, xlab="Average log CPM", ylab="Biological coefficient of variation",
                pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black"))
  dev.off()

  # To find genes that respond differently to treatment vs placebo
  # This tests the interaction between arm and time point
  cat("Fitting model... \n")
  qlf.interaction <- glmQLFTest(fit, contrast = c(0, rep(0, ncol(design)-3), -1, 1))
  counts_norm <- edgeR::cpm(y, log = TRUE, prior.count = 3)
  write_csv(
    rownames_to_column(as.data.frame(counts_norm), var = "gene"),
    paste0(resultspath, "/", cell_type_label, "_countsNorm.csv")
  )
  
  # Lineplot top5 
  cat("Plot top 5 DE genes \n")
  lineplot_top <- function(df, gene_num, path) {
    genes <- rownames(topTags(df))[1:gene_num]
    plot_data <- counts_norm[genes, ] %>%  
      as.data.frame() %>%  
      tibble::rownames_to_column("Gene") %>%  
      pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%   
      separate(Sample, into = c("arm", "patient_id", "visit"), sep = "_")
    
    spag_plots <- lapply(genes, function(x) {
      plot_gene <- plot_data %>% filter(Gene == x)
      p <- ggplot(plot_gene) + 
        aes(x = factor(visit, level = c(baseline_visit, followup_visit)), y = Expression) +
        geom_point(aes(group = patient_id), color = "grey80") +
        geom_line(aes(group = patient_id), color = "grey80") +
        stat_summary(geom = "point") +
        stat_summary(geom = "errorbar", width = 0.15) +
        stat_summary(geom = "line") + 
        facet_wrap(~ arm) +
        theme_bw() +
        ggtitle(x) +
        xlab("Visit")
      return(p)
    })
    
    pdf(file = paste0(resultspath, "/", cell_type_label, path), height = 18, width = 10)
    print(cowplot::plot_grid(plotlist = spag_plots, nrow = gene_num))
    dev.off()
  }
  
  lineplot_top(qlf.interaction, 5, "_tx_top5.pdf")
  
  # Volcano plot
  plot <- function(df, csv_path, pdf_path) {
    de <- as.data.frame(topTags(df, n = nrow(df)))
    de$cell_type <- cell_type
    de$significance <- ifelse(de$logFC > 0 & de$FDR < 0.05, "Upregulated", 
                              ifelse(de$logFC < 0 & de$FDR < 0.05, "Downregulated", "Not Significant"))
    setorder(de, PValue)
    de$gene <- rownames(de)
    de$label <- rownames(de)
    de$model_name <- "Linear Model"
    de$formula <- "~ treatment + end + treatment:end + (1|patient_id)"
    write.csv(de, paste0(resultspath, "/", cell_type_label, csv_path))
    p <- ggplot(de, aes(x=logFC, y=-log10(FDR), color=significance, label=gene)) +
        geom_point(alpha=0.8, size=2) +
        scale_color_manual(values=c(
          "Upregulated"="firebrick",
          "Downregulated"="royalblue",
          "Not Significant"="gray"
        )) +
        geom_text(data=subset(de, -log10(FDR) > 30 & abs(logFC) > 4),
                  aes(label=gene), vjust=1.5, size=3, check_overlap = TRUE) +
        labs(
          x="Effect Size",
          y="-log10(PValue adjusted)",
        ) +
        theme_minimal()

    pdf(file = paste0(resultspath, "/", cell_type_label, pdf_path), height = 12, width = 20)
    print(p)
    dev.off()
  }
  
  cat("Plot Volcano plot \n")
  plot(qlf.interaction, "_tx.csv", "_tx_volcano.pdf")
  
}

# Run edgeR --------------------------------------------------------------------
## input column name 
s$type <- s[[yml$de_annotation]]
cell_types <- unique(s$type)
for (cell_type in cell_types) {
  run_edgeR(seurat.obj = s, cell_type, resultspath = resultspath)
}

# Paste csv file for table -----------------------------------------------------
file_names <- list.files(resultspath, pattern = "\\_tx.csv$", full.names = TRUE)

# Create empty list, save files and merge them
data_frames <- list()
for (file in file_names) {
  data_frames[[file]] <- read.csv(file)
}
merged_data <- dplyr::bind_rows(data_frames)

# Remove unnecessary columns
merged_data <- merged_data[1:7]
colnames(merged_data)[1] <- "gene"
write.csv(merged_data, paste0(resultspath, "/summary_DE_analysis.csv"))
