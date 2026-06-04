#!/usr/bin/env Rscript

###############################################################################
#LOAD LIBRARIES#
###############################################################################

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(arrow))

###############################################################################
#GENERIC HELPER FUNCTIONS#
###############################################################################

#checks if directory ends in a backslash and adds one if not
add_back_slash <- function(directory) {
  return(ifelse(
    endsWith(directory, "/"),
    directory,
    paste0(directory, "/")
  ))
}

#checks if directory exists
#option to exit program when the directory is missing
check_dir_exists <- function(directory, force_stop = FALSE) {
  if (!force_stop) {
    return(dir.exists(directory))
  } else {
    if (!dir.exists(directory)) {
      stop(sprintf("Directory %s is missing! Exiting program!", directory))
    }
  }
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

#saves sessionInfo information to outfile (/path/to/sessionInfo.txt)
print_sessioninfo <- function(outfile) {
  writeLines(capture.output(sessionInfo()), outfile)
}

#current seurat version requires a transcripts.csv.gz file
parquet_to_csv <- function(parquet_file, outdir) {
  parquet_df <- arrow::read_parquet(parquet_file)
  #gets the base of the file name and removes file extension
  outstem <- gsub(".parquet", "", basename(parquet_file))
  write.csv(
    parquet_df,
    file = gzfile(paste0(outdir, outstem, ".csv.gz")),
    row.names = FALSE
  )
}

#checks if directory exists and creates directory if needed
create_directory <- function(newdir) {
  #recursive creates all needed directorys on a path
  #so if only /foo/ exists and you ask for /foo/fum/fee
  #it creates both /foo/fum and then /foo/fum/fee
  if (!check_dir_exists(newdir)) {
    dir.create(newdir, recursive = TRUE)
  }
  check_dir_exists(newdir, force_stop = TRUE)
}

create_symlink <- function(infile, symfile) {
  if (!check_file_exists(symfile)) {
    file.symlink(infile, symfile)
  }
}

n_highest <- function(vec, n) {
  sorted_vec <- sort(vec, decreasing = TRUE)
  return(sorted_vec[n])
}

###############################################################################
#PARSE COMMAND LINE ARGUMENTS FUNCTION#
###############################################################################

#parses command line arguments using argparse library
#returns named list of arguments
parse_command_line <- function() {
  #create the parser object
  parser <- argparser::arg_parser(
    "Performs QC and initial cleaning for one individual sample of Xenium spatial transcriptomics data." # nolint
  )
  #set up the arguments
  parser <- argparser::add_argument(
    parser,
    "--indir",
    type = "character",
    help = "Directory to an individual Xenium sample: /path/to/directory/"
  )
  parser <- argparser::add_argument(
    parser,
    "--outdir",
    type = "character",
    help = "Directory to store output data: /path/to/directory/"
  )
  parser <- argparser::add_argument(
    parser,
    "--sample_name",
    type = "character",
    help = "Sample name"
  )
  #using default of 10 reads per cell to remove dead/lowly expressed cells
  #in practice (n=4) this threshold has worked well but may need more inspection
  parser <- argparser::add_argument(
    parser,
    "--min_reads_per_cell",
    type = "integer",
    default = 10,
    help = "Minimum read count filter to remove low expressing cells. Default = 10" # nolint
  )
  parser <- argparser::add_argument(
    parser,
    "--max_reads_per_cell",
    type = "integer",
    default = 1000000,
    help = "Maximum read count filter to remove high expressing cells. Default: 1000000 (essentially no filter)" # nolint
  )
  # using default of 1% to remove cells with high proportions
  # of control probe and control codeword counts
  #in practice (n=4) this threshold has worked well but may need more inspection
  parser <- argparser::add_argument(
    parser,
    "--max_control_prop_per_cell",
    type = "double",
    default = 1.0,
    help = "Maximum proportion of control probes or control codewords per cell. Default = 1%" # nolint
  )
  parser <- argparser::add_argument(
    parser,
    "--metadata_path",
    type = "character",
    default = NULL,
    help = "Metadata to merge with the by sample QC metrics file - expecting a csv with one row of sample characteristics. Default = NULL" # nolint
  )
  parser <- argparser::add_argument(
    parser, 
    "--seed",
    type = "integer",
    default = 42, 
    help = "Random seed for reproducibility (PCA, UMAP, clustering)"
  )
  # I'm not sure what parameter will work best here
  # so lets make it easy to adjust as needed
  parser <- argparser::add_argument(
    parser,
    "--n_highly_var_genes",
    type = "integer",
    default = 1000,
    help = "Number of highly variable genes to select in FindVariableFeatures Seurat function. Default = 1000" # nolint
  )
  parser <- argparser::add_argument(
    parser,
    "--min_quality_per_gene",
    type = "integer",
    default = 20,
    help = "Average quality value for each gene to remain in the dataset: Default: 20" # nolint
  )
  parser <- argparser::add_argument(
    parser,
    "--genericName",
    type = "logical",
    default = FALSE,
    help = "Use generic Xenium_Seurat_Object.rds filename. Default = FALSE"
  )

  #collect command line arguments
  tryCatch(
    {args <- argparser::parse_args(parser)},
    error = function(e) {
      print(parser)
      quit(status = 0)
    }
  )

  return(args)
}

#cleans up and validates input data
validate_input_data <- function(argparse_args) {

  #make sure the directories ends with a backslash so we can save files
  argparse_args$outdir <- add_back_slash(argparse_args$outdir)
  argparse_args$indir <- add_back_slash(argparse_args$indir)

  #check if directories and files are all there
  invisible(lapply(c(argparse_args$outdir,
                     argparse_args$indir,
                     paste0(argparse_args$indir, "cell_feature_matrix/")),
                   check_dir_exists,
                   force_stop = TRUE))

  # Make transcripts file conversion mandatory
  parquet_to_csv(paste0(argparse_args$indir, "transcripts.parquet"),
                 argparse_args$indir)

  #critical files to process the data
  #lets make sure these are all accounted for before we start
  #trying to process the data
  infiles <- paste0(argparse_args$indir,
                    c("cells.csv.gz",
                      "transcripts.parquet", # "transcripts.csv.gz",
                      "metrics_summary.csv",
                      "analysis_summary.html",
                      "cell_feature_matrix/features.tsv.gz",
                      "cell_feature_matrix/barcodes.tsv.gz",
                      "cell_feature_matrix/matrix.mtx.gz"))

  invisible(lapply(infiles, check_file_exists, force_stop = TRUE))

  if (!is.na(argparse_args$metadata_path)) {
    check_file_exists(argparse_args$metadata_path, force_stop = TRUE)
    #makes sure a csv is input
    if (
      (!endsWith(argparse_args$metadata_path, ".csv")) &&
        (!endsWith(argparse_args$metadata_path, ".csv.gz"))
    ) {
      stop(
        sprintf(
          "Metadata file %s is not detected as a csv - please input a csv!",
          argparse_args$metadata_path
        )
      )
    }
  }

  return(argparse_args)
}

###############################################################################
#PLOTTING FUNCTIONS#
###############################################################################

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
      y = .data[[y_col]],
      color = if (!is.null(color_column)) .data[[color_column]] else NULL,
      size = if (!is.null(marker_size_column)) .data[[marker_size_column]] else NULL # nolint
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
  #I think it will complain if we don"t need a legend (ie no color column)
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
  scatter_plot <- ggExtra::ggMarginal(scatter_plot,
                                      type = "histogram",
                                      bins = hist_bins,
                                      groupColour = !is.null(color_column),
                                      groupFill = !is.null(color_column))

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

waterfall_plot <- function(bctbl,
                           col_read_counts,
                           percentile_v,
                           pts_for_extrema = 1000,
                           title = "",
                           x_ax_title = "Barcode Rank Log10",
                           y1_ax_title = "Read Count Log10",
                           y2_ax_title = "Cumulative Read Count Percentile",
                           verbose = FALSE,
                           y_log = TRUE,
                           savefig = NULL) {

  #copy the dataframe so we don"t keep the plot specific variables
  plot_df <- data.frame(bctbl)
  #sort the dataframe by read counts
  #in descending order (so first cell has most reads)
  plot_df <- plot_df[order(plot_df[, col_read_counts], decreasing = TRUE), ]

  #create ranking and cumulative read count variables for plotting
  plot_df[paste0(col_read_counts, "_log10")] <- log10(plot_df[col_read_counts] + 1) # nolint
  plot_df[paste0(col_read_counts, "_rank")] <- seq(1, nrow(plot_df))
  plot_df[paste0(col_read_counts, "_rank_log10")] <- log10(plot_df[paste0(col_read_counts, "_rank")]) # nolint
  plot_df["cumulative_read_per"] <- 100 * (cumsum(plot_df[col_read_counts]) / sum(plot_df[col_read_counts])) # nolint

  #this gets the next value after each percentile cutoff
  #rank_cutoff is the x-axis coordinate of where the red cumulative
  #read line intersects the y2 axis percentile cutoffs defined in percentile_v
  #so its the log10 read count rank at the nth percentile
  #read_cutoff is the y1-axis coordinate of where the blue cumulative read line
  #intersects x-axis coord defined by the y2 axis percentile cutoffs in
  #percentile_v so its the read count at the nth percentile
  rank_cutoff <- lapply(
    percentile_v,
    function(x) {
      plot_df[
        which(plot_df$cumulative_read_per > x),
        paste0(col_read_counts, "_rank_log10")
      ][1]
    }
  )
  names(rank_cutoff) <- percentile_v
  read_cutoff <- lapply(
    percentile_v,
    function(x) {
      plot_df[which(plot_df$cumulative_read_per > x), col_read_counts][1]
    }
  )
  names(read_cutoff) <- percentile_v

  #ggplot doesn't like secondary y axes so we have to compute an adjustment
  #to get the values on the same scale read count y limits
  if (y_log) {
    ylim1 <- c(
      min(plot_df[paste0(col_read_counts, "_log10")]),
      max(plot_df[paste0(col_read_counts, "_log10")])
    )
  } else {
    ylim1 <- c(
      min(plot_df[col_read_counts]),
      max(plot_df[col_read_counts])
    )
  }
  #cumulative read percentage y limits
  ylim2 <- c(0, 100)
  #get the slope of the adjustment line
  b <- diff(ylim1) / diff(ylim2)
  #and the intercept
  a <- ylim1[1] - b * ylim2[1]

  wfall <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[paste0(col_read_counts, "_rank_log10")]]) # nolint
  ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = if (y_log) .data[[paste0(col_read_counts, "_log10")]]
        else .data[[col_read_counts]]
      ), color = "blue"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = a + b * cumulative_read_per), # nolint
      color = "red"
    ) +
    ggplot2::geom_vline(xintercept = as.numeric(rank_cutoff), color = "black") +
    ggplot2::scale_x_continuous(x_ax_title) +
    ggplot2::scale_y_continuous(
      y1_ax_title,
      sec.axis = ggplot2::sec_axis(~ (. - a) / b, name = y2_ax_title)
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(x_ax_title) +
    ggplot2::ylab(y1_ax_title) +
    ggplot2::theme_minimal()

  if (!is.null(savefig)) {
    ggplot2::ggsave(
      filename = savefig,
      plot = wfall,
      bg = "white"
    )
  }

  return(list(
    "rank_cutoff" = rank_cutoff,
    "read_cutoff" = read_cutoff
  ))
}

#plot the distribution of celltype probabilities
#am interested if there is a cut off where the second
#most likely cell is almost as probable as the first
celltype_prob_waterfalls <- function(bycell_pred_df,
                                     col_max_prob = "max",
                                     title = "",
                                     x_ax_title = "Barcode Rank Log10",
                                     y_ax_title = "Celltype Probability",
                                     savefig = NULL) {

  #lets not propogate these changes - copy the dataframe
  plot_df <- data.frame(bycell_pred_df)
  #sort the dataframe by celltype probability
  #(so first cell has highest probability of belonging to one celltype)
  plot_df <- plot_df[order(plot_df[, col_max_prob], decreasing = TRUE), ]

  #create ranking and cumulative read count variables for plotting
  plot_df[paste0(col_max_prob, "_rank")] <- seq(1, nrow(plot_df))
  plot_df[paste0(col_max_prob, "_rank_log10")] <- log10(plot_df[paste0(col_max_prob, "_rank")]) # nolint

  wfall <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[paste0(col_max_prob, "_rank_log10")]]) # nolint
  ) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[col_max_prob]])) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(x_ax_title) +
    ggplot2::ylab(y_ax_title) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()

  if (!is.null(savefig)) {
    ggplot2::ggsave(
      filename = savefig,
      plot = wfall,
      bg = "white"
    )
  }
}

proportional_barplot <- function(in_table,
                                 x_ax_bars_column,
                                 y_ax_prop_column,
                                 figsize = c(7, 7),
                                 savefig = NULL) {

  #make a copy so we don"t mess up the metadata
  tbl <- data.frame(in_table)

  #dummy numeric column to aggregate on
  tbl$count <- seq(1, nrow(tbl))

  #gets counts of y_ax categories within each x_ax category
  counts_tbl <- aggregate(
    as.formula(paste("count ~", x_ax_bars_column, "+", y_ax_prop_column)),
    tbl,
    length
  )

  # Calculate proportions of y_ax categories within each x_ax category
  #start by getting the denominator (counts of the x_ax categories)
  counts_tbl <- merge(
    counts_tbl,
    aggregate(
      as.formula(paste("count ~", x_ax_bars_column)),
      data = counts_tbl, sum
    ),
    by =  x_ax_bars_column,
    suffixes = c("", "_x_total")
  )
  counts_tbl$prop <- 100 * (counts_tbl$count / counts_tbl$count_x_total)

  # Plot the stacked bar chart colored by each y_ax category
  bar <- ggplot2::ggplot(
    counts_tbl,
    ggplot2::aes(
      x = .data[[x_ax_bars_column]], # nolint
      y = prop, # nolint
      fill = .data[[y_ax_prop_column]]
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = x_ax_bars_column, y = "Proportion (%)") +
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

###############################################################################
#SEURAT FUNCTIONS#
###############################################################################

open_seurat <- function(indir, bybc_qc_df, sample_name) {
  #seurat will alter the name of the fov if it contains hyphens
  #want to add the sample name so we can easily keep track of
  #the images after the merge in the next script
  #add "sample_" in front because must start with a non-numeric character
  fov <- gsub("-", "_", paste0("sample_", sample_name, "_fov"))

  seurat_obj <- Seurat::LoadXenium(indir, fov = fov, assay = "Xenium")

  #add some meta data that doesn't seem to get processed
  #this gets us cell and nucleus area as well as centroid location of each cell
  #make sure rownames are the cell_id
  rownames(bybc_qc_df) <- bybc_qc_df$cell_id
  #rename column so it better matches with seurat naming conventions
  bybc_qc_df$nCount_Total <- bybc_qc_df$total_counts
  #merge with the meta data
  seurat_obj@meta.data <- merge(
    seurat_obj@meta.data,
    bybc_qc_df[
      c("nCount_Total", "x_centroid", "y_centroid", "cell_area", "nucleus_area")
    ],
    by = "row.names"
  )
  rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$Row.names
  seurat_obj@meta.data$Row.names <- NULL

  #create a few columns for QC checks
  seurat_obj$Prop_BlankCodeword <- 100 * (seurat_obj$nCount_BlankCodeword / seurat_obj$nCount_Total) # nolint
  seurat_obj$Prop_ControlCodeword <- 100 * (seurat_obj$nCount_ControlCodeword / seurat_obj$nCount_Total) # nolint
  seurat_obj$Prop_ControlProbe <- 100 * (seurat_obj$nCount_ControlProbe / seurat_obj$nCount_Total) # nolint
  seurat_obj$Nucleus_Ratio <- seurat_obj$nucleus_area / seurat_obj$cell_area
  seurat_obj$nCount_Xenium_log10 <- log10(seurat_obj$nCount_Xenium + 1)
  seurat_obj$nFeature_Xenium_log10 <- log10(seurat_obj$nFeature_Xenium + 1)
  seurat_obj$sample <- sample_name

  return(seurat_obj)
}

#create initial QC plots for inspection
plot_qc <- function(seurat_obj, sample_name, bysamp_qc_df, outdir) {
  #create scatter plots of qc metrics
  qc_cols <- c(
    "nFeature_Xenium",
    "Nucleus_Ratio",
    "Prop_BlankCodeword",
    "Prop_ControlCodeword",
    "Prop_ControlProbe",
    "cell_area"
  )
  for (col in qc_cols) {
    scatter_with_hist(
      seurat_obj[[]],
      "nCount_Xenium",
      col,
      title = sample_name,
      savefig = paste0(
        outdir, "plots/initial_QC/",
        sample_name, "_", col, "_by_nCount_Xenium.png"
      )
    )
  }

  scatter_with_hist(
    seurat_obj[[]],
    "Prop_ControlCodeword",
    "Prop_BlankCodeword",
    title = sample_name,
    savefig = paste0(
      outdir, "plots/initial_QC/",
      sample_name, "_Prop_BlankCodeword_by_Prop_ControlCodeword.png"
    )
  )
  scatter_with_hist(
    seurat_obj[[]],
    "Prop_ControlCodeword",
    "Prop_ControlProbe",
    title = sample_name,
    savefig = paste0(
      outdir, "plots/initial_QC/",
      sample_name, "_Prop_ControlProbe_by_Prop_ControlCodeword.png"
    )
  )

  #create waterfall plots of reads per cell and genes per cell
  #we also get threshold of read counts at the 90th,95th, and 99th percentile
  rpc_thresh <- waterfall_plot(
    seurat_obj[[]],
    "nCount_Xenium",
    c(90, 95, 99),
    title = sample_name,
    savefig = paste0(
      outdir, "plots/initial_QC/",
      sample_name, "_nCount_Xenium_waterfall.png"
    )
  )
  gpc_thresh <- waterfall_plot(
    seurat_obj[[]],
    "nFeature_Xenium",
    c(90, 95, 99),
    title = sample_name,
    y1_ax_title = "Gene Count Log10",
    y2_ax_title = "Cumulative Gene Count Percentile",
    savefig = paste0(
      outdir, "plots/initial_QC/",
      sample_name, "_nFeature_Xenium_waterfall.png"
    )
  )

  #store these in a QC metrics dataframe to look at across samples later
  #use 10** to get the rank (cell count) out of log space
  for (thresh in names(rpc_thresh$read_cutoff)) {
    bysamp_qc_df[sprintf("reads_per_cell_thresh_%s", thresh)] <- as.numeric(rpc_thresh$read_cutoff[thresh]) # nolint
    bysamp_qc_df[sprintf("cell_thresh_reads_%s", thresh)] <- 10**as.numeric(rpc_thresh$rank_cutoff[thresh]) # nolint
  }
  for (thresh in names(gpc_thresh$read_cutoff)) {
    bysamp_qc_df[sprintf("features_per_cell_thresh_%s", thresh)] <- as.numeric(gpc_thresh$read_cutoff[thresh]) # nolint
    bysamp_qc_df[sprintf("cell_thresh_features_%s", thresh)] <- 10**as.numeric(gpc_thresh$rank_cutoff[thresh]) # nolint
  }

  return(bysamp_qc_df)
}

create_bygene_df <- function(seurat_obj,
                             outcol_suffix = "",
                             layer = "counts",
                             assay = "all") {
  #can count all assays or specific assay provided
  #using built-in Assays function had unpredictable behavior
  #Would sometimes return a class instead of a list
  if (assay == "all") {
    count_assays <- names(seurat_obj@assays)
  } else {
    if (!assay %in% names(seurat_obj@assays)) {
      stop(sprintf("Assay %s is not in the seurat object provided!", assay))
    }
    count_assays <- assay
  }

  #get all gene names across assays
  #the do.call concatenates into one list
  genes <- sapply(
    count_assays,
    function(x) { Features(seurat_obj, assay = x) } # nolint
  )

  if (assay == "all") genes <- do.call(c, genes)

  bygene_df <- data.frame(matrix(ncol = 0, nrow = length(genes)))
  rownames(bygene_df) <- genes

  #takes some time to extract expression matrices so let just do it once
  exp_mats <- lapply(
    count_assays,
    function(x) { Matrix::Matrix(seurat_obj[[x]][layer], sparse = TRUE) } # nolint
  )

  #extracts the count matrix and then sums across genes
  cnt <- sapply(exp_mats, rowSums)
  col <- paste0("nCount", outcol_suffix)
  bygene_df[col] <- if (assay != "all") as.numeric(cnt) else do.call(c, cnt)
  bygene_df[paste0(col, "_log10")] <- log10(bygene_df[col] + 1)
  #extracts the count matrix and then counts
  #the number of cells with non-zero gene counts
  cells <- sapply(exp_mats, function(x) rowSums(x != 0))
  col <- paste0("nCell", outcol_suffix)
  bygene_df[col] <- if (assay != "all") as.numeric(cells) else do.call(c, cells)
  bygene_df[paste0(col, "_log10")] <- log10(bygene_df[col] + 1)

  return(bygene_df)
}


#Change the filter_gene function to something more simple

filter_genes <- function(seurat_obj,
                         bygene_df,
                         indir,
                         sample_name,
                         outdir,
                         min_qual_per_gene = 20) {

  #open this within the function because it takes up a lot of memory
  bytrans_df <- arrow::read_parquet(
    paste0(indir, "transcripts.parquet"),
    col_select = c(
      "cell_id",
      "overlaps_nucleus",
      "qv",
      "nucleus_distance",
      "feature_name"
    )
  )
  #some transcripts are not assigned to cells
  bytrans_df <- bytrans_df[bytrans_df$cell_id != "UNASSIGNED", ]

  agg_df <- aggregate(
    bytrans_df[c("overlaps_nucleus", "qv", "nucleus_distance")],
    by = list(bytrans_df$feature_name),
    FUN = mean
  )
  names(agg_df)[names(agg_df) == "qv"] <- "average_qv"
  names(agg_df)[names(agg_df) == "nucleus_distance"] <- "average_nucleus_distance" # nolint
  agg_df$feature_name <- agg_df$Group.1

  rm(bytrans_df)
  invisible(gc())

  #seurat automatically changes the gene names - adjust these to match
  rownames(agg_df) <- agg_df$feature_name

  #mark everything as a "Gene" unless it matches control names below
  bygene_df$transcript_type <- "Gene"

  #merge in the new columns from agg_df
  merged_df <- merge(
    bygene_df,
    agg_df[c(
      "feature_name",
      "overlaps_nucleus",
      "average_qv",
      "average_nucleus_distance"
    )],
    by.x = "row.names",
    by.y = "feature_name",
    all.x = TRUE
  )
  row.names(merged_df) <- merged_df$Row.names
  merged_df$Row.names <- NULL

  merged_df$transcript_type[grepl("^NegControlCodeword", rownames(merged_df))] <- "Negative control codeword" # nolint
  merged_df$transcript_type[grepl("^NegControlProbe", rownames(merged_df))]   <- "Negative control probe" # nolint
  merged_df$transcript_type[grepl("^UnassignedCodeword", rownames(merged_df))] <- "Unassigned codeword" # nolint

  # Only attempt the scatter plot if nCount_log10 and average_qv are present
  if (
    "nCount_log10" %in% colnames(merged_df) &&
      "average_qv" %in% colnames(merged_df)
  ) {
    scatter_with_hist(
      merged_df,
      "nCount_log10",
      "average_qv",
      color_column = "transcript_type",
      hline = 20,
      scatter_marker_size = 1,
      title = sample_name,
      legend_title = "transcript_type",
      savefig = sprintf(
        "%splots/filtering_QC/%s_average_qv_by_nCount_log10_bygene_transcript_type.png", # nolint
        outdir,
        sample_name
      )
    )
  }

  # Only assign filt_low_qv if we have both transcript_type and average_qv
  if (
    "transcript_type" %in% colnames(merged_df) &&
      "average_qv" %in% colnames(merged_df)
  ) {
    merged_df$filt_low_qv <- merged_df$transcript_type == "Gene" &
      merged_df$average_qv < min_qual_per_gene
  } else {
    merged_df$filt_low_qv <- FALSE
  }

  # Filter out low-qv genes (keep everything else)
  keep_rows <- rownames(merged_df[!merged_df$filt_low_qv, ])
  seurat_obj <- subset(seurat_obj, features = keep_rows)

  return(list(
    "seurat_obj" = seurat_obj,
    "bygene_df" = merged_df
  ))
}

#lowered nfeatures in FindVariableFeatures to 1000 since we only have 5k genes
##############################################################################
#CHECK
#leaving as is for now but might want algorithm <- 4
#on FindClusters, right now set to Louvain (1)
##############################################################################
run_seurat <- function(seurat_obj,
                       norm_meth = "LogNormalize",
                       scale_factor = 10000,
                       high_var_select_meth = "vst",
                       n_high_var_genes = 1000,
                       neigh_umap_dims = 1:30,
                       cluster_res = 0.6,
                       cluster_algo = 1,
                       seed = 42) {

  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = norm_meth,
    scale.factor = scale_factor
  )
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    selection.method = high_var_select_meth,
    nfeatures = n_high_var_genes
  )
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(
    seurat_obj,
    features = Seurat::VariableFeatures(object = seurat_obj),
    seed.use = seed
  )
  seurat_obj <- Seurat::FindNeighbors(
    seurat_obj,
    dims = neigh_umap_dims,
    reduction = "pca"
  )
  seurat_obj <- Seurat::RunUMAP(
    seurat_obj, 
    dims = neigh_umap_dims, 
    seed.use = seed
    )
  seurat_obj <- Seurat::FindClusters(
    seurat_obj,
    resolution = cluster_res,
    algorithm = cluster_algo,
    random.seed = seed
  )
  return(seurat_obj)
}

#Change the gene hard coded
plot_processed_data <- function(seurat_obj,
                                col_read_counts,
                                sample_name,
                                outdir) {
  #start by plotting correlations with the PCs and the read counts
  pc_df <- data.frame(Seurat::Embeddings(seurat_obj, reduction = "pca"))
  #add the counts to the pc dataframe for plotting
  pc_df <- merge(
    pc_df,
    seurat_obj@meta.data[c(col_read_counts, paste0(col_read_counts, "_log10"))],
    by = "row.names"
  )
  #run against the top 10 PCs
  for (i in seq(1, 10)) {
    scatter_with_hist(
      pc_df,
      paste0(col_read_counts, "_log10"),
      sprintf("PC_%i", i),
      title = sample_name,
      savefig = sprintf(
        "%splots/PC_plots/%s_PC_%i_by_%s_log10.png",
        outdir, sample_name, i, col_read_counts
      )
    )
  }

  #plot and save the umaps
  umap <- Seurat::DimPlot(
    seurat_obj,
    reduction = "umap",
    label = TRUE
  ) +
    ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    umap,
    filename = sprintf(
      "%splots/seurat_output/%s_umap_cellclusters.png",
      outdir, sample_name
    )
  )
  #could put this in with the bygene plot but don"t want
  #to scale up to the total counts across genes
  readcount_umap <- Seruat::FeaturePlot(
    seurat_obj,
    features = paste0(col_read_counts, "_log10")
  ) + ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    readcount_umap,
    filename = sprintf(
      "%splots/seurat_output/%s_umap_readcounts.png",
      outdir, sample_name
    )
  )

  spat <- Seurat::ImageDimPlot(seurat_obj,
    fov = gsub("-", "_", paste0("sample_", sample_name, "_fov")),
    cols = "polychrome",
  ) + ggplot2::ggtitle(sample_name)
  ggplot2::ggsave(
    spat,
    filename = sprintf(
      "%splots/seurat_output/%s_insitu_cellclusters.png",
      outdir, sample_name
    ),
    width = 21,
    bg = "white"
  )
}

get_cell_markers <- function(seurat_obj,
                             sample_name,
                             outdir,
                             min_pct = .25,
                             logfc_thresh = .25,
                             keep_n_markers = 20,
                             plot_n_markers = 5) {

  markers <- Seurat::FindAllMarkers(seurat_obj,
                                    min.pct = min_pct,
                                    only.pos = TRUE,
                                    logfc.threshold = logfc_thresh)

  #group by each cluster and get the top n markers
  byclust <- split(markers, markers$cluster)
  #sorts so top FC is at the top
  byclust <- lapply(byclust, function(x) {
    x <- x[order(x[, "avg_log2FC"], decreasing = TRUE), ]
    return(x)
  })
  #top n markers for dot plots
  v_smol_markers <- do.call(rbind, lapply(byclust, head, n = plot_n_markers))
  top_genes <- unique(v_smol_markers$gene)

  dot_plot <- Seurat::DotPlot(seurat_obj,
                              features = top_genes) + Seurat::RotatedAxis()
  ggplot2::ggsave(dot_plot,
                  filename = sprintf(
                    "%splots/seurat_output/%s_cluster_markers_dotplot.png",
                    outdir, sample_name
                  ),
                  bg = "white",
                  width = 7 * plot_n_markers)

  #code will fail if not enough top markers in each cluster
  #setting as the min of the input value and the shortest bycluster marker list
  #sometimes when a sample is not clustering well certain clusters will have
  # no diff expr genes in this case, require the keep_n_markers to be at least
  # 1 to keep code from failing
  keep_n_markers <- max(min(keep_n_markers, min(sapply(byclust, nrow))), 1)
  #top n markers to create a wide dataframe
  smol_markers <- do.call(rbind, lapply(byclust, head, n = keep_n_markers))
  rownames(smol_markers) <- NULL

  #make an id variable for reshape
  if (keep_n_markers > 1) {
    smol_markers$rank <- rep(seq(1, keep_n_markers), length(byclust))
  } else {
    smol_markers$rank <- rep(1, nrow(smol_markers))
  }

  #switch to wide format for easier viewing
  markers_wide <- reshape(smol_markers,
                          idvar = "rank",
                          timevar = "cluster",
                          direction = "wide")

  return(markers_wide)
}


###############################################################################
#DRIVER FUNCTION#
###############################################################################

#main driver function
Process_Individual_Xenium_Sample <- function() { # nolint
  #gets the command line arguments into a named list
  args <- parse_command_line()

  #checks input data for issues and checks key files exist
  args <- validate_input_data(args)

  #sets RNG seed for reproducibility across all stochastic steps 
  set.seed(args$seed)

  #create new folders if needed within the out directory to keep things tidy
  folders <- c("plots/PC_plots",
               "plots/initial_QC/",
               "plots/filtering_QC/",
               "plots/seurat_output/",
               "plots/label_transfer_QC/",
               "data",
               "sessionInfo")
  dirs <- paste0(args$outdir,
                 folders)
  #invisible hides the T/F vector output from lapply
  #(ie directory exists or not)
  invisible(lapply(dirs, create_directory))

  #extracts the sample name from the 10x directory
  sample_name <- args$sample_name

  #save sessionInfo before we start
  #would do this sooner but need to parse the outfile first
  print_sessioninfo(sprintf(
    "%ssessionInfo/Process_Individual_Xenium_Sample_%s_sessionInfo_start.txt",
    args$outdir, sample_name
  ))

  #symlink key 10x files so critical data is all together
  data_files_10x <- c("analysis_summary.html")
  files_to_link <- paste0(args$indir, data_files_10x)
  linked_files <- paste0(args$outdir, "data/", sample_name, "_", data_files_10x)
  for (i in seq_along(files_to_link)) {
    create_symlink(files_to_link[i],
                   linked_files[i])
  }

  bysamp_qc <- read.csv(paste0(args$indir, "metrics_summary.csv"))
  bybc_qc <- read.csv(paste0(args$indir, "cells.csv.gz"))
  rownames(bybc_qc) <- bybc_qc$cell_id

  #merges any additional per sample metadata information
  if (!is.na(args$metadata_path)) {
    metadata <- read.csv(args$metadata_path)
    bysamp_qc <- cbind(bysamp_qc, metadata)
  }

  #open seurat object and add some QC and metadata columns
  sample_data <- open_seurat(args$indir,
                             bybc_qc,
                             sample_name)

  #create initial QC plots of the data
  bysamp_qc <- plot_qc(sample_data,
                       sample_name,
                       bysamp_qc,
                       args$outdir)

  #creates a bygene data frame to store gene level info
  #akin to anndata.var
  bygene_counts <- create_bygene_df(sample_data)

  sample_data_bygene_counts <- filter_genes(
    sample_data,
    bygene_counts,
    args$indir,
    sample_name,
    args$outdir,
    min_qual_per_gene = args$min_quality_per_gene
  )
  sample_data <- sample_data_bygene_counts$seurat_obj
  bygene_counts <- sample_data_bygene_counts$bygene_df

  #normalize and create low dimensional embeddings of data
  sample_data <- run_seurat(sample_data,
                            n_high_var_genes = args$n_highly_var_genes,
                            seed = args$seed)


  markers_wide <- get_cell_markers(sample_data,
                                   sample_name,
                                   args$outdir)

  #save to the outdirectory
  #Add an if to change the name of the file to generic of sample name,
  #by default it is the sample name but can be changed to a generic name
  #"Xenium_Seurat_Object.rds"
  if (args$genericName) {
    output_file <- file.path(args$outdir,
                             "Xenium_Seurat_Object_processed.rds")
  } else {
    output_file <- file.path(args$outdir,
                             sprintf("%s_processed.rds", sample_name))
  }
  saveRDS(sample_data, output_file)
  write.csv(sample_data@meta.data,
            file = sprintf("%sdata/%s_metadata.csv", args$outdir, sample_name))
  write.csv(bybc_qc,
            file = gzfile(sprintf("%sdata/%s_cells.csv.gz",
                                  args$outdir, sample_name)),
            row.names = FALSE)
  write.csv(bysamp_qc,
            file = sprintf("%sdata/%s_metrics_summary.csv",
                           args$outdir, sample_name),
            row.names = FALSE)
  write.csv(bygene_counts,
            file = sprintf("%sdata/%s_bygene_tbl.csv",
                           args$outdir, sample_name))
  write.csv(markers_wide,
            file = sprintf("%sdata/%s_cluster_markers_wide.csv",
                           args$outdir, sample_name),
            row.names = FALSE)
  print("Before running Integrate_All_Xenium_Samples, add a column called celltype to the metadata!") # nolint
  print(sprintf("Processing complete for sample %s!", sample_name))

  #save sessionInfo at the end too - will it change? Why would it?
  print_sessioninfo(sprintf(
    "%ssessionInfo/Process_Individual_Xenium_Sample_%s_sessionInfo_end.txt",
    args$outdir, sample_name
  ))

  return(args$outdir)

}

###############################################################################
#RUN DRIVER SCRIPT#
###############################################################################
future_mem_limit <- 8e+09
options(future.globals.maxSize = future_mem_limit)
outdir <- Process_Individual_Xenium_Sample()
