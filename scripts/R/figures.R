#!/usr/bin/env Rscript
#
# Step 9: Nature-style summary figure of one group (183 mm wide, 7 pt text,
# panels a-i). Called by scripts/09_figures.sh.
#
#   a  peaks of each sample                      f  differential peaks
#   b  FRiP of each sample                       g  volcano plots
#   c  TSS enrichment of each sample             h  average signal profiles
#   d  genomic features of consensus peaks       i  signal heatmap
#   e  PCA of samples
#
# Output (Figures/<group>/): Figure_<group>_summary.pdf,
#   signal_profile_<group>.pdf, signal_heatmap_<group>.pdf

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

kFeatureLevels <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "TTS",
                    "Non-Coding", "Intergenic")

#' a. Number of peaks of each sample.
PeakNumberPlot <- function(workdir, meta) {
  counts <- data.frame(sample = meta$sample, condition = meta$condition)
  counts$peaks <- vapply(counts$sample, function(sample) {
    file <- file.path(workdir, "Peak", paste0("peak_", sample, ".final.bed"))
    if (!file.exists(file)) return(NA_integer_)
    return(length(readLines(file)))
  }, integer(1))
  counts$sample <- factor(counts$sample, levels = meta$sample)
  counts$condition <- factor(counts$condition,
                             levels = unique(meta$condition))

  plot <- ggplot2::ggplot(counts,
                          ggplot2::aes(sample, peaks, fill = condition)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = ConditionColors(meta$condition)) +
    ggplot2::labs(x = NULL, y = "Peaks", fill = NULL) +
    ThemeNature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                   legend.position = "none")
  return(plot)
}

#' b, c. One value of each sample (e.g. FRiP), by condition: dots and mean.
#'
#' @param table Data frame with sample, condition and the value column.
#' @param column Value column.
#' @param y_label Axis label.
#' @param conditions Conditions in plot order.
#' @return ggplot object, or NULL if there are no values.
SampleValuePlot <- function(table, column, y_label, conditions) {
  table$value <- suppressWarnings(as.numeric(table[[column]]))
  table <- table[!is.na(table$value), ]
  if (nrow(table) == 0) return(NULL)
  table$condition <- factor(table$condition, levels = conditions)
  plot <- ggplot2::ggplot(table,
                          ggplot2::aes(condition, value, color = condition)) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.5,
                          linewidth = 0.25, color = "black") +
    ggplot2::geom_point(size = 1.2,
                        position = ggplot2::position_jitter(0.1, seed = 1)) +
    ggplot2::scale_color_manual(values = ConditionColors(conditions)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(x = NULL, y = y_label) +
    ThemeNature() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  return(plot)
}

#' TSS enrichment of each sample (step 3).
ReadTssEnrichment <- function(workdir, meta) {
  meta$tss_enrichment <- vapply(meta$sample, function(sample) {
    file <- file.path(workdir, "QC", "TSS",
                      paste0(sample, ".tss_enrichment.txt"))
    if (!file.exists(file)) return(NA_character_)
    return(as.character(data.table::fread(file)$tss_enrichment[1]))
  }, character(1))
  return(meta)
}

#' d. Genomic features of the consensus peaks of each condition (% stacked).
FeaturePlot <- function(consensus_dir, anno, conditions) {
  if (!"Anno" %in% colnames(anno)) return(NULL)
  features <- do.call(rbind, lapply(conditions, function(cond) {
    ids <- readLines(file.path(consensus_dir,
                               paste0("peak_", cond, "_in_union.txt")))
    return(data.frame(condition = cond,
                      Anno = anno$Anno[match(ids, anno$peakN)]))
  }))
  features <- features[!is.na(features$Anno), ]
  if (nrow(features) == 0) return(NULL)
  features$Anno <- factor(features$Anno,
                          levels = intersect(kFeatureLevels, features$Anno))
  features$condition <- factor(features$condition, levels = conditions)

  plot <- ggplot2::ggplot(features, ggplot2::aes(condition, fill = Anno)) +
    ggplot2::geom_bar(position = "fill", width = 0.7) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::labs(x = NULL, y = "Consensus peaks", fill = NULL) +
    ThemeNature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  return(plot)
}

#' g. Volcano plots of all comparisons.
VolcanoPanel <- function(diff_dir, group) {
  files <- list.files(diff_dir, pattern = paste0("^deseq2_All_", group, "_"),
                      full.names = TRUE)
  if (length(files) == 0) return(NULL)
  results <- do.call(rbind, lapply(files, function(file) {
    table <- as.data.frame(data.table::fread(file))
    columns <- c("peakN", "Gene.Name", "log2FoldChange", "padj", "direction",
                 "comparison")
    return(table[, intersect(columns, colnames(table))])
  }))
  results$comparison <- gsub("_vs_", "\nvs ", results$comparison)
  plot <- VolcanoPlot(results, NULL, n_label = 0) +
    ggplot2::facet_wrap(~ comparison, nrow = 1) +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::theme(legend.position = "none")
  return(plot)
}

#' Groups peaks by their direction in any comparison.
#'
#' @return Named vector: class of each peak (Up / Down / Not differential).
PeakClasses <- function(diff_dir, group, peaks) {
  files <- list.files(diff_dir, pattern = paste0("^deseq2_Sig_", group,
                                                 "_.*_vs_.*\\.txt$"),
                      full.names = TRUE)
  classes <- stats::setNames(rep("Not differential", length(peaks)), peaks)
  for (file in files) {
    table <- as.data.frame(data.table::fread(file))
    classes[table$peakN[table$direction == "Down"]] <- "Down"
  }
  for (file in files) {
    table <- as.data.frame(data.table::fread(file))
    classes[table$peakN[table$direction == "Up"]] <- "Up"
  }
  return(factor(classes, levels = c("Up", "Down", "Not differential")))
}

#' h. Average signal around peak centers, by condition and peak class.
ProfilePlot <- function(signal_by_condition, bin_centers, classes) {
  profiles <- do.call(rbind, lapply(names(signal_by_condition), function(cond) {
    values <- signal_by_condition[[cond]]
    by_class <- lapply(levels(classes), function(class) {
      in_class <- names(classes)[classes == class]
      if (length(in_class) == 0) return(NULL)
      return(data.frame(condition = cond, class = class,
                        position = bin_centers,
                        signal = colMeans(values[in_class, , drop = FALSE])))
    })
    return(do.call(rbind, by_class))
  }))
  profiles$condition <- factor(profiles$condition,
                               levels = names(signal_by_condition))
  profiles$class <- factor(profiles$class, levels = levels(classes))

  plot <- ggplot2::ggplot(profiles,
                          ggplot2::aes(position / 1000, signal,
                                       color = condition)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::facet_wrap(~ class, nrow = 1) +
    ggplot2::scale_color_manual(
      values = ConditionColors(names(signal_by_condition))
    ) +
    ggplot2::labs(x = "Distance to peak center (kb)", y = "Mean signal (CPM)",
                  color = NULL) +
    ThemeNature()
  return(plot)
}

#' i. Heatmap of the signal at differential peaks (all peaks if none).
HeatmapPlot <- function(signal_by_condition, bin_centers, classes,
                        max_peaks = 3000) {
  peaks <- names(classes)[classes != "Not differential"]
  if (length(peaks) == 0) peaks <- names(classes)
  total <- Reduce(`+`, lapply(signal_by_condition, function(m) {
    return(rowMeans(m[peaks, , drop = FALSE]))
  }))
  peaks <- utils::head(peaks[order(total, decreasing = TRUE)], max_peaks)

  conditions <- names(signal_by_condition)
  tiles <- do.call(rbind, lapply(conditions, function(cond) {
    values <- signal_by_condition[[cond]][peaks, , drop = FALSE]
    return(data.frame(condition = cond, peak = rep(peaks, ncol(values)),
                      class = rep(as.character(classes[peaks]), ncol(values)),
                      position = rep(bin_centers, each = nrow(values)),
                      signal = as.vector(values)))
  }))
  tiles$peak <- factor(tiles$peak, levels = rev(peaks))
  tiles$condition <- factor(tiles$condition, levels = conditions)
  cap <- stats::quantile(tiles$signal, 0.99)
  tiles$signal <- pmin(tiles$signal, cap)

  plot <- ggplot2::ggplot(tiles,
                          ggplot2::aes(position / 1000, peak, fill = signal)) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(class ~ condition, scales = "free_y",
                        space = "free_y") +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") +
    ggplot2::labs(x = "Distance to peak center (kb)", y = NULL,
                  fill = "Signal\n(CPM)") +
    ThemeNature() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank())
  return(plot)
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--group"),
  optparse::make_option("--workdir"),
  optparse::make_option("--sample-list", dest = "sample_list")
)))
group <- opts$group
workdir <- opts$workdir
out_dir <- file.path(workdir, "Figures", group)
diff_dir <- file.path(workdir, "DiffPeak", group)
count_dir <- file.path(workdir, "ReadCount", group)

meta <- ReadTable(opts$sample_list)
meta <- meta[meta$group == group, ]
conditions <- unique(meta$condition)
sample_condition <- stats::setNames(meta$condition, meta$sample)
anno <- as.data.frame(data.table::fread(file.path(
  workdir, "Annotation", group, paste0("peak_", group, "_union.anno.txt")
)))
vst <- as.data.frame(data.table::fread(file.path(count_dir,
                                                 paste0("VST_", group,
                                                        ".txt"))))
size_factors <- as.data.frame(data.table::fread(
  file.path(count_dir, paste0("sizeFactors_", group, ".txt"))
))

signal <- ReadSignalMatrix(file.path(out_dir,
                                     paste0("signal_", group, ".matrix.gz")))
signal_by_condition <- ConditionSignal(signal, sample_condition)
classes <- PeakClasses(diff_dir, group, signal$regions)

# Panels
panels <- list(
  a = PeakNumberPlot(workdir, meta),
  b = SampleValuePlot(size_factors, "FRiP", "FRiP", conditions),
  c = SampleValuePlot(ReadTssEnrichment(workdir, meta), "tss_enrichment",
                      "TSS enrichment", conditions),
  d = FeaturePlot(file.path(workdir, "Consensus", group), anno, conditions),
  e = PcaPlot(as.matrix(vst[, meta$sample]), meta$condition, NULL)
)
summary_file <- file.path(diff_dir, paste0("deseq2_summary_", group, ".txt"))
if (file.exists(summary_file)) {
  summary <- as.data.frame(data.table::fread(summary_file))
  summary$comparison <- gsub("_vs_", " vs ", summary$comparison)
  panels$f <- DiffSummaryPlot(summary, NULL)
}
volcano <- VolcanoPanel(diff_dir, group)
profile <- ProfilePlot(signal_by_condition, signal$bin_centers, classes)
heatmap <- HeatmapPlot(signal_by_condition, signal$bin_centers, classes)
panels <- panels[!vapply(panels, is.null, logical(1))]

# Layout: 3 small panels per row, then volcano plots, profiles and heatmap
RowOf <- function(widths) {
  names <- intersect(names(widths), names(panels))
  if (length(names) == 0) return(NULL)
  return(patchwork::wrap_plots(panels[names], nrow = 1,
                               widths = widths[names]))
}
rows <- list(RowOf(c(a = 1.4, b = 0.8, c = 0.8)),
             RowOf(c(d = 1, e = 1, f = 0.8)), volcano, profile, heatmap)
heights <- c(1, 1, 0.9, 0.8, 1.6)[!vapply(rows, is.null, logical(1))]
rows <- rows[!vapply(rows, is.null, logical(1))]
figure <- patchwork::wrap_plots(rows, ncol = 1, heights = heights) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste(group, "ATAC-seq summary"),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold")
    )
  )
SaveFigure(figure, file.path(out_dir, paste0("Figure_", group, "_summary.pdf")),
           width = kDoubleColumn, height = 247)
SaveFigure(profile, file.path(out_dir, paste0("signal_profile_", group,
                                              ".pdf")),
           width = kDoubleColumn, height = 55)
SaveFigure(heatmap, file.path(out_dir, paste0("signal_heatmap_", group,
                                              ".pdf")),
           width = kDoubleColumn, height = 120)
message(group, ": summary figure in ", out_dir)
