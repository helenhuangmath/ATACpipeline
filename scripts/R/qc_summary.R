#!/usr/bin/env Rscript
#
# QC summary of all samples: read pairs, alignment rate, mitochondrial reads,
# duplicate rate, nucleosome-free fragments, TSS enrichment, number of peaks,
# FRiP, fragment sizes and TSS profiles.
# Called by scripts/run_all_samples.sh. Output: QC/Summary/

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

#' Reads one small table of each sample, e.g. the alignment stats (step 2).
#'
#' @param samples Sample IDs.
#' @param PathOf Function: sample ID -> file.
#' @return Data frame, one row per sample.
ReadSampleTables <- function(samples, PathOf) {
  tables <- lapply(samples, function(sample) {
    file <- PathOf(sample)
    if (!file.exists(file)) return(data.frame(sample = sample))
    return(as.data.frame(data.table::fread(file, colClasses = "character")))
  })
  return(as.data.frame(data.table::rbindlist(tables, fill = TRUE)))
}

#' Number of cleaned peaks (step 4) of each sample.
#'
#' @param workdir Pipeline WORKDIR.
#' @param samples Sample IDs.
#' @return Integer vector.
CountPeaks <- function(workdir, samples) {
  return(vapply(samples, function(sample) {
    file <- file.path(workdir, "Peak", paste0("peak_", sample, ".final.bed"))
    if (!file.exists(file)) return(NA_integer_)
    return(length(readLines(file)))
  }, integer(1)))
}

#' FRiP of each sample, from the size factor tables of step 7.
#'
#' @param workdir Pipeline WORKDIR.
#' @param groups Groups.
#' @return Data frame: sample, FRiP.
ReadFrip <- function(workdir, groups) {
  tables <- lapply(groups, function(group) {
    file <- file.path(workdir, "ReadCount", group,
                      paste0("sizeFactors_", group, ".txt"))
    if (!file.exists(file)) return(NULL)
    return(as.data.frame(data.table::fread(file))[, c("sample", "FRiP")])
  })
  return(do.call(rbind, tables))
}

#' Bar plot of one QC value of all samples, one panel per group.
#'
#' @param summary QC summary table.
#' @param column Column to plot.
#' @param y_label Axis label.
#' @param file Output PDF.
PlotQcValue <- function(summary, column, y_label, file) {
  if (!column %in% colnames(summary)) return(invisible(NULL))
  summary$value <- suppressWarnings(as.numeric(summary[[column]]))
  summary <- summary[!is.na(summary$value), ]
  if (nrow(summary) == 0) return(invisible(NULL))

  plot <- ggplot2::ggplot(summary,
                          ggplot2::aes(sample, value, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(~ group, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = ConditionColors(summary$condition)) +
    ggplot2::labs(x = NULL, y = y_label, fill = NULL) +
    ThemeNature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  SaveFigure(plot, file, width = min(kDoubleColumn, 30 + 6 * nrow(summary)),
             height = 60)
  return(invisible(file))
}

#' Reads one two-column file (position or size, count) of each sample.
#'
#' @param meta Sample list.
#' @param PathOf Function: sample ID -> file.
#' @param columns Names of the two columns.
#' @param header Whether the files have a header line.
#' @return Data frame with sample, group and condition; NULL if no files.
ReadProfiles <- function(meta, PathOf, columns, header) {
  tables <- lapply(seq_len(nrow(meta)), function(i) {
    file <- PathOf(meta$sample[i])
    if (!file.exists(file) || file.size(file) == 0) return(NULL)
    table <- as.data.frame(data.table::fread(file, header = header))
    colnames(table) <- columns
    table$sample <- meta$sample[i]
    table$group <- meta$group[i]
    table$condition <- meta$condition[i]
    return(table)
  })
  return(do.call(rbind, tables))
}

#' Line plot of one profile per sample (fragment sizes, TSS profiles).
#'
#' @param profiles Data frame with x, y, sample, group, condition.
#' @param x_label,y_label,title Labels.
#' @return ggplot object.
ProfileLinePlot <- function(profiles, x_label, y_label, title) {
  plot <- ggplot2::ggplot(
    profiles, ggplot2::aes(x, y, group = sample, color = condition)
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_wrap(~ group, scales = "free_y") +
    ggplot2::scale_color_manual(values = ConditionColors(profiles$condition)) +
    ggplot2::labs(x = x_label, y = y_label, color = NULL, title = title) +
    ThemeNature()
  return(plot)
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--workdir"),
  optparse::make_option("--sample-list", dest = "sample_list")
)))
workdir <- opts$workdir
out_dir <- file.path(workdir, "QC", "Summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- ReadTable(opts$sample_list)
InWorkdir <- function(...) file.path(workdir, ...)

# 1. QC table of all samples
summary <- meta[, c("sample", "group", "condition", "replicate", "batch")]
align_stats <- ReadSampleTables(meta$sample, function(sample) {
  return(InWorkdir("BAM", paste0(sample, ".align_stats.txt")))
})
tss <- ReadSampleTables(meta$sample, function(sample) {
  return(InWorkdir("QC", "TSS", paste0(sample, ".tss_enrichment.txt")))
})
summary <- merge(summary, align_stats, by = "sample", all.x = TRUE,
                 sort = FALSE)
summary <- merge(summary, tss, by = "sample", all.x = TRUE, sort = FALSE)
summary$n_peaks <- CountPeaks(workdir, summary$sample)
frip <- ReadFrip(workdir, unique(meta$group))
if (!is.null(frip)) {
  summary <- merge(summary, frip, by = "sample", all.x = TRUE, sort = FALSE)
}
summary <- summary[match(meta$sample, summary$sample), ]
WriteTable(summary, file.path(out_dir, "sample_QC_summary.txt"))

# 2. Bar plots
summary$sample <- factor(summary$sample, levels = meta$sample)
qc_values <- c(unique_pairs = "Unique read pairs (fragments)",
               align_rate = "Alignment rate (%)",
               mito_rate = "Mitochondrial read fraction",
               dup_rate = "Duplicate rate",
               nfr_fraction = "Nucleosome-free fragments (< 100 bp)",
               tss_enrichment = "TSS enrichment",
               n_peaks = "Number of peaks",
               FRiP = "FRiP (fraction of fragments in peaks)")
for (column in names(qc_values)) {
  PlotQcValue(summary, column, qc_values[[column]],
              file.path(out_dir, paste0("plot_", column, ".pdf")))
}

# 3. Fragment size distribution: nucleosome-free peak, then one peak per
# nucleosome
fragments <- ReadProfiles(meta, function(sample) {
  return(InWorkdir("BAM", paste0(sample, ".fragment_size.txt")))
}, c("x", "n"), header = FALSE)
if (!is.null(fragments)) {
  fragments$y <- stats::ave(fragments$n, fragments$sample,
                            FUN = function(n) n / sum(n))
  fragments <- fragments[fragments$x <= 1000, ]
  SaveFigure(ProfileLinePlot(fragments, "Fragment size (bp)",
                             "Fraction of fragments", "Fragment size"),
             file.path(out_dir, "plot_fragment_size.pdf"),
             width = kDoubleColumn, height = 60)
}

# 4. TSS profiles: Tn5 cut sites around TSSs, relative to the flanks
profiles <- ReadProfiles(meta, function(sample) {
  return(InWorkdir("QC", "TSS", paste0(sample, ".tss_profile.txt")))
}, c("x", "n"), header = TRUE)
if (!is.null(profiles)) {
  profiles$y <- stats::ave(profiles$n, profiles$sample, FUN = function(n) {
    flank <- mean(c(utils::head(n, 10), utils::tail(n, 10)))
    return(if (flank > 0) n / flank else n * NA)
  })
  SaveFigure(ProfileLinePlot(profiles, "Distance to TSS (bp)",
                             "Cut sites / flanks", "TSS enrichment"),
             file.path(out_dir, "plot_tss_profile.pdf"),
             width = kDoubleColumn, height = 60)
}
message("QC summary: ", file.path(out_dir, "sample_QC_summary.txt"))
