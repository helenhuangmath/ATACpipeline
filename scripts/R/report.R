#!/usr/bin/env Rscript
#
# Step 10: interactive HTML report of all groups. Search a gene (or peak ID or
# locus) to see its peaks, their differential statistics, signal profiles
# and normalized counts. One self-contained file: it opens offline in any
# browser. Called by scripts/10_report.sh; writes the report given by --out
# (Report/ATAC_report.html).

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

kReportBins <- 20  # bins of the signal profile in the report (+-1 kb)

#' Averages neighbouring bins of a signal matrix down to n_bins.
#'
#' @param values Matrix, peaks x bins.
#' @param n_bins Number of bins to keep.
#' @return Matrix, peaks x n_bins.
CoarsenBins <- function(values, n_bins) {
  bin_group <- rep(seq_len(n_bins), each = ncol(values) / n_bins)
  coarse <- sapply(seq_len(n_bins), function(i) {
    return(rowMeans(values[, bin_group == i, drop = FALSE]))
  })
  return(matrix(coarse, nrow = nrow(values)))
}

#' Collects the results of one group for the report.
#'
#' @param workdir Pipeline WORKDIR.
#' @param group Group.
#' @param meta Sample list rows of this group.
#' @return List, converted to JSON for the report.
GroupData <- function(workdir, group, meta) {
  anno <- as.data.frame(data.table::fread(file.path(
    workdir, "Annotation", group, paste0("peak_", group, "_union.anno.txt")
  )))
  peak_ids <- anno$peakN
  Column <- function(name) {
    if (name %in% colnames(anno)) return(anno[[name]])
    return(rep(NA, nrow(anno)))
  }

  # All genes linked to each peak, for searching
  genes <- Column("Gene.Name")
  genes_file <- file.path(workdir, "Annotation", group,
                          paste0("peak_", group, "_union.anno_allgenes.txt"))
  if (file.exists(genes_file)) {
    all_genes <- as.data.frame(data.table::fread(genes_file))
    all_genes <- all_genes[!is.na(all_genes$Gene.Name), ]
    linked <- tapply(all_genes$Gene.Name, all_genes$peakN,
                     function(g) paste(unique(g), collapse = ";"))
    genes <- unname(linked[peak_ids])
  }

  # Normalized counts, samples in sample list order
  norm <- as.data.frame(data.table::fread(file.path(
    workdir, "ReadCount", group, paste0("NormRC_", group, ".txt")
  )))
  norm <- norm[match(peak_ids, norm$peakN), ]
  counts <- round(as.matrix(norm[, meta$sample, drop = FALSE]), 1)

  # Differential statistics of each comparison
  de <- list()
  files <- list.files(file.path(workdir, "DiffPeak", group),
                      pattern = paste0("^deseq2_All_", group, "_.*\\.txt$"),
                      full.names = TRUE)
  for (file in files) {
    result <- as.data.frame(data.table::fread(file))
    result <- result[match(peak_ids, result$peakN), ]
    de[[result$comparison[!is.na(result$comparison)][1]]] <- list(
      lfc = round(result$log2FoldChange, 2), padj = signif(result$padj, 3),
      dir = result$direction
    )
  }

  # Mean signal profile of each condition around each peak center
  conditions <- unique(meta$condition)
  signal <- ReadSignalMatrix(file.path(workdir, "Figures", group,
                                       paste0("signal_", group, ".matrix.gz")))
  by_condition <- ConditionSignal(signal,
                                  stats::setNames(meta$condition, meta$sample))
  profile <- do.call(cbind, lapply(conditions, function(cond) {
    values <- by_condition[[cond]][match(peak_ids, signal$regions), ,
                                   drop = FALSE]
    return(round(CoarsenBins(values, kReportBins), 2))
  }))
  bins <- CoarsenBins(matrix(signal$bin_centers, nrow = 1), kReportBins)

  return(list(
    samples = meta$sample, sampleCondition = meta$condition,
    conditions = conditions, colors = as.list(ConditionColors(conditions)),
    comparisons = names(de), bins = as.vector(bins),
    peaks = list(id = peak_ids, chr = anno$Chr, start = anno$Start,
                 end = anno$End, anno = Column("Anno"),
                 dist = Column("Distance.to.TSS"), gene = Column("Gene.Name"),
                 genes = genes),
    counts = counts, de = de, profile = profile
  ))
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--workdir"),
  optparse::make_option("--sample-list", dest = "sample_list"),
  optparse::make_option("--species", default = "mm10"),
  optparse::make_option("--project", default = "ATAC-seq"),
  optparse::make_option("--version", default = ""),
  optparse::make_option("--out")
)))

meta <- ReadTable(opts$sample_list)
groups <- unique(meta$group)
data <- list(project = opts$project, species = opts$species,
             version = opts$version, date = format(Sys.Date()), groups = list())
for (group in groups) {
  data$groups[[group]] <- GroupData(opts$workdir, group,
                                    meta[meta$group == group, ])
}

json <- as.character(jsonlite::toJSON(data, auto_unbox = TRUE, digits = NA,
                                      na = "null"))
json <- gsub("</", "<\\/", json, fixed = TRUE)  # keep the <script> block intact

template <- paste(readLines(file.path(script_dir, "report_template.html")),
                  collapse = "\n")
parts <- strsplit(template, "__DATA__", fixed = TRUE)[[1]]
html <- paste0(parts[1], json, parts[2])
html <- gsub("__TITLE__", paste(opts$project, "ATAC-seq results"), html,
             fixed = TRUE)
dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
writeLines(html, opts$out)
message("Report: ", opts$out, " (", round(file.size(opts$out) / 1e6, 1),
        " MB)")
