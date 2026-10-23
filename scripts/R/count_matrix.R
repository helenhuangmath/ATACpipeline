#!/usr/bin/env Rscript
#
# Step 7: raw and normalized read count matrices of all samples of one group.
# Called by scripts/07_count_matrix.sh.
#
# Output (in --outdir):
#   RawRC_<group>.txt            raw fragment counts, union peaks x samples
#   NormRC_<group>.txt           counts normalized with DESeq2 size factors
#   NormRC_<group>.CondMean.txt  mean normalized counts of each condition
#   VST_<group>.txt              variance-stabilized counts (log2 scale)
#   sizeFactors_<group>.txt      size factor, fragments in peaks and FRiP of
#                                each sample (step 8 uses these size factors)
#   plot_PCA_<group>.pdf, plot_Corr_<group>.pdf   sample QC (VST counts)

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

#' Reads featureCounts output; sample columns are named by sample ID.
#'
#' @param file featureCounts output file.
#' @return Data frame: Geneid, Chr, Start, End, Strand, Length, samples.
ReadFeatureCounts <- function(file) {
  counts <- as.data.frame(data.table::fread(file, skip = "Geneid"))
  sample_columns <- 7:ncol(counts)
  colnames(counts)[sample_columns] <- sub(
    "\\.map\\.dedup\\.bam$", "", basename(colnames(counts)[sample_columns])
  )
  return(counts)
}

#' Fraction of fragments in peaks (FRiP) from the featureCounts summary.
#'
#' @param file featureCounts summary file.
#' @return Named vector: FRiP of each sample.
ReadFrip <- function(file) {
  summary <- as.data.frame(data.table::fread(file))
  colnames(summary)[-1] <- sub("\\.map\\.dedup\\.bam$", "",
                               basename(colnames(summary)[-1]))
  assigned <- unlist(summary[summary$Status == "Assigned", -1])
  total <- colSums(summary[, -1, drop = FALSE])
  return(round(assigned / total, 4))
}

#' DESeq2 size factors (median of ratios); "poscounts" if every peak has a
#' zero count in some sample.
#'
#' @param dds DESeqDataSet.
#' @return DESeqDataSet with size factors.
EstimateSizeFactors <- function(dds) {
  dds <- tryCatch(
    DESeq2::estimateSizeFactors(dds),
    error = function(e) DESeq2::estimateSizeFactors(dds, type = "poscounts")
  )
  return(dds)
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--group"),
  optparse::make_option("--counts", help = "featureCounts output"),
  optparse::make_option("--anno", default = "", help = "peak annotation"),
  optparse::make_option("--sample-list", dest = "sample_list"),
  optparse::make_option("--outdir")
)))
group <- opts$group
OutFile <- function(name) file.path(opts$outdir, name)

# 1. Samples of this group and their raw counts
meta <- ReadTable(opts$sample_list)
meta <- meta[meta$group == group, ]
feature_counts <- ReadFeatureCounts(opts$counts)

peaks <- data.frame(peakN = feature_counts$Geneid, Chr = feature_counts$Chr,
                    Start = feature_counts$Start - 1, End = feature_counts$End)
peaks$Width <- peaks$End - peaks$Start
if (IsFile(opts$anno)) {
  annotation <- as.data.frame(data.table::fread(opts$anno))
  annotation_columns <- intersect(
    c("peakN", "Anno", "Distance.to.TSS", "Gene.Name"), colnames(annotation)
  )
  peaks <- merge(peaks, annotation[, annotation_columns], by = "peakN",
                 all.x = TRUE, sort = FALSE)
  peaks <- peaks[match(feature_counts$Geneid, peaks$peakN), ]
}

raw_counts <- as.matrix(feature_counts[, meta$sample, drop = FALSE])
rownames(raw_counts) <- feature_counts$Geneid

# 2. Normalization with DESeq2 size factors
col_data <- data.frame(condition = factor(meta$condition),
                       row.names = meta$sample)
dds <- DESeq2::DESeqDataSetFromMatrix(raw_counts, col_data, design = ~ 1)
dds <- EstimateSizeFactors(dds)
norm_counts <- round(DESeq2::counts(dds, normalized = TRUE), 2)

# 3. Write the matrices
WriteTable(cbind(peaks, raw_counts), OutFile(paste0("RawRC_", group, ".txt")))
WriteTable(cbind(peaks, norm_counts),
           OutFile(paste0("NormRC_", group, ".txt")))

conditions <- unique(meta$condition)
condition_means <- sapply(conditions, function(cond) {
  in_cond <- meta$sample[meta$condition == cond]
  return(rowMeans(norm_counts[, in_cond, drop = FALSE]))
})
condition_means <- matrix(round(condition_means, 2), ncol = length(conditions),
                          dimnames = list(NULL, conditions))
WriteTable(cbind(peaks, condition_means),
           OutFile(paste0("NormRC_", group, ".CondMean.txt")))

frip <- ReadFrip(paste0(opts$counts, ".summary"))
size_factors <- data.frame(
  sample = meta$sample, condition = meta$condition,
  size_factor = round(DESeq2::sizeFactors(dds)[meta$sample], 4),
  fragments_in_peaks = colSums(raw_counts)[meta$sample],
  FRiP = frip[meta$sample]
)
WriteTable(size_factors, OutFile(paste0("sizeFactors_", group, ".txt")))

# 4. Variance-stabilized counts (log2 scale, for heatmaps, PCA, clustering)
vst_counts <- VstMatrix(dds)
WriteTable(cbind(peaks, round(vst_counts, 3)),
           OutFile(paste0("VST_", group, ".txt")))

# 5. Sample QC: PCA and correlation of VST counts
SaveFigure(PcaPlot(vst_counts, meta$condition, paste("PCA", group)),
           OutFile(paste0("plot_PCA_", group, ".pdf")))
PlotCorrelation(vst_counts, meta$condition, paste("Correlation", group),
                OutFile(paste0("plot_Corr_", group, ".pdf")))

message(group, ": count matrices of ", ncol(raw_counts), " samples x ",
        nrow(raw_counts), " peaks")
