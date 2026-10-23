#!/bin/bash
#
# Step 7: read count matrix of all samples of one group in its union peaks.
# Counts fragments of every sample (featureCounts), then writes the raw and
# the normalized count matrix, size factors, FRiP and a sample PCA.
#
# Usage: 07_count_matrix.sh <config.sh> <group>
# Output: ReadCount/<group>/RawRC_<group>.txt   raw counts, peaks x samples
#         ReadCount/<group>/NormRC_<group>.txt  DESeq2-normalized counts
#         ReadCount/<group>/sizeFactors_<group>.txt, plot_PCA_<group>.pdf, ...

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local group="$2"
  activate_env tools

  local out_dir="${WORKDIR}/ReadCount/${group}"
  mkdir -p "${out_dir}"

  local bams=()
  local sample
  while read -r sample; do
    bams+=("${WORKDIR}/BAM/${sample}.map.dedup.bam")
  done < <(samples_of "${group}")

  # 7.1 Count fragments (read pairs) of each sample in each union peak
  #   -p --countReadPairs : count fragments, not reads (subread >= 2.0.2)
  featureCounts -p --countReadPairs -F SAF -T "${THREADS}" \
    -a "${WORKDIR}/Consensus/${group}/peak_${group}_union.saf" \
    -o "${out_dir}/featureCounts_${group}.txt" "${bams[@]}" \
    2> "${out_dir}/featureCounts_${group}.log"

  # 7.2 Raw and normalized count matrices, FRiP, PCA
  activate_env R
  Rscript "${R_CODE_DIR}/count_matrix.R" \
    --group="${group}" \
    --counts="${out_dir}/featureCounts_${group}.txt" \
    --anno="${WORKDIR}/Annotation/${group}/peak_${group}_union.anno.txt" \
    --sample-list="${SAMPLE_LIST}" \
    --outdir="${out_dir}"
}

main "$@"
