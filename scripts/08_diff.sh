#!/bin/bash
#
# Step 8: differential peaks of one group with DESeq2, for each comparison in
# comparison_list.txt. Uses the raw count matrix and size factors of step 7.
#
# Usage: 08_diff.sh <config.sh> <group>
# Output: DiffPeak/<group>/deseq2_All_/deseq2_Sig_<group>_<trt>_vs_<ref>.txt,
#         _Up.bed, _Down.bed, plots, deseq2_summary_<group>.txt

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local group="$2"
  activate_env R

  local count_dir="${WORKDIR}/ReadCount/${group}"
  mkdir -p "${WORKDIR}/DiffPeak/${group}"

  Rscript "${R_CODE_DIR}/deseq2_diff.R" \
    --group="${group}" \
    --rawrc="${count_dir}/RawRC_${group}.txt" \
    --size-factors="${count_dir}/sizeFactors_${group}.txt" \
    --sample-list="${SAMPLE_LIST}" \
    --comparison-list="${COMPARISON_LIST}" \
    --design="${DESIGN}" \
    --min-count="${MIN_COUNT}" \
    --padj="${PADJ}" \
    --lfc="${LFC}" \
    --outdir="${WORKDIR}/DiffPeak/${group}"
}

main "$@"
