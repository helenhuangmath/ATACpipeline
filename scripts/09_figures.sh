#!/bin/bash
#
# Step 9: signal around each union peak (deepTools) and a Nature-style
# summary figure of one group.
#
# Usage: 09_figures.sh <config.sh> <group>
# Output: Figures/<group>/signal_<group>.matrix.gz   signal of every sample,
#                                                    +-1 kb of each peak center
#         Figures/<group>/Figure_<group>_summary.pdf  multi-panel figure
#         Figures/<group>/signal_profile_<group>.pdf, signal_heatmap_<group>.pdf

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local group="$2"
  activate_env tools

  local out_dir="${WORKDIR}/Figures/${group}"
  mkdir -p "${out_dir}"

  local bigwigs=() labels=() sample
  while read -r sample; do
    bigwigs+=("${WORKDIR}/Wig/${sample}.cpm.bw")
    labels+=("${sample}")
  done < <(samples_of "${group}")

  # 9.1 Signal of each sample around the center of each union peak
  #   -b / -a 1000 : 1 kb before and after the center, in 20-bp bins
  #   --sortRegions keep : same order as the union peaks
  computeMatrix reference-point --referencePoint center -b 1000 -a 1000 \
    --binSize 20 --missingDataAsZero --sortRegions keep \
    -R "${WORKDIR}/Consensus/${group}/peak_${group}_union.bed" \
    -S "${bigwigs[@]}" --samplesLabel "${labels[@]}" -p "${THREADS}" \
    -o "${out_dir}/signal_${group}.matrix.gz" \
    2> "${out_dir}/computeMatrix_${group}.log"

  # 9.2 Summary figure
  activate_env R
  Rscript "${R_CODE_DIR}/figures.R" --group="${group}" --workdir="${WORKDIR}" \
    --sample-list="${SAMPLE_LIST}"
}

main "$@"
