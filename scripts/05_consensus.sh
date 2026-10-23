#!/bin/bash
#
# Step 5: consensus peaks of each condition and union peaks of one group.
#
# Usage: 05_consensus.sh <config.sh> <group>
# Output: Consensus/<group>/peak_consensus_<condition>.bed
#         Consensus/<group>/peak_<group>_union.bed (Chr Start End peakN Width .)
#         Consensus/<group>/peak_<group>_union.saf (for featureCounts)
#         Consensus/<group>/plot_*.pdf (Venn diagrams, peak width)

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local group="$2"
  activate_env tools

  local out_dir="${WORKDIR}/Consensus/${group}"
  mkdir -p "${out_dir}"
  cd "${out_dir}" || exit 1
  rm -f ./*.bed ./*.saf ./*.txt ./*.pdf  # results of an earlier run

  local conditions=() cond
  while read -r cond; do
    conditions+=("${cond}")
  done < <(conditions_of "${group}")

  # 5.1 Consensus peaks: in >= MIN_REPS replicates of a condition. Replicate
  # peaks closer than CONSENSUS_GAP count as one. Column 4 of
  # rep_support.bed = number of replicates with the peak, column 5 = which.
  local sample reps n_reps min_reps
  for cond in "${conditions[@]}"; do
    reps=()
    while read -r sample; do
      reps+=("${sample}")
    done < <(samples_of "${group}" "${cond}")
    n_reps="${#reps[@]}"
    min_reps=$(( n_reps < MIN_REPS ? n_reps : MIN_REPS ))

    for sample in "${reps[@]}"; do
      awk -v sample="${sample}" 'BEGIN { OFS = "\t" }
          { print $1, $2, $3, sample }' \
        "${WORKDIR}/Peak/peak_${sample}.final.bed"
    done \
      | sort -k1,1 -k2,2n \
      | bedtools merge -d "${CONSENSUS_GAP}" -c 4,4 \
        -o count_distinct,distinct -i - > "peak_${cond}.rep_support.bed"

    awk -v min_reps="${min_reps}" 'BEGIN { OFS = "\t" }
        $4 >= min_reps { print $1, $2, $3 }' \
      "peak_${cond}.rep_support.bed" > "peak_consensus_${cond}.bed"
    log "${group} ${cond}: $(wc -l < "peak_consensus_${cond}.bed")" \
      "consensus peaks (in >= ${min_reps} of ${n_reps} replicates)"
  done

  # 5.2 Union peaks: merge consensus peaks of all conditions, remove narrow
  # peaks, add peak IDs
  local union="peak_${group}_union.bed"
  for cond in "${conditions[@]}"; do
    cat "peak_consensus_${cond}.bed"
  done \
    | sort -k1,1 -k2,2n \
    | bedtools merge -d "${UNION_GAP}" -i - \
    | awk -v min_width="${MIN_WIDTH}" -v group="${group}" '
        BEGIN { OFS = "\t" }
        $3 - $2 >= min_width {
          n++
          print $1, $2, $3, group "_peak" n, $3 - $2, "."
        }' > "${union}"
  [[ -s "${union}" ]] || die "${group}: no union peaks; check peak calling"
  log "${group}: $(wc -l < "${union}") union peaks"

  # SAF format for featureCounts: Start is 1-based
  awk 'BEGIN { OFS = "\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
       { print $4, $1, $2 + 1, $3, "." }' "${union}" > "peak_${group}_union.saf"

  # Union peaks found in each condition (for the Venn diagram)
  for cond in "${conditions[@]}"; do
    bedtools intersect -u -a "${union}" -b "peak_consensus_${cond}.bed" \
      | cut -f 4 > "peak_${cond}_in_union.txt"
  done

  # 5.3 Plots: replicate Venn, condition Venn, peak widths
  activate_env R
  Rscript "${R_CODE_DIR}/peak_summary.R" --dir="${out_dir}" \
    --group="${group}" --conditions="$(IFS=,; echo "${conditions[*]}")" \
    --min-reps="${MIN_REPS}"
}

main "$@"
