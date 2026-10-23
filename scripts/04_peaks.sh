#!/bin/bash
#
# Step 4: peak calling of one sample with MACS2. Weak peaks and blacklist
# peaks are removed.
#
# Usage: 04_peaks.sh <config.sh> <sample_id>
# Output: Peak/MACS2/                  raw MACS2 output
#         Peak/peak_<id>.final.bed     cleaned peaks used for the analysis;
#                                      columns: chr, start, end, name, score,
#                                      strand, fold enrichment, -log10(p),
#                                      -log10(q)

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

#######################################
# Calls peaks with MACS2, on whole fragments or on Tn5 cut sites
# (PEAK_MODE).
# Globals:
#   PEAK_MODE, MACS2_ARGS, REFSIZE
# Arguments:
#   Sample ID.
# Outputs:
#   Path of the MACS2 peak file.
#######################################
call_macs2() {
  local sample_id="$1"
  local bam="BAM/${sample_id}.map.dedup.bam"
  local macs2_args
  read -ra macs2_args <<< "${MACS2_ARGS}"
  mkdir -p Peak/MACS2

  # fragment : -f BAMPE, the real fragment of each read pair
  # cutsite  : Tn5 cut sites, each extended to 150 bp centered on the cut
  #            (--shift -75 --extsize 150), as ENCODE does
  local input=(-t "${bam}" -f BAMPE)
  local cut_sites="Peak/MACS2/peak_${sample_id}.cutsites.bed"
  if [[ "${PEAK_MODE}" == cutsite ]]; then
    tn5_cut_sites "${bam}" > "${cut_sites}"
    input=(-t "${cut_sites}" -f BED --nomodel --shift -75 --extsize 150)
  fi

  # --keep-dup all : duplicates were already handled in step 2
  log "${sample_id}: macs2 (${PEAK_MODE}) ${macs2_args[*]}"
  macs2 callpeak "${input[@]}" -g "${REFSIZE}" --keep-dup all \
    "${macs2_args[@]}" --outdir Peak/MACS2 -n "peak_${sample_id}" \
    2> "Peak/MACS2/peak_${sample_id}.macs2.log"
  rm -f "${cut_sites}"

  local peak_type=narrowPeak
  if [[ " ${macs2_args[*]} " == *" --broad "* ]]; then peak_type=broadPeak; fi
  echo "Peak/MACS2/peak_${sample_id}_peaks.${peak_type}"
}

#######################################
# Removes peaks with low fold enrichment and peaks in blacklist regions.
# Arguments:
#   MACS2 peak file, minimum fold enrichment, output BED.
#######################################
clean_peaks() {
  awk -v min_fold="$2" 'BEGIN { OFS = "\t" }
      !/^(#|track|browser)/ && $7 >= min_fold {
        print $1, $2, $3, $4, $5, ".", $7, $8, $9
      }' "$1" \
    | bedtools intersect -v -a - -b "${BLACKLIST}" \
    | sort -k1,1 -k2,2n > "$3"
}

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p Peak

  local raw_peaks final="Peak/peak_${sample_id}.final.bed"
  raw_peaks="$(call_macs2 "${sample_id}")"
  clean_peaks "${raw_peaks}" "${MIN_FOLD}" "${final}"
  log "${sample_id}: $(awk '!/^#/' "${raw_peaks}" | wc -l) peaks," \
    "$(wc -l < "${final}") after cleaning"
}

main "$@"
