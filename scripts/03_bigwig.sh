#!/bin/bash
#
# Step 3: normalized bigwig track (for IGV / UCSC) and TSS enrichment of one
# sample.
#
# Usage: 03_bigwig.sh <config.sh> <sample_id>
# Output: Wig/<sample_id>.cpm.bw
#         QC/TSS/<sample_id>.tss_profile.txt     Tn5 cut sites around TSSs
#         QC/TSS/<sample_id>.tss_enrichment.txt  TSS enrichment score

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

#######################################
# TSS enrichment: Tn5 cut sites around all TSSs, in 10-bp bins from -2 kb to
# +2 kb (upstream negative). Score = highest bin within 500 bp of the TSS /
# mean of the 10 outer bins on each side (1.9-2 kb away), as by ENCODE.
# Globals:
#   TSS_BED
# Arguments:
#   Sample ID.
#######################################
tss_enrichment() {
  local out="QC/TSS/$1"
  mkdir -p QC/TSS

  # Columns of bedtools window: cut site (1-6), TSS (7-12; 8 = position,
  # 12 = strand)
  tn5_cut_sites "BAM/$1.map.dedup.bam" \
    | bedtools window -w 2000 -a - -b "${TSS_BED}" \
    | awk 'BEGIN { OFS = "\t" }
        {
          d = ($12 == "-" ? $8 - $2 : $2 - $8)
          if (d >= -2000 && d < 2000) n[int((d + 2000) / 10)]++
        }
        END {
          print "position", "cut_sites"
          for (b = 0; b < 400; b++) print b * 10 - 1995, n[b] + 0
        }' > "${out}.tss_profile.txt"

  awk -v sample="$1" 'NR > 1 { pos[NR - 1] = $1; n[NR - 1] = $2 }
      END {
        bins = NR - 1
        for (i = 1; i <= 10; i++) flank += n[i] + n[bins - i + 1]
        flank /= 20
        for (i = 1; i <= bins; i++) {
          if (pos[i] > -500 && pos[i] < 500 && n[i] > top) top = n[i]
        }
        print "sample\ttss_enrichment"
        if (flank > 0) {
          printf "%s\t%.2f\n", sample, top / flank
        } else {
          print sample "\tNA"
        }
      }' "${out}.tss_profile.txt" > "${out}.tss_enrichment.txt"
  log "$1: TSS enrichment $(awk 'NR == 2 { print $2 }' \
    "${out}.tss_enrichment.txt")"
}

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p Wig

  # 3.1 bamCoverage, CPM (counts per million fragments)
  #   --extendReads --samFlagInclude 64 : each fragment counted once (read 1),
  #                                       over the whole fragment
  #   -bl : skip blacklist regions
  local bw_args
  read -ra bw_args <<< "${BW_ARGS}"
  bamCoverage -b "BAM/${sample_id}.map.dedup.bam" \
    -o "Wig/${sample_id}.cpm.bw" --normalizeUsing CPM \
    --extendReads --samFlagInclude 64 -bl "${BLACKLIST}" -p "${THREADS}" \
    "${bw_args[@]}" 2> "Wig/${sample_id}.bamCoverage.log"

  # 3.2 TSS enrichment (quality of the ATAC-seq signal)
  if [[ -f "${TSS_BED}" ]]; then
    tss_enrichment "${sample_id}"
  else
    log "TSS_BED not found: no TSS enrichment"
  fi
}

main "$@"
