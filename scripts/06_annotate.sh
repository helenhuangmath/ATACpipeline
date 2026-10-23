#!/bin/bash
#
# Step 6: annotation of the union peaks of one group: genomic feature,
# nearest gene and genes nearby (ChIPseeker; Homer optional).
#
# Usage: 06_annotate.sh <config.sh> <group>
# Output: Annotation/<group>/peak_<group>_union.anno.txt (one row per peak)
#         Annotation/<group>/peak_<group>_union.anno_allgenes.txt

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local group="$2"

  local out_dir="${WORKDIR}/Annotation/${group}"
  local union="${WORKDIR}/Consensus/${group}/peak_${group}_union.bed"
  mkdir -p "${out_dir}"

  # 6.1 Homer (optional): nearest gene
  local homer_anno=""
  if [[ -n "${HOMER_GENOME}" ]]; then
    activate_env tools
    homer_anno="${out_dir}/peak_${group}_union.homer_anno"
    annotatePeaks.pl "${union}" "${HOMER_GENOME}" > "${homer_anno}" \
      2> "${out_dir}/homer_${group}.log"
  fi

  # 6.2 ChIPseeker: promoter = TSS - TSS_UP to TSS + TSS_DOWN
  activate_env R
  Rscript "${R_CODE_DIR}/annotate_peaks.R" \
    --peaks="${union}" --txdb="${TXDB}" --genes="${GENE_TABLE}" \
    --homer="${homer_anno}" --tss-up="${TSS_UP}" --tss-down="${TSS_DOWN}" \
    --flank="${FLANK}" --out="${out_dir}/peak_${group}_union"
}

main "$@"
