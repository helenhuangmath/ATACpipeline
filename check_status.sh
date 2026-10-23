#!/bin/bash
#
# Shows which steps are finished for each sample and group, and which job
# logs have errors.
#
# Usage: check_status.sh <config.sh>

source "$(dirname "${BASH_SOURCE[0]}")/scripts/common.sh"

#######################################
# Prints "done" if a step is finished, otherwise ".".
# Globals:
#   WORKDIR
# Arguments:
#   Step name.
#######################################
status_of() {
  if [[ -f "${WORKDIR}/Status/$1.done" ]]; then echo "done"; else echo "."; fi
}

main() {
  read_config "$1"
  local sample_id group log_file

  echo
  printf '  %-30s %s\n' "FastQC/MultiQC of raw fastq" "$(status_of qc)"
  echo
  echo "Samples ('.' = not finished)"
  printf '  %-30s %-6s %-6s %-11s %-6s\n' sample trim align bigwig+TSS peaks
  while read -r sample_id; do
    printf '  %-30s %-6s %-6s %-11s %-6s\n' "${sample_id}" \
      "$(status_of "${sample_id}.trim")" "$(status_of "${sample_id}.align")" \
      "$(status_of "${sample_id}.bigwig")" "$(status_of "${sample_id}.peaks")"
  done < <(all_samples)

  echo
  echo "Groups"
  printf '  %-30s %-7s %-5s\n' group matrix diff
  while read -r group; do
    printf '  %-30s %-7s %-5s\n' "${group}" \
      "$(status_of "${group}.matrix")" "$(status_of "${group}.diff")"
  done < <(all_groups)
  printf '  %-30s %-7s\n' "QC summary" "$(status_of qc_summary)"
  printf '  %-30s %-7s\n' "HTML report" "$(status_of report)"

  echo
  while read -r log_file; do
    echo "Error in ${log_file}:"
    grep ERROR "${log_file}" | tail -2 | sed 's/^/    /'
  done < <(grep -l ERROR "${WORKDIR}"/Log/*.err 2> /dev/null || true)
  echo
}

main "$@"
