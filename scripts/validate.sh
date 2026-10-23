#!/bin/bash
#
# Checks config.sh, the sample list, the comparison list and the reference
# files before any job is submitted (run_pipeline.sh runs it first).
#
# Usage: validate.sh <config.sh>
# Exit status: 0 if the inputs are good, 1 if there are errors.

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

readonly HEADER="sample fastq1 fastq2 group condition replicate batch"
n_errors=0

#######################################
# Prints an input error and counts it.
# Globals:
#   n_errors
# Arguments:
#   Message.
#######################################
input_error() {
  err "$*"
  n_errors=$(( n_errors + 1 ))
}

#######################################
# Prints a warning.
# Arguments:
#   Message.
#######################################
warning() {
  log "WARNING: $*"
}

#######################################
# Returns 0 if a name has only letters, numbers, _ . and -.
# Arguments:
#   Name.
#######################################
is_good_name() {
  [[ "$1" =~ ^[A-Za-z0-9_.-]+$ ]]
}

#######################################
# Checks the header, columns and sample IDs of the sample list.
# Globals:
#   SAMPLE_LIST
#######################################
check_sample_list_format() {
  local header
  header="$(read_table "${SAMPLE_LIST}" | head -1 | tr '\t' ' ')"
  if [[ "${header}" != "${HEADER}" ]]; then
    die "The sample list header must be (tab-separated): ${HEADER}"
  fi

  local row
  while read -r row; do
    input_error "Sample ${row} does not have 7 tab-separated columns"
  done < <(read_table "${SAMPLE_LIST}" \
    | awk -F'\t' 'NR > 1 && NF != 7 { print $1 }')

  local sample_id
  while IFS= read -r sample_id; do
    is_good_name "${sample_id}" \
      || input_error "Sample ID '${sample_id}' may only have letters," \
        "numbers, _ . -"
  done < <(read_table "${SAMPLE_LIST}" | awk -F'\t' 'NR > 1 { print $1 }')

  while read -r sample_id; do
    input_error "Sample ID ${sample_id} is used more than once"
  done < <(all_samples | sort | uniq -d)

  (( n_errors == 0 )) || die "Please fix the sample list first."
}

#######################################
# Checks names and fastq files of every sample.
# Globals:
#   SAMPLE_LIST, DATADIR
#######################################
check_samples() {
  local samples=() sample_id column value read fastq
  while read -r sample_id; do samples+=("${sample_id}"); done < <(all_samples)

  for sample_id in "${samples[@]}"; do
    for column in group condition; do
      value="$(sample_info "${sample_id}" "${column}")"
      is_good_name "${value}" || input_error "${sample_id}: ${column}" \
        "'${value}' may only have letters, numbers, _ . -"
    done

    for read in 1 2; do
      local fastqs
      IFS=',' read -ra fastqs <<< "$(sample_info "${sample_id}" "fastq${read}")"
      for fastq in "${fastqs[@]}"; do
        [[ "${fastq}" == *.gz ]] \
          || input_error "${sample_id}: fastq must be gzipped: ${fastq}"
        [[ -f "$(full_path "${fastq}")" ]] || input_error "${sample_id}:" \
          "fastq${read} not found: $(full_path "${fastq}")"
      done
    done
  done
}

#######################################
# Checks the replicates of each group; prints a summary.
# Globals:
#   MIN_REPS
#######################################
check_groups() {
  local group cond n_reps summary
  while read -r group; do
    summary=""
    while read -r cond; do
      n_reps="$(samples_of "${group}" "${cond}" | wc -l | tr -d ' ')"
      summary+=" ${cond}=${n_reps}"
      if (( n_reps < MIN_REPS )); then
        warning "${group} ${cond} has ${n_reps} replicate(s);" \
          "its consensus peaks use all of them"
      fi
    done < <(conditions_of "${group}")
    log "${group}:${summary}"
  done < <(all_groups)
}

#######################################
# Checks that both conditions of each comparison exist for the group.
# Globals:
#   COMPARISON_LIST
#######################################
check_comparisons() {
  if [[ -z "${COMPARISON_LIST}" ]]; then
    warning "COMPARISON_LIST is empty: no differential peaks"
    return 0
  fi
  if [[ ! -f "${COMPARISON_LIST}" ]]; then
    input_error "Comparison list not found: ${COMPARISON_LIST}"
    return 0
  fi

  local comparison_group treatment reference group conditions
  while read -r comparison_group treatment reference; do
    local groups=("${comparison_group}")
    if [[ "${comparison_group}" == "*" ]]; then
      groups=()
      while read -r group; do groups+=("${group}"); done < <(all_groups)
    fi
    for group in "${groups[@]}"; do
      conditions=" $(conditions_of "${group}" | tr '\n' ' ')"
      if [[ "${conditions}" != *" ${treatment} "* \
        || "${conditions}" != *" ${reference} "* ]]; then
        warning "${group}: ${treatment} vs ${reference} is skipped," \
          "a condition is not found"
      fi
    done
  done < <(read_table "${COMPARISON_LIST}" \
    | awk -F'\t' 'NR > 1 { print $1, $2, $3 }')
}

#######################################
# Checks the peak calling settings.
# Globals:
#   PEAK_MODE
#######################################
check_peak_settings() {
  case "${PEAK_MODE}" in
    fragment | cutsite) ;;
    *) input_error "PEAK_MODE must be fragment or cutsite" ;;
  esac
}

#######################################
# Checks the reference files and the scheduler.
# Globals:
#   BOWTIE_IDX, GNM_SIZES, BLACKLIST, CHR1TOY_BED, TSS_BED, TXDB, SCHEDULER
#######################################
check_reference_and_cluster() {
  [[ -f "${BOWTIE_IDX}.1.bt2" || -f "${BOWTIE_IDX}.1.bt2l" ]] \
    || input_error "bowtie2 index not found: ${BOWTIE_IDX}"
  [[ -f "${GNM_SIZES}" ]] || input_error "GNM_SIZES not found: ${GNM_SIZES}"
  [[ -f "${BLACKLIST}" ]] || input_error "BLACKLIST not found: ${BLACKLIST}"
  [[ -f "${CHR1TOY_BED}" ]] \
    || input_error "CHR1TOY_BED not found: ${CHR1TOY_BED}"
  [[ -f "${TSS_BED}" ]] || warning "TSS_BED not found: no TSS enrichment" \
    "(run reference/build_reference.sh to make it)"
  [[ -f "${TXDB}" ]] || warning "TXDB not found: no ChIPseeker annotation"

  case "${SCHEDULER}" in
    lsf) command -v bsub > /dev/null || input_error "bsub not found" ;;
    slurm) command -v sbatch > /dev/null || input_error "sbatch not found" ;;
    local) ;;
    *) input_error "SCHEDULER must be lsf, slurm or local" ;;
  esac
}

main() {
  read_config "$1"
  [[ -f "${SAMPLE_LIST}" ]] || die "Sample list not found: ${SAMPLE_LIST}"

  check_sample_list_format
  check_samples
  check_groups
  check_comparisons
  check_peak_settings
  check_reference_and_cluster

  (( n_errors == 0 )) || die "${n_errors} error(s). Please fix them."
  log "Input check passed: $(all_samples | wc -l | tr -d ' ') samples"
}

main "$@"
