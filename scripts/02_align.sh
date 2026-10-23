#!/bin/bash
#
# Step 2: alignment of one sample: bowtie2 -> filter -> remove mitochondrial
# reads -> remove duplicates.
#
# Usage: 02_align.sh <config.sh> <sample_id>
# Output: BAM/<sample_id>.map.dedup.bam (+ .bai), .bowtie2.log, .idxstats.txt,
#         .align_stats.txt, .fragment_size.txt

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p BAM

  local trim1="Fastq/${sample_id}_R1.trim.fastq.gz"
  local trim2="Fastq/${sample_id}_R2.trim.fastq.gz"
  local sort_bam="BAM/${sample_id}.map.sort.bam"
  local dedup_bam="BAM/${sample_id}.map.dedup.bam"
  local bowtie2_args
  read -ra bowtie2_args <<< "${BOWTIE2_ARGS}"

  # 2.1 Align (ENCODE ATAC-seq settings) and keep good read pairs
  #   --very-sensitive           : end-to-end alignment (adapters are trimmed)
  #   --no-mixed --no-discordant : both reads of a pair align together
  #   -X 2000                    : fragments up to 2 kb (several nucleosomes)
  #   -F 1804 : remove unmapped, mate unmapped, secondary, QC-fail, duplicate
  #   -f 2    : keep properly paired reads
  #   -q      : minimum mapping quality, removes multi-mapping reads
  #   fixmate : flags reads whose mate was removed; the 2nd filter drops them
  bowtie2 "${bowtie2_args[@]}" -p "${THREADS}" -x "${BOWTIE_IDX}" \
    -1 "${trim1}" -2 "${trim2}" 2> "BAM/${sample_id}.bowtie2.log" \
    | samtools view -u -F 1804 -f 2 -q "${MIN_MAPQ}" - \
    | samtools fixmate -m - - \
    | samtools view -u -F 1804 -f 2 - \
    | samtools sort -@ "${THREADS}" -m 1G -T "BAM/${sample_id}.tmp" \
      -o "${sort_bam}" -
  samtools index "${sort_bam}"

  # 2.2 Reads on each chromosome. ATAC-seq has many mitochondrial reads (Tn5
  # cuts the open mitochondrial DNA); they are counted, then removed.
  local idxstats="BAM/${sample_id}.idxstats.txt"
  samtools idxstats "${sort_bam}" > "${idxstats}"
  local good_pairs mito_pairs filtered_pairs
  good_pairs="$(awk '{ n += $3 } END { print n / 2 }' "${idxstats}")"
  mito_pairs="$(awk -v chr="${MITO_CHR}" '$1 == chr { n += $3 }
      END { print n / 2 }' "${idxstats}")"
  filtered_pairs="$(awk 'NR == FNR { keep[$1] = 1; next }
      $1 in keep { n += $3 } END { print n / 2 }' \
    "${CHR1TOY_BED}" "${idxstats}")"

  # 2.3 Keep chr1 to chrY (removes chrM, random and chrUn), then remove PCR
  # duplicates (REMOVE_DUPLICATES=false: only mark them). The BAM is read
  # as a stream (-), which keeps it sorted for markdup.
  local dedup_args=(-s -@ "${THREADS}")
  if [[ "${REMOVE_DUPLICATES}" == true ]]; then dedup_args+=(-r); fi
  samtools view -u -L "${CHR1TOY_BED}" - < "${sort_bam}" \
    | samtools markdup "${dedup_args[@]}" - "${dedup_bam}" \
      2> "BAM/${sample_id}_dedup.log"
  samtools index "${dedup_bam}"
  samtools flagstat "${dedup_bam}" > "${dedup_bam}_flagstat"

  # 2.4 Fragment sizes (insert size of read 1): columns size, count.
  # ATAC-seq shows a nucleosome-free peak (< 100 bp) and peaks every ~200 bp
  # (1, 2, 3 nucleosomes).
  samtools view -f 66 "${dedup_bam}" \
    | awk '{ s = ($9 < 0 ? -$9 : $9); n[s]++ }
           END { for (s in n) print s "\t" n[s] }' \
    | sort -k1,1n > "BAM/${sample_id}.fragment_size.txt"

  # 2.5 Read pairs at each step (read 1 only, -f 64)
  #   input_pairs    : read pairs given to bowtie2
  #   mito_rate      : fraction of good pairs (2.1) on MITO_CHR
  #   filtered_pairs : good pairs on chr1-chrY
  #   unique_pairs   : filtered pairs that are not PCR duplicates
  #   nfr_fraction   : fraction of unique fragments < 100 bp (nucleosome-free)
  local input_pairs align_rate unique_pairs mito_rate dup_rate nfr_fraction
  input_pairs="$(awk '/reads; of these:/ { print $1 }' \
    "BAM/${sample_id}.bowtie2.log")"
  align_rate="$(awk '/overall alignment rate/ { sub("%", "", $1); print $1 }' \
    "BAM/${sample_id}.bowtie2.log")"
  unique_pairs="$(samtools view -c -f 64 -F 1024 "${dedup_bam}")"
  mito_rate="$(awk -v m="${mito_pairs}" -v g="${good_pairs}" \
    'BEGIN { printf "%.4f", (g > 0 ? m / g : 0) }')"
  dup_rate="$(awk -v u="${unique_pairs}" -v f="${filtered_pairs}" \
    'BEGIN { printf "%.4f", (f > 0 ? 1 - u / f : 0) }')"
  nfr_fraction="$(awk '{ all += $2; if ($1 < 100) nfr += $2 }
      END { printf "%.4f", (all > 0 ? nfr / all : 0) }' \
    "BAM/${sample_id}.fragment_size.txt")"

  {
    printf 'sample\tinput_pairs\talign_rate\tmito_rate\tfiltered_pairs'
    printf '\tunique_pairs\tdup_rate\tnfr_fraction\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${sample_id}" "${input_pairs}" \
      "${align_rate}" "${mito_rate}" "${filtered_pairs}" "${unique_pairs}" \
      "${dup_rate}" "${nfr_fraction}"
  } > "BAM/${sample_id}.align_stats.txt"
  log "${sample_id}: ${input_pairs} read pairs, ${align_rate}% aligned," \
    "mitochondrial ${mito_rate}, ${unique_pairs} unique pairs" \
    "(duplicate rate ${dup_rate})"

  # Remove intermediate files to save disk space
  rm -f "${sort_bam}" "${sort_bam}.bai" "${trim1}" "${trim2}" \
    "Fastq/${sample_id}_R1.merged.fastq.gz" \
    "Fastq/${sample_id}_R2.merged.fastq.gz"
}

main "$@"
