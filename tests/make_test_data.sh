#!/bin/bash
#
# Makes a small simulated ATAC-seq data set to test the whole pipeline.
#
# Usage: make_test_data.sh <output_dir>
#
# Genome : chr1, chr2 (500 kb each, random sequence) and chrM
# Genes  : 20 genes; the TSS of each is at an open site (TSS enrichment)
# Sites  : 100 open sites in Ctrl, KO1 and KO2 (group Tcell, 2 replicates)
#            every 7th site  4x lower in Ctrl -> Up vs Ctrl
#            every 10th site 4x higher in KO2 -> Up in KO2
# Reads  : paired-end 50 bp; fragments of 60-100 bp (nucleosome-free),
#          ~210 bp (1 nucleosome) and ~370 bp (2 nucleosomes) centered on the
#          sites; background reads; reads on chrM (to be removed)
# Lanes  : Ctrl_R1 is split into 2 lanes (tests lane merging)

set -eo pipefail

#######################################
# Writes the genome, blacklist and genes. build_reference.sh uses these files
# instead of downloading them.
# Arguments:
#   Genome folder.
#######################################
make_reference() {
  local gnm_idx="$1"
  mkdir -p "${gnm_idx}"

  awk 'BEGIN {
      srand(1); split("A C G T", base, " ")
      n = split("chr1 500000 chr2 500000 chrM 16000", genome, " ")
      for (c = 1; c < n; c += 2) {
        print ">" genome[c]
        for (i = 0; i < genome[c + 1] / 1000; i++) {
          line = ""
          for (j = 0; j < 1000; j++) line = line base[int(rand() * 4) + 1]
          print line
        }
      }
    }' | gzip > "${gnm_idx}/mm10.fa.gz"

  # One blacklist region, between two sites
  printf 'chr1\t207500\t208500\n' | gzip > "${gnm_idx}/mm10-blacklist.v2.bed.gz"

  # A gene every 45 kb, starting at an open site (every 5th site of the
  # chromosome), with 2 exons
  awk 'BEGIN {
      OFS = "\t"
      for (c = 1; c <= 2; c++) for (s = 5000; s < 450000; s += 45000) {
        g++; e = s + 20000
        id = "gene_id \"ENSMUSG" g ".1\"; gene_name \"Gene" g "\"; " \
          "gene_type \"protein_coding\";"
        tx = id " transcript_id \"ENSMUST" g ".1\";"
        print "chr" c, "TEST", "gene", s, e, ".", "+", ".", id
        print "chr" c, "TEST", "transcript", s, e, ".", "+", ".", tx
        print "chr" c, "TEST", "exon", s, s + 2000, ".", "+", ".", tx
        print "chr" c, "TEST", "exon", e - 3000, e, ".", "+", ".", tx
      }
    }' | gzip > "${gnm_idx}/genes.gtf.gz"
}

#######################################
# Writes the sample list and the comparison list.
# Arguments:
#   Output folder.
#######################################
make_lists() {
  local ctrl1="Ctrl_R1_L1_R1.fastq.gz,Ctrl_R1_L2_R1.fastq.gz"
  local ctrl2="Ctrl_R1_L1_R2.fastq.gz,Ctrl_R1_L2_R2.fastq.gz"
  {
    printf 'sample\tfastq1\tfastq2\tgroup\tcondition\treplicate\tbatch\n'
    printf 'Ctrl_R1\t%s\t%s\tTcell\tCtrl\t1\tb1\n' "${ctrl1}" "${ctrl2}"
    local row sample cond rep batch
    for row in Ctrl_R2:Ctrl:2:b2 KO1_R1:KO1:1:b1 KO1_R2:KO1:2:b2 \
      KO2_R1:KO2:1:b1 KO2_R2:KO2:2:b2; do
      IFS=: read -r sample cond rep batch <<< "${row}"
      printf '%s\t%s_R1.fastq.gz\t%s_R2.fastq.gz\tTcell\t%s\t%s\t%s\n' \
        "${sample}" "${sample}" "${sample}" "${cond}" "${rep}" "${batch}"
    done
  } > "$1/sample_list.txt"

  printf 'group\ttreatment\treference\n*\tKO1\tCtrl\n*\tKO2\tCtrl\n' \
    > "$1/comparison_list.txt"
  printf '*\tKO2\tKO1\n' >> "$1/comparison_list.txt"
}

#######################################
# Writes the paired-end reads of one sample.
# Arguments:
#   Genome fasta.gz, sample ID, random seed, output prefix.
#######################################
make_reads() {
  gzip -dc "$1" \
    | awk -v sample="$2" -v seed="$3" -v r1="$4_R1.fastq" -v r2="$4_R2.fastq" '
      function revcomp(s,   r, i) {
        r = ""
        for (i = length(s); i > 0; i--) r = r comp[substr(s, i, 1)]
        return r
      }
      function fragment(chr, start, len,   seq) {
        seq = substr(genome[chr], start, len); n++
        print "@" sample "_" n "/1\n" substr(seq, 1, 50) "\n+\n" qual > r1
        print "@" sample "_" n "/2\n" revcomp(substr(seq, len - 49, 50)) \
          "\n+\n" qual > r2
      }
      function fragment_length(   r) {
        r = rand()
        if (r < 0.6) return 60 + int(rand() * 40)
        if (r < 0.9) return 180 + int(rand() * 60)
        return 340 + int(rand() * 60)
      }
      /^>/ { chr = substr($1, 2); next }
      { genome[chr] = genome[chr] $0 }
      END {
        srand(seed)
        comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"
        qual = sprintf("%50s", ""); gsub(/ /, "I", qual)

        for (k = 1; k <= 100; k++) {
          depth = 40
          if (sample ~ /Ctrl/ && k % 7 == 0) depth = 10
          if (sample ~ /KO2/ && k % 10 == 0) depth = 160
          chr = (k % 2 ? "chr1" : "chr2"); site = 5000 + int((k - 1) / 2) * 9000
          for (j = 0; j < depth; j++) {
            len = fragment_length()
            fragment(chr, site - int(len / 2) + int(rand() * 100) - 50, len)
          }
        }
        for (j = 0; j < 3000; j++) {
          start = 1 + int(rand() * 498000)
          fragment("chr" (1 + j % 2), start, fragment_length())
        }
        for (j = 0; j < 800; j++) {
          fragment("chrM", 1 + int(rand() * 15500), fragment_length())
        }
      }'
  gzip -f "$4_R1.fastq" "$4_R2.fastq"
}

#######################################
# Splits the reads of one sample into 2 lanes.
# Arguments:
#   Raw folder, sample ID.
#######################################
split_lanes() {
  local read
  for read in R1 R2; do
    gzip -dc "$1/$2_${read}.fastq.gz" > "$1/lanes.tmp"
    head -n 8000 "$1/lanes.tmp" | gzip > "$1/$2_L1_${read}.fastq.gz"
    tail -n +8001 "$1/lanes.tmp" | gzip > "$1/$2_L2_${read}.fastq.gz"
    rm "$1/lanes.tmp" "$1/$2_${read}.fastq.gz"
  done
}

main() {
  local out_dir="$1"
  if [[ -z "${out_dir}" ]]; then
    echo "Usage: make_test_data.sh <output_dir>" >&2
    exit 1
  fi
  local gnm_idx="${out_dir}/indexes/mm10"
  mkdir -p "${out_dir}/Raw"

  make_reference "${gnm_idx}"
  make_lists "${out_dir}"

  local seed=0 sample
  for sample in Ctrl_R1 Ctrl_R2 KO1_R1 KO1_R2 KO2_R1 KO2_R2; do
    seed=$(( seed + 1 ))
    make_reads "${gnm_idx}/mm10.fa.gz" "${sample}" "${seed}" \
      "${out_dir}/Raw/${sample}"
  done
  split_lanes "${out_dir}/Raw" Ctrl_R1
  echo "Test data ready in ${out_dir}"
}

main "$@"
