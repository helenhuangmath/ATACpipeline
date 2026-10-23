#!/bin/bash
#
# Settings of one ATAC-seq project.
# Copy this file, sample_list.txt and comparison_list.txt into your project
# folder, edit them, then run:
#   bash ~/ATACpipeline/run_pipeline.sh config.sh
# Use full paths. No spaces around "=".

# 1. Project
PROJECT_NAME="MyATAC"                          # short name for cluster jobs
WORKDIR="/home/you/MyATAC/Data"                # all results go here
DATADIR="/home/you/MyATAC/Raw"                 # raw fastq files
SAMPLE_LIST="/home/you/MyATAC/sample_list.txt"
COMPARISON_LIST="/home/you/MyATAC/comparison_list.txt"  # empty: no DESeq2

# 2. Reference genome (made by reference/build_reference.sh; a reference
#    made by CUTRUNpipeline works too after running build_reference.sh once)
SPECIES="mm10"
LIBRARIES="/home/you/indexes"                  # one folder per genome
GNM_IDX="${LIBRARIES}/${SPECIES}"
BOWTIE_IDX="${GNM_IDX}/bowtie2/${SPECIES}"     # without .1.bt2
GNM_SIZES="${GNM_IDX}/${SPECIES}.genomesize"
BLACKLIST="${GNM_IDX}/${SPECIES}-blacklist.clean.bed"
CHR1TOY_BED="${GNM_IDX}/${SPECIES}.chr1toY.bed"  # keep reads on chr1-chrY
TSS_BED="${GNM_IDX}/${SPECIES}.tss.bed"        # TSSs for TSS enrichment
TXDB="${GNM_IDX}/genes.txdb.sqlite"            # genes for ChIPseeker
GENE_TABLE="${GNM_IDX}/genes.txt"              # gene ID -> gene name
REFSIZE="1.87e9"                               # mouse 1.87e9, human 2.7e9
HOMER_GENOME=""                                # e.g. mm10 (empty: no Homer)

# 3. Software: conda environments from envs/*.yml (empty: tools on PATH)
CONDA_ENV="atac"
CONDA_ENV_R="atac_r"
MODULES=""                                     # or modules, e.g. "bowtie2"
MODULES_R=""                                   # e.g. "R/4.3"

# 4. Cluster
SCHEDULER="lsf"                                # lsf | slurm | local
QUEUE=""                                       # empty: default queue
THREADS=8                                      # CPUs per job
MEMORY=32000                                   # memory per job (MB)
EXTRA_SUBMIT_ARGS=""                           # e.g. "-W 24:00"

# 5. Trimming (Nextera adapters are found from the overlap of R1 and R2)
RUN_FASTQC=true
RUN_TRIMMING=true                              # fastp adapter trimming
FASTP_ARGS="--detect_adapter_for_pe --length_required 20"

# 6. Alignment (ENCODE ATAC-seq settings: fragments up to 2 kb)
BOWTIE2_ARGS="--very-sensitive --no-unal --no-mixed --no-discordant \
--phred33 -X 2000"
MIN_MAPQ=30                                    # removes multi-mapping reads
REMOVE_DUPLICATES=true                         # false: only mark them
MITO_CHR="chrM"                                # mitochondria (removed)

# 7. Bigwig (CPM: counts per million fragments)
BW_ARGS="--binSize 10"                         # more bamCoverage options

# 8. Peaks (MACS2)
PEAK_MODE="fragment"                           # fragment: whole fragments
                                               #   (-f BAMPE)
                                               # cutsite: Tn5 cut sites,
                                               #   smoothed 150 bp (ENCODE)
MACS2_ARGS="-q 0.01"                           # more options, e.g. --broad
MIN_FOLD=1                                     # remove peaks with fold
                                               #   enrichment below this
MIN_REPS=2                                     # consensus: in >= 2 reps
CONSENSUS_GAP=0                                # replicate peaks closer than
                                               #   this (bp) are one peak
UNION_GAP=0                                    # merge consensus peaks of
                                               #   conditions closer than this
MIN_WIDTH=50                                   # remove narrower union peaks

# 9. Annotation
TSS_UP=2000                                    # promoter: TSS - 2 kb
TSS_DOWN=500                                   #       to: TSS + 500 bp
FLANK=10000                                    # genes within 10 kb

# 10. Differential peaks (DESeq2)
DESIGN="~ condition"                           # or "~ batch + condition"
MIN_COUNT=20                                   # max raw count of a peak
PADJ=0.05                                      # FDR cutoff
LFC=0.585                                      # |log2FC| (0.585: 1.5 fold)

# true: redo all steps, also finished ones
FORCE=false
