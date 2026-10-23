# Changelog

All notable changes are listed here. Version numbers follow [Semantic Versioning](https://semver.org):
MAJOR (results or inputs change in a way that is not backward compatible), MINOR (new features),
PATCH (bug fixes).

## [0.1.0] - 2026-10-23

First version, rebuilt from the lab's ATAC-seq scripts (2020-2026) into one pipeline, with the same
design as CUTRUNpipeline.

### Added
- One command for the whole analysis (`run_pipeline.sh`), driven by `config.sh`, `sample_list.txt`
  and `comparison_list.txt`; jobs for LSF, SLURM or a single computer, with job dependencies.
- Input checks before any job is submitted (`scripts/validate.sh`).
- Trimming with fastp; FastQC and MultiQC reports of raw and trimmed reads.
- Alignment with ENCODE ATAC-seq settings (bowtie2 `--very-sensitive -X 2000`, MAPQ ≥ 30);
  mitochondrial reads counted and removed; duplicates removed with samtools markdup.
- ATAC-seq QC: mitochondrial read fraction, nucleosome-free fraction, fragment size distribution,
  TSS enrichment from Tn5 cut sites (+4/−5 bp shift) and TSS profiles.
- MACS2 peak calling on whole fragments (BAMPE) or Tn5 cut sites (`PEAK_MODE`).
- Groups of samples (e.g. cell types) analysed separately in one run (`group` column).
- Rerun after a failure: finished steps are skipped (`Status/*.done`); `check_status.sh`.
- Consensus peaks for any number of replicates; union peaks of all conditions.
- Count matrix step (`07_count_matrix.sh`): raw, DESeq2-normalized and VST matrices of all samples;
  the DESeq2 step reuses its size factors.
- Generic DESeq2 analysis for any conditions and comparisons, with batch correction.
- QC summary of all samples (`QC/Summary/`).
- Nature-style plots (7 pt, journal column widths, colour-blind-safe colours) and a multi-panel
  summary figure of each group with signal profiles and heatmaps (`09_figures.sh`).
- Searchable HTML report of all groups (`10_report.sh`): find a gene and see its peaks, statistics,
  signal profile and counts.
- Reproducibility records: settings of every run, software versions, R sessionInfo.
- Reference builder for mm10, hg38 and hg19 (`reference/build_reference.sh`); reuses a reference made
  by CUTRUNpipeline and only adds the TSS file.
- Conda environments (`envs/`), automated tests on simulated data (`tests/`), GitHub Actions CI.
- Code style: Google Shell and R style guides, checked with ShellCheck and lintr in CI.
- MIT license; contributor guide, contributors list, issue and pull request templates.
