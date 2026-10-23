# ATACpipeline

An ATAC-seq analysis pipeline, from paired-end FASTQ files to peaks, ATAC quality control (mitochondrial
reads, fragment sizes, TSS enrichment), gene annotation, count matrices, differential accessibility,
publication-style figures and a searchable HTML report.
Written in plain bash and R; runs on an HPC cluster (LSF or SLURM) or on one computer.

Fill in **a sample list** and **a config file**, and run **one command**.

```
FASTQ ─▶ FastQC + MultiQC ─▶ fastp ─▶ bowtie2 ─▶ filter, remove chrM + duplicates ─▶ BAM
                                                                                     │
                                  bigwig + TSS enrichment ─▶ peaks of each sample (MACS2)
                                                                                     │
                               consensus peaks of each condition ─▶ union peaks of each group
                                                                                     │
                 annotation ─▶ raw / normalized / VST count matrices of all samples ─▶ QC summary
                                                                                     │
              DESeq2 differential peaks ─▶ Nature-style figures ─▶ HTML report (search a gene)
```

Each step, what it does and what to check: [PIPELINE_STEPS.txt](PIPELINE_STEPS.txt).

---

## 1. Install (once)

Use conda (Miniconda or Miniforge). Create the two software environments:

```bash
git clone https://github.com/helenhuangmath/ATACpipeline.git ~/ATACpipeline
cd ~/ATACpipeline
conda env create -f envs/atac.yml      # tools: fastp, bowtie2, samtools, MACS2, deepTools, ...
conda env create -f envs/atac_r.yml    # R: DESeq2, ChIPseeker, ggplot2, patchwork, ...
```

Software versions are fixed in `envs/*.yml`, so every user gets the same versions. If your cluster uses
modules instead, leave `CONDA_ENV` empty in the config and list the modules in `MODULES`.

## 2. Reference genome (once per genome)

```bash
bsub -n 8 -M 32000 -o ref_mm10.log bash ~/ATACpipeline/reference/build_reference.sh mm10 ~/indexes 8
```

Supported: `mm10`, `hg38` and `hg19`. It downloads the genome (UCSC), the blacklist (ENCODE v2) and genes
(GENCODE), and makes the bowtie2 index (1-3 hours), `mm10.genomesize`, `mm10.chr1toY.bed`,
`mm10-blacklist.clean.bed`, the gene annotation and `mm10.tss.bed` (TSSs for the TSS enrichment) in
`~/indexes/mm10/`.

## 3. Set up a project

```bash
mkdir -p ~/MyProject && cd ~/MyProject
cp ~/ATACpipeline/config/{config.sh,sample_list.txt,comparison_list.txt} .
```

**sample_list.txt**: one row per sample, tab-separated (in Excel: save as *Tab delimited text*).

```
sample          fastq1                       fastq2                       group  condition  replicate  batch
CD8_Ctrl_R1     CD8_Ctrl_R1_R1.fastq.gz      CD8_Ctrl_R1_R2.fastq.gz      CD8    Ctrl       1          b1
CD8_Ctrl_R2     CD8_Ctrl_R2_R1.fastq.gz      CD8_Ctrl_R2_R2.fastq.gz      CD8    Ctrl       2          b2
CD8_Arid2KO_R1  CD8_Arid2KO_R1_R1.fastq.gz   CD8_Arid2KO_R1_R2.fastq.gz   CD8    Arid2KO    1          b1
```

- `fastq1/2` are inside `DATADIR`; several lanes: `L1_R1.fastq.gz,L2_R1.fastq.gz`
- `group`: samples analysed together, e.g. a cell type. Each group gets its own consensus peaks, count
  matrices and DESeq2 model. If all samples belong together, use one value (e.g. `ATAC`)
- `condition`: replicates of one condition are combined into consensus peaks and compared by DESeq2

**comparison_list.txt**: comparisons for DESeq2 (`*` = every group; log2FC > 0 = more open in treatment).

```
group  treatment  reference
*      Arid2KO    Ctrl
*      Pbrm1KO    Ctrl
```

**config.sh**: fill in sections 1-4 (project folders, reference, software, cluster). The other settings
have defaults; all are explained in the file. Peak calling:

| setting | meaning |
|---|---|
| `PEAK_MODE="fragment"` | MACS2 on whole fragments (`-f BAMPE`); `cutsite`: on Tn5 cut sites, as ENCODE |
| `MACS2_ARGS="-q 0.01"` | more MACS2 options, e.g. `-p 1e-5`, `--call-summits` |
| `MIN_REPS=2` | a consensus peak must be in at least 2 replicates of a condition |

## 4. Run

```bash
bash ~/ATACpipeline/run_pipeline.sh config.sh
```

First the inputs are checked (nothing is submitted if something is wrong). Then the jobs are submitted;
each starts when the jobs it needs are finished:

| part | jobs | what |
|---|---|---|
| `qc` | 1 | FastQC and MultiQC of all raw fastq files (runs next to the alignment) |
| `align` | 1 per sample | fastp trimming, FastQC of trimmed reads, bowtie2, filtering, duplicates, bigwig, TSS enrichment |
| `peaks` | 1 per sample | MACS2 |
| `matrix` | 1 | consensus and union peaks, annotation, count matrices, QC summary |
| `diff` | 1 | DESeq2, summary figures, HTML report |

Run one part only with `run_pipeline.sh config.sh qc|align|peaks|matrix|diff`, e.g. `qc` first to check
the data, or `diff` again after changing the comparisons or `DESIGN`.

## 5. Check progress and rerun

```bash
bash ~/ATACpipeline/check_status.sh config.sh
```

- **A job failed**: fix the problem and run `run_pipeline.sh config.sh` again; finished steps are skipped.
  On LSF, first `bkill` the jobs waiting for the failed one.
- **Changed peak settings**: delete `Data/Status/*.peaks.done` and run again.
- **Redo everything**: `FORCE=true` in config.sh.

## 6. Results (in WORKDIR)

| folder | files |
|---|---|
| `Report/` | **`ATAC_report.html`**: search a gene, peak or locus; see its peaks, differential statistics, signal profile and normalized counts; links to UCSC and IGV. One file, opens offline |
| `Figures/<group>/` | **`Figure_<group>_summary.pdf`**: multi-panel figure (183 mm, 7 pt): peaks, FRiP, TSS enrichment, genomic features, PCA, differential peaks, volcano plots, signal profiles, signal heatmap |
| `ReadCount/<group>/` | **count matrices of all samples** (union peaks x samples): `RawRC_<group>.txt` (raw), `NormRC_<group>.txt` (DESeq2-normalized), `VST_<group>.txt` (log2 scale), `NormRC_<group>.CondMean.txt`, `sizeFactors_<group>.txt` (size factors, FRiP), PCA, correlation |
| `DiffPeak/<group>/` | per comparison: `deseq2_All_...`, `deseq2_Sig_...`, `_Up.bed`, `_Down.bed`, MA, volcano and heatmap plots; `deseq2_summary_<group>.txt` |
| `QC/` | `MultiQC/multiqc_raw_fastq.html`, `MultiQC/multiqc_report.html` (all steps), `FastQC_raw/`, `FastQC_trimmed/`, fastp reports, `TSS/` (TSS profiles and scores), `Summary/sample_QC_summary.txt` + plots |
| `BAM/` | `<sample>.map.dedup.bam` (filtered, chrM and duplicates removed), alignment stats, reads per chromosome, fragment sizes |
| `Wig/` | `<sample>.cpm.bw` tracks for IGV / UCSC |
| `Peak/` | `MACS2/` (raw output), `peak_<sample>.final.bed` (cleaned peaks used for the analysis) |
| `Consensus/<group>/` | `peak_consensus_<condition>.bed`, `peak_<group>_union.bed`, Venn diagrams, peak width |
| `Annotation/<group>/` | `peak_<group>_union.anno.txt` (nearest gene, promoter/intron/...), `.anno_allgenes.txt` |
| `Log/` | job logs, settings of every run (`run_<date>/`), `software_versions.txt` |

**ATAC-seq quality** (`QC/Summary/sample_QC_summary.txt`; ENCODE guidelines in brackets):

| column | meaning | good |
|---|---|---|
| `align_rate` | % of read pairs aligned | > 80% (> 95% ideal) |
| `mito_rate` | fraction of good read pairs on chrM | low; < 0.2 is typical for Omni-ATAC |
| `unique_pairs` | fragments after filtering and duplicate removal | ≥ 25 M per replicate |
| `nfr_fraction` | fraction of fragments < 100 bp (nucleosome-free) | fragment plot shows the nucleosome-free peak and a peak every ~200 bp |
| `tss_enrichment` | Tn5 cut sites at TSSs / 1.9-2 kb away | hg38 > 7, mm10 > 10 (computed slightly differently) |
| `FRiP` | fraction of fragments in union peaks | > 0.3 (> 0.2 acceptable) |

## 7. Tests

```bash
TEST_CONDA_ENV=atac TEST_CONDA_ENV_R=atac_r bash tests/run_tests.sh
```

Makes a small simulated data set with known open sites, TSSs, mitochondrial reads and known differences
between conditions ([tests/make_test_data.sh](tests/make_test_data.sh)), runs the whole pipeline on it
(a few minutes), and checks 99 results: input checks catch common mistakes, chrM reads are counted and
removed, lanes are merged, the fragment sizes and TSS enrichment are right, MACS2 finds the simulated
sites in both peak modes, the expected differential peaks are found, the figures and the HTML report are
made, and a rerun skips finished steps.

Code follows the [Google Shell Style Guide](https://google.github.io/styleguide/shellguide.html) and the
[Google R Style Guide](https://google.github.io/styleguide/Rguide.html), checked with ShellCheck
([.shellcheckrc](.shellcheckrc)) and lintr ([.lintr](.lintr)).

## 8. Troubleshooting

| message | what to do |
|---|---|
| `The sample list header must be ...` | save the sample list as tab-separated text with this header |
| `fastq1 not found` | check `DATADIR` and the file names |
| `conda: command not found` | load conda (`module load miniconda`), or use `MODULES` |
| `TSS_BED not found` | run `reference/build_reference.sh` once more on the reference folder |
| high `mito_rate` | normal for standard ATAC (up to ~50%); chrM reads are removed, but fewer fragments are left |
| very few peaks | check `QC/Summary` (unique fragments, TSS enrichment, FRiP); loosen the cutoff in `MACS2_ARGS`, e.g. `-p 1e-3` |
| `not enough replicates` | each condition of a comparison needs at least 2 replicates |

## 9. Methods text

> Raw reads were checked with FastQC and MultiQC, trimmed with fastp, and aligned to the genome with
> bowtie2 (`--very-sensitive --no-mixed --no-discordant -X 2000`). Properly paired reads with mapping
> quality ≥ 30 on chr1-chrY were kept (mitochondrial reads removed), and PCR duplicates were removed with
> samtools markdup. CPM-normalized coverage tracks were made with deepTools bamCoverage. TSS enrichment
> was computed from Tn5 insertion sites (reads shifted +4/−5 bp) around GENCODE protein-coding TSSs.
> Peaks were called for each sample with MACS2 (BAMPE mode, q < 0.01), and peaks in ENCODE blacklist
> regions were removed. Consensus peaks were peaks found in at least two replicates of a condition;
> consensus peaks of all conditions were merged into union peaks. Fragments in union peaks were counted
> with featureCounts, normalized with DESeq2 size factors, and differentially accessible peaks were found
> with DESeq2 (FDR < 0.05, fold change > 1.5). Peaks were annotated with ChIPseeker (promoter: −2 kb to
> +500 bp of the TSS) using GENCODE genes.

Exact versions are in `Log/software_versions.txt` and `DiffPeak/<group>/sessionInfo_<group>.txt`.

Downstream (not automatic): motifs (Homer `findMotifsGenome.pl`) and pathways (GREAT, clusterProfiler) of
the `_Up.bed` / `_Down.bed` peaks.

## 10. Contributing, citation, license

- **Contributing**: questions, bug reports and changes are welcome. Open an issue (bug report or feature
  request template) or a pull request; changes are listed in [CHANGELOG.md](CHANGELOG.md).
- **Citation**: please cite this pipeline ([CITATION.cff](CITATION.cff)) and the tools it uses: FastQC,
  MultiQC, fastp, Bowtie2, SAMtools, deepTools, MACS2, BEDTools, Subread/featureCounts, DESeq2,
  ChIPseeker.
- **License**: [MIT](LICENSE). You may use, change and share the pipeline, also commercially, if you keep
  the copyright and license notice.
