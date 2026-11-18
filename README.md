# mTracePatho

**Version:** 1.0.1
**Web platform:** https://mngs.dtsunlab.com/

## Introduction

Clinical metagenomic next-generation sequencing (mNGS) enables culture-independent pathogen detection, but the shallow coverage typical of routine workflows usually yields only tens to thousands of pathogen-specific reads—enough for species identification, yet long considered insufficient for reliable strain-level analysis. Strain resolution, however, is essential for precision infectious-disease care: within-species diversity shapes virulence, antimicrobial resistance, and transmissibility, as illustrated by hypervirulent *Klebsiella pneumoniae*.

**mTracePatho** is a reference-guided framework for rapid, nearest-reference strain placement from ultra-low-coverage clinical mNGS data. It addresses the extreme imbalance between sparse sample SNP sets and comprehensive reference profiles (|S| ≪ |G|) through a two-stage algorithm: (i) asymmetry-aware containment scoring with a size-penalty term to preselect candidate genomes, and (ii) iterative phylogenetic single-nucleotide polymorphism (SNP) assignment via lowest-common-ancestor (LCA) refinement on a BIONJ tree. The companion module `mTracePathoMakeDB.py` builds species-specific SNP–phylogeny databases from quality-filtered assemblies using MinHash clustering for representative selection. Version 1.0 focuses on dominant-strain placement within the reference tree and does not estimate strain abundance or fully deconvolute mixed infections.

mTracePatho is available as a command-line tool and an interactive web platform for upload-based tracing with phylogenetic visualization—supporting precision diagnostics and real-time molecular epidemiological surveillance from routine low-coverage mNGS.

---

## Overview

```
FASTQ
  │
  ├─ Step 1  Reads preprocessing
  │            fastp (QC, optional)
  │            → kraken2 host removal (homo_kraken2)
  │            → kraken2 plasmid filter (PLSDB_kraken2_index)
  │            → kraken2 pathogen classification
  │            → pullfq (extract target taxid reads)
  │            → kz (complexity filter)
  │
  ├─ Step 2  Alignment & variant detection
  │            bwa mem → samtools → BAM filter
  │            → bcftools mpileup|call → VCF filter
  │            → SNP → node mapping (.filt.node)
  │
  ├─ Step 3  Source tracing
  │            mTracePathoAlgorithm (LCA + containment scoring)
  │            → Trace.result.tsv  + figures
  │
  └─ Step 4  Interactive HTML viewer (optional, --taxname)
               Generate_trace_viewer.v2.py → mTracePatho_viewer.html
```

### Core scripts

| File | Role |
|------|------|
| `mTracePathoRun.py` | End-to-end tracing pipeline (CLI entry point) |
| `mTracePathoMakeDB.py` | Build species-specific tracing database from genomes |
| `mTracePathoAlgorithm.py` | LCA tracing algorithm & scoring |
| `Generate_trace_viewer.v2.py` | Interactive phylogenetic HTML viewer |

---

## Requirements

### System

- Linux x86_64
- Python ≥ 3.8
- Sufficient RAM for Kraken2 indexes (typically ≥ 16 GB recommended)
- Disk: Kraken2 indexes ~8 GB + pathogen DB + working space for BAM/VCF

### Python packages

```bash
pip install numpy scipy matplotlib pysam psutil
```

### External tools

| Tool | Purpose |
|------|---------|
| [fastp](https://github.com/OpenGene/fastp) | Read QC / adapter trimming |
| [Kraken2](https://github.com/DerrickWood/kraken2) | Host / plasmid / pathogen classification |
| [BWA](https://github.com/lh3/bwa) | Read alignment |
| [samtools](https://github.com/samtools/samtools) | BAM processing |
| [bcftools](https://github.com/samtools/bcftools) | Variant calling |
| [bamdst](https://github.com/shiquan/bamdst) | Coverage statistics |
| [kz](https://github.com/kloetzl/kz) | Low-complexity read filter |

A `soft/` directory with common binaries may be shipped with the release; otherwise install the tools yourself and point to them in `tools.ini` (see below).

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/GoGoGao/mTracePatho.git
cd mTracePatho

# 2. Create / activate a conda environment (recommended)
conda create -n mtracepatho python=3.10 -y
conda activate mtracepatho

# 3. Install Python dependencies
pip install numpy scipy matplotlib pysam psutil

# 4. Install bioinformatics tools (example with conda)
conda install -c bioconda fastp kraken2 bwa samtools bcftools -y
# Install bamdst and kz separately (see soft/ or upstream repos)

# 5. Make scripts executable (optional)
chmod +x mTracePathoRun.py mTracePathoMakeDB.py
```

### Configure tool & database paths

Generate a config template:

```bash
python mTracePathoRun.py --print_config > tools.ini
# or: python mTracePathoMakeDB.py --print_config
```

Edit `tools.ini` to match your environment:

```ini
# Search order: --config ARG > ./tools.ini > ~/.mtracepatho/tools.ini > built-in defaults

[tools]
fastp     = /path/to/fastp
kraken2   = /path/to/kraken2
bwa       = /path/to/bwa
samtools  = /path/to/samtools
bcftools  = /path/to/bcftools
bamdst    = /path/to/bamdst
kz        = /path/to/kz

[databases]
kraken2_pathogen_db = /path/to/kraken2_pathogen_db
kraken2_plsdb       = /path/to/PLSDB_kraken2_index
kraken2_homo_db     = /path/to/homo_kraken2
```

You can also place the file at `~/.mtracepatho/tools.ini`, or pass `--config /path/to/tools.ini` on each run.

---

## Databases

mTracePatho needs **two kinds** of databases:

1. **Kraken2 classification indexes** — used in Step 1 (host / plasmid / pathogen filtering)
2. **Species tracing database (`--db_dir`)** — built by `mTracePathoMakeDB.py`, used in Steps 2–3

### 1. Kraken2 indexes

Configure these under `[databases]` in `tools.ini`.

#### `homo_kraken2/` — host (human) removal

Used to remove human host reads before plasmid/pathogen classification (`--skip_rmhost` to disable).

```
homo_kraken2/
├── hash.k2d    # ~4.3 GB
├── opts.k2d
└── taxo.k2d
```

| File | Description |
|------|-------------|
| `hash.k2d` | Kraken2 k-mer hash table |
| `opts.k2d` | Index build options |
| `taxo.k2d` | Taxonomy tree |

Typical size: **~4.3 GB** total.

#### `PLSDB_kraken2_index/` — plasmid filter

Used to remove plasmid-derived reads (PLSDB) so chromosomal pathogen SNPs are not confounded.

```
PLSDB_kraken2_index/
├── hash.k2d    # ~3.3 GB
├── opts.k2d
└── taxo.k2d
```

Typical size: **~3.3 GB** total.

```
Files shared via cloud storage: mTracePatho_Database Baidu Netdisk link
Link: https://pan.baidu.com/s/1-znpClAIiAJxRZXrVTfKNA Extraction code: a2n4
```

### 2. Species tracing database (`--db_dir`)

Build with `mTracePathoMakeDB.py` from reference genomes. Database files are written to `<outdir>/02_database/` by default.

```bash
python mTracePathoMakeDB.py \
  --genome_list genomes.tsv \
  --reference ref.fa \
  --outdir Klebsiella_build/ \
  --taxid 1125630 \
  --threads 16
```

`genomes.tsv` format (tab-separated):

```text
genome1	/path/to/genome1.fa
genome2	/path/to/genome2.fa
```

Required files in `--db_dir` (e.g. `02_database/`):

| File | Description |
|------|-------------|
| `total.var.tsv` | SNP presence/absence matrix |
| `Jaccard.fastme.input.nwk` | Phylogenetic tree (Newick) |
| `node.tsv` | Tree node annotations |
| `numbered_tree.tree` | Numbered tree for LCA |
| `snp_lingage.tsv` | SNP–lineage mapping |
| `lingage.stat.tsv` | Lineage statistics |
| `LCA.output.tsv` | Precomputed LCA table |

---

## Usage examples

### Recommended: batch tracing (*Klebsiella pneumoniae*)

Prepare `in.list` (tab-separated: sample name, FASTQ path):

```text
sample1	/data/fastq/sample1.fq.gz
sample2	/data/fastq/sample2.fq.gz
```

```bash
python mTracePathoRun.py \
  --fastq_list in.list \
  --db_dir 02_database \
  --outdir ./result-containment_weighted \
  --taxid "1125630" \
  --threads 32 \
  --min_snp 5 \
  --top_n_cutoff 10 \
  --taxname "Klebsiella pneumoniae" \
  --skip_pathogen
  # Optional (uncomment as needed):
  # --skip_fastp --skip_rmhost
  # --clean_temp
```

Notes for this command:

| Flag | Meaning |
|------|---------|
| `--db_dir 02_database` | Tracing DB from `mTracePathoMakeDB.py` output |
| `--taxid "1125630"` | Target NCBI taxid |
| `--min_snp 5` | Minimum SNPs required to run tracing |
| `--top_n_cutoff 10` | Keep top 10 distinct containment score levels |
| `--taxname "Klebsiella pneumoniae"` | Enable interactive HTML viewer |
| `--skip_pathogen` | Skip pathogen Kraken2 + pullfq (use when reads are already species-filtered); still runs PLSDB / host steps unless also skipped |
| `--skip_fastp` / `--skip_rmhost` | Optional: skip QC or host removal |
| `--clean_temp` | Optional: delete intermediates after tracing |

### Single sample

```bash
python mTracePathoRun.py \
  --fastq sample.fq.gz \
  --db_dir 02_database \
  --taxid 1125630 \
  --outdir results/ \
  --threads 8 \
  --taxname "Klebsiella pneumoniae"
```

### Multiple taxids

```bash
python mTracePathoRun.py \
  --fastq sample.fq.gz \
  --db_dir 02_database \
  --taxid "470-1234-5678" \
  --outdir results/
```

### Skip entire Step 1 preprocessing

```bash
python mTracePathoRun.py \
  --fastq sample.fq.gz \
  --db_dir 02_database \
  --taxid 1125630 \
  --outdir results/ \
  --skip_step1
```

### Re-trace from existing results

If Steps 1–2 already finished and you only want to re-run tracing:

```bash
python mTracePathoRun.py \
  --result_dir ./result-containment_weighted \
  --fastq_list in.list \
  --db_dir 02_database \
  --outdir results_retrace/
```

### SGE / qsub cluster mode

```bash
python mTracePathoRun.py \
  --fastq_list in.list \
  --db_dir 02_database \
  --taxid 1125630 \
  --outdir results/ \
  --use_sge \
  --sge_nodes "R7525|big2|big3|big4" \
  --sge_mem 1G \
  --sge_slots 1 \
  --sge_queue all.q
```

---

## Main arguments (`mTracePathoRun.py`)

### Required

| Argument | Description |
|----------|-------------|
| `--fastq` **or** `--fastq_list` | Single FASTQ, or TSV (`sample_name`, `fastq_path`) |
| `--db_dir` | Species tracing database directory (e.g. `02_database`) |
| `--taxid` | Target NCBI taxonomy ID(s); join multiple with `-` |
| `--outdir` | Output directory |

### Common options

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | 4 | Local workers / tool threads |
| `--config` | auto | Path to `tools.ini` |
| `--min_snp` | 5 | Minimum SNPs required for tracing |
| `--top_n_cutoff` | 99999 | Top N distinct containment score levels |
| `--filter_identity` | 0.95 | Min alignment identity for BAM filter |
| `--filter_coverage` | 0.95 | Min coverage fraction for BAM filter |
| `--skip_step1` | off | Skip all Step 1 preprocessing |
| `--skip_fastp` | off | Skip fastp QC |
| `--skip_rmhost` | off | Skip human host removal |
| `--skip_pathogen` | off | Skip pathogen Kraken2 + pullfq |
| `--force` | off | Ignore `.done` checkpoints and re-run |
| `--jaccard_percent` | 0.8 | Jaccard threshold for tracing |
| `--use_containment` | on | Containment-based genome pre-selection |
| `--taxname` | — | Latin name; enables HTML viewer |
| `--no_html_viewer` | off | Disable HTML viewer |
| `--clean_temp` | off | Delete intermediate files after tracing |
| `--print_config` | — | Print `tools.ini` template and exit |

Run `python mTracePathoRun.py -h` for the full list (scoring modes, SGE options, etc.).

### Build database (`mTracePathoMakeDB.py`)

| Argument | Description |
|----------|-------------|
| `--genome_list` | TSV: sample name + genome FASTA path |
| `--reference` | Reference FASTA |
| `--outdir` | Build output directory (`02_database/` is created inside) |
| `--taxid` | Pathogen NCBI taxid(s) |
| `--threads` | Parallel workers (default: 4) |
| `--print_config` | Print tools.ini template and exit |

---

## Output structure

```text
{outdir}/
├── logs/
│   ├── pipeline.log          # Run log
│   └── time.RSS.TSV          # Per-step timing / memory
├── {sample}/
│   ├── {sample}.homo.kraken
│   ├── {sample}.PLSDB.kraken
│   ├── {sample}.kraken
│   ├── {sample}.kz.fq.gz
│   ├── {sample}.sorted.bam
│   ├── {sample}.vcf
│   ├── {sample}.filt.node    # SNP→node mapping for tracing
│   └── ...
└── results/
    ├── {sample}/
    │   ├── Trace.result.tsv
    │   ├── iteration_statistics.tsv
    │   ├── tracing_report.txt
    │   ├── figures/                    # PNG visualizations
    │   ├── mTracePatho_viewer.html     # if --taxname
    │   └── mTracePatho_viewer.data.json
    └── batch_summary.tsv
```

Primary result file: **`results/{sample}/Trace.result.tsv`** (source genomes, paths, confidence scores).

---

## Checkpoint & resume

Each sample writes `.step1.done` / `.step2.done` markers. Re-running the same `--outdir` resumes from unfinished steps. Use `--force` to ignore checkpoints and recompute everything.

---

## Citation

If you use mTracePatho in your research, please cite the forthcoming paper and this repository:

> mTracePatho enables rapid reference-guided strain-level pathogen tracing from low-coverage clinical metagenomic sequencing data.  
> Source code: https://github.com/GoGoGao/mTracePatho  
> Web platform: https://mngs.dtsunlab.com/

---

## License

[Choose a license, e.g. MIT / GPL-3.0]

---

## Contact

- Issues: https://github.com/GoGoGao/mTracePatho/issues
- Authors: gjpbioinfo@163.com; linwenchao@yeah.net



👋 Hi，想聊点什么～
输入您的问题... (Enter 发送，Shift+Enter 换行)

内容由AI生成，仅供参考
