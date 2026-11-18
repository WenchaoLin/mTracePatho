#!/usr/bin/env python3
"""
mTracePathoRun.py — Run mTracePatho pathogen source tracing from FASTQ inputs.

Given one or more clinical/environmental FASTQ files, this program executes a
full reads-to-tracing pipeline:

  Step 1 — Reads preprocessing
            kraken2 (PLSDB plasmid filter) → kraken2 (pathogen classification)
            → pullfq (extract target taxid reads) → kz (complexity filter)

  Step 2 — Alignment & variant detection
            bwa mem → samtools sort/index → filter_bam (identity/coverage)
            → bcftools mpileup|call → filter_vcf
            → samtools depth → SNP→node mapping → .filt.node

  Step 3 — Pathogen source tracing
            Collect per-sample SNP sets → mTracePatho tracing algorithm
            → Trace.result.tsv + confidence scores + visualizations

  Step 4a — Interactive HTML Viewer (optional, requires --taxname)
             Generate per-sample interactive phylogenetic tree viewer
             using Generate_trace_viewer_v2.py → mTracePatho_viewer.html

Version:   1.0.0
Authors:   [Author Name]
Citation:  [To be filled after publication]
GitHub:    https://github.com/[To be filled]

Requires mTracePathoMakeDB.py output database (--db_dir):
    total.var.tsv, Jaccard.fastme.input.nwk, node.tsv,
    numbered_tree.tree, snp_lingage.tsv, lingage.stat.tsv, LCA.output.tsv

Tracing algorithm provided by mTracePathoAlgorithm.py (must be in same directory
or on PYTHONPATH).

Usage:
    # Single FASTQ sample:
    mTracePathoRun.py --fastq sample.fq.gz \\
                      --db_dir dbdir/ --taxid 470 \\
                      --outdir results/

    # Batch mode (TSV list):
    mTracePathoRun.py --fastq_list samples.tsv \\
                      --db_dir dbdir/ --taxid 470 \\
                      --outdir results/ --threads 16

    # SGE cluster mode:
    mTracePathoRun.py --fastq_list samples.tsv \\
                      --db_dir dbdir/ --taxid 470 \\
                      --outdir results/ \\
                      --use_sge --sge_nodes "R7525|big2" --sge_mem 1G --sge_slots 1

    # Skip preprocessing (FASTQ already pre-filtered):
    mTracePathoRun.py --fastq sample.fq.gz \\
                      --db_dir dbdir/ --taxid 470 \\
                      --outdir results/ --skip_steps step1

    # Skip Step 1 directly (shortcut flag):
    mTracePathoRun.py --fastq sample.fq.gz \\
                      --db_dir dbdir/ --taxid 470 \\
                      --outdir results/ --skip_step1

    # Print config template:
    mTracePathoRun.py --print_config

Required arguments:
    --fastq             Single FASTQ input (.fq / .fq.gz)
    --fastq_list        TSV: col1=sample_name, col2=fastq_path (batch mode)
    --db_dir            Database directory (mTracePathoMakeDB.py output)
    --taxid             Target taxid(s), multiple joined by '-' (e.g. "470-1234")
    --outdir            Output directory

Optional arguments:
    --threads           Local parallel workers and tool threads (default: 4)
    --filter_identity   BAM filter: minimum identity fraction (default: 0.95)
    --filter_coverage   BAM filter: minimum coverage fraction (default: 0.95)
    --skip_steps        Steps to skip: step1,step2 (comma-separated)
    --skip_step1        Skip entire Step 1 preprocessing and start from Step 2
                        (equivalent to --skip_steps step1)
    --force             Ignore .done checkpoints, re-run all steps
    --config            Path to tools.ini config file
    --print_config      Print config template and exit

Tracing algorithm parameters (passed to mTracePathoAlgorithm core):
    --jaccard_percent   Jaccard similarity threshold (default: 0.8)
    --max_iterations    Maximum LCA iteration rounds (default: 10)
    --use_containment   Enable containment-based genome pre-selection
    --min_containment   Minimum containment score threshold (default: 0.0)
    --top_n_cutoff      Top N distinct containment score levels (default: 99999 = all)
    --min_candidates    Minimum candidate genomes to retain (default: 3)
    --no_visualization  Skip visualization output (faster)
    --keep_temp         [DEPRECATED] Ignored; intermediate files kept by default
    --clean_temp        Delete ALL intermediate files after tracing (saves disk)

HTML Viewer options:
    --taxname           Species Latin name for the HTML viewer
                        (e.g. "Klebsiella pneumoniae"); required to enable viewer
    --meta              Genome metadata TSV (e.g. Klebsiella_pneumoniae-gene.tsv);
                        auto-detected from --db_dir if not provided
    --no_html_viewer    Explicitly skip HTML viewer generation

SGE options (Step 1+2 only; disabled by default):
    --use_sge           Enable SGE cluster submission via qsub
    --sge_nodes         Pipe-separated node list, e.g. "R7525|big2|big3|big4"
    --sge_mem           Memory per job, vf= value (default: 1G)
    --sge_slots         CPU slots per job, num_proc= (default: 1)
    --sge_threads       Tool threads inside each job (default: = sge_slots)
    --sge_queue         SGE queue name (default: all.q)
    --sge_max_jobs      Max concurrent SGE jobs (default: 500)
    --sge_poll_sec      qstat polling interval seconds (default: 30)
    --sge_extra         Extra qsub options appended verbatim

Output structure:
    {outdir}/
      logs/
        pipeline.log      Structured run log
        time.RSS.TSV      Per-step timing and memory records
        sge_scripts/      Per-sample qsub shell scripts (SGE mode)
        sge_logs/         qsub stdout/stderr (SGE mode)
      {sample}/
        .step1.done / .step1.done.time
        .step2.done / .step2.done.time
        {sample}.PLSDB.kraken
        {sample}.kraken
        {sample}.filter.fq
        {sample}.kz.fq.gz
        {sample}.sorted.bam / .bai
        {sample}.sorted.filter.bam
        {sample}.raw.vcf  →  {sample}.vcf  →  {sample}.filt.node
        coverage_report/  (bamdst output)
        filtered.bam.log
        {sample}.depth
      results/
        {sample}/
          Trace.result.tsv
          iteration_statistics.tsv
          tracing_report.txt
          figures/        (PNG visualizations)
          mTracePatho_viewer.html   (interactive tree viewer, if --taxname)
          mTracePatho_viewer.data.json
        batch_summary.tsv
"""

# ── Standard library ───────────────────────────────────────────────────────────
import os
import sys
import re
import csv
import gzip
import time
import math
import shutil
import logging
import hashlib
import argparse
import tempfile
import traceback
import subprocess
import configparser
from pathlib import Path
from multiprocessing import Pool
from typing import List, Dict, Set, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

# ── Third-party (required) ─────────────────────────────────────────────────────
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from datetime import datetime
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import warnings

# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

VERSION = '1.0.0'

DEFAULT_TOOL_PATHS: Dict[str, str] = {
    'fastp':    '/data/home/gaojianpeng/miniconda3/bin/fastp',
    'kraken2':  '/data/home/gaojianpeng/mTracePatho/software/kraken2-2.1.2/kraken2',
    'bwa':      '/data/home/gaojianpeng/mTracePatho/software/bwa',
    'samtools': '/data/home/gaojianpeng/miniconda3/bin/samtools',
    'bcftools': '/data/home/gaojianpeng/miniconda3/bin/bcftools',
    'bamdst':   '/data/home/gaojianpeng/mTracePatho/software/bamdst/bamdst',
    'kz':       '/data/home/gaojianpeng/mTracePatho/software/kz',
}

DEFAULT_DB_PATHS: Dict[str, str] = {
    'kraken2_pathogen_db': '/ssd/kraken2_db1',
    'kraken2_plsdb':       '/data/home/gaojianpeng/mTracePatho/PLSDB_kraken2_index',
    'kraken2_homo_db':     '/data/home/gaojianpeng/mTracePatho/homo_kraken2/',
}

CONFIG_SEARCH_PATHS = [
    './tools.ini',
    os.path.expanduser('~/.mtracepatho/tools.ini'),
]

# Database files expected in --db_dir (note: total.var.tsv, not .xls)
DB_FILES_REQUIRED = [
    'total.var.tsv',
    'Jaccard.fastme.input.nwk',
    'node.tsv',
    'numbered_tree.tree',
    'snp_lingage.tsv',
    'lingage.stat.tsv',
    'LCA.output.tsv',
]

CONFIG_TEMPLATE = """\
# mTracePathoRun tools configuration
# Edit paths to match your environment and save as tools.ini
# Search order: --config arg > ./tools.ini > ~/.mtracepatho/tools.ini > built-in defaults

[tools]
fastp     = /data/home/gaojianpeng/miniconda3/bin/fastp
kraken2   = /root/mTracePatho/software/kraken2-2.1.2/kraken2
bwa       = /software/bin/bwa
samtools  = /software/samtools-1.16/samtools
bcftools  = /software/bcftools-1.16/bcftools
bamdst    = /gfs/users/gaojianpeng/software/bamdst/bamdst
kz        = /newVol/gaojianpeng/software/kz

[databases]
kraken2_pathogen_db = /ssd/kraken2_db1
kraken2_plsdb       = /ssd/PLSDB_kraken2_index
kraken2_homo_db     = /ssd/homo_kraken2
"""


# ═══════════════════════════════════════════════════════════════════════════════
# INFRASTRUCTURE: TimeRSSTracker
# ═══════════════════════════════════════════════════════════════════════════════

class TimeRSSTracker:
    """
    Write per-step timing and RSS memory records to logs/time.RSS.TSV (append).

    TSV columns: Timestamp / Module / Submodule / Sample /
                 ElapsedSec / RSS_MB / Status / Notes
    """
    HEADER = "Timestamp\tModule\tSubmodule\tSample\tElapsedSec\tRSS_MB\tStatus\tNotes\n"

    def __init__(self, tsv_path: str):
        self.tsv_path = tsv_path
        os.makedirs(os.path.dirname(os.path.abspath(tsv_path)), exist_ok=True)
        if not os.path.isfile(tsv_path):
            with open(tsv_path, 'w', encoding='utf-8') as f:
                f.write(self.HEADER)

    def record(self, module: str, submodule: str, sample: str,
               elapsed_sec: float, rss_mb: float,
               status: str = 'DONE', notes: str = ''):
        ts = time.strftime('%Y-%m-%dT%H:%M:%S')
        notes = notes[:120]
        row = (f"{ts}\t{module}\t{submodule}\t{sample}\t"
               f"{elapsed_sec:.3f}\t{rss_mb:.1f}\t{status}\t{notes}\n")
        with open(self.tsv_path, 'a', encoding='utf-8') as f:
            f.write(row)


# ═══════════════════════════════════════════════════════════════════════════════
# INFRASTRUCTURE: PipelineLogger
# ═══════════════════════════════════════════════════════════════════════════════

class PipelineLogger:
    """Structured logger with step timing, memory tracking, and TSV recording."""

    def __init__(self, log_file: str, tsv_file: str = None,
                 level=logging.INFO):
        self.log_file = log_file
        self._step_starts: Dict[str, float] = {}
        self._step_rss_start: Dict[str, float] = {}
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)

        fmt = logging.Formatter(
            '[%(asctime)s] [%(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger('mTracePathoRun')
        self.logger.setLevel(level)
        self.logger.handlers.clear()

        fh = logging.FileHandler(log_file, encoding='utf-8')
        fh.setFormatter(fmt)
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(fmt)
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

        self.tracker: Optional[TimeRSSTracker] = (
            TimeRSSTracker(tsv_file) if tsv_file else None
        )

    def set_tracker(self, tracker: 'TimeRSSTracker'):
        self.tracker = tracker

    def info(self, msg):    self.logger.info(msg)
    def warning(self, msg): self.logger.warning(msg)
    def error(self, msg):   self.logger.error(msg)
    def debug(self, msg):   self.logger.debug(msg)

    def step_start(self, step_name: str):
        self._step_starts[step_name] = time.time()
        self._step_rss_start[step_name] = self._rss_mb()
        mem = self._step_rss_start[step_name]
        self.info('=' * 60)
        self.info(f"START  {step_name}  [RSS: {mem:.1f} MB]")

    def step_end(self, step_name: str, success: bool = True,
                 module: str = '', submodule: str = '',
                 sample: str = '-', notes: str = ''):
        elapsed = time.time() - self._step_starts.get(step_name, time.time())
        mem = self._rss_mb()
        status_str = "DONE " if success else "FAILED"
        self.info(f"{status_str} {step_name}  [{elapsed:.1f}s | RSS: {mem:.1f} MB]")
        self.info('=' * 60)
        if self.tracker:
            mod = module or step_name.split(':')[0].strip().replace(' ', '_').lower()
            sub = submodule or step_name
            self.tracker.record(
                module=mod, submodule=sub, sample=sample,
                elapsed_sec=elapsed, rss_mb=mem,
                status='DONE' if success else 'FAILED',
                notes=notes
            )

    def record_submodule(self, module: str, submodule: str,
                         elapsed_sec: float, rss_mb: float,
                         sample: str = '-', status: str = 'DONE',
                         notes: str = ''):
        if self.tracker:
            self.tracker.record(module=module, submodule=submodule,
                                sample=sample, elapsed_sec=elapsed_sec,
                                rss_mb=rss_mb, status=status, notes=notes)

    def record_sample(self, sample: str, elapsed_sec: float,
                      rss_mb: float, status: str = 'DONE', notes: str = ''):
        if self.tracker:
            self.tracker.record(module='step12', submodule='sample',
                                sample=sample, elapsed_sec=elapsed_sec,
                                rss_mb=rss_mb, status=status, notes=notes)

    @staticmethod
    def _rss_mb() -> float:
        try:
            import psutil
            return psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024
        except Exception:
            return 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# INFRASTRUCTURE: StepManager
# ═══════════════════════════════════════════════════════════════════════════════

class StepManager:
    """Manage .done checkpoint files for resume-on-restart capability."""

    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)

    def _done_file(self, step: str) -> Path:
        return self.base_dir / f".{step}.done"

    def _time_file(self, step: str) -> Path:
        return self.base_dir / f".{step}.done.time"

    def _fail_file(self, step: str) -> Path:
        return self.base_dir / f".{step}.fail"

    def is_done(self, step: str) -> bool:
        return self._done_file(step).exists()

    def mark_done(self, step: str, elapsed: float = 0.0):
        self._done_file(step).touch()
        self._time_file(step).write_text(f"{elapsed:.3f}")

    def mark_failed(self, step: str, reason: str):
        self._fail_file(step).write_text(reason[:2000])

    def clear_done(self, step: str):
        for f in [self._done_file(step), self._time_file(step),
                  self._fail_file(step)]:
            if f.exists():
                f.unlink()

    def read_elapsed(self, step: str) -> float:
        try:
            return float(self._time_file(step).read_text().strip())
        except Exception:
            return 0.0

    def read_fail(self, step: str) -> str:
        try:
            return self._fail_file(step).read_text()[:400].strip()
        except Exception:
            return 'unknown'


# ═══════════════════════════════════════════════════════════════════════════════
# INFRASTRUCTURE: ExternalTool
# ═══════════════════════════════════════════════════════════════════════════════

class ExternalTool:
    """Run external commands with logging."""

    def __init__(self, logger: PipelineLogger):
        self.logger = logger

    def run(self, cmd: List[str] = None, shell_pipe: str = None,
            cwd: str = None, log_file: str = None,
            step_label: str = '') -> int:
        if shell_pipe:
            cmd_str = shell_pipe
            t0  = time.time()
            rss0 = self.logger._rss_mb()
            if step_label:
                self.logger.info(f"  {step_label} START  [RSS: {rss0:.1f} MB]")
            self.logger.info(f"  [CMD] {cmd_str}")
            result = subprocess.run(cmd_str, shell=True, cwd=cwd,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            cmd_str = ' '.join(str(c) for c in cmd)
            t0  = time.time()
            rss0 = self.logger._rss_mb()
            if step_label:
                self.logger.info(f"  {step_label} START  [RSS: {rss0:.1f} MB]")
            self.logger.info(f"  [CMD] {cmd_str}")
            result = subprocess.run(cmd, cwd=cwd,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        elapsed = time.time() - t0
        rss1    = self.logger._rss_mb()
        stdout_t = result.stdout.decode('utf-8', errors='replace')
        stderr_t = result.stderr.decode('utf-8', errors='replace')

        if step_label:
            self.logger.info(
                f"  {step_label} DONE   [{elapsed:.1f}s | RSS: {rss1:.1f} MB "
                f"(Δ{rss1-rss0:+.1f} MB)]")

        if log_file:
            with open(log_file, 'a', encoding='utf-8') as lf:
                lf.write(f"CMD: {cmd_str}\n")
                if stderr_t:
                    lf.write(f"STDERR:\n{stderr_t}\n")
                lf.write(f"RETURNCODE: {result.returncode}\n\n")

        if result.returncode != 0:
            self.logger.error(f"Command failed (rc={result.returncode}): {cmd_str}")
            if stderr_t:
                self.logger.error(f"STDERR: {stderr_t[:500]}")
            raise RuntimeError(f"External tool failed: {cmd_str[:200]}")
        return result.returncode


# ═══════════════════════════════════════════════════════════════════════════════
# INFRASTRUCTURE: PipelineConfig
# ═══════════════════════════════════════════════════════════════════════════════

class PipelineConfig:
    """
    Load tool and database paths.
    Priority: --config arg > ./tools.ini > ~/.mtracepatho/tools.ini > built-in defaults
    """
    def __init__(self, config_path: Optional[str] = None,
                 logger: 'PipelineLogger' = None):
        self.tools = dict(DEFAULT_TOOL_PATHS)
        self.databases = dict(DEFAULT_DB_PATHS)
        self._logger = logger
        self._source = 'built-in defaults'
        self._load(config_path)

    def _load(self, config_path: Optional[str]):
        search = ([config_path] if config_path else []) + CONFIG_SEARCH_PATHS
        for path in search:
            if path and os.path.isfile(path):
                cfg = configparser.ConfigParser()
                cfg.read(path, encoding='utf-8')
                if 'tools' in cfg:
                    self.tools.update(dict(cfg['tools']))
                if 'databases' in cfg:
                    self.databases.update(dict(cfg['databases']))
                self._source = path
                break
        if self._logger:
            msg = f"Config loaded from: {self._source}"
            if self._source == 'built-in defaults':
                self._logger.warning(msg + " (recommend creating tools.ini, see --print_config)")
            else:
                self._logger.info(msg)

    def get_tool(self, name: str) -> str:
        return self.tools.get(name, name)

    def get_db(self, name: str) -> str:
        return self.databases.get(name, '')


# ═══════════════════════════════════════════════════════════════════════════════
# EMBEDDED TOOL FUNCTIONS (replacing external Perl/Python scripts)
# ═══════════════════════════════════════════════════════════════════════════════

# ── pullfq (from pullfq.pl) ───────────────────────────────────────────────────
def pullfq(conf_txt: str, fastq_path: str, taxid_str: str, output_fq: str):
    """
    Extract reads from fastq_path whose read_id appears in the kraken2
    --output file (conf_txt) classified to any taxid in taxid_str, plus
    all unclassified reads (flag 'U').

    conf_txt:    kraken2 --output file (tab-sep: C/U, read_id, taxid, ...)
    fastq_path:  FASTQ input, plain or gzip
    taxid_str:   target taxids joined by '-', e.g. "470" or "470-1234"
    output_fq:   output FASTQ path (plain text)

    Equivalent to pullfq.pl.
    """
    taxid_set: Set[str] = set(taxid_str.split('-'))

    keep_reads: Set[str] = set()
    if os.path.isfile(conf_txt) and os.path.getsize(conf_txt) > 0:
        with open(conf_txt, 'r') as fin:
            for line in fin:
                line = line.rstrip('\n\r')
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                flag, read_id, taxid = parts[0], parts[1], parts[2]
                if flag == 'U' or taxid in taxid_set:
                    keep_reads.add(read_id)

    opener = gzip.open if fastq_path.endswith('.gz') else open
    mode   = 'rt' if fastq_path.endswith('.gz') else 'r'
    with opener(fastq_path, mode) as fin, open(output_fq, 'w') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq   = fin.readline()
            plus  = fin.readline()
            qual  = fin.readline()
            if not qual:
                break
            # strip '@' and take first whitespace-separated token
            read_id = header.rstrip('\n\r')[1:].split()[0] if header[1:].strip() else ''
            if read_id in keep_reads:
                fout.write(header if header.endswith('\n') else header + '\n')
                fout.write(seq    if seq.endswith('\n')    else seq    + '\n')
                fout.write(plus   if plus.endswith('\n')   else plus   + '\n')
                fout.write(qual   if qual.endswith('\n')   else qual   + '\n')


# ── filter_bam (from filter_sam.py) ──────────────────────────────────────────
def filter_bam(input_bam: str, output_bam: str, log_file: str,
               identity_thr: float = 0.95, coverage_thr: float = 0.95):
    """
    Filter BAM: retain alignments where identity >= identity_thr AND
    coverage >= coverage_thr.

    identity = (match_length − NM) / read_length
    coverage = (match_length + deletion_length) / read_length

    Equivalent to filter_sam.py with configurable thresholds.
    """
    import pysam

    def _calc(cigar, nm, read_length):
        match_len = del_len = 0
        for op, length in cigar:
            if op == 0:
                match_len += length   # M
            elif op == 2:
                del_len   += length   # D
        identity = (match_len - nm) / read_length if read_length > 0 else 0.0
        coverage = (match_len + del_len) / read_length if read_length > 0 else 0.0
        return identity, coverage

    bam_in  = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)
    with open(log_file, 'w', encoding='utf-8') as log:
        log.write("Filtered Reads Log\n")
        log.write("ReadName\tIdentity\tCoverage\tCIGAR\tNM\tReason\n")
        for read in bam_in.fetch():
            if read.is_unmapped or not read.cigar:
                continue
            try:
                nm = read.get_tag("NM")
            except KeyError:
                nm = 0
            rl = read.query_length
            if not rl:
                continue
            identity, coverage = _calc(read.cigar, nm, rl)
            if identity >= identity_thr and coverage >= coverage_thr:
                bam_out.write(read)
            else:
                reasons = []
                if identity < identity_thr:
                    reasons.append(f"Identity({identity:.2f})<{identity_thr}")
                if coverage < coverage_thr:
                    reasons.append(f"Coverage({coverage:.2f})<{coverage_thr}")
                cigar_s = "".join(f"{l}{'MIDNSHP=XB'[op]}"
                                  for op, l in read.cigar)
                log.write(f"{read.query_name}\t{identity:.2f}\t{coverage:.2f}"
                          f"\t{cigar_s}\t{nm}\t{'; '.join(reasons)}\n")
    bam_in.close()
    bam_out.close()


# ── filter_vcf_py (from filter_vcf.pl) ───────────────────────────────────────
def filter_vcf_py(input_vcf: str, output_vcf: str):
    """
    Filter VCF: keep comment lines; for data lines retain when
      VAF = (alt_fwd + alt_rev) / DP > 0.1
      AND DP > 0
      AND INFO does NOT contain ';AC=1;'

    Equivalent to filter_vcf.pl.
    """
    _dp4_re = re.compile(r'DP4=([^;]+)')
    _dp_re  = re.compile(r'(?:^|;)DP=(\d+)')

    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            stripped = line.rstrip('\n\r')
            if not stripped:
                continue
            fields = stripped.split('\t')
            if len(fields) < 8:
                continue
            info = fields[7]

            dp4_m = _dp4_re.search(info)
            dp_m  = _dp_re.search(info)
            if not dp4_m or not dp_m:
                continue
            dp4_vals = dp4_m.group(1).split(',')
            if len(dp4_vals) < 4:
                continue
            try:
                af = int(dp4_vals[2])
                ar = int(dp4_vals[3])
                dp = int(dp_m.group(1))
            except ValueError:
                continue
            if dp <= 0 or ';AC=1;' in info:
                continue
            if (af + ar) / dp > 0.1:
                fout.write(line)


# ── snp2node (from SNP2node.pl) ──────────────────────────────────────────────
def snp2node(vcf_path: str, lca_output_tsv: str) -> List[Tuple[int, int, List[str]]]:
    """
    Map VCF variant calls to phylogenetic node IDs using the LCA database.

    Reads:
        vcf_path        — filtered VCF file
        lca_output_tsv  — LCA.output.tsv (col0=SNP, col1=LCA_NodeID)

    Returns list of (node_id, count, [snp_strings]) for nodes != 1.

    Variant string format: "{CHROM}:{POS}({REF}->{ALT})"
    Equivalent to SNP2node.pl.
    """
    # Load LCA database
    snp_to_node: Dict[str, str] = {}
    liansuo:     Dict[str, str] = {}

    with open(lca_output_tsv, 'r') as fin:
        reader = csv.reader(fin, delimiter='\t')
        next(reader, None)    # skip header
        for row in reader:
            if len(row) < 2:
                continue
            snp, node_id = row[0], row[1]
            if node_id == 'None':
                continue
            # Strip "ent-" prefix used by some database versions
            if snp.startswith('ent-'):
                snp = snp[4:]
            snp_to_node[snp] = node_id
            if '+' in snp:
                liansuo[snp] = node_id

    # Read VCF and collect all variant strings
    all_vars: List[str] = []
    node_snps: Dict[str, List[str]] = defaultdict(list)

    with open(vcf_path, 'r') as fin:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            var = f"{cols[0]}:{cols[1]}({cols[3]}->{cols[4]})"
            all_vars.append(var)
            if var in snp_to_node:
                node_snps[snp_to_node[var]].append(var)

    # Handle linked (combined) SNP pairs
    if liansuo:
        all_vars_set = set(all_vars)
        for combined, nid in liansuo.items():
            parts = combined.split('+')
            if len(parts) == 2:
                v1, v2 = parts[0], parts[1]
                if v1 in all_vars_set and v2 in all_vars_set:
                    node_snps[nid].append(combined)

    # Build output list (skip node 1, the root)
    results: List[Tuple[int, int, List[str]]] = []
    for nid_str in sorted(node_snps.keys(), key=lambda x: int(x) if x.isdigit() else 0):
        try:
            nid = int(nid_str)
        except ValueError:
            continue
        if nid == 1:
            continue
        snps = node_snps[nid_str]
        results.append((nid, len(snps), snps))
    return results


# ── filt_absent (from filt_absent.pl) ────────────────────────────────────────
def filt_absent(node_results: List[Tuple[int, int, List[str]]]) -> List[Tuple[int, int, List[str]]]:
    """
    Filter out node records whose SNP list contains any item with 'None'.
    All other records are kept as-is.
    Equivalent to filt_absent.pl (most of the Perl code was commented out).
    """
    return [
        (nid, cnt, snps)
        for nid, cnt, snps in node_results
        if not any('None' in s for s in snps)
    ]


# ── write_filt_node ───────────────────────────────────────────────────────────
def write_filt_node(filtered: List[Tuple[int, int, List[str]]], output_path: str):
    """Write .filt.node in the canonical format: node_id,count \\t snp1,snp2,..."""
    with open(output_path, 'w') as fout:
        for nid, cnt, snps in filtered:
            fout.write(f"{nid},{cnt}\t{','.join(snps)}\n")


# ── make_bed (from make_bed.pl) ───────────────────────────────────────────────
def make_bed(input_fa: str, output_bed: str):
    """
    Generate 3-column BED from FASTA: seqid \t 1 \t length
    Start position is 1 (1-based), matching original make_bed.pl behaviour.
    Used internally to produce ref.bed from --ref FASTA automatically.
    """
    with open(input_fa, 'r') as fin, open(output_bed, 'w') as fout:
        seq_id    = None
        seq_parts: List[str] = []
        for line in fin:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if seq_id is not None:
                    fout.write(f"{seq_id}\t1\t{sum(len(p) for p in seq_parts)}\n")
                seq_id    = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if seq_id is not None:
            fout.write(f"{seq_id}\t1\t{sum(len(p) for p in seq_parts)}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# INTERMEDIATE FILE CLEANUP HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def _cleanup_sample_intermediate_files(
        sample: str,
        sample_dir: str,
        files: dict,
        log_fn=None) -> None:
    """
    Remove intermediate files generated by Steps 1+2 after tracing completes.

    Called only when --clean_temp is set.  Preserved files (not in ``files``):
        *.filt.node, *.vcf, *.fastp.json, *.fastp.html,
        *.kraken, *.PLSDB.kraken, *.homo.kraken,
        coverage_report/, TraceResult/,
        .step12.done, .step12.done.time

    Args:
        sample:     Sample name (for log messages).
        sample_dir: Per-sample working directory.
        files:      {label: absolute_path} mapping of files to attempt deletion.
        log_fn:     Optional callable(msg) for per-sample logging.
    """
    def _log(msg):
        if log_fn:
            log_fn(msg)

    removed, skipped, errors = [], [], []
    for label, path in files.items():
        if not path:
            continue
        try:
            if os.path.isfile(path):
                os.remove(path)
                removed.append(os.path.basename(path))
            else:
                skipped.append(os.path.basename(path))
        except OSError as e:
            errors.append(f"{os.path.basename(path)}: {e}")

    if removed:
        _log(f"[{sample}] --clean_temp: removed {len(removed)} intermediate file(s): "
             f"{', '.join(removed)}")
    if errors:
        _log(f"[{sample}] --clean_temp WARNING: failed to remove {len(errors)} file(s): "
             f"{'; '.join(errors)}")


# ═══════════════════════════════════════════════════════════════════════════════
# PER-SAMPLE WORKER FUNCTION (runs in multiprocessing.Pool or SGE job)
# ═══════════════════════════════════════════════════════════════════════════════

def _process_single_sample_run(args: dict) -> Tuple[str, bool, str, float]:
    """
    Process one FASTQ sample through Steps 1, 2, and 3 (per-sample tracing):
      Step 1: kraken2 PLSDB → kraken2 pathogen → pullfq → kz
      Step 2: bwa mem → samtools sort/index → filter_bam → bcftools
              → filter_vcf → SNP2node → filt_absent → .filt.node
      Step 3: build per-sample in.tsv → mTracePathoAlgorithm tracing
              → results written to {sample_dir}/TraceResult/

    Returns (sample_name, success, filt_node_path_or_error, elapsed_seconds).
    Designed to run inside multiprocessing.Pool or an SGE qsub shell script.
    """
    sample        = args['sample']
    fastq_path    = args['fastq_path']
    ref_fa        = args['ref_fa']
    ref_bed       = args['ref_bed']
    db_dir        = args['db_dir']
    sample_dir    = args['sample_dir']
    taxid_str     = args['taxid_str']
    threads       = args['threads']
    tools         = args['tools']
    databases     = args['databases']
    force         = args['force']
    log_path      = args['log_path']
    identity_thr  = args.get('identity_thr', 0.95)
    coverage_thr  = args.get('coverage_thr', 0.95)
    skip_step1    = args.get('skip_step1', False)
    skip_fastp    = args.get('skip_fastp', False)
    skip_rmhost   = args.get('skip_rmhost', False)
    skip_pathogen = args.get('skip_pathogen', False)
    clean_temp    = args.get('clean_temp', False)
    # Step 3 tracing parameters (optional — if absent, Step 3 is skipped)
    tracing_args  = args.get('tracing_args', None)

    os.makedirs(sample_dir, exist_ok=True)
    _sample_t0 = time.time()

    # Per-sample file logger (avoids stdout collision in multiprocess mode)
    _logger = logging.getLogger(f'sample.{sample}')
    _logger.setLevel(logging.DEBUG)
    if not _logger.handlers:
        fh = logging.FileHandler(log_path, encoding='utf-8')
        fh.setFormatter(logging.Formatter('[%(asctime)s] %(message)s',
                                          datefmt='%Y-%m-%d %H:%M:%S'))
        _logger.addHandler(fh)

    def log(msg: str):
        _logger.info(msg)

    def run(cmd_list: list = None, shell_str: str = None):
        # Shell pipes must use bash with pipefail so failures in any pipe
        # stage are propagated correctly (default /bin/sh has no pipefail)
        if shell_str:
            r = subprocess.run(
                ['bash', '-o', 'pipefail', '-c', shell_str],
                cwd=sample_dir,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        else:
            r = subprocess.run(cmd_list, cwd=sample_dir,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stderr = r.stderr.decode('utf-8', errors='replace').strip()
        stdout = r.stdout.decode('utf-8', errors='replace').strip()
        cmd_s  = shell_str or ' '.join(str(c) for c in cmd_list)
        with open(log_path, 'a') as lf:
            lf.write(f"CMD: {cmd_s}\n")
            if stdout:
                lf.write(f"STDOUT: {stdout[:500]}\n")
            if stderr:
                lf.write(f"STDERR: {stderr[:2000]}\n")
            lf.write(f"RC: {r.returncode}\n\n")
        if r.returncode != 0:
            # Include stderr in the exception so .fail files are informative
            stderr_hint = f"\nSTDERR: {stderr[:800]}" if stderr else ""
            raise RuntimeError(
                f"[{sample}] Command failed (rc={r.returncode}): {cmd_s}"
                f"{stderr_hint}")
        return r

    # ── Sentinel files ─────────────────────────────────────────────────────────
    done_file      = os.path.join(sample_dir, '.step12.done')
    filt_node_path = os.path.join(sample_dir, f"{sample}.filt.node")

    if not force and os.path.isfile(done_file) and os.path.isfile(filt_node_path):
        log(f"[{sample}] Already done, skipping")
        return sample, True, filt_node_path, time.time() - _sample_t0

    # ── Intermediate paths ─────────────────────────────────────────────────────
    fastp_fq     = os.path.join(sample_dir, f"{sample}.qc.fq.gz")
    fastp_json   = os.path.join(sample_dir, f"{sample}.fastp.json")
    fastp_html   = os.path.join(sample_dir, f"{sample}.fastp.html")
    homo_rep     = os.path.join(sample_dir, f"{sample}.homo.kraken")
    homo_cls     = os.path.join(sample_dir, "homo.classified.fq")
    homo_unc     = os.path.join(sample_dir, "homo.unclassified.fq")
    homo_conf    = os.path.join(sample_dir, "homo.conf.txt")
    plsdb_rep    = os.path.join(sample_dir, f"{sample}.PLSDB.kraken")
    plsdb_cls    = os.path.join(sample_dir, "PLSDB.classified.fq")
    plsdb_unc    = os.path.join(sample_dir, "PLSDB.unclassified.fq")
    plsdb_conf   = os.path.join(sample_dir, "PLSDB.conf.txt")
    kraken_rep   = os.path.join(sample_dir, f"{sample}.kraken")
    conf_txt     = os.path.join(sample_dir, "conf.txt")
    cls_fq       = os.path.join(sample_dir, "classified.fq")
    uncls_fq     = os.path.join(sample_dir, "unclassified.fq")
    filter_fq    = os.path.join(sample_dir, f"{sample}.filter.fq")
    kz_fq_gz     = os.path.join(sample_dir, f"{sample}.kz.fq.gz")
    sorted_bam   = os.path.join(sample_dir, f"{sample}.sorted.bam")
    filter_bam_f = os.path.join(sample_dir, f"{sample}.sorted.filter.bam")
    filter_log   = os.path.join(sample_dir, "filtered.bam.log")
    raw_vcf      = os.path.join(sample_dir, f"{sample}.raw.vcf")
    vcf_final    = os.path.join(sample_dir, f"{sample}.vcf")
    depth_file   = os.path.join(sample_dir, f"{sample}.depth")
    cov_dir      = os.path.join(sample_dir, "coverage_report")
    os.makedirs(cov_dir, exist_ok=True)

    # Tool shortcuts
    fastp = tools.get('fastp',    '/software/fastp-0.23.4/fastp')
    kr2   = tools.get('kraken2',  'kraken2')
    bwa   = tools.get('bwa',      'bwa')
    samt  = tools.get('samtools', 'samtools')
    bcf   = tools.get('bcftools', 'bcftools')
    bam_t = tools.get('bamdst',   'bamdst')
    kz_t  = tools.get('kz',       'kz')
    plsdb  = databases.get('kraken2_plsdb',       '')
    pathdb = databases.get('kraken2_pathogen_db', '')
    homodb = databases.get('kraken2_homo_db',     '')

    lca_tsv = os.path.join(db_dir, 'LCA.output.tsv')

    try:
        log(f"[{sample}] === Starting sample processing ===")

        # ── STEP 1: Reads preprocessing ────────────────────────────────────────
        if skip_step1:
            log(f"[{sample}] Step 1 skipped; using input FASTQ directly")
            # Use input fastq as kz_fq_gz (may not be kz-filtered)
            kz_fq_gz = fastq_path
        else:
            # Track the current input fastq as we pass through optional steps
            current_fq = fastq_path

            # ── STEP 1-0a: fastp QC (optional) ───────────────────────────────
            if skip_fastp:
                log(f"[{sample}] Step 1-0a: fastp QC skipped (--skip_fastp)")
            else:
                log(f"[{sample}] Step 1-0a: fastp QC")
                run([fastp,
                     '-q', '15', '-u', '10', '--dedup',
                     '-w', str(threads), '-e', '30', '--length_required', '40',
                     '-i', current_fq,
                     '-j', fastp_json,
                     '-h', fastp_html,
                     '-o', fastp_fq])
                if not os.path.isfile(fastp_fq) or os.path.getsize(fastp_fq) == 0:
                    log(f"[{sample}] WARNING: fastp output empty, creating empty .filt.node.")
                    write_filt_node([], filt_node_path)
                    open(done_file, 'w').close()
                    return sample, True, filt_node_path, time.time() - _sample_t0
                current_fq = fastp_fq

            # ── STEP 1-0b: kraken2 remove host reads (optional) ──────────────
            if skip_rmhost:
                log(f"[{sample}] Step 1-0b: host removal skipped (--skip_rmhost)")
            else:
                log(f"[{sample}] Step 1-0b: kraken2 host removal ({homodb})")
                run([kr2, '--db', homodb, '--threads', str(threads),
                     '--report', homo_rep,
                     '--classified-out', homo_cls,
                     '--unclassified-out', homo_unc,
                     '--output', homo_conf,
                     '--memory-mapping', current_fq])
                if not os.path.isfile(homo_unc) or os.path.getsize(homo_unc) == 0:
                    log(f"[{sample}] WARNING: host-removal unclassified file absent/empty, "
                        f"creating empty .filt.node.")
                    write_filt_node([], filt_node_path)
                    open(done_file, 'w').close()
                    return sample, True, filt_node_path, time.time() - _sample_t0
                current_fq = homo_unc  # uncompressed fastq from kraken2 --unclassified-out
            log(f"[{sample}] Step 1-1: kraken2 PLSDB plasmid filter")
            run([kr2, '--db', plsdb, '--threads', str(threads),
                 '--report', plsdb_rep,
                 '--classified-out', plsdb_cls,
                 '--unclassified-out', plsdb_unc,
                 '--output', plsdb_conf,
                 '--memory-mapping', current_fq])

            if not os.path.isfile(plsdb_unc) or os.path.getsize(plsdb_unc) == 0:
                log(f"[{sample}] WARNING: PLSDB unclassified file absent/empty, creating empty")
                open(plsdb_unc, 'w').close()

            if skip_pathogen:
                log(f"[{sample}] Step 1-2/1-3 skipped: pathogen Kraken2 classification "
                    f"and target-taxid pullfq (--skip_pathogen)")
                # PLSDB unclassified reads are already the output of the preceding
                # plasmid-filtering step. Use them directly for complexity filtering.
                shutil.copyfile(plsdb_unc, filter_fq)
            else:
                log(f"[{sample}] Step 1-2: kraken2 pathogen classification")
                run([kr2, '--db', pathdb, '--threads', str(threads),
                     '--report', kraken_rep,
                     '--classified-out', cls_fq,
                     '--unclassified-out', uncls_fq,
                     '--output', conf_txt,
                     '--memory-mapping', plsdb_unc])

                if not os.path.isfile(conf_txt):
                    log(f"[{sample}] WARNING: conf.txt absent, creating empty")
                    open(conf_txt, 'w').close()

                log(f"[{sample}] Step 1-3: pullfq — extract target taxid reads")
                # conf_txt was generated from plsdb_unc, so extract from the same FASTQ.
                pullfq(conf_txt, plsdb_unc, taxid_str, filter_fq)

            if not os.path.isfile(filter_fq) or os.path.getsize(filter_fq) == 0:
                log(f"[{sample}] WARNING: No reads retained after taxid filtering "
                    f"(taxid={taxid_str}). Creating empty .filt.node.")
                write_filt_node([], filt_node_path)
                open(done_file, 'w').close()
                return sample, True, filt_node_path, time.time() - _sample_t0

            log(f"[{sample}] Step 1-4: kz complexity filter")
            run(shell_str=f"{kz_t} -t 0.6 -F < {filter_fq} | gzip > {kz_fq_gz}")

            if not os.path.isfile(kz_fq_gz) or os.path.getsize(kz_fq_gz) == 0:
                log(f"[{sample}] WARNING: kz output is empty, creating empty .filt.node.")
                write_filt_node([], filt_node_path)
                open(done_file, 'w').close()
                return sample, True, filt_node_path, time.time() - _sample_t0

        # ── STEP 2: Alignment and variant detection ─────────────────────────────
        # Auto-build BWA index if any index file is missing.
        # Uses a lock file to prevent concurrent builds when running in parallel
        # (all workers share the same ref_fa; only one should run bwa index).
        _bwa_idx_exts = ['.bwt', '.pac', '.ann', '.amb', '.sa']
        _bwa_lock     = ref_fa + '.bwa_index.lock'
        _missing_idx  = [ext for ext in _bwa_idx_exts
                         if not os.path.isfile(ref_fa + ext)]
        if _missing_idx:
            log(f"[{sample}] BWA index missing ({', '.join(_missing_idx)}), "
                f"building index for {ref_fa} ...")
            # Try to acquire an exclusive lock via atomic file creation
            _lock_acquired = False
            _lock_wait_secs = 0
            while True:
                try:
                    _lf = open(_bwa_lock, 'x')   # atomic; fails if file exists
                    _lf.close()
                    _lock_acquired = True
                    break
                except FileExistsError:
                    # Another worker holds the lock; wait until index is built
                    if _lock_wait_secs == 0:
                        log(f"[{sample}] Another worker is building BWA index, "
                            f"waiting...")
                    import time as _time
                    _time.sleep(5)
                    _lock_wait_secs += 5
                    # Re-check: if index appeared while we waited, stop waiting
                    _still_missing = [ext for ext in _bwa_idx_exts
                                      if not os.path.isfile(ref_fa + ext)]
                    if not _still_missing:
                        log(f"[{sample}] BWA index now present after waiting "
                            f"{_lock_wait_secs}s, continuing.")
                        break
                    if _lock_wait_secs > 3600:
                        raise RuntimeError(
                            f"Waited >1 h for BWA index lock on {ref_fa}; "
                            f"lock file: {_bwa_lock}")

            if _lock_acquired:
                try:
                    # Re-check inside lock: another process may have just finished
                    _still_missing = [ext for ext in _bwa_idx_exts
                                      if not os.path.isfile(ref_fa + ext)]
                    if _still_missing:
                        log(f"[{sample}] Running: {bwa} index {ref_fa}")
                        _idx_r = subprocess.run(
                            [bwa, 'index', ref_fa],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE
                        )
                        _idx_stderr = _idx_r.stderr.decode('utf-8',
                                                           errors='replace').strip()
                        with open(log_path, 'a') as _lf2:
                            _lf2.write(f"CMD: {bwa} index {ref_fa}\n")
                            if _idx_stderr:
                                _lf2.write(f"STDERR: {_idx_stderr[:2000]}\n")
                            _lf2.write(f"RC: {_idx_r.returncode}\n\n")
                        if _idx_r.returncode != 0:
                            raise RuntimeError(
                                f"bwa index failed (rc={_idx_r.returncode})"
                                + (f"\nSTDERR: {_idx_stderr[:800]}"
                                   if _idx_stderr else ""))
                        log(f"[{sample}] BWA index built successfully.")
                    else:
                        log(f"[{sample}] BWA index already built by another "
                            f"worker, skipping.")
                finally:
                    # Always release the lock
                    try:
                        os.remove(_bwa_lock)
                    except OSError:
                        pass

        log(f"[{sample}] Step 2-1: bwa mem + samtools sort + index")
        # Use bash -o pipefail (handled by run()) so bwa failures are caught
        run(shell_str=(
            f"{bwa} mem -t {threads} {ref_fa} {kz_fq_gz} "
            f"| {samt} view -b -@ {threads} "
            f"| {samt} sort -o {sorted_bam} -@ {threads} "
            f"&& {samt} index {sorted_bam} -@ {threads}"
        ))

        log(f"[{sample}] Step 2-2: filter_bam (identity/coverage)")
        filter_bam(sorted_bam, filter_bam_f, filter_log,
                   identity_thr=identity_thr, coverage_thr=coverage_thr)

        log(f"[{sample}] Step 2-3: samtools depth")
        with open(depth_file, 'w') as df_out:
            r = subprocess.run([samt, 'depth', '-aa', filter_bam_f],
                               stdout=df_out, stderr=subprocess.PIPE,
                               cwd=sample_dir)
            if r.returncode != 0:
                _se = r.stderr.decode('utf-8', errors='replace').strip()
                raise RuntimeError(
                    f"samtools depth failed (rc={r.returncode})"
                    + (f"\nSTDERR: {_se[:800]}" if _se else ""))

        log(f"[{sample}] Step 2-4: bamdst coverage statistics")
        run(shell_str=(
            f"{samt} view {filter_bam_f} -u "
            f"| {bam_t} -p {ref_bed} -o {cov_dir} -"
        ))

        log(f"[{sample}] Step 2-5: bcftools mpileup | bcftools call")
        run(shell_str=(
            f"{bcf} mpileup -Ou -f {ref_fa} {filter_bam_f} "
            f"| {bcf} call -mv -Ov -P 0.1 -o {raw_vcf}"
        ))

        log(f"[{sample}] Step 2-6: filter_vcf")
        filter_vcf_py(raw_vcf, vcf_final)

        log(f"[{sample}] Step 2-7: SNP→node mapping + filt_absent")
        node_results = snp2node(vcf_final, lca_tsv)
        filtered     = filt_absent(node_results)
        write_filt_node(filtered, filt_node_path)

        # ── STEP 3: Per-sample tracing ─────────────────────────────────────────
        if tracing_args is not None:
            log(f"[{sample}] Step 3: Per-sample tracing → TraceResult/")
            _run_per_sample_tracing(
                sample=sample,
                filt_node_path=filt_node_path,
                sample_dir=sample_dir,
                db_dir=db_dir,
                tracing_args=tracing_args,
                log_fn=log,
                force=force,
            )

        # ── Clean intermediate files (--clean_temp) ────────────────────────────
        if clean_temp:
            _cleanup_sample_intermediate_files(sample, sample_dir, {
                'fastp_fq':     fastp_fq,
                'homo_cls':     homo_cls,
                'homo_unc':     homo_unc,
                'homo_conf':    homo_conf,
                'plsdb_cls':    plsdb_cls,
                'plsdb_unc':    plsdb_unc,
                'plsdb_conf':   plsdb_conf,
                'cls_fq':       cls_fq,
                'uncls_fq':     uncls_fq,
                'conf_txt':     conf_txt,
                'filter_fq':    filter_fq,
                'kz_fq_gz':     kz_fq_gz,
                'sorted_bam':   sorted_bam,
                'sorted_bai':   sorted_bam + '.bai',
                'filter_bam':   filter_bam_f,
                'filter_bai':   filter_bam_f + '.bai',
                'raw_vcf':      raw_vcf,
                'depth_file':   depth_file,
            }, log_fn=log)

        log(f"[{sample}] Step 2-8: Marking done")
        open(done_file, 'w').close()
        elapsed = time.time() - _sample_t0
        with open(os.path.join(sample_dir, '.step12.done.time'), 'w') as tf:
            tf.write(f"{elapsed:.3f}")
        log(f"[{sample}] === Completed successfully [{elapsed:.1f}s] ===")
        return sample, True, filt_node_path, elapsed

    except Exception as exc:
        elapsed = time.time() - _sample_t0
        err_msg = traceback.format_exc()
        log(f"[{sample}] FAILED after {elapsed:.1f}s:\n{err_msg}")
        fail_f = os.path.join(sample_dir, '.step12.fail')
        with open(fail_f, 'w') as ff:
            ff.write(err_msg)
        return sample, False, str(exc), elapsed


# ═══════════════════════════════════════════════════════════════════════════════
# SGE CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

class SGEConfig:
    """SGE/qsub submission parameters."""

    def __init__(self, nodes: str = '', mem: str = '4G', slots: int = 10,
                 queue: str = 'all.q', max_jobs: int = 200,
                 poll_sec: int = 30, extra: str = ''):
        self.nodes    = nodes.strip()
        self.mem      = mem.strip()
        self.slots    = slots
        self.queue    = queue.strip()
        self.max_jobs = max_jobs
        self.poll_sec = poll_sec
        self.extra    = extra.strip()

    def build_resource_str(self) -> str:
        parts = [f"vf={self.mem}", f"num_proc={self.slots}"]
        if self.nodes:
            parts.append(f"h='{self.nodes}'")
        return ",".join(parts)

    def build_qsub_cmd(self, job_name: str, out_log: str,
                       err_log: str, script_path: str) -> List[str]:
        cmd = ['qsub', '-cwd',
               '-N', job_name,
               '-o', out_log,
               '-e', err_log,
               '-l', self.build_resource_str(),
               '-q', self.queue]
        if self.extra:
            cmd.extend(self.extra.split())
        cmd.append(script_path)
        return cmd

    def __repr__(self) -> str:
        return (f"SGEConfig(nodes={self.nodes!r}, mem={self.mem!r}, "
                f"slots={self.slots}, queue={self.queue!r}, "
                f"max_jobs={self.max_jobs}, poll_sec={self.poll_sec}s)")


# ═══════════════════════════════════════════════════════════════════════════════
# LOCAL PARALLEL PROCESSOR (imap_unordered with real-time progress)
# ═══════════════════════════════════════════════════════════════════════════════

class SampleRunProcessor:
    """
    Run per-sample Steps 1+2+3 in local parallel with real-time progress output.
    Uses Pool.imap_unordered so progress is shown as each sample finishes.
    """

    def __init__(self, sample_list: List[Tuple[str, str]],
                 ref_fa: str, ref_bed: str, db_dir: str, outdir: str,
                 taxid_str: str, config: PipelineConfig,
                 logger: PipelineLogger, threads: int, force: bool,
                 log_dir: str, identity_thr: float, coverage_thr: float,
                 tracing_args: dict = None,
                 skip_step1: bool = False,
                 skip_fastp: bool = False,
                 skip_rmhost: bool = False,
                 skip_pathogen: bool = False,
                 clean_temp: bool = False):
        self.sample_list  = sample_list   # [(sample_name, fastq_path), ...]
        self.ref_fa       = ref_fa
        self.ref_bed      = ref_bed
        self.db_dir       = db_dir
        self.outdir       = outdir
        self.taxid_str    = taxid_str
        self.config       = config
        self.logger       = logger
        self.threads      = threads
        self.force        = force
        self.log_dir      = log_dir
        self.identity_thr = identity_thr
        self.coverage_thr = coverage_thr
        self.tracing_args = tracing_args   # dict of Step 3 params; None = skip
        self.skip_step1   = skip_step1
        self.skip_fastp   = skip_fastp
        self.skip_rmhost  = skip_rmhost
        self.skip_pathogen = skip_pathogen
        self.clean_temp   = clean_temp

    @staticmethod
    def _eta_str(elapsed: float, done: int, total: int) -> str:
        if done <= 0 or elapsed <= 0:
            return "calculating..."
        remaining = max(0, (total - done) / (done / elapsed))
        if remaining < 120:   return f"{remaining:.0f}s"
        if remaining < 7200:  return f"{remaining/60:.0f}min"
        return f"{remaining/3600:.1f}h"

    def run(self) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
        """
        Returns (succeeded_list, failed_list) where each entry is
        (sample_name, filt_node_path_or_error).
        """
        n = len(self.sample_list)
        self.logger.step_start(f"Steps 1+2: Per-sample processing ({n} samples, local)")

        tools     = dict(self.config.tools)
        databases = dict(self.config.databases)

        # Partition: cached vs. needs running
        cached:   List[Tuple[str, str]] = []
        job_args: List[dict] = []

        threads_per_job = max(1, self.threads // min(self.threads, n))

        for sample, fastq_path in self.sample_list:
            sample_dir     = os.path.join(self.outdir, sample)
            done_f         = os.path.join(sample_dir, '.step12.done')
            filt_node_path = os.path.join(sample_dir, f"{sample}.filt.node")
            if not self.force and os.path.isfile(done_f) and os.path.isfile(filt_node_path):
                cached.append((sample, filt_node_path))
                self.logger.debug(f"  [CACHED] {sample}")
                continue
            log_path = os.path.join(self.log_dir, f"{sample}.log")
            job_args.append({
                'sample':       sample,
                'fastq_path':   fastq_path,
                'ref_fa':       self.ref_fa,
                'ref_bed':      self.ref_bed,
                'db_dir':       self.db_dir,
                'sample_dir':   sample_dir,
                'taxid_str':    self.taxid_str,
                'threads':      threads_per_job,
                'tools':        tools,
                'databases':    databases,
                'force':        self.force,
                'log_path':     log_path,
                'identity_thr':  self.identity_thr,
                'coverage_thr':  self.coverage_thr,
                'skip_step1':    self.skip_step1,
                'skip_fastp':    self.skip_fastp,
                'skip_rmhost':   self.skip_rmhost,
                'skip_pathogen': self.skip_pathogen,
                'clean_temp':    self.clean_temp,
                'tracing_args':  self.tracing_args,
            })

        if cached:
            self.logger.info(f"  {len(cached)}/{n} already done (checkpoint), skipping.")

        n_cached  = len(cached)
        results   = list(cached)
        failed:   List[Tuple[str, str]] = []
        n_pending = len(job_args)

        if n_pending == 0:
            self.logger.info("  All samples already complete — nothing to run.")
        else:
            pool_size = min(self.threads, n_pending)
            self.logger.info(
                f"  Parallel workers: {pool_size}  threads/job: {threads_per_job}  "
                f"pending: {n_pending}")

            t_start    = time.time()
            n_finished = n_cached
            n_ok       = n_cached
            n_fail     = 0

            with Pool(pool_size) as pool:
                for ret in pool.imap_unordered(_process_single_sample_run, job_args):
                    sample, success, info, s_elapsed = ret
                    n_finished += 1
                    rss  = self.logger._rss_mb()
                    pct  = n_finished / n * 100
                    wall = time.time() - t_start
                    eta  = self._eta_str(wall, n_finished - n_cached, n_pending)

                    if success:
                        n_ok += 1
                        results.append((sample, info))
                        self.logger.info(
                            f"  [OK]   {sample}  [{s_elapsed:.1f}s | RSS:{rss:.0f}MB]\n"
                            f"  Progress: {n_finished}/{n} ({pct:.1f}%)  "
                            f"done={n_ok} fail={n_fail}  ETA≈{eta}")
                        self.logger.record_sample(
                            sample=sample, elapsed_sec=s_elapsed,
                            rss_mb=rss, status='DONE')
                    else:
                        n_fail += 1
                        failed.append((sample, str(info)))
                        self.logger.error(
                            f"  [FAIL] {sample}  [{s_elapsed:.1f}s]\n"
                            f"  Progress: {n_finished}/{n} ({pct:.1f}%)  "
                            f"done={n_ok} fail={n_fail}  ETA≈{eta}\n"
                            f"  reason: {str(info)[:80]}")
                        self.logger.record_sample(
                            sample=sample, elapsed_sec=s_elapsed,
                            rss_mb=rss, status='FAILED', notes=str(info)[:120])

        self.logger.step_end(
            "Steps 1+2: Per-sample processing",
            success=(len(failed) == 0),
            module='step12', submodule='step12_all_samples',
            notes=f"total={n} ok={len(results)} fail={len(failed)} cached={n_cached}")
        self.logger.info(
            f"  mTracePathoRun 完成: total={n} ok={len(results)} "
            f"fail={len(failed)} cached={n_cached}")
        return results, failed


# ═══════════════════════════════════════════════════════════════════════════════
# SGE JOB SUBMITTER
# ═══════════════════════════════════════════════════════════════════════════════

class SGERunSubmitter:
    """
    Orchestrates per-sample Steps 1+2+3 across an SGE cluster.

    Generates self-contained bash scripts that re-import _process_single_sample_run
    from this module and execute it (includes per-sample Step 3 tracing).
    Sentinel files:
        .step12.done          — Steps 1+2 (+3) succeeded
        .step12.done.time     — wall-clock seconds (float)
        .step12.fail          — failure traceback
    TraceResult per sample:
        {sample_dir}/TraceResult/.trace.done   — tracing succeeded
        {sample_dir}/TraceResult/Trace.result.tsv
    """
    DONE_SENTINEL = '.step12.done'
    FAIL_SENTINEL = '.step12.fail'
    TIME_SENTINEL = '.step12.done.time'

    def __init__(self, sample_list: List[Tuple[str, str]],
                 ref_fa: str, ref_bed: str, db_dir: str, outdir: str,
                 taxid_str: str, config: PipelineConfig,
                 sge_cfg: SGEConfig, logger: PipelineLogger,
                 threads_per_job: int, force: bool, log_dir: str,
                 identity_thr: float, coverage_thr: float,
                 tracing_args: dict = None,
                 skip_step1: bool = False,
                 skip_fastp: bool = False,
                 skip_rmhost: bool = False,
                 skip_pathogen: bool = False,
                 clean_temp: bool = False,
                 python_exe: str = None,
                 sge_qsubname: str = 'mTracePatho'):
        self.sample_list     = sample_list
        self.sge_qsubname    = sge_qsubname
        self._job_counter    = 0  # sequential job number for qsub job names
        self.ref_fa          = ref_fa
        self.ref_bed         = ref_bed
        self.db_dir          = db_dir
        self.outdir          = outdir
        self.taxid_str       = taxid_str
        self.config          = config
        self.sge_cfg         = sge_cfg
        self.logger          = logger
        self.threads_per_job = threads_per_job
        self.force           = force
        self.log_dir         = log_dir
        self.identity_thr    = identity_thr
        self.coverage_thr    = coverage_thr
        self.tracing_args    = tracing_args   # dict of Step 3 params; None = skip
        self.skip_step1      = skip_step1
        self.skip_fastp      = skip_fastp
        self.skip_rmhost     = skip_rmhost
        self.skip_pathogen   = skip_pathogen
        self.clean_temp      = clean_temp
        self.python_exe      = python_exe or sys.executable
        self.module_path     = os.path.abspath(__file__)
        self.module_dir      = os.path.dirname(self.module_path)
        self.script_dir      = os.path.join(log_dir, 'sge_scripts')
        self.sge_log_dir     = os.path.join(log_dir, 'sge_logs')
        os.makedirs(self.script_dir,  exist_ok=True)
        os.makedirs(self.sge_log_dir, exist_ok=True)

    def _sample_dir(self, sample: str) -> str:
        return os.path.join(self.outdir, sample)

    def _sentinel(self, sample: str, name: str) -> str:
        return os.path.join(self._sample_dir(sample), name)

    def _is_done(self, sample: str) -> bool:
        return os.path.isfile(self._sentinel(sample, self.DONE_SENTINEL))

    def _is_failed(self, sample: str) -> bool:
        return os.path.isfile(self._sentinel(sample, self.FAIL_SENTINEL))

    def _fail_reason(self, sample: str) -> str:
        try:
            return open(self._sentinel(sample, self.FAIL_SENTINEL)).read(3000).strip()
        except Exception:
            return 'unknown'

    def _read_elapsed(self, sample: str) -> float:
        try:
            return float(open(self._sentinel(sample, self.TIME_SENTINEL)).read().strip())
        except Exception:
            return 0.0

    def _write_job_script(self, sample: str, fastq_path: str) -> str:
        sample_dir  = self._sample_dir(sample)
        os.makedirs(sample_dir, exist_ok=True)
        script_path = os.path.join(self.script_dir, f"{sample}.sh")
        log_path    = os.path.join(self.log_dir, f"{sample}.log")

        tools_repr = repr(dict(self.config.tools))
        dbs_repr   = repr(dict(self.config.databases))

        tracing_repr = repr(self.tracing_args)

        py_code = (
            "import sys, os, time\n"
            f"sys.path.insert(0, {repr(self.module_dir)})\n"
            "from mTracePathoRun import _process_single_sample_run, "
            "_run_per_sample_tracing, "
            "pullfq, filter_bam, filter_vcf_py, snp2node, filt_absent, "
            "write_filt_node, make_bed\n"
            "\n"
            f"args = dict(\n"
            f"    sample        = {repr(sample)},\n"
            f"    fastq_path    = {repr(fastq_path)},\n"
            f"    ref_fa        = {repr(self.ref_fa)},\n"
            f"    ref_bed       = {repr(self.ref_bed)},\n"
            f"    db_dir        = {repr(self.db_dir)},\n"
            f"    sample_dir    = {repr(sample_dir)},\n"
            f"    taxid_str     = {repr(self.taxid_str)},\n"
            f"    threads       = {self.threads_per_job},\n"
            f"    tools         = {tools_repr},\n"
            f"    databases     = {dbs_repr},\n"
            f"    force         = {repr(self.force)},\n"
            f"    log_path      = {repr(log_path)},\n"
            f"    identity_thr  = {self.identity_thr},\n"
            f"    coverage_thr  = {self.coverage_thr},\n"
            f"    skip_step1    = {repr(self.skip_step1)},\n"
            f"    skip_fastp    = {repr(self.skip_fastp)},\n"
            f"    skip_rmhost   = {repr(self.skip_rmhost)},\n"
            f"    skip_pathogen = {repr(self.skip_pathogen)},\n"
            f"    clean_temp    = {repr(self.clean_temp)},\n"
            f"    tracing_args  = {tracing_repr},\n"
            ")\n"
            "\n"
            f"done_f = {repr(self._sentinel(sample, self.DONE_SENTINEL))}\n"
            f"time_f = {repr(self._sentinel(sample, self.TIME_SENTINEL))}\n"
            f"fail_f = {repr(self._sentinel(sample, self.FAIL_SENTINEL))}\n"
            "\n"
            "try:\n"
            "    _s, success, info, elapsed = _process_single_sample_run(args)\n"
            "    if success:\n"
            "        open(time_f, 'w').write(f'{elapsed:.3f}')\n"
            "        sys.exit(0)\n"
            "    else:\n"
            "        open(fail_f, 'w').write(str(info))\n"
            "        sys.exit(1)\n"
            "except Exception as _ex:\n"
            "    import traceback\n"
            "    open(fail_f, 'w').write(traceback.format_exc())\n"
            "    sys.exit(1)\n"
        )

        script = (
            "#!/bin/bash\n"
            "#$ -S /bin/bash\n"
            "#$ -V\n"
            f"# mTracePathoRun SGE job — sample: {sample}\n"
            f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            "\n"
            "set -euo pipefail\n"
            f'echo "[$(date \'+%Y-%m-%d %H:%M:%S\')] SGE job START: {sample}  host=$(hostname)"\n'
            "\n"
            f"{self.python_exe} << 'MTRACEPATHORUN_PYEOF'\n"
            f"{py_code}"
            "MTRACEPATHORUN_PYEOF\n"
            "\n"
            "_exit=$?\n"
            f'echo "[$(date \'+%Y-%m-%d %H:%M:%S\')] SGE job END: {sample}  exit=$_exit"\n'
            "exit $_exit\n"
        )

        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        return script_path

    def _submit(self, sample: str, script_path: str) -> Optional[str]:
        self._job_counter += 1
        prefix   = getattr(self, 'sge_qsubname', 'mTracePatho')
        job_name = f"{prefix}_{self._job_counter}"[:64]
        out_log  = os.path.join(self.sge_log_dir, f"{sample}.o")
        err_log  = os.path.join(self.sge_log_dir, f"{sample}.e")
        cmd      = self.sge_cfg.build_qsub_cmd(job_name, out_log, err_log, script_path)
        cmd_str  = ' '.join(str(c) for c in cmd)
        self.logger.info(f"  [qsub] {cmd_str}")
        try:
            res    = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout = res.stdout.decode('utf-8', errors='replace').strip()
            stderr = res.stderr.decode('utf-8', errors='replace').strip()
            if res.returncode != 0:
                self.logger.error(f"  [qsub FAIL] {sample}: {stderr[:200]}")
                return None
            m = re.search(r'Your job(?:-array)?\s+(\d+)', stdout)
            if m:
                jid = m.group(1)
                self.logger.debug(f"  [qsub OK]  {sample} → job {jid}")
                return jid
            self.logger.warning(f"  [qsub] Cannot parse job ID from: {stdout!r}")
            return None
        except FileNotFoundError:
            self.logger.error("  [SGE] 'qsub' not found in PATH.")
            return None

    @staticmethod
    def _active_job_ids() -> Set[str]:
        try:
            user = os.environ.get('USER', os.environ.get('LOGNAME', ''))
            cmd  = ['qstat'] + (['-u', user] if user else [])
            res  = subprocess.run(cmd, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
            if res.returncode != 0:
                return set()
            ids = set()
            for line in res.stdout.decode('utf-8', errors='replace').splitlines():
                parts = line.split()
                if parts and parts[0].isdigit():
                    ids.add(parts[0])
            return ids
        except FileNotFoundError:
            return set()

    @staticmethod
    def _eta_str(elapsed: float, done: int, total: int) -> str:
        if done <= 0 or elapsed <= 0:
            return "calculating..."
        remaining = max(0, (total - done) / (done / elapsed))
        if remaining < 120:   return f"{remaining:.0f}s"
        if remaining < 7200:  return f"{remaining/60:.0f}min"
        return f"{remaining/3600:.1f}h"

    def _log_progress(self, done: int, failed: int, total: int,
                      t_start: float, n_active: int,
                      n_submitted: int = 0, n_pending_submit: int = 0):
        finished = done + failed
        pct      = finished / total * 100 if total else 0.0
        elapsed  = time.time() - t_start
        eta      = self._eta_str(elapsed, finished, total)
        avg_t    = elapsed / finished if finished > 0 else 0
        wall_str = (f"{elapsed/60:.1f}min" if elapsed < 7200
                    else f"{elapsed/3600:.1f}h")
        parts = [
            f"[SGE Progress] {finished}/{total} ({pct:.1f}%)",
            f"done={done}", f"fail={failed}", f"active_jobs={n_active}",
        ]
        if n_submitted:
            parts.append(f"submitted={n_submitted}")
        if n_pending_submit:
            parts.append(f"pending_submit={n_pending_submit}")
        parts += [f"avg={avg_t:.0f}s/sample", f"ETA≈{eta}", f"wall={wall_str}"]
        self.logger.info("  " + "  ".join(parts))

    def run(self) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
        n = len(self.sample_list)
        self.logger.step_start(
            f"Steps 1+2 [SGE mode]: {n} samples  {self.sge_cfg!r}")
        self.logger.info(
            f"  qsub resource: -l {self.sge_cfg.build_resource_str()}"
            f"  -q {self.sge_cfg.queue}")

        already_done: List[str] = []
        pending: List[Tuple[str, str]] = []
        for sample, fastq_path in self.sample_list:
            if not self.force and self._is_done(sample):
                already_done.append(sample)
                self.logger.debug(f"  [CACHED] {sample}")
            else:
                fail_f = self._sentinel(sample, self.FAIL_SENTINEL)
                if os.path.isfile(fail_f):
                    os.remove(fail_f)
                pending.append((sample, fastq_path))

        if already_done:
            self.logger.info(
                f"  {len(already_done)}/{n} already done (checkpoint), skipping.")
        if not pending:
            self.logger.info("  All samples already done — nothing to submit.")
            succeeded = [(s, os.path.join(self._sample_dir(s), f"{s}.filt.node"))
                         for s in already_done]
            self.logger.step_end("Steps 1+2 [SGE]: All cached", success=True,
                                 module='step12', submodule='step12_sge_all')
            return succeeded, []

        # Write job scripts
        self.logger.info(f"  Writing {len(pending)} job scripts → {self.script_dir}/")
        scripts: Dict[str, str] = {}
        for sample, fastq_path in pending:
            scripts[sample] = self._write_job_script(sample, fastq_path)
        self.logger.info(f"  {len(scripts)} scripts written. Starting submission...")

        # Submission + monitoring loop
        t_start          = time.time()
        submit_queue     = list(pending)
        submitted_jids:  Dict[str, str] = {}
        failed_submit:   List[str]      = []
        n_submitted_total = 0

        while submit_queue:
            active_ids = self._active_job_ids()
            slots_free = max(0, self.sge_cfg.max_jobs - len(active_ids))
            if slots_free == 0:
                self.logger.info(
                    f"  [SGE] Cluster at capacity ({len(active_ids)} active). "
                    f"Waiting {self.sge_cfg.poll_sec}s...")
                time.sleep(self.sge_cfg.poll_sec)
                continue

            batch        = submit_queue[:slots_free]
            submit_queue = submit_queue[slots_free:]
            for sample, fastq_path in batch:
                jid = self._submit(sample, scripts[sample])
                if jid:
                    submitted_jids[sample] = jid
                    n_submitted_total += 1
                else:
                    failed_submit.append(sample)
                    ff = self._sentinel(sample, self.FAIL_SENTINEL)
                    if not os.path.isfile(ff):
                        open(ff, 'w').write("qsub submission failed")

            all_sub   = list(submitted_jids.keys())
            done_now  = sum(1 for s in all_sub if self._is_done(s))
            fail_now  = sum(1 for s in all_sub if self._is_failed(s))
            self._log_progress(done_now + len(already_done),
                               fail_now + len(failed_submit), n,
                               t_start, len(self._active_job_ids()),
                               n_submitted_total, len(submit_queue))
            if submit_queue:
                time.sleep(self.sge_cfg.poll_sec)

        # Wait for all submitted jobs to finish
        self.logger.info(
            f"  All {len(submitted_jids)} jobs submitted. Waiting for completion...")
        while True:
            all_sub = list(submitted_jids.keys())
            unresolved = [s for s in all_sub
                          if not self._is_done(s) and not self._is_failed(s)]
            if not unresolved:
                break
            active_ids = self._active_job_ids()
            my_active  = {s: jid for s, jid in submitted_jids.items()
                          if jid in active_ids}
            done_now   = sum(1 for s in all_sub if self._is_done(s))
            fail_now   = sum(1 for s in all_sub if self._is_failed(s))
            self._log_progress(done_now + len(already_done),
                               fail_now + len(failed_submit), n,
                               t_start, len(my_active), n_submitted_total)
            if not my_active and unresolved:
                self.logger.warning(
                    f"  {len(unresolved)} job(s) vanished from qstat "
                    f"without sentinel files — marking FAILED.")
                for s in unresolved:
                    ff = self._sentinel(s, self.FAIL_SENTINEL)
                    if not os.path.isfile(ff):
                        open(ff, 'w').write(
                            "Job vanished from qstat without completion sentinel")
                break
            time.sleep(self.sge_cfg.poll_sec)

        # Collect final results
        succeeded: List[Tuple[str, str]] = []
        failed:    List[Tuple[str, str]] = []
        for sample, _ in self.sample_list:
            filt_node = os.path.join(self._sample_dir(sample), f"{sample}.filt.node")
            if self._is_done(sample):
                s_elapsed = self._read_elapsed(sample)
                succeeded.append((sample, filt_node))
                self.logger.info(f"  [OK]   {sample}  [{s_elapsed:.1f}s]")
                self.logger.record_sample(sample=sample, elapsed_sec=s_elapsed,
                                          rss_mb=0.0, status='DONE', notes='sge')
            else:
                reason = self._fail_reason(sample)
                failed.append((sample, reason))
                # Log full reason (truncated only for TSV notes column)
                self.logger.error(f"  [FAIL] {sample}:")
                for _line in reason.splitlines():
                    self.logger.error(f"    {_line}")
                self.logger.record_sample(sample=sample, elapsed_sec=0.0,
                                          rss_mb=0.0, status='FAILED',
                                          notes=f"sge:{reason[:100]}")

        elapsed_total = time.time() - t_start
        self.logger.step_end(
            "Steps 1+2 [SGE mode]", success=(len(failed) == 0),
            module='step12', submodule='step12_sge_all',
            notes=(f"total={n} ok={len(succeeded)} fail={len(failed)} "
                   f"cached={len(already_done)} wall={elapsed_total:.0f}s"))
        return succeeded, failed



# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: TRACING ALGORITHM — imported from mTracePathoAlgorithm
# ═══════════════════════════════════════════════════════════════════════════════

from mTracePathoAlgorithm import (
    calculate_confidence_score, load_lca_output, load_snp_lineage,
    load_snp_lineage_with_genome_map, load_snp_matrix_for_containment,
    snp_to_node_mapping, calculate_cumulative_scores,
    find_group_lca, find_lca_two_nodes,
    merge_genome_info_with_globalpath,
    process_single_sample, save_trace_result, compute_global_path_scores,
    select_candidates_by_containment,
    _worker_init,          # multiprocessing initializer for genome_snps_map
)


# ═══════════════════════════════════════════════════════════════════════════════
# QC METRIC EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════════

def _update_trace_qc_metrics(sample: str, sample_dir: str,
                              filt_node_path: str, trace_dir: str,
                              log_fn=None) -> None:
    """
    After tracing completes, fill QC placeholder columns in Trace.result.tsv:
      Reads, Coverage, Depth, FractionMappedReads, TotalSNP, TracingSnpUsed

    Reads:
      {sample_dir}/{sample}.kz.fq.gz         → Reads (skip if >200 MB)
      {sample_dir}/{sample}.vcf               → TotalSNP
      {sample_dir}/{sample}.filt.node         → TracingSnpUsed count+ratio
      {sample_dir}/coverage_report/coverage.report → Coverage, Depth, FractionMappedReads
    """
    import csv as _csv
    import gzip as _gzip

    result_tsv = os.path.join(trace_dir, 'Trace.result.tsv')
    if not os.path.isfile(result_tsv):
        return

    # ── KzReads ──────────────────────────────────────────────────────────────
    kz_fq_gz = os.path.join(sample_dir, f"{sample}.kz.fq.gz")
    kz_reads = 'N/A'
    try:
        if os.path.isfile(kz_fq_gz):
            size_mb = os.path.getsize(kz_fq_gz) / (1024 * 1024)
            if size_mb > 200:
                kz_reads = '/'
            else:
                line_count = 0
                with _gzip.open(kz_fq_gz, 'rb') as _gz:
                    for _ in _gz:
                        line_count += 1
                kz_reads = str(line_count // 4)  # FASTQ: 4 lines per read
    except Exception:
        kz_reads = 'N/A'

    # ── TotalSNP (from .vcf) ─────────────────────────────────────────────────
    vcf_path  = os.path.join(sample_dir, f"{sample}.vcf")
    total_snp = 'N/A'
    vcf_count = 0
    try:
        if os.path.isfile(vcf_path):
            with open(vcf_path, 'r') as _v:
                for _line in _v:
                    if not _line.startswith('#'):
                        vcf_count += 1
            total_snp = str(vcf_count)
    except Exception:
        total_snp = 'N/A'

    # ── TracingSnpUsed (from .filt.node) ──────────────────────────────────────
    tracing_snp_used = 'N/A'
    try:
        if os.path.isfile(filt_node_path):
            filt_snp_set: set = set()
            with open(filt_node_path, 'r') as _fn:
                for _line in _fn:
                    _parts = _line.rstrip('\n\r').split('\t')
                    if len(_parts) >= 2:
                        for _s in _parts[1].split(','):
                            _s = _s.strip()
                            if _s:
                                filt_snp_set.add(_s)
            filt_count = len(filt_snp_set)
            if vcf_count > 0:
                ratio = round(filt_count / vcf_count, 4)
                tracing_snp_used = f"{filt_count}|{ratio}"
            else:
                tracing_snp_used = str(filt_count)
    except Exception:
        tracing_snp_used = 'N/A'

    # ── Coverage, Depth, FractionMappedReads (from coverage_report) ──────────
    cov_report = os.path.join(sample_dir, 'coverage_report', 'coverage.report')
    coverage_val    = 'N/A'
    depth_val       = 'N/A'
    fraction_mapped = 'N/A'
    try:
        if os.path.isfile(cov_report):
            with open(cov_report, 'r') as _cr:
                for _line in _cr:
                    # coverage.report uses \t as field separator (bamdst format)
                    if '\t' not in _line:
                        continue
                    _label, _sep, _value = _line.rpartition('\t')
                    _label  = _label.strip()
                    _value  = _value.strip()
                    if '[flank] Coverage (>0x)' in _label:
                        coverage_val = _value
                    elif '[Target] Average depth' in _label:
                        depth_val = _value
                    elif '[Total] Fraction of Mapped Reads' in _label:
                        fraction_mapped = _value
    except Exception as e:
        log_fn(f"[{sample}] WARNING: failed to parse coverage.report: {type(e).__name__}: {e}")

    # ── Update Trace.result.tsv ───────────────────────────────────────────────
    try:
        with open(result_tsv, 'r', newline='') as _f:
            reader = _csv.DictReader(_f, delimiter='\t')
            rows   = list(reader)
            fieldnames = reader.fieldnames or []

        update_map = {
            'Reads':                kz_reads,
            'Coverage':             coverage_val,
            'Depth':                depth_val,
            'FractionMappedReads':  fraction_mapped,
            'TotalSNP':             total_snp,
            'TracingSnpUsed':       tracing_snp_used,
        }

        # Columns that should always be overwritten (not just when N/A),
        # because the algorithm pre-fills them with proxy values (filt.node count).
        FORCE_OVERWRITE = {'TotalSNP', 'TracingSnpUsed'}
        updated = False
        for row in rows:
            for col, val in update_map.items():
                if col not in row:
                    continue
                if col in FORCE_OVERWRITE:
                    if val not in ('N/A', None) and val != '':
                        row[col] = val
                        updated = True
                elif row[col] in ('N/A', '', None):
                    row[col] = val
                    updated = True

        if updated:
            with open(result_tsv, 'w', newline='') as _f:
                writer = _csv.DictWriter(_f, delimiter='\t',
                                         fieldnames=fieldnames,
                                         extrasaction='ignore')
                writer.writeheader()
                writer.writerows(rows)
            if log_fn:
                log_fn(f"[{sample}] QC metrics updated in Trace.result.tsv")
    except Exception as _e:
        if log_fn:
            log_fn(f"[{sample}] WARNING: failed to update QC metrics: {_e}")

def _run_per_sample_tracing(sample: str, filt_node_path: str, sample_dir: str,
                             db_dir: str, tracing_args: dict,
                             log_fn, force: bool = False):
    """
    Step 3 for a single sample: build in.tsv → run mTracePathoAlgorithm tracing.

    Outputs are written to  {sample_dir}/TraceResult/
    Sentinel files:
        {sample_dir}/TraceResult/.trace.done   — tracing succeeded
        {sample_dir}/TraceResult/.trace.fail   — failure message

    tracing_args keys (all from argparse namespace converted to dict):
        jaccard_percent, max_iterations, use_containment, min_containment,
        top_n_cutoff, min_candidates, no_visualization, keep_temp
    """
    import csv as _csv

    trace_dir  = os.path.join(sample_dir, 'TraceResult')
    done_file  = os.path.join(trace_dir, '.trace.done')
    fail_file  = os.path.join(trace_dir, '.trace.fail')
    os.makedirs(trace_dir, exist_ok=True)

    # Checkpoint: skip if already done
    if not force and os.path.isfile(done_file):
        log_fn(f"[{sample}] Step 3: TraceResult/.trace.done found — skipping")
        return

    # Remove stale fail file if rerunning
    if os.path.isfile(fail_file):
        os.remove(fail_file)

    try:
        # Functions are embedded directly in this module (no import needed)
        # ── Build per-sample in.tsv ──────────────────────────────────────────
        # Store in.tsv at trace_dir level (sentinels also here)
        in_tsv  = os.path.join(trace_dir, 'in.tsv')
        snp_set: Set[str] = set()
        if os.path.isfile(filt_node_path):
            with open(filt_node_path, 'r') as _fin:
                for _line in _fin:
                    _line = _line.rstrip('\n\r')
                    if not _line:
                        continue
                    _parts = _line.split('\t')
                    if len(_parts) >= 2:
                        for _snp in _parts[1].split(','):
                            _s = _snp.strip()
                            if _s:
                                snp_set.add(_s)

        with open(in_tsv, 'w') as _fout:
            if snp_set:
                _fout.write(f"{sample}\t{','.join(sorted(snp_set))}\n")

        if not snp_set:
            log_fn(f"[{sample}] Step 3: No SNPs in .filt.node — writing empty TraceResult")
            open(done_file, 'w').close()
            return

        log_fn(f"[{sample}] Step 3: {len(snp_set)} SNPs → in.tsv → tracing ...")

        # ── Load database ────────────────────────────────────────────────────
        _, original_genome_to_lineage, _, global_max_depth = \
            load_snp_lineage_with_genome_map(db_dir)

        genome_snps_map = None
        if tracing_args.get('use_containment', False):
            genome_snps_map, _, _ = load_snp_matrix_for_containment(db_dir)

        total_var   = os.path.join(db_dir, 'total.var.tsv')
        jaccard_nwk = os.path.join(db_dir, 'Jaccard.fastme.input.nwk')

        # ── Run tracing via temp parent dir ──────────────────────────────────
        # mTracePathoAlgorithm.process_single_sample() internally creates:
        #   {output_dir}/result/{sample_name}/
        # To get output directly in trace_dir (without extra layers), we use a
        # temporary parent so that: {_algo_parent}/result/{sample} == trace_dir,
        # then move the algo's output into trace_dir and clean up.
        _algo_parent = os.path.join(sample_dir, '.trace_work')
        os.makedirs(_algo_parent, exist_ok=True)
        try:
            result = process_single_sample(
                sample, snp_set, db_dir, _algo_parent,
                tracing_args.get('jaccard_percent',   0.8),
                tracing_args.get('max_iterations',    10),
                tracing_args.get('keep_temp',         True),
                tracing_args.get('no_visualization',  False),
                total_var, jaccard_nwk,
                original_genome_to_lineage, global_max_depth,
                tracing_args.get('use_containment',   False),
                tracing_args.get('min_containment',   0.0),
                tracing_args.get('top_n_cutoff',      99999),
                tracing_args.get('min_candidates',    3),
                genome_snps_map,
                tracing_args.get('phylo_gap_ratio',   999999999999),
                tracing_args.get('phylo_share_thresh', 0.7),
                tracing_args.get('scoring_mode',      'containment_weighted'),
                tracing_args.get('cont_alpha',        0.5),
                tracing_args.get('cont_ref_size',     None),
                tracing_args.get('score_cutoff_ratio', 0.5),
                tracing_args.get('min_snp', 5),
                tracing_args.get('min_unique_snp_n_cutoff', 1),
            )
            # Algorithm writes to: {_algo_parent}/result/{sample}/
            _algo_out = os.path.join(_algo_parent, 'result', sample)
            if os.path.isdir(_algo_out):
                # Merge into trace_dir (already exists; merge rather than replace)
                for _item in os.listdir(_algo_out):
                    _src = os.path.join(_algo_out, _item)
                    _dst = os.path.join(trace_dir, _item)
                    if os.path.isdir(_src):
                        if os.path.isdir(_dst):
                            shutil.rmtree(_dst)
                        shutil.move(_src, _dst)
                    else:
                        shutil.move(_src, _dst)
        finally:
            # Clean up the temp work directory only when --clean_temp is active
            # (keep_temp=False means clean_temp was set)
            if not tracing_args.get('keep_temp', True):
                shutil.rmtree(_algo_parent, ignore_errors=True)

        # ── Write Trace.result.tsv summary entry ─────────────────────────────
        # The algorithm already writes its own Trace.result.tsv inside trace_dir.
        # We additionally write a clean single-row summary for our batch collector.
        result_tsv = os.path.join(trace_dir, 'Trace.result.tsv')
        if not os.path.isfile(result_tsv):
            # Algorithm didn't write it; write from return dict
            fieldnames = ['sample_name', 'iterations', 'final_genomes', 'target',
                          'trace_level', 'confidence_score', 'confidence_level', 'time']
            if tracing_args.get('use_containment', False):
                fieldnames.append('containment_candidates')
            with open(result_tsv, 'w', newline='') as _rt:
                writer = _csv.DictWriter(_rt, delimiter='\t', fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow({k: result.get(k, 'N/A') for k in fieldnames})

        conf_str = (f"{result.get('confidence_score', 'N/A')}"
                    f"({result.get('confidence_level', 'N/A')})")
        log_fn(f"[{sample}] Step 3: DONE  target={result.get('target','?')}  "
               f"level={result.get('trace_level','?')}  conf={conf_str}")

        # ── Populate QC metrics in Trace.result.tsv ──────────────────────────
        _update_trace_qc_metrics(sample, sample_dir, filt_node_path,
                                  trace_dir, log_fn=log_fn)

        open(done_file, 'w').close()

    except Exception as _te:
        import traceback as _tb
        _msg = _tb.format_exc()
        log_fn(f"[{sample}] Step 3: FAILED\n{_msg}")
        with open(fail_file, 'w') as _ff:
            _ff.write(_msg)
        # Step 3 failure is non-fatal for Step 1+2 sentinel; re-raise so caller
        # can decide.  Wrapped to distinguish from Step 1/2 errors.
        raise RuntimeError(f"[{sample}] Step 3 (tracing) failed: {_te}") from _te


# ═══════════════════════════════════════════════════════════════════════════════
# REUSE-MODE STEP 3 WORKER  (module-level so it is picklable by multiprocessing)
# ═══════════════════════════════════════════════════════════════════════════════

def _reuse_trace_worker(args: dict) -> Tuple[str, bool, str, float]:
    """
    Module-level (picklable) worker for --result_dir reuse mode Step 3.

    Mirrors _process_single_sample_run for Steps 1+2 but only performs
    Step 3 tracing on a pre-existing .filt.node file.

    Returns (sample_name, success:bool, filt_node_or_error:str, elapsed:float).
    """
    sample        = args['sample']
    filt_node     = args['filt_node']
    sample_dir    = args['sample_dir']
    db_dir        = args['db_dir']
    tracing_args  = args['tracing_args']
    force         = args.get('force', False)
    log_path      = args.get('log_path', None)

    t0 = time.time()

    # Per-sample file logger (same pattern as _process_single_sample_run)
    _logger = logging.getLogger(f'reuse.{sample}')
    _logger.setLevel(logging.DEBUG)
    if log_path:
        os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
        _fh = logging.FileHandler(log_path, mode='a')
        _fh.setFormatter(logging.Formatter('[%(asctime)s] %(message)s',
                                           datefmt='%H:%M:%S'))
        _logger.addHandler(_fh)

    def _log(msg):
        _logger.info(msg)

    try:
        _run_per_sample_tracing(
            sample         = sample,
            filt_node_path = filt_node,
            sample_dir     = sample_dir,
            db_dir         = db_dir,
            tracing_args   = tracing_args,
            log_fn         = _log,
            force          = force,
        )
        return sample, True, filt_node, time.time() - t0
    except Exception as exc:
        return sample, False, str(exc), time.time() - t0


def _build_in_tsv(succeeded: List[Tuple[str, str]], tsv_path: str):
    """
    Convert per-sample .filt.node files to in.tsv for tracing.

    .filt.node format:  node_id,count \\t snp1,snp2,...
    in.tsv format:      sample_name \\t snp1,snp2,...

    All unique SNP variant strings from a sample's .filt.node are merged
    into one tab-separated record.
    """
    with open(tsv_path, 'w') as fout:
        for sample, filt_node_path in succeeded:
            snp_set: Set[str] = set()
            if os.path.isfile(filt_node_path):
                with open(filt_node_path, 'r') as fin:
                    for line in fin:
                        line = line.rstrip('\n\r')
                        if not line:
                            continue
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            for snp in parts[1].split(','):
                                s = snp.strip()
                                if s:
                                    snp_set.add(s)
            if snp_set:
                fout.write(f"{sample}\t{','.join(sorted(snp_set))}\n")


def run_step3_tracing(succeeded: List[Tuple[str, str]],
                      db_dir: str, outdir: str, args,
                      logger: PipelineLogger):
    """
    Step 3: Load tracing database and run the mTracePatho tracing algorithm
    for all samples.
    """
    # Functions are embedded directly in this module (no import needed)
    result_dir = os.path.join(outdir, 'results')
    os.makedirs(result_dir, exist_ok=True)

    # Build in.tsv
    in_tsv = os.path.join(result_dir, 'in.tsv')
    logger.step_start("Step 3: Build in.tsv from .filt.node files")
    _build_in_tsv(succeeded, in_tsv)
    logger.info(f"  Written: {in_tsv}")
    logger.step_end("Step 3: Build in.tsv", success=True,
                    module='step3', submodule='3-0_build_in_tsv')

    # Verify database files (note: total.var.tsv, not .xls)
    for fname in DB_FILES_REQUIRED:
        fpath = os.path.join(db_dir, fname)
        if not os.path.isfile(fpath):
            logger.error(f"Missing database file: {fpath}")
            raise FileNotFoundError(f"Missing database file: {fpath}")

    total_var  = os.path.join(db_dir, 'total.var.tsv')
    jaccard_nwk = os.path.join(db_dir, 'Jaccard.fastme.input.nwk')

    logger.step_start("Step 3: Load tracing database")
    _, original_genome_to_lineage, _, global_max_depth = \
        load_snp_lineage_with_genome_map(db_dir)
    logger.info(
        f"  Loaded {len(original_genome_to_lineage)} genome lineages, "
        f"max depth: {global_max_depth}")

    genome_snps_map = None
    if args.use_containment:
        logger.info("  Loading SNP matrix for containment...")
        genome_snps_map, _, _ = load_snp_matrix_for_containment(db_dir)
        logger.info(f"  Loaded SNP data for {len(genome_snps_map)} genomes")
    logger.step_end("Step 3: Load database", success=True,
                    module='step3', submodule='3-1_load_db')

    # Load samples
    import csv as _csv
    samples = []
    with open(in_tsv, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                sname = parts[0].strip()
                snps  = set(s.strip() for s in parts[1].split(',') if s.strip())
                if sname and snps:
                    samples.append((sname, snps))

    if not samples:
        logger.warning("  No valid samples in in.tsv — Step 3 skipped.")
        return []

    logger.info(f"  Processing {len(samples)} sample(s) with {args.threads} worker(s)...")
    logger.step_start("Step 3: Tracing analysis")

    batch_results = []

    if args.threads == 1:
        for sample_name, sample_snps in samples:
            t0 = time.time()
            result = process_single_sample(
                sample_name, sample_snps, db_dir, result_dir,
                args.jaccard_percent, args.max_iterations,
                not args.clean_temp, args.no_visualization,
                total_var, jaccard_nwk,
                original_genome_to_lineage, global_max_depth,
                args.use_containment, args.min_containment,
                args.top_n_cutoff, args.min_candidates,
                genome_snps_map,
                args.phylo_gap_ratio, args.phylo_share_thresh,
                args.scoring_mode, args.cont_alpha, args.cont_ref_size,
                args.score_cutoff_ratio, args.min_snp,
                args.min_unique_snp_n_cutoff,
            )
            batch_results.append(result)
            rss = logger._rss_mb()
            conf_str = f"{result['confidence_score']:.2f}({result['confidence_level']})"
            elapsed  = time.time() - t0
            logger.info(
                f"  [OK]   {result['sample_name']}  [{elapsed:.1f}s | RSS:{rss:.0f}MB]"
                f"  Conf:{conf_str}  Level:{result['trace_level']}")
            logger.record_submodule(
                module='step3', submodule='sample',
                sample=result['sample_name'],
                elapsed_sec=elapsed, rss_mb=rss)
    else:
        # Pass None for genome_snps_map; workers load it via _worker_init() to
        # avoid pickling the full map N times (same fix as Algorithm BUG#11).
        sample_args = [
            (sample_name, sample_snps, db_dir, result_dir,
             args.jaccard_percent, args.max_iterations,
             not args.clean_temp, args.no_visualization,
             total_var, jaccard_nwk,
             original_genome_to_lineage, global_max_depth,
             args.use_containment, args.min_containment,
             args.top_n_cutoff, args.min_candidates,
             None,                              # genome_snps_map loaded by _worker_init
             args.phylo_gap_ratio, args.phylo_share_thresh,
             args.scoring_mode, args.cont_alpha, args.cont_ref_size,
             args.score_cutoff_ratio, args.min_snp,
             args.min_unique_snp_n_cutoff)
            for sample_name, sample_snps in samples
        ]
        t_start    = time.time()
        n          = len(samples)
        n_finished = 0
        n_ok       = 0
        n_fail     = 0

        _init_args = (db_dir,) if args.use_containment else None
        with ProcessPoolExecutor(
                max_workers=args.threads,
                initializer=_worker_init if args.use_containment else None,
                initargs=_init_args) as executor:
            futures = {executor.submit(process_sample_wrapper, arg): arg[0]
                       for arg in sample_args}
            for future in as_completed(futures):
                result = future.result()
                batch_results.append(result)
                n_finished += 1
                rss  = logger._rss_mb()
                pct  = n_finished / n * 100
                wall = time.time() - t_start
                remaining = max(0, (n - n_finished) / (n_finished / wall)) if wall > 0 and n_finished > 0 else 0
                if remaining < 120:   eta_s = f"{remaining:.0f}s"
                elif remaining < 7200: eta_s = f"{remaining/60:.0f}min"
                else:                 eta_s = f"{remaining/3600:.1f}h"

                conf_str = f"{result['confidence_score']:.2f}({result['confidence_level']})"
                if result['trace_level'] != 'Error':
                    n_ok += 1
                    logger.info(
                        f"  [OK]   {result['sample_name']}  "
                        f"[RSS:{rss:.0f}MB]  Conf:{conf_str}\n"
                        f"  Progress: {n_finished}/{n} ({pct:.1f}%)  "
                        f"done={n_ok} fail={n_fail}  ETA≈{eta_s}")
                else:
                    n_fail += 1
                    logger.error(
                        f"  [FAIL] {result['sample_name']}  "
                        f"reason: {result.get('target','')[:80]}\n"
                        f"  Progress: {n_finished}/{n} ({pct:.1f}%)  "
                        f"done={n_ok} fail={n_fail}  ETA≈{eta_s}")
                logger.record_submodule(
                    module='step3', submodule='sample',
                    sample=result['sample_name'],
                    elapsed_sec=float(result.get('time', '0').rstrip('s') or 0),
                    rss_mb=rss,
                    status='DONE' if result['trace_level'] != 'Error' else 'FAILED')

    # Sort results to match input order
    sample_order = {name: i for i, (name, _) in enumerate(samples)}
    batch_results.sort(key=lambda x: sample_order.get(x['sample_name'], len(samples)))

    # Write batch summary
    summary_path = os.path.join(result_dir, 'batch_summary.tsv')
    with open(summary_path, 'w', newline='') as f:
        fieldnames = ['sample_name', 'iterations', 'final_genomes', 'target',
                      'trace_level', 'confidence_score', 'confidence_level', 'time']
        if args.use_containment:
            fieldnames.append('containment_candidates')
        writer = _csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for r in batch_results:
            writer.writerow({k: r.get(k, 'N/A') for k in fieldnames})
    logger.info(f"  Batch summary: {summary_path}")

    logger.step_end("Step 3: Tracing analysis", success=True,
                    module='step3', submodule='step3_tracing_all',
                    notes=f"samples={len(samples)}")
    return batch_results



# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4a: PER-SAMPLE HTML VIEWER GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

def _find_meta_tsv(db_dir: str) -> Optional[str]:
    """Locate the *-gene.tsv metadata file inside the database directory.

    Convention: the file is named ``<SpeciesName>-gene.tsv`` (e.g.
    ``Klebsiella_pneumoniae-gene.tsv``).  There should be exactly one such
    file; if several exist the first one (alphabetically) is returned.

    Returns
    -------
    str or None
        Absolute path to the meta TSV, or *None* if not found.
    """
    candidates = sorted(
        p for p in os.listdir(db_dir)
        if p.endswith('-gene.tsv') and not p.startswith('.')
    )
    if candidates:
        return os.path.join(db_dir, candidates[0])
    return None


def _extract_trace_genomes(trace_result_tsv: str) -> Optional[str]:
    """Read Trace.result.tsv and return semicolon-joined genome IDs
    from the row(s) where FinalResult == 'Yes'.

    Returns
    -------
    str or None
        e.g. ``"GCF_900451285.1;GCF_019796015.1"`` or *None* if no
        FinalResult=Yes row exists.
    """
    import csv as _csv
    genomes_parts = []
    if not os.path.isfile(trace_result_tsv):
        return None
    with open(trace_result_tsv, 'r', newline='') as fh:
        reader = _csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row.get('FinalResult', '') == 'Yes':
                g = row.get('Genomes', '').strip()
                if g:
                    genomes_parts.append(g)
    if not genomes_parts:
        return None
    # Genomes column may already be semicolon-separated; combine all rows
    return ';'.join(genomes_parts)


def run_html_viewer_for_samples(succeeded, outdir, db_dir, taxname, logger,
                                meta_tsv=None):
    """Generate an interactive HTML viewer for each successfully traced sample.

    This is **Step 4a** in the pipeline — executed after per-sample tracing
    (Step 3) completes and before the final summary.

    For each sample whose ``TraceResult/Trace.result.tsv`` contains at least
    one ``FinalResult=Yes`` row, we invoke
    ``Generate_trace_viewer_v2.generate_viewer_from_data()`` to produce:

        ``{outdir}/{sample}/TraceResult/mTracePatho_viewer.html``
        ``{outdir}/{sample}/TraceResult/mTracePatho_viewer.data.json``

    Parameters
    ----------
    succeeded : list of (sample_name, filt_node_path)
    outdir : str
    db_dir : str
    taxname : str
        Species Latin name (e.g. "Klebsiella pneumoniae").
    logger : PipelineLogger
    meta_tsv : str or None
        Path to *-gene.tsv metadata file. If None, auto-detected from db_dir.
    """
    logger.step_start("Step 4a: HTML Viewer generation")

    # ── Locate database files ─────────────────────────────────────────────
    tree_nwk = os.path.join(os.path.abspath(db_dir), 'Jaccard.fastme.input.nwk')
    if not os.path.isfile(tree_nwk):
        logger.warning(f"  [SKIP] Jaccard.fastme.input.nwk not found in {db_dir}")
        logger.step_end("Step 4a: HTML Viewer generation", success=False)
        return

    # ── Locate meta file: --meta > auto-detect in db_dir ──────────────────
    if meta_tsv:
        meta_tsv = os.path.abspath(meta_tsv)
        if not os.path.isfile(meta_tsv):
            logger.warning(f"  [SKIP] --meta file not found: {meta_tsv}")
            logger.step_end("Step 4a: HTML Viewer generation", success=False)
            return
    else:
        meta_tsv = _find_meta_tsv(os.path.abspath(db_dir))
        if not meta_tsv:
            logger.warning(f"  [SKIP] No *-gene.tsv metadata file found in {db_dir}. "
                           f"Please provide --meta explicitly.")
            logger.step_end("Step 4a: HTML Viewer generation", success=False)
            return

    logger.info(f"  Tree:   {tree_nwk}")
    logger.info(f"  Meta:   {meta_tsv}")
    logger.info(f"  Taxname: {taxname}")

    # ── Import the viewer generator ───────────────────────────────────────
    try:
        # Try importing as a module first (if installed or on PYTHONPATH)
        from Generate_trace_viewer_v2 import generate_viewer_from_data
    except ImportError:
        # Fall back: look for the script next to this file
        _this_dir = os.path.dirname(os.path.abspath(__file__))
        _viewer_script = os.path.join(_this_dir, 'Generate_trace_viewer.v2.py')
        if not os.path.isfile(_viewer_script):
            # Also try the name variant with dots
            _viewer_script = os.path.join(_this_dir, 'Generate_trace_viewer.v2.py')
        if os.path.isfile(_viewer_script):
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                'Generate_trace_viewer_v2', _viewer_script)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            generate_viewer_from_data = mod.generate_viewer_from_data
        else:
            logger.warning(
                "  [SKIP] Generate_trace_viewer_v2.py not found. "
                "Place it in the same directory as mTracePathoRun.py or on PYTHONPATH.")
            logger.step_end("Step 4a: HTML Viewer generation", success=False)
            return

    # ── Generate viewer for each sample ───────────────────────────────────
    n_ok, n_skip = 0, 0
    for sample, _filt_node in succeeded:
        trace_dir = os.path.join(outdir, sample, 'TraceResult')
        trace_result_tsv = os.path.join(trace_dir, 'Trace.result.tsv')

        trace_genomes_str = _extract_trace_genomes(trace_result_tsv)
        if not trace_genomes_str:
            logger.info(f"  [{sample}] No FinalResult=Yes genomes — viewer skipped")
            n_skip += 1
            continue

        # Also write the trace_result.tsv expected by the viewer
        # (single line of semicolon-separated genome IDs)
        trace_input_file = os.path.join(trace_dir, 'trace_result.tsv')
        with open(trace_input_file, 'w') as _tf:
            _tf.write(trace_genomes_str + '\n')

        output_html = os.path.join(trace_dir, 'mTracePatho_viewer.html')
        ok = generate_viewer_from_data(
            tree_path=tree_nwk,
            meta_path=meta_tsv,
            trace_genomes_str=trace_genomes_str,
            taxname=taxname,
            output_html=output_html,
            dist_unit='Jaccard distance',
            log_fn=lambda msg, s=sample: logger.info(f"  [{s}] {msg}"),
        )
        if ok:
            n_ok += 1
        else:
            n_skip += 1

    logger.info(f"  Viewer generated: {n_ok} sample(s), skipped: {n_skip}")
    logger.step_end("Step 4a: HTML Viewer generation", success=True,
                    notes=f"ok={n_ok} skip={n_skip}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: BATCH SUMMARY STATISTICS & VISUALIZATION (CNS-level output)
# ═══════════════════════════════════════════════════════════════════════════════

def run_batch_statistics(summary_tsv: str, output_dir: str, logger: 'PipelineLogger'):
    """
    Generate statistical tables and publication-quality figures from
    result_comb/batch_summary.tsv (column schema mirrors Trace.result.tsv,
    FinalResult=Yes rows only).

    Key algorithm-performance columns used:
        Sample, TraceLevel, ConfidenceScore, ConfidenceLevel,
        LCA_Round, GenomesNumber, TotalSNP, TracingSnpUsed,
        FractionMappedReads, ZeroRatio, SnpRank

    Outputs written to {output_dir}/statistics/:
        Table1_overview.tsv
        Table2_by_trace_level.tsv
        Table3_by_confidence_level.tsv
        Fig1_confidence_dist.pdf/jpg
        Fig2_trace_level_bar.pdf/jpg
        Fig3_confidence_by_level.pdf/jpg
        Fig4_score_vs_reads.pdf/jpg
        Fig5_heatmap_level_conf.pdf/jpg
        Fig6_iterations_dist.pdf/jpg
    """
    import csv as _csv

    logger.step_start("Step 4: Batch statistics & visualization")
    stat_dir = os.path.join(output_dir, 'statistics')
    os.makedirs(stat_dir, exist_ok=True)

    # ── Load data ────────────────────────────────────────────────────────────
    all_rows = []
    with open(summary_tsv, 'r', newline='') as _f:
        reader = _csv.DictReader(_f, delimiter='\t')
        for row in reader:
            all_rows.append(row)

    if not all_rows:
        logger.warning("  No rows in batch_summary.tsv — skipping statistics")
        logger.step_end("Step 4", success=False, module='step4', submodule='statistics')
        return

    # batch_summary already contains only FinalResult=Yes rows; accept all.
    # Optionally still filter to be safe.
    rows = [r for r in all_rows if r.get('FinalResult', 'Yes') in ('Yes', 'N/A')]
    if not rows:
        rows = all_rows
        logger.warning("  No FinalResult='Yes' rows found — using all rows for stats")

    logger.info(f"  Rows for statistics: {len(rows)} (from {len(all_rows)} total)")

    # ── Helper parsers ───────────────────────────────────────────────────────
    def _f(v, default=float('nan')):
        try: return float(v)
        except (ValueError, TypeError): return default

    def _i(v, default=0):
        try: return int(v)
        except (ValueError, TypeError): return default

    # ── Parse algorithm-performance columns ──────────────────────────────────
    samples        = [r.get('Sample', '') for r in rows]
    conf_scores    = np.array([_f(r.get('ConfidenceScore', 'nan')) for r in rows])
    conf_levels    = [r.get('ConfidenceLevel', 'N/A') for r in rows]
    trace_levels   = [r.get('TraceLevel', 'N/A') for r in rows]
    lca_rounds     = np.array([_i(r.get('LCA_Round', '0')) for r in rows])
    genomes_num    = np.array([_i(r.get('GenomesNumber', '0')) for r in rows])
    total_snps     = np.array([_i(r.get('TotalSNP', '0')) for r in rows])
    tracing_snps   = np.array([_i(r.get('TracingSnpUsed', '0')) for r in rows])
    frac_mapped    = np.array([_f(r.get('FractionMappedReads', 'nan')) for r in rows])
    zero_ratio     = np.array([_f(r.get('ZeroRatio', 'nan')) for r in rows])
    snp_rank       = np.array([_i(r.get('SnpRank', '0')) for r in rows])

    valid_mask = ~np.isnan(conf_scores)
    cs_valid   = conf_scores[valid_mask]
    n_total    = len(rows)
    n_valid    = int(valid_mask.sum())

    # ── CNS-level style ──────────────────────────────────────────────────────
    _STYLE = {
        'font.family': 'sans-serif', 'font.size': 9,
        'axes.titlesize': 10, 'axes.labelsize': 9,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 8, 'figure.dpi': 300, 'savefig.dpi': 300,
        'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
    }
    plt.rcParams.update(_STYLE)

    _CONF_COL  = {'High': '#2166ac', 'Medium': '#74add1', 'Low': '#d73027', 'N/A': '#969696'}
    _LEVEL_COL = {'Genome': '#1a9641', 'Node': '#fdae61', 'Error': '#d73027', 'N/A': '#969696'}

    def _save(fig, stem):
        for ext, kw in [('.pdf', {'bbox_inches': 'tight'}),
                        ('.jpg', {'bbox_inches': 'tight', 'dpi': 300})]:
            fig.savefig(os.path.join(stat_dir, stem + ext), **kw)
        plt.close(fig)

    def _pct(n, d): return f"{100*n/d:.1f}" if d > 0 else "N/A"
    def _fmt(v): return f"{v:.4f}" if not np.isnan(v) else "N/A"

    # ════════════════════════════════════════════════════════════════════════
    # TABLES
    # ════════════════════════════════════════════════════════════════════════

    # Table 1: Overview
    tl_counts = {k: trace_levels.count(k) for k in sorted(set(trace_levels))}
    cl_counts = {k: conf_levels.count(k) for k in ['High', 'Medium', 'Low']}
    mean_snp  = np.nanmean(total_snps) if n_total > 0 else float('nan')
    mean_tracing_snp = np.nanmean(tracing_snps) if n_total > 0 else float('nan')
    mean_frac_mapped = np.nanmean(frac_mapped[~np.isnan(frac_mapped)]) if n_total > 0 else float('nan')

    t1 = [
        ('total_samples',              str(n_total),
         'Total samples analysed (FinalResult=Yes)'),
        ('samples_with_score',         str(n_valid),
         'Samples with valid ConfidenceScore'),
        ('mean_confidence_score',      _fmt(np.nanmean(cs_valid))   if n_valid>0 else 'N/A',
         'Mean ConfidenceScore (0–1)'),
        ('median_confidence_score',    _fmt(np.nanmedian(cs_valid)) if n_valid>0 else 'N/A',
         'Median ConfidenceScore'),
        ('sd_confidence_score',        _fmt(np.nanstd(cs_valid,ddof=1)) if n_valid>1 else 'N/A',
         'SD of ConfidenceScore'),
        ('min_confidence_score',       _fmt(np.nanmin(cs_valid))   if n_valid>0 else 'N/A',
         'Minimum ConfidenceScore'),
        ('max_confidence_score',       _fmt(np.nanmax(cs_valid))   if n_valid>0 else 'N/A',
         'Maximum ConfidenceScore'),
        ('pct_high_confidence',        _pct(cl_counts.get('High',0), n_total)+'%',
         'Percent High confidence (≥0.8)'),
        ('pct_medium_confidence',      _pct(cl_counts.get('Medium',0), n_total)+'%',
         'Percent Medium confidence (0.5–0.8)'),
        ('pct_low_confidence',         _pct(cl_counts.get('Low',0), n_total)+'%',
         'Percent Low confidence (<0.5)'),
        ('pct_genome_level',           _pct(tl_counts.get('Genome',0), n_total)+'%',
         'Percent Genome-level resolution'),
        ('pct_node_level',             _pct(tl_counts.get('Node',0), n_total)+'%',
         'Percent Node-level resolution'),
        ('mean_lca_rounds',            f"{np.nanmean(lca_rounds):.2f}",
         'Mean LCA_Round (iteration count)'),
        ('mean_genomes_number',        f"{np.nanmean(genomes_num):.2f}",
         'Mean GenomesNumber (final candidate count)'),
        ('mean_total_snp',             f"{mean_snp:.1f}" if not np.isnan(mean_snp) else 'N/A',
         'Mean TotalSNP per sample'),
        ('mean_tracing_snp_used',      f"{mean_tracing_snp:.1f}" if not np.isnan(mean_tracing_snp) else 'N/A',
         'Mean TracingSnpUsed per sample'),
        ('mean_fraction_mapped_reads', f"{mean_frac_mapped:.4f}" if not np.isnan(mean_frac_mapped) else 'N/A',
         'Mean FractionMappedReads'),
    ]
    t1_path = os.path.join(stat_dir, 'Table1_overview.tsv')
    with open(t1_path, 'w', newline='') as _f:
        w = _csv.writer(_f, delimiter='\t')
        w.writerow(['metric', 'value', 'description'])
        w.writerows(t1)
    logger.info(f"  Table1: {t1_path}")

    # Table 2: By trace level
    t2 = []
    for lv in sorted(set(trace_levels)):
        idx = [i for i, t in enumerate(trace_levels) if t == lv]
        lv_cs  = conf_scores[idx]; lv_cs  = lv_cs[~np.isnan(lv_cs)]
        lv_it  = lca_rounds[idx]
        lv_snp = total_snps[idx]
        lv_tracing = tracing_snps[idx]
        lv_frac = frac_mapped[idx]; lv_frac = lv_frac[~np.isnan(lv_frac)]
        t2.append({
            'TraceLevel':              lv,
            'n_samples':               len(idx),
            'pct':                     _pct(len(idx), n_total)+'%',
            'mean_ConfidenceScore':    _fmt(np.nanmean(lv_cs))    if len(lv_cs)>0 else 'N/A',
            'median_ConfidenceScore':  _fmt(np.nanmedian(lv_cs))  if len(lv_cs)>0 else 'N/A',
            'sd_ConfidenceScore':      _fmt(np.nanstd(lv_cs, ddof=1)) if len(lv_cs)>1 else 'N/A',
            'mean_LCA_Round':          f"{np.nanmean(lv_it):.2f}",
            'mean_GenomesNumber':      f"{np.nanmean(genomes_num[idx]):.2f}",
            'mean_TotalSNP':           f"{np.nanmean(lv_snp):.1f}",
            'mean_TracingSnpUsed':     f"{np.nanmean(lv_tracing):.1f}",
            'mean_FractionMappedReads': _fmt(np.nanmean(lv_frac)) if len(lv_frac)>0 else 'N/A',
        })
    if t2:
        t2_path = os.path.join(stat_dir, 'Table2_by_trace_level.tsv')
        with open(t2_path, 'w', newline='') as _f:
            w = _csv.DictWriter(_f, delimiter='\t', fieldnames=list(t2[0].keys()))
            w.writeheader(); w.writerows(t2)
        logger.info(f"  Table2: {t2_path}")

    # Table 3: By confidence level
    t3 = []
    for cl in ['High', 'Medium', 'Low', 'N/A']:
        idx = [i for i, c in enumerate(conf_levels) if c == cl]
        if not idx:
            continue
        cl_cs  = conf_scores[idx]; cl_cs = cl_cs[~np.isnan(cl_cs)]
        cl_it  = lca_rounds[idx]
        cl_tl  = [trace_levels[i] for i in idx]
        cl_snp = total_snps[idx]
        cl_tracing = tracing_snps[idx]
        cl_frac = frac_mapped[idx]; cl_frac = cl_frac[~np.isnan(cl_frac)]
        cl_zr   = zero_ratio[idx];  cl_zr   = cl_zr[~np.isnan(cl_zr)]
        t3.append({
            'ConfidenceLevel':         cl,
            'n_samples':               len(idx),
            'pct':                     _pct(len(idx), n_total)+'%',
            'mean_ConfidenceScore':    _fmt(np.nanmean(cl_cs))    if len(cl_cs)>0 else 'N/A',
            'median_ConfidenceScore':  _fmt(np.nanmedian(cl_cs))  if len(cl_cs)>0 else 'N/A',
            'sd_ConfidenceScore':      _fmt(np.nanstd(cl_cs, ddof=1)) if len(cl_cs)>1 else 'N/A',
            'pct_Genome':              _pct(cl_tl.count('Genome'), len(idx))+'%',
            'pct_Node':                _pct(cl_tl.count('Node'),   len(idx))+'%',
            'mean_LCA_Round':          f"{np.nanmean(cl_it):.2f}",
            'mean_TotalSNP':           f"{np.nanmean(cl_snp):.1f}",
            'mean_TracingSnpUsed':     f"{np.nanmean(cl_tracing):.1f}",
            'mean_FractionMappedReads': _fmt(np.nanmean(cl_frac)) if len(cl_frac)>0 else 'N/A',
            'mean_ZeroRatio':          _fmt(np.nanmean(cl_zr))    if len(cl_zr)>0 else 'N/A',
        })
    if t3:
        t3_path = os.path.join(stat_dir, 'Table3_by_confidence_level.tsv')
        with open(t3_path, 'w', newline='') as _f:
            w = _csv.DictWriter(_f, delimiter='\t', fieldnames=list(t3[0].keys()))
            w.writeheader(); w.writerows(t3)
        logger.info(f"  Table3: {t3_path}")

    # ════════════════════════════════════════════════════════════════════════
    # FIGURES
    # ════════════════════════════════════════════════════════════════════════
    if n_valid == 0:
        logger.warning("  No valid confidence scores — skipping figures")
        logger.step_end("Step 4", success=True, module='step4', submodule='statistics')
        return

    # Fig 1: Confidence score distribution
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    valid_pairs = [(conf_scores[i], conf_levels[i]) for i in range(n_total) if not np.isnan(conf_scores[i])]
    for cl, col in _CONF_COL.items():
        vals = [v for v, c in valid_pairs if c == cl]
        if vals:
            ax.hist(vals, bins=20, range=(0, 1), alpha=0.75, color=col,
                    label=cl, edgecolor='white', linewidth=0.4)
    ax.axvline(np.nanmean(cs_valid), color='#333', lw=1.2, ls='--',
               label=f'Mean ({np.nanmean(cs_valid):.3f})', zorder=5)
    ax.set_xlabel('ConfidenceScore'); ax.set_ylabel('Number of Samples')
    ax.set_title('Confidence Score Distribution', fontweight='bold', pad=6)
    ax.set_xlim(0, 1); ax.legend(frameon=False, loc='upper left', ncol=2)
    plt.tight_layout(); _save(fig, 'Fig1_confidence_dist')
    logger.info("  Fig1: confidence distribution")

    # Fig 2: Trace level + confidence level bars (2-panel)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.0))
    tl_labs = sorted(set(trace_levels), key=lambda x: -trace_levels.count(x))
    bars1 = ax1.bar(tl_labs, [trace_levels.count(l) for l in tl_labs],
                    color=[_LEVEL_COL.get(l, '#969696') for l in tl_labs],
                    edgecolor='white', linewidth=0.6, width=0.55)
    for b in bars1:
        ax1.text(b.get_x()+b.get_width()/2, b.get_height()+0.2, str(int(b.get_height())),
                 ha='center', va='bottom', fontsize=7.5)
    ax1.set_ylabel('Count'); ax1.set_title('A  TraceLevel', fontweight='bold', pad=5, loc='left')
    cl_ord  = ['High', 'Medium', 'Low']
    cl_vals = [conf_levels.count(l) for l in cl_ord]
    bars2 = ax2.bar(cl_ord, cl_vals, color=[_CONF_COL[l] for l in cl_ord],
                    edgecolor='white', linewidth=0.6, width=0.55)
    for b in bars2:
        ax2.text(b.get_x()+b.get_width()/2, b.get_height()+0.2, str(int(b.get_height())),
                 ha='center', va='bottom', fontsize=7.5)
    ax2.set_ylabel('Count'); ax2.set_title('B  ConfidenceLevel', fontweight='bold', pad=5, loc='left')
    plt.tight_layout(w_pad=2.5); _save(fig, 'Fig2_trace_level_bar')
    logger.info("  Fig2: trace/confidence level bars")

    # Fig 3: Violin plot — ConfidenceScore by ConfidenceLevel
    cl_groups = ['High', 'Medium', 'Low']
    cl_data = [conf_scores[[i for i, c in enumerate(conf_levels) if c == cl and not np.isnan(conf_scores[i])]]
               for cl in cl_groups]
    cl_data_nz = [(d, cl) for d, cl in zip(cl_data, cl_groups) if len(d) > 0]
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    if cl_data_nz:
        vp = ax.violinplot([d for d, _ in cl_data_nz],
                           positions=range(1, len(cl_data_nz)+1),
                           showmedians=True, widths=0.65)
        for i, (pc, (_, cl)) in enumerate(zip(vp['bodies'], cl_data_nz)):
            pc.set_facecolor(_CONF_COL.get(cl, '#969696')); pc.set_alpha(0.7)
        for part in ('cmedians', 'cmins', 'cmaxes', 'cbars'):
            if part in vp:
                vp[part].set_edgecolor('#333'); vp[part].set_linewidth(1.0)
        for i, (d, cl) in enumerate(cl_data_nz):
            jit = np.random.default_rng(42).uniform(-0.12, 0.12, len(d))
            ax.scatter([i+1+j for j in jit], d, s=6, alpha=0.45,
                       color=_CONF_COL.get(cl, '#969696'), zorder=3)
        ax.set_xticks(range(1, len(cl_data_nz)+1))
        ax.set_xticklabels([cl for _, cl in cl_data_nz])
    ax.set_ylabel('ConfidenceScore'); ax.set_ylim(-0.05, 1.05)
    ax.set_title('ConfidenceScore by ConfidenceLevel', fontweight='bold', pad=6)
    ax.axhline(0.5, color='#aaa', lw=0.8, ls=':')
    plt.tight_layout(); _save(fig, 'Fig3_confidence_by_level')
    logger.info("  Fig3: violin plot")

    # Fig 4: ConfidenceScore vs TotalSNP (scatter)
    snp_avail = [(total_snps[i], conf_scores[i], conf_levels[i])
                 for i in range(n_total) if total_snps[i] > 0 and not np.isnan(conf_scores[i])]
    if len(snp_avail) >= 3:
        snp_x  = np.array([r[0] for r in snp_avail])
        snp_y  = np.array([r[1] for r in snp_avail])
        snp_cl = [r[2] for r in snp_avail]
        fig, ax = plt.subplots(figsize=(4.8, 3.5))
        for cl in ['Low', 'Medium', 'High']:
            mask = [c == cl for c in snp_cl]
            if any(mask):
                ax.scatter(snp_x[mask], snp_y[mask], s=18, alpha=0.7,
                           color=_CONF_COL.get(cl, '#969696'), label=cl,
                           zorder=3, edgecolors='none')
        log_x = np.log10(snp_x + 1)
        coef  = np.polyfit(log_x, snp_y, 1)
        x_lin = np.linspace(log_x.min(), log_x.max(), 200)
        ax.plot(10**x_lin, np.polyval(coef, x_lin), '--', color='#333',
                lw=1.2, alpha=0.8, label='Trend (log-linear)', zorder=4)
        ax.set_xscale('log')
        ax.set_xlabel('TotalSNP'); ax.set_ylabel('ConfidenceScore')
        ax.set_title('ConfidenceScore vs TotalSNP', fontweight='bold', pad=6)
        ax.set_ylim(-0.05, 1.05); ax.legend(frameon=False, loc='lower right')
        ax.axhline(0.5, color='#aaa', lw=0.8, ls=':')
        plt.tight_layout(); _save(fig, 'Fig4_score_vs_snp')
        logger.info("  Fig4: score vs TotalSNP")

    # Fig 5: Heatmap TraceLevel × ConfidenceLevel
    tl_uniq = sorted(set(trace_levels) - {'N/A', 'Error', 'Unknown'})
    if tl_uniq:
        cl_uniq = ['High','Medium','Low']
        matrix  = np.zeros((len(cl_uniq), len(tl_uniq)), dtype=int)
        for tl, cl in zip(trace_levels, conf_levels):
            if tl in tl_uniq and cl in cl_uniq:
                matrix[cl_uniq.index(cl), tl_uniq.index(tl)] += 1
        fig, ax = plt.subplots(figsize=(3.5+0.6*len(tl_uniq), 2.8))
        im = ax.imshow(matrix, cmap='Blues', aspect='auto', vmin=0)
        ax.set_xticks(range(len(tl_uniq))); ax.set_xticklabels(tl_uniq)
        ax.set_yticks(range(len(cl_uniq))); ax.set_yticklabels(cl_uniq)
        ax.set_xlabel('Trace Level'); ax.set_ylabel('Confidence Level')
        ax.set_title('Samples: Trace Level × Confidence Level', fontweight='bold', pad=6)
        fig.colorbar(im, ax=ax, shrink=0.8, label='Count')
        for ci in range(len(cl_uniq)):
            for ti in range(len(tl_uniq)):
                v = matrix[ci,ti]
                ax.text(ti, ci, str(v), ha='center', va='center',
                        fontsize=8, fontweight='bold',
                        color='white' if v > matrix.max()*0.6 else '#222')
        plt.tight_layout(); _save(fig, 'Fig5_heatmap_level_conf')
        logger.info("  Fig5: heatmap")

    # Fig 6: LCA rounds distribution
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    rvals, rcounts = np.unique(lca_rounds, return_counts=True)
    ax.bar(rvals, rcounts, color='#4393c3', edgecolor='white', linewidth=0.5, width=0.6)
    for v,c in zip(rvals, rcounts):
        ax.text(v, c+0.2, str(c), ha='center', va='bottom', fontsize=7.5)
    ax.set_xlabel('Number of LCA Rounds'); ax.set_ylabel('Number of Samples')
    ax.set_title('LCA Iteration Count Distribution', fontweight='bold', pad=6)
    plt.tight_layout(); _save(fig, 'Fig6_iterations_dist')
    logger.info("  Fig6: LCA rounds distribution")

    # Fig 7: SNP count distribution (new column available)
    if np.any(total_snps > 0):
        fig, ax = plt.subplots(figsize=(4.5, 3.2))
        ax.hist(total_snps[total_snps > 0], bins=30, color='#5aae61',
                edgecolor='white', linewidth=0.4, log=True)
        ax.axvline(np.nanmedian(total_snps[total_snps > 0]), color='#333', lw=1.2,
                   ls='--', label=f'Median ({np.nanmedian(total_snps[total_snps>0]):.0f})')
        ax.set_xlabel('Total SNP Count per Sample'); ax.set_ylabel('Number of Samples (log)')
        ax.set_title('SNP Count Distribution', fontweight='bold', pad=6)
        ax.legend(frameon=False)
        plt.tight_layout(); _save(fig, 'Fig7_snp_count_dist')
        logger.info("  Fig7: SNP count distribution")

    logger.step_end("Step 4: Batch statistics & visualization", success=True,
                    module='step4', submodule='statistics',
                    notes=f"n={n_total} final_rows={len(rows)}")
    logger.info(f"  Statistics output: {stat_dir}")


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def load_fastq_list(fastq_list_path: str) -> List[Tuple[str, str]]:
    """Load batch sample list TSV: col1=sample_name, col2=fastq_path."""
    entries = []
    with open(fastq_list_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            entries.append((parts[0].strip(), parts[1].strip()))
    return entries


def check_dependencies(config: PipelineConfig, logger: PipelineLogger) -> bool:
    """Check required Python packages and external tools."""
    ok = True
    required_py = {
        'pysam':      'pip install pysam',
        'numpy':      'pip install numpy',
        'scipy':      'pip install scipy',
        'matplotlib': 'pip install matplotlib',
        'psutil':     'pip install psutil',
    }
    for pkg, install_cmd in required_py.items():
        try:
            __import__(pkg)
            logger.debug(f"  [OK] Python: {pkg}")
        except ImportError:
            logger.error(f"  [MISSING] Python: {pkg}  →  {install_cmd}")
            ok = False

    required_tools = ['fastp', 'kraken2', 'bwa', 'samtools', 'bcftools', 'bamdst', 'kz']
    for tool in required_tools:
        path = config.get_tool(tool)
        if (os.path.isfile(path) and os.access(path, os.X_OK)) or shutil.which(tool):
            logger.debug(f"  [OK] Tool: {tool}")
        else:
            logger.warning(f"  [WARN] Tool: {tool} not found at '{path}'")

    for db_name in ['kraken2_pathogen_db', 'kraken2_plsdb', 'kraken2_homo_db']:
        db_path = config.get_db(db_name)
        if os.path.isdir(db_path):
            logger.debug(f"  [OK] Database: {db_name}")
        else:
            logger.warning(f"  [WARN] Database: {db_name} not found at '{db_path}'")

    return ok


def verify_db_dir(db_dir: str, logger: PipelineLogger) -> bool:
    """Check all required database files are present."""
    ok = True
    for fname in DB_FILES_REQUIRED:
        fpath = os.path.join(db_dir, fname)
        if os.path.isfile(fpath):
            logger.debug(f"  [OK] DB file: {fname}")
        else:
            logger.error(f"  [MISSING] DB file: {fpath}")
            ok = False
    return ok


# ═══════════════════════════════════════════════════════════════════════════════
# ARGUMENT PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='mTracePathoRun.py',
        description=(
            "mTracePathoRun v{ver} — Pathogen source tracing from FASTQ inputs.\n\n"
            "Runs a full reads→tracing pipeline:\n"
            "  Step 1: Reads pre-processing (kraken2 filter → pullfq → kz)\n"
            "  Step 2: Alignment + variant calling + SNP→node mapping\n"
            "  Step 3: Source tracing (mTracePatho2 LCA algorithm)\n\n"
            "Requires a reference database built with mTracePathoMakeDB.py."
        ).format(ver=VERSION),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Single sample:\n"
            "  mTracePathoRun.py --fastq sample.fq.gz \\\n"
            "    --db_dir dbdir/ --taxid 470 --outdir results/\n\n"
            "  # Batch with 16 threads:\n"
            "  mTracePathoRun.py --fastq_list samples.tsv \\\n"
            "    --db_dir dbdir/ --taxid 470 \\\n"
            "    --outdir results/ --threads 16\n\n"
            "  # Multiple taxids:\n"
            "  mTracePathoRun.py --fastq sample.fq.gz \\\n"
            "    --db_dir dbdir/ \\\n"
            "    --taxid '470-1234-5678' --outdir results/\n\n"
            "  # Skip Step 1 (reads already pre-filtered):\n"
            "  mTracePathoRun.py --fastq sample.fq.gz \\\n"
            "    --db_dir dbdir/ --taxid 470 \\\n"
            "    --outdir results/ --skip_steps step1\n\n"
            "  # SGE cluster mode:\n"
            "  mTracePathoRun.py --fastq_list samples.tsv \\\n"
            "    --db_dir dbdir/ --taxid 470 --outdir results/ \\\n"
            "    --use_sge --sge_nodes 'R7525|big2|big3|big4' \\\n"
            "    --sge_mem 1G --sge_slots 1 --sge_queue all.q\n\n"
            "  # Print config template:\n"
            "  mTracePathoRun.py --print_config\n"
        )
    )

    # ── Input (mutually exclusive) ─────────────────────────────────────────────
    inp = parser.add_argument_group(
        'Input',
        '--fastq and --fastq_list are mutually exclusive. '
        '--result_dir can be combined with --fastq_list (for sample name validation).')
    inp_ex = inp.add_mutually_exclusive_group()
    inp_ex.add_argument('--fastq', metavar='FQ',
                        help='Single FASTQ file (.fq or .fq.gz)')
    inp_ex.add_argument('--fastq_list', metavar='TSV',
                        help='Batch mode TSV: col1=sample_name, col2=fastq_path')
    # --result_dir is intentionally NOT in the mutually_exclusive_group so that
    # it can be freely combined with --fastq_list (which then acts only as a
    # sample-name validation list, not as a source of FASTQ paths).
    inp.add_argument('--result_dir', metavar='DIR', default=None,
                     help='Re-use an existing mTracePathoRun output directory. '
                          'Reads result/{sample}/{sample}.filt.node directly '
                          '(skips Steps 1+2). Combine with --fastq_list to validate '
                          'expected sample names. --ref and --taxid are not required '
                          'in this mode.')

    # ── Required ───────────────────────────────────────────────────────────────
    req = parser.add_argument_group('Required arguments')
    req.add_argument('--db_dir', metavar='DIR', default=None,
                     help='Database directory (mTracePathoMakeDB.py output)')
    req.add_argument('--taxid', metavar='TAXID', default=None,
                     help='Target taxid(s); multiple taxids joined by "-" '
                          '(e.g. "470" or "470-1234-5678")')
    req.add_argument('--outdir', metavar='DIR', default=None,
                     help='Output directory (created if not exists)')

    # ── Pipeline options ───────────────────────────────────────────────────────
    pipe = parser.add_argument_group('Pipeline options')
    pipe.add_argument('--threads', type=int, default=4, metavar='N',
                      help='Parallel workers (local) and tool threads per job (default: 4)')
    pipe.add_argument('--filter_identity', type=float, default=0.95, metavar='F',
                      help='BAM filter: minimum alignment identity fraction (default: 0.95)')
    pipe.add_argument('--filter_coverage', type=float, default=0.95, metavar='F',
                      help='BAM filter: minimum alignment coverage fraction (default: 0.95)')
    pipe.add_argument('--skip_steps', metavar='STEPS',
                      help='Comma-separated steps to skip: step1,step2 '
                           '(step3 cannot be skipped independently)')
    pipe.add_argument('--skip_fastp', action='store_true',
                      help='Skip fastp QC step (default: run fastp)')
    pipe.add_argument('--skip_rmhost', action='store_true',
                      help='Skip kraken2 host read removal step (default: run host removal)')
    pipe.add_argument('--skip_pathogen', action='store_true',
                      help='Skip pathogen Kraken2 classification and target-taxid read '
                           'extraction; pass PLSDB-unclassified reads directly to kz/BWA')
    pipe.add_argument('--skip_step1', action='store_true',
                      help='Skip entire Step 1 preprocessing (fastp QC, host removal, '
                           'kraken2 PLSDB, kraken2 pathogen, pullfq, kz filter) and '
                           'jump directly to Step 2 (bwa mem alignment). '
                           'Equivalent to --skip_steps step1')
    pipe.add_argument('--force', action='store_true',
                      help='Ignore .done checkpoints and re-run all steps')
    pipe.add_argument('--config', metavar='INI',
                      help=f'Path to tools.ini (default search: '
                           f'{", ".join(CONFIG_SEARCH_PATHS)})')

    # ── Tracing algorithm ──────────────────────────────────────────────────────
    tr = parser.add_argument_group(
        'Tracing algorithm parameters',
        'Passed directly to the mTracePathoAlgorithm tracing core.')
    tr.add_argument('--jaccard_percent', type=float, default=0.8, metavar='F',
                    help='Jaccard similarity threshold for candidate selection (default: 0.8)')
    tr.add_argument('--max_iterations', type=int, default=10, metavar='N',
                    help='Maximum LCA iteration rounds (default: 10)')
    tr.add_argument('--use_containment', action='store_true', default=True,
                    help='Enable containment-based genome pre-selection (default: enabled; pass --no_use_containment to disable)')
    tr.add_argument('--no_use_containment', dest='use_containment', action='store_false',
                    help='Disable containment-based genome pre-selection')
    tr.add_argument('--min_containment', type=float, default=0.0, metavar='F',
                    help='Minimum containment score threshold (default: 0.0)')
    tr.add_argument('--top_n_cutoff', type=int, default=99999, metavar='N',
                    help='Top N distinct containment score levels to select (default: 99999 = all levels)')
    tr.add_argument('--min_candidates', type=int, default=3, metavar='N',
                    help='Minimum candidate genomes to guarantee (default: 3)')
    tr.add_argument('--phylo_gap_ratio', type=float, default=999999999999, metavar='F',
                    help='Phylo pre-filter activation: top1_cont / level2_cont must be >= this '
                         'value to trigger pre-filtering (default: 1.4). '
                         'Set >99 to effectively disable phylo pre-filter.')
    tr.add_argument('--phylo_share_thresh', type=float, default=0.7, metavar='F',
                    help='Phylo pre-filter removal threshold: candidates sharing < this '
                         'fraction of top-1 lineage depth are removed as distant relatives '
                         '(default: 0.7).')
    tr.add_argument('--scoring_mode', type=str, default='containment_weighted',
                    choices=['containment', 'containment_weighted', 'jaccard', 'cosine',
                             'cosine_weighted', 'jaccard_weighted'],
                    help='Scoring mode: containment (original), '
                         'containment_weighted (default, containment × genome-size penalty), jaccard (symmetric), '
                         'cosine (cosine similarity), cosine_weighted (cosine × genome-size penalty), '
                         'jaccard_weighted (jaccard × genome-size penalty). '
                         'Use containment_weighted/jaccard/cosine_weighted to reduce large-genome bias.')
    tr.add_argument('--cont_alpha', type=float, default=0.5, metavar='F',
                    help='Exponent for weighted scoring: score = base_metric × '
                         '(ref_size/GenomeSNPs)^alpha. '
                         '0=no penalty, 0.5=sqrt (default), 1=linear. '
                         'Used when --scoring_mode is containment_weighted/cosine_weighted/jaccard_weighted.')
    tr.add_argument('--cont_ref_size', type=float, default=None, metavar='F',
                    help='Reference genome size for weighted scoring. '
                         'Auto-computed as median(GenomeSNPs) when not specified.')
    tr.add_argument('--score_cutoff_ratio', type=float, default=0.5, metavar='F',
                    help='Hard cutoff: all candidates must have score >= top1_score × ratio. '
                         'Candidates below this threshold are removed. If only 1 genome survives, '
                         'it is returned directly as the tracing result. (default: 0.5)')
    tr.add_argument('--min_unique_snp_n_cutoff', type=int, default=1, metavar='N',
                    help='Remove candidate genomes whose terminal-node unique SNPs hit <N times in sample (default: 1; 0=disable)')
    tr.add_argument('--min_snp', type=int, default=5, metavar='N',
                    help='Minimum non-root SNPs required; results with fewer are marked SNP<=N (default: 5)')
    tr.add_argument('--no_visualization', action='store_true',
                    help='Skip visualization outputs (faster, less disk)')
    tr.add_argument('--keep_temp', action='store_true',
                    help='[DEPRECATED] Ignored. Intermediate files are kept by default; '
                         'use --clean_temp to delete them.')
    tr.add_argument('--clean_temp', action='store_true',
                    help='Delete Step 1+2 intermediate files (*.qc.fq.gz, *.filter.fq, '
                         '*.kz.fq.gz, *.bam, *.depth, conf.txt, classified.fq, etc.) '
                         'after each sample completes tracing (saves disk space)')

    # ── HTML Viewer options ───────────────────────────────────────────────────
    viewer = parser.add_argument_group(
        'HTML Viewer options',
        'Generate an interactive HTML phylogenetic tree viewer for each sample '
        'tracing result. Requires Generate_trace_viewer_v2.py in the same directory '
        'or on PYTHONPATH.')
    viewer.add_argument('--taxname', metavar='NAME', default=None,
                        help='Species Latin name for the HTML viewer title, '
                             'e.g. "Klebsiella pneumoniae". Underscores are '
                             'automatically converted to spaces. If not provided, '
                             'the viewer step is skipped.')
    viewer.add_argument('--meta', metavar='TSV', default=None,
                        help='Genome metadata TSV file for the HTML viewer '
                             '(e.g. Klebsiella_pneumoniae-gene.tsv). '
                             'Each species has its own *-gene.tsv with different '
                             'columns; the viewer auto-classifies typing vs gene '
                             'columns. If not provided, the program searches '
                             '--db_dir for a *-gene.tsv file automatically.')
    viewer.add_argument('--no_html_viewer', action='store_true',
                        help='Explicitly skip HTML viewer generation even if '
                             '--taxname is provided')

    # ── SGE options ────────────────────────────────────────────────────────────
    sge = parser.add_argument_group(
        'SGE options (Steps 1+2 only)',
        'Submit per-sample preprocessing+variant jobs to SGE cluster via qsub. '
        'Disabled by default; local Pool is used instead.')
    sge.add_argument('--use_sge', action='store_true',
                     help='Enable SGE mode (default: local multiprocessing)')
    sge.add_argument('--sge_qsubname', metavar='NAME', default='mTracePatho',
                     help='Base name for SGE jobs (default: mTracePatho)')
    sge.add_argument('--sge_nodes', metavar='NODES', default='',
                     help='Pipe-separated target nodes: "R7525|big2|big3|big4"  '
                          '(default: any node)')
    sge.add_argument('--sge_mem', metavar='MEM', default='1G',
                     help='Memory per job, SGE vf= value (default: 1G)')
    sge.add_argument('--sge_slots', type=int, metavar='N', default=1,
                     help='CPU slots per job, SGE num_proc= (default: 1)')
    sge.add_argument('--sge_threads', type=int, metavar='N', default=1,
                     help='Tool threads inside each SGE job (default: 1)')
    sge.add_argument('--sge_queue', metavar='QUEUE', default='all.q',
                     help='SGE queue name (default: all.q)')
    sge.add_argument('--sge_max_jobs', type=int, metavar='N', default=500,
                     help='Maximum concurrent SGE jobs (default: 500)')
    sge.add_argument('--sge_poll_sec', type=int, metavar='SEC', default=30,
                     help='qstat polling interval in seconds (default: 30)')
    sge.add_argument('--sge_extra', metavar='OPTS', default='',
                     help='Extra qsub options appended verbatim (default: none)')

    # ── Utility ────────────────────────────────────────────────────────────────
    util = parser.add_argument_group('Utility')
    util.add_argument('--print_config', action='store_true',
                      help='Print tools.ini template and exit')

    return parser.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY PRINTER
# ═══════════════════════════════════════════════════════════════════════════════

def print_summary(sample_list, succeeded, failed, batch_results,
                  outdir, elapsed, logger):
    sep = '=' * 70
    logger.info(sep)
    logger.info("  mTracePathoRun — Pipeline Complete")
    logger.info(sep)
    logger.info(f"  Total samples  : {len(sample_list)}")
    logger.info(f"  Step1+2 OK     : {len(succeeded)}")
    logger.info(f"  Step1+2 Failed : {len(failed)}")
    if failed:
        logger.warning("  FAILED SAMPLES (Steps 1+2):")
        for s, reason in failed:
            logger.warning(f"    - {s}: {reason[:80]}")

    if batch_results:
        high  = sum(1 for r in batch_results if r.get('confidence_level') == 'High')
        med   = sum(1 for r in batch_results if r.get('confidence_level') == 'Medium')
        low   = sum(1 for r in batch_results if r.get('confidence_level') == 'Low')
        g_lvl = sum(1 for r in batch_results if r.get('trace_level') == 'Genome')
        n_lvl = sum(1 for r in batch_results if r.get('trace_level') == 'Node')
        logger.info(f"  Tracing confidence: High={high} Medium={med} Low={low}")
        logger.info(f"  Trace level:        Genome={g_lvl} Node={n_lvl}")

    m, s = divmod(elapsed, 60)
    h, m = divmod(m, 60)
    if h:
        t_str = f"{h:.0f}h {m:.0f}min {s:.0f}s"
    elif m:
        t_str = f"{m:.0f}min {s:.0f}s"
    else:
        t_str = f"{s:.1f}s"
    logger.info(f"  Total wall time: {t_str}")
    logger.info(f"  Results in: {os.path.join(outdir, 'results/')}")
    logger.info(sep)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()

    # ── --print_config ─────────────────────────────────────────────────────────
    if args.print_config:
        print(CONFIG_TEMPLATE)
        sys.exit(0)

    # ── Required arg validation ────────────────────────────────────────────────
    missing = [f"--{a}" for a in ("db_dir", "outdir")
               if not getattr(args, a, None)]
    if missing:
        print(f"Error: the following arguments are required: {chr(39)}{chr(39).join(missing)}{chr(39)}",
              file=sys.stderr)
        sys.exit(1)
    if not args.fastq and not args.fastq_list and not args.result_dir:
        print("Error: one of --fastq, --fastq_list, or --result_dir is required.", file=sys.stderr)
        sys.exit(1)
    # ref + taxid are not needed when reusing existing .filt.node files
    _reuse_mode = bool(args.result_dir)
    if not _reuse_mode:
        missing_rt = [f"--{a}" for a in ("taxid",) if not getattr(args, a, None)]
        if missing_rt:
            print(f"Error: {chr(39)}{chr(39).join(missing_rt)}{chr(39)} required when not using --result_dir.",
                  file=sys.stderr)
            sys.exit(1)

    # ── Directory layout ───────────────────────────────────────────────────────
    outdir  = os.path.abspath(args.outdir)
    log_dir = os.path.join(outdir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    logger = PipelineLogger(
        log_file=os.path.join(log_dir, 'pipeline.log'),
        tsv_file=os.path.join(log_dir, 'time.RSS.TSV'),
    )

    logger.info("=" * 70)
    logger.info(f"  mTracePathoRun v{VERSION} — Pathogen Source Tracing from FASTQ")
    logger.info("=" * 70)
    logger.info(f"  Command: {' '.join(sys.argv)}")
    logger.info(f"  outdir : {outdir}")
    logger.info(f"  db_dir : {os.path.abspath(args.db_dir)}")
    if not _reuse_mode:
        logger.info(f"  ref    : {os.path.join(os.path.abspath(args.db_dir), 'ref.fa')}")
        logger.info(f"  taxid  : {args.taxid}")
    else:
        logger.info(f"  result_dir (reuse): {os.path.abspath(args.result_dir)}")
    logger.info(f"  threads: {args.threads}")

    # ── Load config ────────────────────────────────────────────────────────────
    config = PipelineConfig(config_path=args.config, logger=logger)
    logger.info(f"  Config source: {config._source}")

    # ── Dependency and database check ──────────────────────────────────────────
    if not _reuse_mode:
        logger.info("Checking dependencies...")
        check_dependencies(config, logger)
    logger.info(f"Checking database files in {args.db_dir} ...")
    if not verify_db_dir(args.db_dir, logger):
        logger.error("Missing required database files. Exiting.")
        sys.exit(1)

    # ── ref.bed: auto-generate from FASTA if not already present ─────────────
    if not _reuse_mode:
        ref_fa  = os.path.join(os.path.abspath(args.db_dir), 'ref.fa')
        if not os.path.isfile(ref_fa):
            logger.error(f"Reference FASTA not found: {ref_fa}  (expected ref.fa inside --db_dir)")
            sys.exit(1)
        ref_bed = os.path.splitext(ref_fa)[0] + '.bed'
        if os.path.isfile(ref_bed):
            logger.info(f"  ref.bed found: {ref_bed}")
        else:
            logger.info(f"  ref.bed not found at {ref_bed} — generating from {ref_fa} ...")
            try:
                make_bed(ref_fa, ref_bed)
                logger.info(f"  ref.bed generated: {ref_bed}")
            except Exception as e:
                logger.error(f"  Failed to generate ref.bed: {e}")
                sys.exit(1)

    # ── Build sample list ──────────────────────────────────────────────────────
    if _reuse_mode:
        # Discover samples from existing .filt.node files inside result_dir
        # Support both flat layout ({result_dir}/{sample}/) and
        # nested layout ({result_dir}/result/{sample}/).
        _rd_abs = os.path.abspath(args.result_dir)
        _nested = os.path.join(_rd_abs, 'result')
        _src_result = _nested if os.path.isdir(_nested) else _rd_abs
        logger.info(f"  Scanning for .filt.node in {_src_result}")

        # Optional: validate against --fastq_list if supplied together
        _fastq_list_names: set = set()
        if args.fastq_list:
            _fastq_list_names = {s for s, _ in load_fastq_list(args.fastq_list)}

        sample_list = []
        for _sname in sorted(os.listdir(_src_result)):
            _sdir = os.path.join(_src_result, _sname)
            if not os.path.isdir(_sdir):
                continue
            _filt_node = os.path.join(_sdir, f'{_sname}.filt.node')
            if not os.path.isfile(_filt_node) or os.path.getsize(_filt_node) == 0:
                continue
            if _fastq_list_names and _sname not in _fastq_list_names:
                continue
            # Use the .filt.node path as the "fastq" placeholder
            sample_list.append((_sname, _filt_node))

        logger.info(f"  {len(sample_list)} sample(s) found in {_src_result}")
        if not sample_list:
            logger.error("No valid .filt.node files found in result_dir. Exiting.")
            sys.exit(1)

    elif args.fastq:
        # Single sample: derive name from filename
        sample_name = os.path.basename(args.fastq)
        for ext in ('.fq.gz', '.fastq.gz', '.fq', '.fastq'):
            if sample_name.endswith(ext):
                sample_name = sample_name[:-len(ext)]
                break
        sample_list = [(sample_name, os.path.abspath(args.fastq))]
    else:
        sample_list = [(s, os.path.abspath(p))
                       for s, p in load_fastq_list(args.fastq_list)]
    logger.info(f"  {len(sample_list)} sample(s) loaded")
    if not sample_list:
        logger.error("Sample list is empty. Exiting.")
        sys.exit(1)

    # Create per-sample output directories
    for sample, _ in sample_list:
        os.makedirs(os.path.join(outdir, sample), exist_ok=True)

    # ── Skip steps logic ───────────────────────────────────────────────────────
    skip_steps: Set[str] = set()
    if _reuse_mode:
        # In reuse mode, steps 1+2 are always skipped
        skip_steps = {'step1', 'step2'}
        logger.info("  Reuse mode: Steps 1+2 skipped (using existing .filt.node files)")
    elif args.skip_steps:
        skip_steps = {s.strip().lower() for s in args.skip_steps.split(',')}
    # --skip_step1 flag is a shortcut for --skip_steps step1
    if getattr(args, 'skip_step1', False):
        skip_steps.add('step1')
    skip_step1 = 'step1' in skip_steps
    skip_step2 = 'step2' in skip_steps
    logger.info(f"  Skip steps: {skip_steps if skip_steps else 'none'}")
    logger.info(f"  Skip fastp (QC):    {args.skip_fastp}")
    logger.info(f"  Skip rmhost (homo): {args.skip_rmhost}")
    logger.info(f"  Skip pathogen DB:   {args.skip_pathogen}")
    logger.info(f"  Force rerun: {args.force}")

    # ── Build tracing_args dict (passed into each worker) ────────────────────
    tracing_args = dict(
        jaccard_percent   = args.jaccard_percent,
        max_iterations    = args.max_iterations,
        use_containment   = args.use_containment,
        min_containment   = args.min_containment,
        top_n_cutoff      = args.top_n_cutoff,
        min_candidates    = args.min_candidates,
        phylo_gap_ratio   = args.phylo_gap_ratio,
        phylo_share_thresh= args.phylo_share_thresh,
        scoring_mode      = args.scoring_mode,
        cont_alpha        = args.cont_alpha,
        cont_ref_size     = args.cont_ref_size,
        score_cutoff_ratio        = args.score_cutoff_ratio,
        min_unique_snp_n_cutoff   = args.min_unique_snp_n_cutoff,
        min_snp                   = args.min_snp,
        no_visualization  = args.no_visualization,
        keep_temp         = not args.clean_temp,   # keep tracing temps unless --clean_temp
    )
    logger.info("  Tracing params: "
                f"jaccard={tracing_args['jaccard_percent']}  "
                f"max_iter={tracing_args['max_iterations']}  "
                f"containment={tracing_args['use_containment']}  "
                f"scoring_mode={tracing_args['scoring_mode']}  "
                f"phylo_gap_ratio={tracing_args['phylo_gap_ratio']}  "
                f"phylo_share_thresh={tracing_args['phylo_share_thresh']}")

    # ── HTML Viewer config ─────────────────────────────────────────────────
    if args.taxname and not args.no_html_viewer:
        _meta_src = args.meta if args.meta else "(auto-detect from db_dir)"
        logger.info(f"  HTML Viewer: enabled  taxname=\"{args.taxname}\"  meta={_meta_src}")
    else:
        _reason = "--no_html_viewer" if args.no_html_viewer else "--taxname not set"
        logger.info(f"  HTML Viewer: disabled ({_reason})")

    # ── Pipeline timing ────────────────────────────────────────────────────────
    t_start   = time.time()
    succeeded: List[Tuple[str, str]] = []
    failed:    List[Tuple[str, str]] = []

    try:
        # ── STEPS 1 + 2 ───────────────────────────────────────────────────────
        if _reuse_mode or skip_step2:
            if _reuse_mode:
                logger.info("Reuse mode: loading .filt.node files from --result_dir")
                _reuse_pending = []
                for sample, filt_node_src in sample_list:
                    # Copy .filt.node into outdir so all tracing outputs land there
                    _dst_sample_dir = os.path.join(outdir, sample)
                    os.makedirs(_dst_sample_dir, exist_ok=True)
                    _dst_filt_node = os.path.join(_dst_sample_dir, f"{sample}.filt.node")
                    if not os.path.isfile(_dst_filt_node) or args.force:
                        shutil.copy2(filt_node_src, _dst_filt_node)
                    _reuse_pending.append((sample, _dst_filt_node))

                # ── Run Step 3 (tracing) for each sample ─────────────────────
                # Use multiprocessing.Pool (same as SampleRunProcessor) so that:
                #   1. Each worker is a separate process → GIL is not shared
                #      → real CPU parallelism for the tracing algorithm
                #   2. Each process loads the database independently
                #      → no shared-memory race conditions
                #   3. Module-level _reuse_trace_worker is picklable by Pool
                if tracing_args is not None:
                    _n_workers = min(args.threads, len(_reuse_pending))
                    logger.info(
                        f"  Running Step 3 tracing for {len(_reuse_pending)} samples "
                        f"({_n_workers} parallel workers) ...")

                    _job_args = [
                        {
                            'sample':       _s,
                            'filt_node':    _fn,
                            'sample_dir':   os.path.join(outdir, _s),
                            'db_dir':       os.path.abspath(args.db_dir),
                            'tracing_args': tracing_args,
                            'force':        args.force,
                            'log_path':     os.path.join(log_dir, f"{_s}.trace.log"),
                        }
                        for _s, _fn in _reuse_pending
                    ]

                    _n_total   = len(_job_args)
                    _n_done    = 0
                    _n_ok      = 0
                    _n_fail    = 0
                    _t_pool    = time.time()

                    with Pool(_n_workers) as _pool:
                        for _s, _ok, _info, _elapsed in _pool.imap_unordered(
                                _reuse_trace_worker, _job_args):
                            _n_done += 1
                            _rss = logger._rss_mb()
                            _pct = _n_done / _n_total * 100
                            _wall = time.time() - _t_pool
                            if _ok:
                                _n_ok += 1
                                succeeded.append((_s, _info))
                                logger.info(
                                    f"  [OK]   {_s}  [{_elapsed:.1f}s | RSS:{_rss:.0f}MB]\n"
                                    f"  Progress: {_n_done}/{_n_total} ({_pct:.1f}%)  "
                                    f"ok={_n_ok} fail={_n_fail}")
                            else:
                                _n_fail += 1
                                failed.append((_s, _info))
                                logger.error(
                                    f"  [FAIL] {_s}  [{_elapsed:.1f}s]\n"
                                    f"  Progress: {_n_done}/{_n_total} ({_pct:.1f}%)  "
                                    f"ok={_n_ok} fail={_n_fail}\n"
                                    f"  reason: {str(_info)[:120]}")
                else:
                    # tracing disabled — just mark as succeeded
                    succeeded.extend(_reuse_pending)

            else:
                logger.info("Skipping Steps 1+2 (--skip_steps step2)")
                # Reconstruct from existing .filt.node
                for sample, _ in sample_list:
                    filt_node = os.path.join(outdir, sample, f"{sample}.filt.node")
                    if os.path.isfile(filt_node):
                        succeeded.append((sample, filt_node))
                    else:
                        logger.warning(f"  [WARN] {sample}: .filt.node not found at {filt_node}")

        elif args.use_sge:
            eff_threads = args.sge_threads if args.sge_threads > 0 else args.sge_slots  # default: 1
            sge_cfg = SGEConfig(
                nodes    = args.sge_nodes,
                mem      = args.sge_mem,
                slots    = args.sge_slots,
                queue    = args.sge_queue,
                max_jobs = args.sge_max_jobs,
                poll_sec = args.sge_poll_sec,
                extra    = args.sge_extra,
            )
            logger.info(f"SGE mode enabled: {sge_cfg!r}")
            logger.info(
                f"  SGE resource request: -l {sge_cfg.build_resource_str()}"
                f"  -q {sge_cfg.queue}")
            logger.info(
                f"  Tool threads/job    : {eff_threads} "
                f"({'= sge_slots' if eff_threads == args.sge_slots else '--sge_threads override'})")
            logger.info(f"  Max concurrent jobs : {sge_cfg.max_jobs}")

            # Log per-node capacity
            _node_cores = {
                'R7525': 256, 'R940xa': 160,
                'big2': 128, 'big3': 128, 'big4': 128,
                'Dell7920': 80, 'node10': 80, 'node11': 80, 'node12': 80,
            }
            target_nodes = (args.sge_nodes.split('|')
                            if args.sge_nodes else list(_node_cores.keys()))
            for nd in target_nodes:
                cores = _node_cores.get(nd, '?')
                cap   = (cores // args.sge_slots
                         if isinstance(cores, int) else '?')
                logger.info(f"    {nd:<12} {cores} cores → max {cap} concurrent jobs")

            submitter = SGERunSubmitter(
                sample_list     = sample_list,
                ref_fa          = ref_fa,
                ref_bed         = ref_bed,
                db_dir          = os.path.abspath(args.db_dir),
                outdir          = outdir,
                taxid_str       = args.taxid,
                config          = config,
                sge_cfg         = sge_cfg,
                logger          = logger,
                threads_per_job = eff_threads,
                force           = args.force,
                log_dir         = log_dir,
                identity_thr    = args.filter_identity,
                coverage_thr    = args.filter_coverage,
                tracing_args    = tracing_args,
                skip_step1      = skip_step1,
                skip_fastp      = args.skip_fastp,
                skip_rmhost     = args.skip_rmhost,
                skip_pathogen   = args.skip_pathogen,
                clean_temp      = args.clean_temp,
                sge_qsubname    = args.sge_qsubname,
            )
            succeeded, failed = submitter.run()

        else:
            processor = SampleRunProcessor(
                sample_list  = sample_list,
                ref_fa       = ref_fa,
                ref_bed      = ref_bed,
                db_dir       = os.path.abspath(args.db_dir),
                outdir       = outdir,
                taxid_str    = args.taxid,
                config       = config,
                logger       = logger,
                threads      = args.threads,
                force        = args.force,
                log_dir      = log_dir,
                identity_thr = args.filter_identity,
                coverage_thr = args.filter_coverage,
                tracing_args = tracing_args,
                skip_step1   = skip_step1,
                skip_fastp   = args.skip_fastp,
                skip_rmhost  = args.skip_rmhost,
                skip_pathogen= args.skip_pathogen,
                clean_temp   = args.clean_temp,
            )
            succeeded, failed = processor.run()

        # ── STEP 3: Per-sample tracing results collection ─────────────────────
        # Tracing is now done inside each worker (Step 3 embedded per sample).
        # Here we just collect TraceResult/ entries for the final summary.
        if not succeeded:
            logger.error("No successful samples completed Steps 1+2. Exiting.")
            sys.exit(1)

        batch_results = []
        # Canonical column order matching Trace.result.tsv (FinalResult=Yes rows only)
        _BATCH_COLS = [
            'Sample', 'Genomes', 'GenomesNumber', 'FinalResult', 'TraceLevel',
            'ConfidenceScore', 'ConfidenceLevel',
            'ContainmentRankScore', 'Reads', 'Coverage', 'Depth',
            'FractionMappedReads', 'TotalSNP', 'TracingSnpUsed',
            'SnpRank', 'GlobalPath', 'GlobalPathScores', 'GlobalNodeSNPCounts',
            'LCA_Round', 'SubPath', 'SubPathScores', 'NodeSNPCounts',
            'ZeroRatio', 'JaccardClusterID',
        ]
        if tracing_args is not None:
            import csv as _csv
            logger.info(f"Collecting per-sample tracing results ...")
            for _sample, _filt_node in succeeded:
                _trace_dir  = os.path.join(outdir, _sample, 'TraceResult')
                _result_tsv = os.path.join(_trace_dir, 'Trace.result.tsv')
                _done_f     = os.path.join(_trace_dir, '.trace.done')
                _fail_f     = os.path.join(_trace_dir, '.trace.fail')
                if os.path.isfile(_done_f) and os.path.isfile(_result_tsv):
                    # Read FinalResult=Yes row directly — no column remapping
                    _sample_row = None
                    with open(_result_tsv, 'r', newline='') as _rt:
                        _reader = _csv.DictReader(_rt, delimiter='\t')
                        for _row in _reader:
                            if _row.get('FinalResult', '') == 'Yes':
                                _sample_row = _row
                                break
                    if _sample_row:
                        # Copy all columns from Trace.result.tsv directly
                        batch_results.append(
                            {col: _sample_row.get(col, 'N/A') for col in _BATCH_COLS}
                        )
                    else:
                        # No FinalResult=Yes row (all paths failed / empty)
                        _empty = {col: 'N/A' for col in _BATCH_COLS}
                        _empty.update({'Sample': _sample, 'FinalResult': 'No',
                                       'TraceLevel': 'Unknown', 'ConfidenceScore': '0',
                                       'ConfidenceLevel': 'Low', 'GenomesNumber': '0'})
                        batch_results.append(_empty)
                    logger.info(f"  [OK]   {_sample}: TraceResult/Trace.result.tsv")
                elif os.path.isfile(_fail_f):
                    with open(_fail_f) as _ff:
                        _reason = _ff.read(500).strip()
                    logger.error(f"  [FAIL] {_sample}: Step 3 failed — {_reason[:200]}")
                    _fail_row = {col: 'N/A' for col in _BATCH_COLS}
                    _fail_row.update({'Sample': _sample, 'FinalResult': 'Error',
                                      'TraceLevel': 'Error', 'ConfidenceScore': '0',
                                      'ConfidenceLevel': 'Low', 'GenomesNumber': '0',
                                      'Genomes': f'FAIL:{_reason[:80]}'})
                    batch_results.append(_fail_row)
                else:
                    logger.warning(f"  [WARN] {_sample}: TraceResult not found — "
                                   f"tracing may have been skipped or failed silently")
                    _missing = {col: 'N/A' for col in _BATCH_COLS}
                    _missing.update({'Sample': _sample, 'FinalResult': 'Missing',
                                     'TraceLevel': 'Unknown', 'ConfidenceScore': '0',
                                     'ConfidenceLevel': 'Low', 'GenomesNumber': '0'})
                    batch_results.append(_missing)

            # Write aggregate batch_summary.tsv into outdir/result_comb/
            _res_dir = os.path.join(outdir, 'result_comb')
            os.makedirs(_res_dir, exist_ok=True)
            _summary = os.path.join(_res_dir, 'batch_summary.tsv')
            if batch_results:
                with open(_summary, 'w', newline='') as _sf:
                    writer = _csv.DictWriter(_sf, delimiter='\t', fieldnames=_BATCH_COLS)
                    writer.writeheader()
                    for _r in batch_results:
                        writer.writerow({k: _r.get(k, 'N/A') for k in _BATCH_COLS})
                logger.info(f"  Batch summary: {_summary}")
                # Run statistics + visualization on the combined batch summary
                try:
                    run_batch_statistics(_summary, _res_dir, logger)
                except Exception as _stat_exc:
                    logger.warning(f"  [WARN] Batch statistics failed (non-fatal): {_stat_exc}")

        # ── STEP 4a: Per-sample HTML Viewer ───────────────────────────────
        if (args.taxname and not args.no_html_viewer and succeeded):
            try:
                run_html_viewer_for_samples(
                    succeeded=succeeded,
                    outdir=outdir,
                    db_dir=os.path.abspath(args.db_dir),
                    taxname=args.taxname,
                    logger=logger,
                    meta_tsv=args.meta,
                )
            except Exception as _viewer_exc:
                logger.warning(f"  [WARN] HTML viewer generation failed (non-fatal): {_viewer_exc}")
        elif not args.taxname and not args.no_html_viewer:
            logger.info("  HTML viewer skipped (--taxname not provided)")

    except KeyboardInterrupt:
        logger.error("Pipeline interrupted by user (Ctrl-C)")
        sys.exit(130)
    except Exception as exc:
        logger.error(f"Pipeline failed: {exc}")
        logger.error(traceback.format_exc())
        sys.exit(1)

    # ── Final summary ──────────────────────────────────────────────────────────
    elapsed = time.time() - t_start
    logger.record_submodule(
        module='pipeline', submodule='total',
        elapsed_sec=elapsed, rss_mb=logger._rss_mb(),
        notes=f"samples={len(sample_list)} ok={len(succeeded)} fail={len(failed)}")
    print_summary(sample_list, succeeded, failed, batch_results,
                  outdir, elapsed, logger)
    sys.exit(0 if not failed else 1)


if __name__ == '__main__':
    main()
