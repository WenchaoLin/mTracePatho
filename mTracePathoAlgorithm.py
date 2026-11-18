#!/usr/bin/env python3
"""
mTracePatho v1.0 - Pathogen Source Tracing Tool

Traces the source of a pathogen sample by mapping its SNPs through a
phylogenetic database using LCA (Lowest Common Ancestor) algorithm with
iterative sub-database refinement.

Key features:
  - Containment-based genome pre-selection (several scoring modes)
  - Iterative LCA tracing with sub-database refinement
  - Confidence scoring (SNP mapping rate, path dominance, consistency, uniqueness)
  - Tie-breaking via containment scores + LCA of tied genomes
  - Hard score cutoff ratio filter
  - Minimum SNP warning flag (--min_snp)

Confidence score thresholds:
  High >= 0.8  |  Medium 0.5–0.8  |  Low < 0.5

Version: 1.0
"""

import os
import sys
import argparse
import csv
import math
import time
import shutil
import warnings
import traceback
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from collections import defaultdict
from datetime import datetime
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

sys.setrecursionlimit(10000)

# =============================================================================
# Worker initializer for multiprocessing (BUG#11: avoid N×pickle of genome_snps_map)
# =============================================================================

_worker_genome_snps_map = None

def _worker_init(db_dir):
    """Load genome_snps_map once per worker process instead of pickling per task."""
    global _worker_genome_snps_map
    try:
        _worker_genome_snps_map, _, _ = load_snp_matrix_for_containment(db_dir)
    except Exception as e:
        _worker_genome_snps_map = None
        print(f"  [WARNING] Worker init: failed to load SNP matrix: {e}")

# =============================================================================
# Confidence Score Calculation (Accuracy-focused)
# =============================================================================

def calculate_confidence_score(mapped_snps, total_snps, best_path_info, all_paths_info, 
                                genome_count, max_depth, path_depth, trace_level):
    """
    Calculate confidence score for tracing result - focused on ACCURACY factors.
    
    Components (all relate to accuracy/reliability of result):
    1. SNP Mapping Rate (30%): Proportion of sample SNPs successfully mapped
       - Higher = more information available for tracing
    
    2. Path Dominance (35%): How much better is the best path vs alternatives
       - Higher = clearer distinction, more confident in the choice
    
    3. SNP Consistency (20%): How consistently SNPs support the best path
       - Lower zero_ratio = SNPs distributed along path, not just at end
       - SNPs concentrated at terminal nodes = stronger evidence
    
    4. Path Uniqueness (15%): How many competing paths exist
       - Fewer paths = less ambiguity
    
    Returns:
        (total_score, component_scores_dict)
    """
    components = {}
    
    # ==========================================================================
    # 1. SNP Mapping Rate (30%) - 信息量指标
    # ==========================================================================
    # 样本SNP能够映射到数据库的比例，映射率越高，溯源依据越充分
    if total_snps > 0:
        snp_mapping_rate = mapped_snps / total_snps
    else:
        snp_mapping_rate = 0.0
    components['snp_mapping_rate'] = snp_mapping_rate
    
    # ==========================================================================
    # 2. Path Dominance (35%) - 区分度指标 (最重要)
    # ==========================================================================
    # 最佳路径相对其他路径的优势程度，优势越明显结果越可靠
    if best_path_info and all_paths_info and len(all_paths_info) > 1:
        # 按累计分数排序
        sorted_paths = sorted(all_paths_info, key=lambda x: x.get('total_snp', 0), reverse=True)
        best_snp = sorted_paths[0].get('total_snp', 0)
        second_snp = sorted_paths[1].get('total_snp', 0) if len(sorted_paths) > 1 else 0
        
        # 方法1: SNP数量差异比
        if best_snp > 0:
            if second_snp > 0:
                # (best - second) / best: 差异越大越好
                snp_diff_ratio = (best_snp - second_snp) / best_snp
            else:
                snp_diff_ratio = 1.0  # 次优路径无SNP
        else:
            snp_diff_ratio = 0.0
        
        # 方法2: 累计分数优势
        # 比较路径末端的累计分数
        best_scores = [int(s) for s in best_path_info.get('scores_str', '0').split(',')]
        best_terminal_score = best_scores[-1] if best_scores else 0
        
        second_path = sorted_paths[1] if len(sorted_paths) > 1 else None
        if second_path:
            second_scores = [int(s) for s in second_path.get('scores_str', '0').split(',')]
            second_terminal_score = second_scores[-1] if second_scores else 0
        else:
            second_terminal_score = 0
        
        if best_terminal_score > 0:
            if second_terminal_score > 0:
                score_diff_ratio = (best_terminal_score - second_terminal_score) / best_terminal_score
            else:
                score_diff_ratio = 1.0
        else:
            score_diff_ratio = 0.0
        
        # 方法3: 考虑所有路径的分散度
        all_total_snps = [p.get('total_snp', 0) for p in all_paths_info]
        if sum(all_total_snps) > 0:
            # 最佳路径占总SNP的比例
            best_proportion = best_snp / sum(all_total_snps)
        else:
            best_proportion = 1.0
        
        # 综合三个指标
        path_dominance = (snp_diff_ratio * 0.4 + score_diff_ratio * 0.3 + best_proportion * 0.3)
        
    else:
        path_dominance = 1.0  # 只有一条路径，完全确定
    
    components['path_dominance'] = path_dominance
    
    # ==========================================================================
    # 3. SNP Consistency (20%) - 证据一致性指标
    # ==========================================================================
    # SNP沿路径的分布一致性，反映证据的可靠性
    if best_path_info:
        snp_str = best_path_info.get('snp_str', '')
        if snp_str:
            snp_counts = [int(x) for x in snp_str.split(',')]
            path_length = len(snp_counts)
            
            if path_length > 1:
                # 1) Zero ratio: 路径上有SNP支持的节点比例
                # 排除根节点(通常为0)
                branch_snps = snp_counts[1:]  # 去掉根节点
                if branch_snps:
                    non_zero_count = sum(1 for s in branch_snps if s > 0)
                    non_zero_ratio = non_zero_count / len(branch_snps)
                else:
                    non_zero_ratio = 0.0
                
                # 2) Terminal concentration: SNP是否集中在末端
                # 末端1/3路径的SNP占比
                terminal_portion = max(1, path_length // 3)
                terminal_snps = sum(snp_counts[-terminal_portion:])
                total_path_snps = sum(snp_counts)
                if total_path_snps > 0:
                    terminal_ratio = terminal_snps / total_path_snps
                else:
                    terminal_ratio = 0.0
                
                # 综合: 既要有分布(non_zero_ratio)又要末端集中(terminal_ratio)
                # 理想情况: SNP沿路径都有支持，且末端有独特SNP
                snp_consistency = (non_zero_ratio * 0.4 + terminal_ratio * 0.6)
            else:
                snp_consistency = 0.5  # 路径太短
        else:
            snp_consistency = 0.0
    else:
        snp_consistency = 0.0
    
    components['snp_consistency'] = snp_consistency
    
    # ==========================================================================
    # 4. Path Uniqueness (15%) - 路径唯一性指标
    # ==========================================================================
    # 候选路径数量越少，结果歧义性越低
    if all_paths_info:
        n_paths = len(all_paths_info)
        
        # 同时考虑不同Jaccard cluster的数量
        unique_clusters = len(set(p.get('jaccard_cluster', 0) for p in all_paths_info))
        
        # 路径数量评分
        if n_paths == 1:
            path_count_score = 1.0
        elif n_paths == 2:
            path_count_score = 0.9
        elif n_paths <= 5:
            path_count_score = 0.7
        elif n_paths <= 10:
            path_count_score = 0.5
        elif n_paths <= 20:
            path_count_score = 0.3
        else:
            path_count_score = 0.1
        
        # Cluster数量评分
        if unique_clusters == 1:
            cluster_score = 1.0  # 所有路径属于同一cluster，本质上是一个方向
        elif unique_clusters == 2:
            cluster_score = 0.7
        elif unique_clusters <= 5:
            cluster_score = 0.4
        else:
            cluster_score = 0.2
        
        path_uniqueness = (path_count_score * 0.5 + cluster_score * 0.5)
    else:
        path_uniqueness = 0.0
    
    components['path_uniqueness'] = path_uniqueness
    
    # ==========================================================================
    # Calculate weighted total
    # ==========================================================================
    weights = {
        'snp_mapping_rate': 0.30,   # 30% - 信息充分性
        'path_dominance': 0.35,     # 35% - 区分度(最重要)
        'snp_consistency': 0.20,    # 20% - 证据一致性
        'path_uniqueness': 0.15     # 15% - 歧义性
    }
    
    total_score = sum(components[k] * weights[k] for k in weights)
    
    return total_score, components


def get_confidence_level(score):
    """Get confidence level label based on score."""
    if score >= 0.8:
        return 'High'
    elif score >= 0.5:
        return 'Medium'
    else:
        return 'Low'


# =============================================================================
# Timer and Logger Classes
# =============================================================================

class Timer:
    def __init__(self):
        self.start_time = time.time()
    
    def elapsed(self):
        return time.time() - self.start_time
    
    def format_time(self, seconds):
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            return f"{seconds/60:.1f}min"
        else:
            return f"{seconds/3600:.1f}h"


class TimingLogger:
    def __init__(self, log_file):
        self.log_file = log_file
        self.file_handle = open(log_file, 'w')
        self.file_handle.write(f"mTracePatho v1.0 - Timing Log\n")
        self.file_handle.write(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.file_handle.write("=" * 70 + "\n\n")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False  # 不吞掉异常

    def log(self, message):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.file_handle.write(f"[{timestamp}] {message}\n")
        self.file_handle.flush()
    
    def log_sample(self, sample_name, iterations, genomes, target, elapsed, confidence='N/A'):
        self.log(f"Sample: {sample_name}")
        self.log(f"  Iterations: {iterations}")
        self.log(f"  Final genomes: {genomes}")
        self.log(f"  Target: {target}")
        self.log(f"  Confidence: {confidence}")
        self.log(f"  Time: {elapsed}")
        self.file_handle.write("\n")
        self.file_handle.flush()
    
    def close(self):
        if not self.file_handle.closed:
            self.file_handle.write("=" * 70 + "\n")
            self.file_handle.write(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            self.file_handle.close()


# =============================================================================
# Containment Calculation and Elbow Detection
# =============================================================================

def load_snp_matrix_for_containment(db_dir):
    filepath = os.path.join(db_dir, 'total.var.tsv')
    genome_snps = defaultdict(set)
    all_snps = set()
    
    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        genomes = header[1:]
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            snp_name = fields[0]
            all_snps.add(snp_name)
            for i, val in enumerate(fields[1:]):
                if i < len(genomes):
                    try:
                        if float(val) == 1:
                            genome_snps[genomes[i]].add(snp_name)
                    except ValueError:
                        pass
    return genome_snps, all_snps, genomes


def calculate_containment(sample_snps, genome_snps,
                           scoring_mode='containment_weighted',
                           cont_alpha=0.5, cont_ref_size=None):
    """
    Compute containment (or alternative scoring) for all candidate genomes.

    scoring_mode options
    --------------------
    'containment'  (default / original)
        score = SharedSNPs / SampleSNPs
        Pure sample-side coverage; ignores genome size entirely.

    'containment_weighted'
        score = containment × (ref_size / GenomeSNPs) ^ cont_alpha
        Penalises large genomes; ref_size is the median GenomeSNPs across all
        candidates (auto-computed when cont_ref_size=None) or a user-specified
        constant.  alpha=0 → identical to containment; alpha=1 → linear size
        penalty (≈ Jaccard behaviour in shallow-sequencing limit).

    'jaccard'
        score = SharedSNPs / (SampleSNPs + GenomeSNPs - SharedSNPs)
        Symmetric, parameter-free.  In shallow-sequencing (SampleSNPs << Genome
        SNPs) the denominator is dominated by GenomeSNPs, which provides a
        strong and stable penalty for large genomes.

    'cosine'
        score = SharedSNPs / (sqrt(SampleSNPs) * sqrt(GenomeSNPs))
        Cosine similarity for binary SNP vectors.  Provides moderate genome-size
        penalty (geometric mean normalisation).

    'cosine_weighted'
        score = cosine_sim × (ref_size / GenomeSNPs) ^ cont_alpha
        Cosine similarity with additional genome-size penalty, analogous to
        'containment_weighted' mode but based on cosine distance instead of containment.

    'jaccard_weighted'
        score = jaccard × (ref_size / GenomeSNPs) ^ cont_alpha
        Jaccard similarity with additional genome-size penalty.

    Returns
    -------
    List of dicts sorted by 'score' descending.  Every dict contains:
        genome        : genome identifier
        containment   : raw SharedSNPs / SampleSNPs  (always, for display)
        score         : effective ranking score (mode-dependent)
        shared_snps   : |intersection|
        genome_snps   : |genome SNP set|
        scoring_mode  : the mode used
        (weighted/cosine_weighted/jaccard_weighted only) ref_size, penalty
    """
    results = []
    sample_size = len(sample_snps)
    if sample_size == 0:
        return results

    # Pre-compute ref_size for weighted-family modes (median of all genome sizes)
    needs_ref = scoring_mode in ('containment_weighted', 'cosine_weighted', 'jaccard_weighted')
    if needs_ref and cont_ref_size is None:
        all_sizes = [len(g) for g in genome_snps.values()]
        cont_ref_size = float(np.median(all_sizes)) if all_sizes else 1.0

    import math as _math

    for genome, g_snps in genome_snps.items():
        shared    = sample_snps & g_snps
        n_shared  = len(shared)
        n_genome  = len(g_snps)

        raw_cont = n_shared / sample_size

        if scoring_mode == 'containment_weighted':
            penalty = (cont_ref_size / n_genome) ** cont_alpha if n_genome > 0 else 1.0
            score   = raw_cont * penalty
            results.append({
                'genome': genome, 'containment': raw_cont, 'score': score,
                'base_metric': raw_cont, 'base_metric_name': 'Containment',
                'shared_snps': n_shared, 'genome_snps': n_genome,
                'scoring_mode': scoring_mode,
                'ref_size': cont_ref_size, 'penalty': round(penalty, 6),
            })
        elif scoring_mode == 'jaccard':
            denom = sample_size + n_genome - n_shared
            score = n_shared / denom if denom > 0 else 0.0
            results.append({
                'genome': genome, 'containment': raw_cont, 'score': score,
                'base_metric': score, 'base_metric_name': 'Jaccard',
                'shared_snps': n_shared, 'genome_snps': n_genome,
                'scoring_mode': scoring_mode,
            })
        elif scoring_mode == 'cosine':
            # Cosine similarity for binary vectors:
            # cos_sim = |A ∩ B| / (sqrt(|A|) * sqrt(|B|))
            denom = _math.sqrt(sample_size) * _math.sqrt(n_genome) if n_genome > 0 else 0.0
            score = n_shared / denom if denom > 0 else 0.0
            results.append({
                'genome': genome, 'containment': raw_cont, 'score': score,
                'base_metric': score, 'base_metric_name': 'Cosine',
                'shared_snps': n_shared, 'genome_snps': n_genome,
                'scoring_mode': scoring_mode,
            })
        elif scoring_mode == 'cosine_weighted':
            # Cosine similarity × genome-size penalty
            denom = _math.sqrt(sample_size) * _math.sqrt(n_genome) if n_genome > 0 else 0.0
            cos_sim = n_shared / denom if denom > 0 else 0.0
            penalty = (cont_ref_size / n_genome) ** cont_alpha if n_genome > 0 else 1.0
            score   = cos_sim * penalty
            results.append({
                'genome': genome, 'containment': raw_cont, 'score': score,
                'base_metric': cos_sim, 'base_metric_name': 'Cosine',
                'shared_snps': n_shared, 'genome_snps': n_genome,
                'scoring_mode': scoring_mode,
                'ref_size': cont_ref_size, 'penalty': round(penalty, 6),
                'cosine_sim': round(cos_sim, 6),
            })
        elif scoring_mode == 'jaccard_weighted':
            # Jaccard × genome-size penalty
            j_denom = sample_size + n_genome - n_shared
            jacc = n_shared / j_denom if j_denom > 0 else 0.0
            penalty = (cont_ref_size / n_genome) ** cont_alpha if n_genome > 0 else 1.0
            score   = jacc * penalty
            results.append({
                'genome': genome, 'containment': raw_cont, 'score': score,
                'base_metric': jacc, 'base_metric_name': 'Jaccard',
                'shared_snps': n_shared, 'genome_snps': n_genome,
                'scoring_mode': scoring_mode,
                'ref_size': cont_ref_size, 'penalty': round(penalty, 6),
                'jaccard_sim': round(jacc, 6),
            })
        else:  # 'containment' — original behaviour
            results.append({
                'genome': genome, 'containment': raw_cont, 'score': raw_cont,
                'base_metric': raw_cont, 'base_metric_name': 'Containment',
                'shared_snps': n_shared, 'genome_snps': n_genome,
                'scoring_mode': scoring_mode,
            })

    results.sort(key=lambda x: x['score'], reverse=True)
    return results





def build_node_to_snps(snp_to_node: dict) -> dict:
    """Build reverse index: node_id → set of SNPs assigned to that node."""
    from collections import defaultdict as _dd
    node_to_snps = _dd(set)
    for snp, node_id in snp_to_node.items():
        node_to_snps[node_id].add(snp)
    return dict(node_to_snps)


def filter_by_unique_snp_presence(selected_genomes: list,
                                   sample_snps: set,
                                   snp_to_node: dict,
                                   genome_to_node: dict,
                                   min_unique_snp_n_cutoff: int = 1) -> tuple:
    """
    Hard-filter candidate genomes by terminal-node (leaf) unique SNP presence.

    A genome is removed when:
      - It has ≥1 SNP uniquely assigned to its terminal node in the database, AND
      - The number of those unique SNPs found in sample_snps < min_unique_snp_n_cutoff.

    Genomes whose terminal node has NO assigned SNPs (i.e. no node-unique SNPs
    in the database) are kept regardless, to avoid false exclusions.

    Args:
        selected_genomes:        List of genome names after score_cutoff_ratio filter.
        sample_snps:             Set of SNP strings detected in the sample.
        snp_to_node:             {snp_str: node_id} from LCA.output.tsv.
        genome_to_node:          {genome_name: terminal_node_id} from snp_lingage.tsv.
        min_unique_snp_n_cutoff: Minimum number of a genome's unique SNPs that must
                                 appear in the sample to keep that genome (default 1).

    Returns:
        (kept, info_dict)
        info_dict keys:
            n_before, n_after, removed_genomes (list of (genome, n_unique, n_hit))
    """
    if not selected_genomes or min_unique_snp_n_cutoff <= 0:
        return selected_genomes, {
            'n_before': len(selected_genomes), 'n_after': len(selected_genomes),
            'removed_genomes': [], 'activated': False
        }

    # Build reverse index once (cheap for typical db sizes)
    node_to_snps = build_node_to_snps(snp_to_node)

    kept = []
    removed = []
    for genome in selected_genomes:
        terminal_node = genome_to_node.get(genome)
        if terminal_node is None:
            # Unknown genome — keep conservatively
            kept.append(genome)
            continue
        unique_snps = node_to_snps.get(terminal_node, set())
        if not unique_snps:
            # No unique SNPs defined for this terminal node → keep
            kept.append(genome)
            continue
        n_hit = len(unique_snps & sample_snps)
        if n_hit >= min_unique_snp_n_cutoff:
            kept.append(genome)
        else:
            removed.append((genome, len(unique_snps), n_hit))

    # Safety: if ALL candidates would be removed, restore the original list
    if len(kept) == 0 and removed:
        print(f"  [UniqueSNPFilter] WARNING: filter would remove ALL "
              f"{len(selected_genomes)} candidates; filter skipped.")
        return selected_genomes, {
            'n_before': len(selected_genomes), 'n_after': len(selected_genomes),
            'removed_genomes': removed, 'activated': True, 'skipped': True
        }

    return kept, {
        'n_before':        len(selected_genomes),
        'n_after':         len(kept),
        'removed_genomes': removed,   # [(genome, n_unique_in_db, n_hit_in_sample)]
        'activated':       True,
        'skipped':         False,
    }


def select_candidates_by_containment(containment_results, min_containment=0.0, top_n_cutoff=99999, min_candidates=3,
                                      score_cutoff_ratio=0.5):
    """
    Select candidate genomes based on containment score using top N distinct value levels.

    Uses the 'score' field (set by calculate_containment according to scoring_mode)
    for all cutoff comparisons.  Falls back to 'containment' for backward compatibility
    if 'score' is absent.

    Hard filter (score_cutoff_ratio):
        After top_n_cutoff / min_candidates selection, ALL remaining candidates
        must satisfy:  score >= top1_score × score_cutoff_ratio
        Candidates below this threshold are removed regardless of other criteria.
        If only 1 genome survives, it will be returned directly as tracing result.

    Args:
        containment_results: List of containment results sorted by score (descending)
        min_containment: Minimum score threshold (default 0.0)
        top_n_cutoff: Select genomes within top N distinct score levels (default 99999 = all levels)
        min_candidates: Minimum number of candidates to select (default 5)
        score_cutoff_ratio: Hard cutoff: candidates must have score >= top1 * ratio (default 0.5)

    Returns:
        (selected_genomes, selection_info)
    """
    if not containment_results:
        return [], {'method': 'none', 'cutoff': 0, 'n_selected': 0}

    # Precision for floating point comparison (6 decimal places)
    PRECISION = 6

    # Use 'score' field (mode-aware); fall back to 'containment' for compat
    def _score(r):
        return r.get('score', r['containment'])

    # Get all score values with rounding to handle float precision
    all_values = [round(_score(r), PRECISION) for r in containment_results]
    values_array = np.array(all_values)
    non_zero_values = values_array[values_array > 0]

    if len(non_zero_values) == 0:
        return [], {'method': 'no_match', 'cutoff': 0, 'n_selected': 0}

    # Get DISTINCT non-zero values sorted in descending order
    distinct_values = sorted(set(v for v in all_values if v > 0), reverse=True)

    # Start with top_n_cutoff distinct levels
    effective_n_cutoff = top_n_cutoff
    selected = []
    actual_cutoff = 0

    # Iteratively expand if we don't have enough candidates
    while effective_n_cutoff <= len(distinct_values):
        cutoff_idx   = min(effective_n_cutoff - 1, len(distinct_values) - 1)
        cutoff_value = distinct_values[cutoff_idx]

        # Ensure cutoff is at least min_containment
        actual_cutoff = max(cutoff_value, min_containment)

        # Select all genomes with score >= cutoff (using rounded comparison)
        selected = [r for r in containment_results
                    if round(_score(r), PRECISION) >= actual_cutoff]

        # Check if we have enough candidates
        if len(selected) >= min_candidates:
            break

        # Expand to include more distinct levels
        effective_n_cutoff += 1

    # If still not enough after all levels, use all non-zero
    if len(selected) < min_candidates:
        selected      = [r for r in containment_results if _score(r) > 0]
        actual_cutoff = min(round(_score(r), PRECISION) for r in selected) if selected else 0
        effective_n_cutoff = len(distinct_values)

    # =====================================================================
    # Hard filter: score_cutoff_ratio — all candidates must be >= top1 × ratio
    # =====================================================================
    n_before_hard_filter = len(selected)
    hard_cutoff_value = 0.0
    hard_filtered_count = 0
    if selected and score_cutoff_ratio > 0:
        top1_score = round(_score(selected[0]), PRECISION)  # selected is sorted desc
        hard_cutoff_value = round(top1_score * score_cutoff_ratio, PRECISION)
        before_count = len(selected)
        # Save pre-filter list in case hard filter removes everything
        pre_hard_filter_selected = selected[:]
        selected = [r for r in selected
                    if round(_score(r), PRECISION) >= hard_cutoff_value]
        hard_filtered_count = before_count - len(selected)
        # Issue 5 fix: if hard filter removes ALL candidates, restore pre-filter list
        if len(selected) == 0 and pre_hard_filter_selected:
            selected = pre_hard_filter_selected
            hard_filtered_count = 0
            print(f"  [HardFilter] WARNING: score_cutoff_ratio={score_cutoff_ratio} removed "
                  f"all {before_count} candidates; restored pre-filter selection.")


    # Calculate statistics for visualization
    max_value = _score(containment_results[0]) if containment_results else 0

    # Calculate top1 statistics (genomes at 1st distinct level)
    top1_value = round(max_value, PRECISION)
    top1_count = sum(1 for r in containment_results
                     if round(_score(r), PRECISION) == top1_value)
    top1_ratio = top1_count / len(containment_results) if containment_results else 0

    # Calculate topN statistics (genomes at cutoff value or above)
    topn_count = len(selected)
    topn_ratio = topn_count / len(containment_results) if containment_results else 0

    # Calculate other statistics
    median_value = np.median(non_zero_values) if len(non_zero_values) > 0 else 0
    mean_value   = np.mean(non_zero_values)   if len(non_zero_values) > 0 else 0
    p99_value    = np.percentile(non_zero_values, 99) if len(non_zero_values) > 0 else 0

    scoring_mode = containment_results[0].get('scoring_mode', 'containment') if containment_results else 'containment'

    selection_info = {
        'method': f'top{effective_n_cutoff}_distinct_levels',
        'cutoff': actual_cutoff,
        'cutoff_rank': effective_n_cutoff,
        'original_top_n': top_n_cutoff,
        'effective_top_n': effective_n_cutoff,
        'n_selected': len(selected),
        'total_genomes': len(containment_results),
        'non_zero_count': len(non_zero_values),
        'distinct_levels': len(distinct_values),
        'min_candidates': min_candidates,
        'scoring_mode': scoring_mode,
        'score_cutoff_ratio': score_cutoff_ratio,
        'hard_cutoff_value': hard_cutoff_value,
        'hard_filtered_count': hard_filtered_count,
        'n_before_hard_filter': n_before_hard_filter,
        'statistics': {
            'max': max_value,
            'median': median_value,
            'mean': mean_value,
            'p99': p99_value,
            'top1_value': top1_value,
            'top1_count': top1_count,
            'top1_ratio': top1_ratio,
            'topn_count': topn_count,
            'topn_ratio': topn_ratio,
            'top_n_cutoff': effective_n_cutoff,
            'distinct_values_top5': distinct_values[:5]
        }
    }

    return [r['genome'] for r in selected], selection_info


def find_elbow_for_visualization(values):
    """Find elbow point for visualization only (not used for filtering)."""
    n = len(values)
    if n <= 5:
        return n - 1, values[-1] if values else 0
    
    values = np.array(values)
    
    # Simple gradient-based elbow detection
    drops = []
    for i in range(1, n):
        if values[i-1] > 0:
            relative_drop = (values[i-1] - values[i]) / values[i-1]
            drops.append((i, relative_drop, values[i]))
    
    # Find first significant drop
    for idx, drop, val in drops:
        if drop > 0.2 and idx >= 3:
            return idx, val
    
    # Default: use 20% of data
    default_idx = min(n // 5, n - 1)
    return default_idx, values[default_idx] if default_idx < n else 0


def plot_containment_curve(containment_results, selection_info, output_path, sample_name):
    """Plot containment curve with statistics markers (single panel)."""
    if not containment_results:
        return
    
    fig, ax = plt.subplots(figsize=(12, 7))
    fig.suptitle(f'Containment Analysis: {sample_name}', fontsize=14, fontweight='bold')
    
    n = len(containment_results)
    x = range(1, n + 1)
    # BUG#B fix: use 'score' (mode-aware) so curve aligns with cutoff lines;
    # fall back to 'containment' for backward compatibility.
    y = [r.get('score', r['containment']) for r in containment_results]
    
    # Plot containment curve
    ax.plot(x, y, 'b-', linewidth=2, label='Containment')
    ax.fill_between(x, y, alpha=0.3)
    
    # Get statistics
    stats = selection_info.get('statistics', {})
    cutoff = selection_info.get('cutoff', 0)
    n_selected = selection_info.get('n_selected', 0)
    top_n_cutoff = stats.get('top_n_cutoff', 2)
    
    # Plot cutoff line (actual selection threshold)
    if cutoff > 0:
        ax.axhline(y=cutoff, color='red', linestyle='-', linewidth=2.5, 
                   label=f'Cutoff (Top{top_n_cutoff}={cutoff:.4f}, n={n_selected})')
    
    # Plot statistics lines
    if stats:
        # Median
        if stats.get('median', 0) > 0:
            ax.axhline(y=stats['median'], color='green', linestyle='--', linewidth=1.5, 
                       label=f'Median ({stats["median"]:.4f})')
        # Mean
        if stats.get('mean', 0) > 0:
            ax.axhline(y=stats['mean'], color='orange', linestyle='--', linewidth=1.5, 
                       label=f'Mean ({stats["mean"]:.4f})')
        # 99% percentile
        if stats.get('p99', 0) > 0:
            ax.axhline(y=stats['p99'], color='purple', linestyle=':', linewidth=1.5, 
                       label=f'99% ({stats["p99"]:.4f})')
    
    # Find and mark elbow for visualization
    elbow_idx, elbow_val = find_elbow_for_visualization(y)
    if 0 < elbow_idx < n:
        # BUG#12: x-axis is 1-based (range(1, n+1)), y is 0-based list;
        # use y[elbow_idx] so circle and text label are on the same point.
        elbow_y = y[elbow_idx] if elbow_idx < len(y) else y[-1]
        ax.axvline(x=elbow_idx + 1, color='gray', linestyle='--', linewidth=1.5, alpha=0.8)
        ax.plot(elbow_idx + 1, elbow_y, 'o', color='gray', markersize=10, zorder=5)
        ax.annotate(f'Elbow\n(n={elbow_idx + 1}, v={elbow_y:.4f})',
                   xy=(elbow_idx + 1, elbow_y),
                   xytext=(elbow_idx + 1 + max(n * 0.05, 2), elbow_y),
                   fontsize=9, color='gray',
                   arrowprops=dict(arrowstyle='->', color='gray', lw=1))
    
    ax.set_xlabel('Genome Rank', fontsize=11)
    ax.set_ylabel('Containment Score', fontsize=11)
    ax.set_title('Containment Curve with Selection Threshold', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Add detailed info text box
    info_lines = [
        f"Total Genomes: {n}",
        f"Non-zero: {selection_info.get('non_zero_count', 0)}",
        f"",
        f"Selection (Top{top_n_cutoff} Cutoff):",
        f"  Cutoff Value: {cutoff:.6f}",
        f"  Selected: {n_selected} ({n_selected/n*100:.1f}%)",
        f"",
        f"Top1 Statistics:",
        f"  Value: {stats.get('top1_value', 0):.6f}",
        f"  Count: {stats.get('top1_count', 0)} ({stats.get('top1_ratio', 0)*100:.1f}%)",
    ]
    info_text = '\n'.join(info_lines)
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(output_path.replace('.png', '.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def save_containment_results(containment_results, selection_info, output_path, sample_name):
    """Save containment results with selection info.

    Column order follows the scoring logic:
        BaseMetric → Penalty → Score → SharedSNPs → GenomeSNPs → Selected
    Header names adapt to the actual scoring_mode used.
    """
    PRECISION = 6
    cutoff       = selection_info.get('cutoff', 0) if selection_info else 0
    scoring_mode = selection_info.get('scoring_mode',
                       containment_results[0].get('scoring_mode', 'containment')
                       if containment_results else 'containment')

    # Determine base metric column name from scoring_mode
    if scoring_mode in ('jaccard', 'jaccard_weighted'):
        base_col = 'Jaccard'
    elif scoring_mode in ('cosine', 'cosine_weighted'):
        base_col = 'Cosine'
    else:
        base_col = 'Containment'

    has_penalty = scoring_mode in ('containment_weighted', 'cosine_weighted', 'jaccard_weighted')
    # For non-'containment' modes, score differs from base metric
    has_score   = scoring_mode != 'containment'

    def _score(r):
        return r.get('score', r['containment'])

    def _base_metric(r):
        return r.get('base_metric', r['containment'])

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Build header: Rank, Genome, BaseMetric, [Penalty, Score,] SharedSNPs, GenomeSNPs, Selected
        header = ['Rank', 'Genome', base_col]
        if has_penalty:
            header.append('Penalty')
        if has_score:
            header.append('Score')
        header.extend(['SharedSNPs', 'GenomeSNPs', 'Selected'])
        writer.writerow(header)

        for i, r in enumerate(containment_results):
            score  = _score(r)
            is_sel = 'Yes' if round(score, PRECISION) >= cutoff and score > 0 else 'No'
            row = [i + 1, r['genome'], f"{_base_metric(r):.6f}"]
            if has_penalty:
                row.append(f"{r.get('penalty', 'N/A'):.6f}" if isinstance(r.get('penalty'), (int, float)) else 'N/A')
            if has_score:
                row.append(f"{score:.6f}")
            row.extend([r['shared_snps'], r['genome_snps'], is_sel])
            writer.writerow(row)


def save_containment_filter_log(containment_results_all, phylo_info, selection_info,
                                 output_path, sample_name):
    """
    Write a human-readable filter log tracing every genome through two stages.

    Stage 1 — PhyloFilter  (phylo_prefilter_candidates)
    Stage 2 — TopN cutoff  (select_candidates_by_containment)

    Final verdict:
        PASS   → survived both stages, entered LCA/PTC
        FAIL-P → removed by PhyloFilter
        FAIL-T → survived PhyloFilter but below TopN score cutoff
        ZERO   → score = 0, excluded automatically
    """
    PRECISION = 6
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Scoring mode info
    scoring_mode = (selection_info.get('scoring_mode')
                    or (containment_results_all[0].get('scoring_mode', 'containment')
                        if containment_results_all else 'containment'))
    show_score   = scoring_mode not in ('containment',)

    # Dynamic base metric column name
    if scoring_mode in ('jaccard', 'jaccard_weighted'):
        base_col = 'Jaccard'
    elif scoring_mode in ('cosine', 'cosine_weighted'):
        base_col = 'Cosine'
    else:
        base_col = 'Containment'

    def _score(r):
        return r.get('score', r['containment'])

    # Build PhyloFilter removal lookup
    phylo_removed_set = {
        r['genome']: r['shared_ratio']
        for r in phylo_info.get('removed_genomes', [])
    } if phylo_info.get('activated') else {}

    topn_cutoff = selection_info.get('cutoff', 0) if selection_info else 0

    with open(output_path, 'w', encoding='utf-8') as f:
        # ── Header ────────────────────────────────────────────────────────────
        f.write("=" * 72 + "\n")
        f.write(f"  Containment Filter Log\n")
        f.write(f"  Sample     : {sample_name}\n")
        f.write(f"  Generated  : {now}\n")
        f.write(f"  ScoringMode: {scoring_mode}")
        if scoring_mode in ('containment_weighted', 'cosine_weighted', 'jaccard_weighted') and containment_results_all:
            ref_size = containment_results_all[0].get('ref_size', 'N/A')
            alpha    = selection_info.get('cont_alpha', 'N/A') if selection_info else 'N/A'
            f.write(f"  (alpha={alpha}, ref_size={ref_size})")
        f.write("\n")
        f.write("=" * 72 + "\n\n")

        # ── Stage 1: PhyloFilter summary ──────────────────────────────────────
        f.write("━" * 72 + "\n")
        f.write("STAGE 1 — PhyloFilter (phylo_prefilter_candidates)\n")
        f.write("━" * 72 + "\n")

        if not phylo_info.get('activated'):
            f.write(f"  Status : NOT ACTIVATED\n")
            f.write(f"  Reason : {phylo_info.get('reason', 'unknown')}\n")
            if 'gap_ratio' in phylo_info:
                f.write(f"  gap_ratio computed : {phylo_info['gap_ratio']:.4f}  "
                        f"(threshold={phylo_info.get('phylo_gap_ratio', 'N/A')})\n")
        else:
            f.write(f"  Status             : ACTIVATED\n")
            f.write(f"  top1_genome        : {phylo_info.get('top1_genome', 'N/A')}\n")
            f.write(f"  top1_cont          : {phylo_info.get('top1_cont', 0):.6f}\n")
            f.write(f"  level2_cont        : {phylo_info.get('level2_cont', 0):.6f}\n")
            f.write(f"  gap_ratio          : {phylo_info.get('gap_ratio', 0):.4f}  "
                    f">= threshold {phylo_info.get('phylo_gap_ratio', 1.4)} → activate\n")
            f.write(f"  top1_lineage_depth : {phylo_info.get('top1_depth', 0)}\n")
            f.write(f"  phylo_share_thresh : {phylo_info.get('phylo_share_thresh', 0.7)}\n")
            f.write(f"  Genomes before     : {phylo_info.get('n_before', 0)}\n")
            f.write(f"  Genomes removed    : {phylo_info.get('n_removed', 0)}\n")
            f.write(f"  Genomes after      : {phylo_info.get('n_after', 0)}\n")

            if phylo_removed_set:
                f.write(f"\n  Removed genomes detail:\n")
                f.write(f"  {'Genome':<30}  {'SharedRatio':>12}  {'Threshold':>10}  Decision\n")
                f.write(f"  {'-'*30}  {'-'*12}  {'-'*10}  --------\n")
                thresh = phylo_info.get('phylo_share_thresh', 0.7)
                for genome, ratio in sorted(phylo_removed_set.items(),
                                            key=lambda x: x[1], reverse=True):
                    f.write(f"  {genome:<30}  {ratio:>12.4f}  {thresh:>10.2f}  "
                            f"REMOVE (ratio < thresh)\n")
        f.write("\n")

        # ── Stage 2: TopN selection summary ───────────────────────────────────
        f.write("━" * 72 + "\n")
        f.write("STAGE 2 — TopN Cutoff (select_candidates_by_containment)\n")
        f.write("━" * 72 + "\n")
        if selection_info:
            stats = selection_info.get('statistics', {})
            orig_n     = selection_info.get('original_top_n', stats.get('top_n_cutoff', 'N/A'))
            eff_n      = selection_info.get('effective_top_n', stats.get('top_n_cutoff', 'N/A'))
            expanded   = (orig_n != eff_n)
            f.write(f"  scoring_mode     : {scoring_mode}\n")
            f.write(f"  top_n_cutoff     : {orig_n}"
                    + (f" → effective: {eff_n} (expanded by min_candidates)" if expanded else "") + "\n")
            f.write(f"  score_cutoff     : {topn_cutoff:.6f}"
                    + (f"  [NOTE: cutoff is on {scoring_mode} score, not raw containment]"
                       if show_score else "") + "\n")
            f.write(f"  min_candidates   : {selection_info.get('min_candidates', 'N/A')}\n")
            n_input = phylo_info.get('n_after', len(containment_results_all))
            f.write(f"  Input genomes    : {n_input}  (post-PhyloFilter)\n")
            f.write(f"  Selected (PASS)  : {selection_info.get('n_selected', 0)}\n")
            # Hard filter info
            scr = selection_info.get('score_cutoff_ratio', 0.5)
            hcv = selection_info.get('hard_cutoff_value', 0)
            hfc = selection_info.get('hard_filtered_count', 0)
            nbh = selection_info.get('n_before_hard_filter', 0)
            f.write(f"\n  Hard Score Cutoff (score_cutoff_ratio={scr}):\n")
            f.write(f"    top1_score × ratio = {hcv:.6f}\n")
            f.write(f"    Before hard filter : {nbh}\n")
            f.write(f"    Removed by hard flt: {hfc}\n")
            f.write(f"    After hard filter  : {selection_info.get('n_selected', 0)}\n")
        f.write("\n")

        # ── Per-genome decision table ──────────────────────────────────────────
        f.write("━" * 72 + "\n")
        f.write("PER-GENOME DECISION TABLE  (all original genomes)\n")
        f.write("━" * 72 + "\n")

        col_genome = 30
        col_cont   = 12
        col_score  = 12
        col_shared = 10
        col_gsnp   = 10
        col_s1     = 8
        col_s2     = 8
        col_final  = 8

        if show_score:
            hdr = (f"  {'Rank':>4}  {'Genome':<{col_genome}}  {base_col:>{col_cont}}"
                   f"  {'Score':>{col_score}}  {'SharedSNPs':>{col_shared}}  {'GenomeSNPs':>{col_gsnp}}"
                   f"  {'Stage1':^{col_s1}}  {'Stage2':^{col_s2}}  {'Final':^{col_final}}  Reason")
        else:
            hdr = (f"  {'Rank':>4}  {'Genome':<{col_genome}}  {base_col:>{col_cont}}"
                   f"  {'SharedSNPs':>{col_shared}}  {'GenomeSNPs':>{col_gsnp}}"
                   f"  {'Stage1':^{col_s1}}  {'Stage2':^{col_s2}}  {'Final':^{col_final}}  Reason")
        f.write(hdr + "\n")
        f.write("  " + "-" * (len(hdr) - 2) + "\n")

        phylo_active = phylo_info.get('activated', False)
        thresh_share = phylo_info.get('phylo_share_thresh', 0.7)

        for i, r in enumerate(containment_results_all):
            genome = r['genome']
            bm     = r.get('base_metric', r['containment'])
            score  = _score(r)
            sv     = round(score, PRECISION)

            if score == 0 and bm == 0:
                s1, s2, final = "ZERO", "N/A", "ZERO"
                reason = "containment=0 / score=0, excluded"
            elif score == 0:
                # score=0 can happen if penalty collapses it (e.g. huge genome)
                s1, s2, final = "PASS", "FAIL", "FAIL-T"
                reason = f"score=0 after {scoring_mode} scoring; below cutoff"
            elif not phylo_active:
                s1 = "PASS"
                if sv >= topn_cutoff:
                    s2, final = "PASS", "PASS"
                    reason = f"PhyloFilter not activated; score={score:.6f} >= cutoff={topn_cutoff:.6f}"
                else:
                    s2, final = "FAIL", "FAIL-T"
                    reason = f"PhyloFilter not activated; score={score:.6f} < cutoff={topn_cutoff:.6f}"
            elif genome in phylo_removed_set:
                ratio = phylo_removed_set[genome]
                s1, s2, final = "FAIL", "N/A", "FAIL-P"
                reason = f"PhyloFilter removed: shared_ratio={ratio:.4f} < thresh={thresh_share:.2f}"
            else:
                s1 = "PASS"
                if sv >= topn_cutoff:
                    s2, final = "PASS", "PASS"
                    reason = (f"PhyloFilter kept (ratio>={thresh_share:.2f}); "
                              f"score={score:.6f} >= cutoff={topn_cutoff:.6f}")
                else:
                    s2, final = "FAIL", "FAIL-T"
                    reason = (f"PhyloFilter kept; "
                              f"score={score:.6f} < cutoff={topn_cutoff:.6f}")

            if show_score:
                f.write(f"  {i+1:>4}  {genome:<{col_genome}}  {bm:>{col_cont}.6f}"
                        f"  {score:>{col_score}.6f}"
                        f"  {r['shared_snps']:>{col_shared}}  {r['genome_snps']:>{col_gsnp}}"
                        f"  {s1:^{col_s1}}  {s2:^{col_s2}}  {final:^{col_final}}  {reason}\n")
            else:
                f.write(f"  {i+1:>4}  {genome:<{col_genome}}  {bm:>{col_cont}.6f}"
                        f"  {r['shared_snps']:>{col_shared}}  {r['genome_snps']:>{col_gsnp}}"
                        f"  {s1:^{col_s1}}  {s2:^{col_s2}}  {final:^{col_final}}  {reason}\n")

        # ── Summary counts ─────────────────────────────────────────────────────
        finals = []
        for r in containment_results_all:
            genome = r['genome']
            score  = _score(r)
            sv     = round(score, PRECISION)
            if score == 0 and r.get('base_metric', r['containment']) == 0:
                finals.append("ZERO")
            elif score == 0:
                finals.append("FAIL-T")
            elif not phylo_active:
                finals.append("PASS" if sv >= topn_cutoff else "FAIL-T")
            elif genome in phylo_removed_set:
                finals.append("FAIL-P")
            else:
                finals.append("PASS" if sv >= topn_cutoff else "FAIL-T")

        f.write("\n")
        f.write("━" * 72 + "\n")
        f.write("SUMMARY\n")
        f.write("━" * 72 + "\n")
        f.write(f"  Total genomes evaluated : {len(containment_results_all)}\n")
        f.write(f"  PASS  (→ LCA/PTC input) : {finals.count('PASS')}\n")
        f.write(f"  FAIL-P (PhyloFilter)    : {finals.count('FAIL-P')}\n")
        f.write(f"  FAIL-T (TopN cutoff)    : {finals.count('FAIL-T')}\n")
        f.write(f"  ZERO  (score=0)         : {finals.count('ZERO')}\n")
        f.write("=" * 72 + "\n")


# =============================================================================
# Phylogenetic Pre-filter for Containment Candidates
# =============================================================================

def phylo_prefilter_candidates(containment_results, original_genome_to_lineage,
                                phylo_gap_ratio=999999999999, phylo_share_thresh=0.7):
    """
    Phylogenetic pre-filter applied to the full containment-sorted list BEFORE
    top_n_cutoff selection.  Addresses Type-3 errors where wrong genomes from a
    distant branch inflate PTC scores through their large number of private nodes.

    Activation condition (gap check):
        top1_cont / level2_cont >= phylo_gap_ratio   (default 1.4)
    Meaning: top-1 genome's containment is at least 40 % higher than the 2nd
    distinct level, suggesting it is genuinely closer to the sample than all
    others → use its lineage as the phylogenetic reference.

    Filter condition (per candidate):
        shared_prefix_len(candidate_lineage, top1_lineage) / top1_depth
            < phylo_share_thresh   (default 0.7)
    Candidates that diverge from top-1 before reaching 70 % of top-1's depth
    are considered distant relatives and removed.

    Top-1-level genomes (same containment as top-1) are always kept.
    Genomes with no lineage information are kept (fail-safe).

    Args:
        containment_results : list of dicts sorted by containment descending
        original_genome_to_lineage : {genome_id: [node_id, ...]} lineage lists
        phylo_gap_ratio  : activation threshold for top1/level2 ratio (default 1.4)
        phylo_share_thresh : minimum shared-lineage fraction to keep (default 0.7)

    Returns:
        (filtered_results, prefilter_info)
        prefilter_info keys:
            activated     : bool
            reason        : str (explanation when not activated)
            gap_ratio     : float
            top1_genome   : str
            top1_depth    : int
            n_before      : int
            n_after       : int
            n_removed     : int
            removed_genomes : list of {genome, shared_ratio}
    """
    if len(containment_results) < 2:
        return containment_results, {
            'activated': False,
            'reason': 'fewer than 2 candidates, skip phylo prefilter'
        }

    PRECISION = 6

    # ── 统一使用 score 字段（containment_results 已按 score 降序排列）──
    def _score(r):
        return r.get('score', r['containment'])

    # Collect distinct score levels (positive values only)
    rounded_vals = [round(_score(r), PRECISION) for r in containment_results]
    positive_distinct = sorted(set(v for v in rounded_vals if v > 0), reverse=True)

    if len(positive_distinct) < 2:
        return containment_results, {
            'activated': False,
            'reason': 'only one distinct score level, skip phylo prefilter'
        }

    top1_score   = positive_distinct[0]
    level2_score = positive_distinct[1]

    gap_ratio = top1_score / level2_score  # level2_score > 0 guaranteed

    if gap_ratio < phylo_gap_ratio:
        return containment_results, {
            'activated': False,
            'reason': (f'gap_ratio={gap_ratio:.4f} < phylo_gap_ratio={phylo_gap_ratio}, '
                       'top-1 advantage not significant enough'),
            'gap_ratio': gap_ratio
        }

    # --- Activation confirmed ---
    top1_genome  = containment_results[0]['genome']
    top1_lineage = original_genome_to_lineage.get(top1_genome, [])
    top1_depth   = len(top1_lineage)

    if top1_depth == 0:
        return containment_results, {
            'activated': False,
            'reason': f'top-1 genome {top1_genome} has no lineage entry, skip phylo prefilter',
            'gap_ratio': gap_ratio
        }

    def shared_prefix_len(lin_a, lin_b):
        """Length of the common prefix (root-to-leaf overlap) between two lineage lists."""
        count = 0
        for a, b in zip(lin_a, lin_b):
            if a == b:
                count += 1
            else:
                break
        return count

    filtered = []
    removed  = []

    for r in containment_results:
        genome  = r['genome']
        r_val   = round(_score(r), PRECISION)

        # Always keep top-1 level genomes (by score, not raw containment)
        if r_val == top1_score:
            filtered.append(r)
            continue

        lineage = original_genome_to_lineage.get(genome, [])
        if not lineage:
            # No lineage info — keep to be safe
            filtered.append(r)
            continue

        sp_len = shared_prefix_len(top1_lineage, lineage)
        ratio  = sp_len / top1_depth

        if ratio >= phylo_share_thresh:
            filtered.append(r)
        else:
            removed.append({'genome': genome, 'shared_ratio': round(ratio, 4)})

    prefilter_info = {
        'activated':          True,
        'gap_ratio':          round(gap_ratio, 4),
        'top1_genome':        top1_genome,
        'top1_depth':         top1_depth,
        'top1_score':         top1_score,    # score-based (was top1_cont)
        'level2_score':       level2_score,  # score-based (was level2_cont)
        # backward-compat aliases kept for downstream log readers
        'top1_cont':          top1_score,
        'level2_cont':        level2_score,
        'phylo_gap_ratio':    phylo_gap_ratio,
        'phylo_share_thresh': phylo_share_thresh,
        'n_before':           len(containment_results),
        'n_after':            len(filtered),
        'n_removed':          len(removed),
        'removed_genomes':    removed,
    }

    return filtered, prefilter_info


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_samples_from_tsv(filepath):
    samples = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                sample_name = parts[0].strip()
                snps = set(s.strip() for s in parts[1].split(',') if s.strip())
                if sample_name and snps:
                    samples.append((sample_name, snps))
    return samples


def load_lca_output(db_dir):
    filepath = os.path.join(db_dir, 'LCA.output.tsv')
    snp_to_node = {}
    combined_snps = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            snp = row['SNP']
            node_id_str = row['LCA_NodeID']
            if node_id_str == 'None':
                continue
            node_id = int(node_id_str)
            if snp.startswith('ent-'):
                snp = snp[4:]
            snp_to_node[snp] = node_id
            if '+' in snp:
                combined_snps[snp] = node_id
    return snp_to_node, combined_snps


def load_snp_lineage(db_dir):
    filepath = os.path.join(db_dir, 'snp_lingage.tsv')
    node_data = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            node_id = int(row['NodeID'])
            node_data[node_id] = {
                'genome': row['Genome'],
                'lineage': [int(x) for x in row['Lineage'].split(';')],
                'depth': int(row['TreeDepth'])
            }
    return node_data


def load_snp_lineage_with_genome_map(db_dir):
    filepath = os.path.join(db_dir, 'snp_lingage.tsv')
    node_data = {}
    genome_to_lineage = {}
    genome_to_node = {}
    max_depth = 0
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            node_id = int(row['NodeID'])
            genome = row['Genome']
            lineage = [int(x) for x in row['Lineage'].split(';')]
            depth = int(row['TreeDepth'])
            max_depth = max(max_depth, depth)
            node_data[node_id] = {'genome': genome, 'lineage': lineage, 'depth': depth}
            if genome.startswith(("GCF", "GCA")):
                genome_to_lineage[genome] = lineage
                genome_to_node[genome] = node_id
            genome_to_lineage[str(node_id)] = lineage
    
    return node_data, genome_to_lineage, genome_to_node, max_depth


def load_node_dmp(db_dir):
    filepath = os.path.join(db_dir, 'node.tsv')
    node_dmp = {}
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            if len(row) >= 4:
                node_id = int(row[0])
                parent_id = None if row[1] == 'NULL' else int(row[1])
                node_dmp[node_id] = {'parent': parent_id, 'depth': int(row[2]), 'name': row[3]}
    return node_dmp


def load_lineage_stat(db_dir):
    filepath = os.path.join(db_dir, 'lingage.stat.tsv')
    lineage_stat = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            node_id = row['NodeID']
            count = int(row['GenomeCount'])
            genomes = row['GenomeNames'].strip(';').split(';') if row['GenomeNames'] else []
            lineage_stat[node_id] = {'count': count, 'genomes': genomes}
    return lineage_stat


def load_tree_branch_lengths(db_dir, node_dmp):
    try:
        import dendropy
        filepath = os.path.join(db_dir, 'numbered_tree.tree')
        if not os.path.exists(filepath):
            return {}
        tree = dendropy.Tree.get(path=filepath, schema="newick")
        branch_lengths = {}
        for node in tree.postorder_node_iter():
            if node.label:
                bl = node.edge_length if node.edge_length is not None else 0.0
                branch_lengths[node.label] = bl if bl > 0 else 1e-6
        return branch_lengths
    except Exception as e:
        warnings.warn(f"load_tree_branch_lengths: {type(e).__name__}: {e}")
        return {}


# =============================================================================
# Sub-database Construction (unchanged)
# =============================================================================

def filter_snp_matrix(total_var_file, genome_list_file, output_file, sample_snps=None):
    with open(genome_list_file, 'r') as f:
        keep_genomes = set(line.strip() for line in f if line.strip())
    included_snps = set()
    with open(total_var_file, 'r') as fin:
        header_line = fin.readline().strip()
        headers = header_line.split('\t')
        keep_indices = [0]
        kept_headers = [headers[0]]
        for i, col in enumerate(headers[1:], 1):
            if col in keep_genomes:
                keep_indices.append(i)
                kept_headers.append(col)
        if len(keep_indices) <= 1:
            raise ValueError("No specified samples found")
        with open(output_file, 'w') as fout:
            fout.write('\t'.join(kept_headers) + '\n')
            for line in fin:
                fields = line.strip().split('\t')
                if len(fields) < max(keep_indices) + 1:
                    continue
                snp_name = fields[0]
                has_nonzero = any(float(fields[idx]) != 0 for idx in keep_indices[1:] if fields[idx].replace('.','',1).replace('-','',1).isdigit())
                if has_nonzero or (sample_snps and snp_name in sample_snps):
                    filtered = [fields[i] for i in keep_indices]
                    fout.write('\t'.join(filtered) + '\n')
                    included_snps.add(snp_name)
    return included_snps


def prune_tree_to_samples(input_tree_path, samples_path, output_tree_path):
    import dendropy
    with open(samples_path, 'r') as f:
        retain_labels = {line.strip().replace(' ', '_') for line in f if line.strip()}
    tree = dendropy.Tree.get(path=input_tree_path, schema="newick")
    current_labels = {t.label.replace(' ', '_') for t in tree.taxon_namespace if t.label}
    to_retain = retain_labels & current_labels
    if not to_retain:
        raise ValueError("No matching taxa")
    norm_to_original = {t.label.replace(' ', '_'): t.label for t in tree.taxon_namespace if t.label}
    original_to_retain = {norm_to_original[n] for n in to_retain}
    tree.retain_taxa_with_labels(original_to_retain)
    tree.write(path=output_tree_path, schema="newick")


def number_tree_nodes(input_tree_file, output_tree_file, output_node_file):
    """
    Read a Newick tree, assign sequential integer IDs to every node
    (root first, DFS pre-order), write the re-labelled tree and a node
    table.  Uses dendropy exclusively (BUG#9: ete3 removed).
    """
    import dendropy
    tree = dendropy.Tree.get(path=input_tree_file, schema="newick",
                             preserve_underscores=True)

    node_counter = 0
    node_data = []   # (node_id, parent_id, depth, original_name)

    # DFS pre-order — mirrors the old ete3 stack-based traversal
    for node in tree.preorder_node_iter():
        node_counter += 1
        original_name = node.label if node.label else (
            node.taxon.label if node.taxon else f"Internal_{node_counter}"
        )
        # Store id on the node object so children can retrieve parent id
        node._num_id = node_counter
        parent_id = node.parent_node._num_id if node.parent_node is not None else None
        depth = 0
        p = node.parent_node
        while p is not None:
            depth += 1
            p = p.parent_node
        node_data.append((node_counter, parent_id, depth, original_name))
        # Relabel: internal nodes use node.label, leaves use taxon.label
        if node.taxon is not None:
            node.taxon.label = str(node_counter)
        else:
            node.label = str(node_counter)

    tree.write(path=output_tree_file, schema="newick")

    with open(output_node_file, 'w') as f:
        f.write("NodeID\tParentID\tDepth\tOriginalName\n")
        for node_id, parent_id, depth, original_name in node_data:
            parent = str(parent_id) if parent_id is not None else "NULL"
            f.write(f"{node_id}\t{parent}\t{depth}\t{original_name}\n")


def make_lineage(node_dmp_file, output_file):
    nodes = {}
    with open(node_dmp_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                nodes[parts[0]] = {'parent': parts[1] if parts[1] != 'NULL' else 'root', 'name': parts[3]}
    def get_lineage_path(node_id):
        path = []
        current = node_id
        while current != 'root' and current in nodes:
            path.append(current)
            current = nodes[current]['parent']
        path.reverse()
        return path
    output_data = []
    for node_id, info in nodes.items():
        if not str(node_id).isdigit():
            continue
        path = get_lineage_path(node_id)
        genome = info['name'] if info['name'].startswith(("GCF", "GCA")) else node_id
        output_data.append((int(node_id), genome, ';'.join(path), len(path)))
    output_data.sort(key=lambda x: x[0])
    with open(output_file, 'w') as f:
        f.write("NodeID\tGenome\tLineage\tTreeDepth\n")
        for node_id, genome, lineage, depth in output_data:
            f.write(f"{node_id}\t{genome}\t{lineage}\t{depth}\n")


def make_lineage_stat(snp_lineage_file, output_file):
    data = []
    with open(snp_lineage_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            data.append({'NodeID': row['NodeID'], 'Genome': row['Genome'], 'Lineage': row['Lineage'], 'TreeDepth': int(row['TreeDepth'])})
    genomes = [d for d in data if d['Genome'].startswith(("GCF", "GCA"))]
    total_genomes = len(genomes)
    results = []
    for item in data:
        node_lineage = item['Lineage']
        lineage_prefix = node_lineage + ';'
        node_genomes = [g for g in genomes if g['Lineage'] == node_lineage or g['Lineage'].startswith(lineage_prefix)]
        count = len(node_genomes)
        ratio = count / total_genomes if total_genomes > 0 else 0
        names = ';'.join([g['Genome'] for g in node_genomes])
        results.append({'NodeID': item['NodeID'], 'GenomeCount': count, 'GenomeRatio': ratio, 'GenomeNames': names})
    with open(output_file, 'w') as f:
        f.write("NodeID\tGenomeCount\tGenomeRatio\tGenomeNames\n")
        for r in results:
            f.write(f"{r['NodeID']}\t{r['GenomeCount']}\t{r['GenomeRatio']}\t{r['GenomeNames']}\n")


# =============================================================================
# SNP LCA Analysis (unchanged)
# =============================================================================

class TreeNode:
    def __init__(self, nodeID, parentID, depth, name=None):
        self.nodeID = nodeID
        self.parentID = parentID
        self.depth = depth
        self.name = name


class LineageInfo:
    def __init__(self, nodeID, lineage, tree_depth):
        self.nodeID = nodeID
        self.lineage = lineage
        self.tree_depth = tree_depth


def load_node_file_for_lca(node_file):
    nodes = {}
    nodeID_to_name = {}
    with open(node_file, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                nodeID = int(parts[0])
                parentID = 0 if parts[1] == "NULL" else int(parts[1])
                nodes[nodeID] = TreeNode(nodeID, parentID, int(parts[2]), parts[3])
                nodeID_to_name[nodeID] = parts[3]
    return nodes, nodeID_to_name


def load_lineage_file_for_lca(lineage_file, nodes, nodeID_to_name):
    lineages = {}
    genome_to_node = {}
    node_to_genome = {}
    with open(lineage_file, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                nodeID = int(parts[0])
                genome = parts[1]
                lineage_str = parts[2]
                tree_depth = int(parts[3])
                lineages[nodeID] = LineageInfo(nodeID, lineage_str, tree_depth)
                node_to_genome[nodeID] = genome
                lineage_list = [int(x) for x in lineage_str.split(';')]
                if nodeID not in nodes:
                    parentID = 0 if len(lineage_list) == 1 else lineage_list[-2]
                    name = genome if genome.startswith(("GCF", "GCA")) else str(nodeID)
                    nodes[nodeID] = TreeNode(nodeID, parentID, len(lineage_list) - 1, name)
                    nodeID_to_name[nodeID] = name
                if genome.startswith(("GCF", "GCA")):
                    genome_to_node[genome] = nodeID
                    nodeID_to_name[nodeID] = genome
    return lineages, genome_to_node, node_to_genome


def precompute_binary_lifting(nodes, max_log=20):
    max_id = max(nodes.keys() or [0]) + 1
    parent = [[0] * max_id for _ in range(max_log)]
    depth = [0] * max_id
    for node_id in nodes:
        parent[0][node_id] = nodes[node_id].parentID
        depth[node_id] = nodes[node_id].depth
    for k in range(1, max_log):
        for node_id in nodes:
            if parent[k-1][node_id]:
                parent[k][node_id] = parent[k-1][parent[k-1][node_id]]
    return parent, depth


def find_lca_two_nodes(nodeA, nodeB, parent, depth, max_log=20):
    if depth[nodeA] > depth[nodeB]:
        nodeA, nodeB = nodeB, nodeA
    diff = depth[nodeB] - depth[nodeA]
    k = 0
    while diff:
        if diff & 1:
            nodeB = parent[k][nodeB]
        diff >>= 1
        k += 1
        if k >= max_log:
            return None
    if nodeA == nodeB:
        return nodeA
    for k in range(max_log - 1, -1, -1):
        if parent[k][nodeA] != parent[k][nodeB]:
            nodeA = parent[k][nodeA]
            nodeB = parent[k][nodeB]
    return parent[0][nodeA]


def find_group_lca(genome_nodes, parent, depth, max_log):
    if not genome_nodes:
        return None
    lca = genome_nodes[0]
    for nodeID in genome_nodes[1:]:
        lca = find_lca_two_nodes(lca, nodeID, parent, depth, max_log)
        if lca is None:
            break
    return lca


def precompute_total_genomes(root, children, total_genomes, node_to_genome, nodes):
    if root not in nodes:
        return
    stack = [(root, 0)]
    while stack:
        node, child_idx = stack[-1]
        if child_idx == len(children[node]):
            t = 0
            if children[node]:
                for child in children[node]:
                    t += total_genomes.get(child, 0)
            else:
                genome = node_to_genome.get(node, None)
                if genome and genome.startswith(("GCF", "GCA")):
                    t = 1
            total_genomes[node] = t
            stack.pop()
            continue
        child = children[node][child_idx]
        stack[-1] = (node, child_idx + 1)
        stack.append((child, 0))


def compute_snp_subtree_counts(node_id, children, snp_genomes, input_genomes, node_to_genome, nodes):
    if node_id not in nodes:
        return
    stack = [(node_id, 0)]
    while stack:
        node, child_idx = stack[-1]
        if child_idx == len(children[node]):
            s = 0
            if children[node]:
                for child in children[node]:
                    s += snp_genomes.get(child, 0)
            else:
                genome = node_to_genome.get(node, None)
                if genome and genome.startswith(("GCF", "GCA")):
                    s = 1 if genome in input_genomes else 0
            snp_genomes[node] = s
            stack.pop()
            continue
        child = children[node][child_idx]
        stack[-1] = (node, child_idx + 1)
        stack.append((child, 0))


def correct_lca(lca_node, children, total_genomes, snp_genomes, frequence_cutoff, frequency_cutoff, difference_multiplier, nodes):
    current = lca_node
    while True:
        childs = children.get(current, [])
        if not childs:
            break
        child_data = []
        for child in childs:
            freq = snp_genomes.get(child, 0)
            total = total_genomes.get(child, 0)
            rate = freq / total if total > 0 else 0.0
            child_data.append((child, freq, total, rate))
        child_data.sort(key=lambda x: x[3], reverse=True)
        if len(child_data) == 1:
            top_child, top_freq, top_total, top_rate = child_data[0]
            if top_freq >= frequence_cutoff and top_rate >= frequency_cutoff:
                current = top_child
                continue
            else:
                break
        else:
            top_child, top_freq, top_total, top_rate = child_data[0]
            all_meet_diff = True
            for other_child, other_freq, other_total, other_rate in child_data[1:]:
                rate_diff = top_rate / other_rate if other_rate > 0 else float('inf') if top_rate > 0 else 0.0
                freq_diff = top_freq / other_freq if other_freq > 0 else float('inf') if top_freq > 0 else 0.0
                if rate_diff < difference_multiplier or freq_diff < difference_multiplier:
                    all_meet_diff = False
                    break
            if all_meet_diff and top_freq >= frequence_cutoff and top_rate >= frequency_cutoff:
                current = top_child
                continue
            break
    return current


def run_snp_lca_analysis(matrix_file, node_dmp_file, snp_lineage_file, output_file, difference_multiplier=5000000):
    nodes, nodeID_to_name = load_node_file_for_lca(node_dmp_file)
    lineages, genome_to_node, node_to_genome = load_lineage_file_for_lca(snp_lineage_file, nodes, nodeID_to_name)
    children = defaultdict(list)
    for node_id, node in nodes.items():
        if node.parentID != 0:
            children[node.parentID].append(node_id)
    max_depth = max([node.depth for node in nodes.values()] if nodes else [0])
    max_log = int(math.log2(max_depth)) + 2 if max_depth > 0 else 1
    parent_table, depth = precompute_binary_lifting(nodes, max_log=max_log)
    total_genomes_map = {}
    precompute_total_genomes(1, children, total_genomes_map, node_to_genome, nodes)
    with open(output_file, "w") as output_f:
        output_f.write("SNP\tLCA_NodeID\tGenome\tLineage\tTreeDepth\n")
        with open(matrix_file, "r") as matrix_f:
            reader = csv.reader(matrix_f, delimiter="\t")
            headers = next(reader)
            genomes = headers[1:]
            for row in reader:
                if not row:
                    continue
                snp = row[0]
                input_genomes = set()
                for i in range(1, len(row)):
                    if row[i] == "1":
                        input_genomes.add(genomes[i - 1])
                genome_nodes = [genome_to_node[g] for g in input_genomes if g in genome_to_node]
                lca_node = find_group_lca(genome_nodes, parent_table, depth, max_log)
                corrected_lca = lca_node
                if lca_node is not None and genome_nodes:
                    snp_genomes = {}
                    compute_snp_subtree_counts(lca_node, children, snp_genomes, input_genomes, node_to_genome, nodes)
                    corrected_lca = correct_lca(lca_node, children, total_genomes_map, snp_genomes, 0, 0, difference_multiplier, nodes)
                if corrected_lca is None:
                    output_f.write(f"{snp}\tNone\tNone\tNone\tNone\n")
                else:
                    genome_name = node_to_genome.get(corrected_lca, "Unknown")
                    lineage_info = lineages.get(corrected_lca)
                    lineage_path = lineage_info.lineage if lineage_info else "Unknown"
                    tree_depth = lineage_info.tree_depth if lineage_info else "Unknown"
                    output_f.write(f"{snp}\t{corrected_lca}\t{genome_name}\t{lineage_path}\t{tree_depth}\n")


def build_sub_database(genome_list_file, total_var_file, jaccard_nwk_file, output_dir, sample_snps=None):
    os.makedirs(output_dir, exist_ok=True)
    filtered_matrix = os.path.join(output_dir, 'total.var.select.xls')
    included_snps = filter_snp_matrix(total_var_file, genome_list_file, filtered_matrix, sample_snps)
    pruned_tree = os.path.join(output_dir, 'Jaccard.fastme.input.select.nwk')
    prune_tree_to_samples(jaccard_nwk_file, genome_list_file, pruned_tree)
    numbered_tree = os.path.join(output_dir, 'numbered_tree.tree')
    node_dmp = os.path.join(output_dir, 'node.tsv')
    number_tree_nodes(pruned_tree, numbered_tree, node_dmp)
    snp_lineage = os.path.join(output_dir, 'snp_lingage.tsv')
    make_lineage(node_dmp, snp_lineage)
    lineage_stat = os.path.join(output_dir, 'lingage.stat.tsv')
    make_lineage_stat(snp_lineage, lineage_stat)
    lca_output = os.path.join(output_dir, 'LCA.output.tsv')
    run_snp_lca_analysis(filtered_matrix, node_dmp, snp_lineage, lca_output)
    return included_snps


def get_database_stats(db_dir):
    stats = {'total_genomes': 0, 'total_nodes': 0, 'total_snps': 0, 'max_depth': 0}
    try:
        lineage_stat_file = os.path.join(db_dir, 'lingage.stat.tsv')
        if os.path.exists(lineage_stat_file):
            with open(lineage_stat_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    stats['total_nodes'] += 1
                    if row['GenomeNames'] and row['GenomeNames'].startswith(("GCF", "GCA")) and int(row['GenomeCount']) == 1:
                        stats['total_genomes'] += 1
        lca_file = os.path.join(db_dir, 'LCA.output.tsv')
        if os.path.exists(lca_file):
            with open(lca_file, 'r') as f:
                stats['total_snps'] = sum(1 for _ in f) - 1
        snp_lineage_file = os.path.join(db_dir, 'snp_lingage.tsv')
        if os.path.exists(snp_lineage_file):
            with open(snp_lineage_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    stats['max_depth'] = max(stats['max_depth'], int(row['TreeDepth']))
    except Exception as e:
        print(f"  [WARNING] get_database_stats: {type(e).__name__}: {e}")
    return stats


# =============================================================================
# Core Tracing Algorithm
# =============================================================================

def snp_to_node_mapping(sample_snps, snp_to_node, combined_snps):
    node_snp_data = defaultdict(lambda: {'count': 0, 'snps': []})
    mapped_snps = set()
    snp_detail = []
    for snp in sample_snps:
        if snp in snp_to_node:
            node_id = snp_to_node[snp]
            node_snp_data[node_id]['count'] += 1
            node_snp_data[node_id]['snps'].append(snp)
            mapped_snps.add(snp)
            snp_detail.append({'snp': snp, 'node_id': node_id, 'type': 'single'})
    if combined_snps:
        # O(N) reverse-index: single SNP → set of combo keys that contain it
        # Built once per call; far cheaper than O(N²) pairwise enumeration
        snp_to_combos = defaultdict(set)
        for combo_key in combined_snps:
            for part in combo_key.split('+'):
                snp_to_combos[part].add(combo_key)

        checked = set()
        for snp in sample_snps:
            for combo in snp_to_combos.get(snp, ()):
                if combo in checked:
                    continue
                checked.add(combo)
                # Keep only combos whose every component is in the sample
                if set(combo.split('+')).issubset(sample_snps):
                    node_id = combined_snps[combo]
                    node_snp_data[node_id]['count'] += 1
                    node_snp_data[node_id]['snps'].append(combo)
                    mapped_snps.add(combo)
                    snp_detail.append({'snp': combo, 'node_id': node_id, 'type': 'combined'})
    return {k: v for k, v in node_snp_data.items() if k != 1 and v['count'] > 0}, mapped_snps, snp_detail


def calculate_cumulative_scores(node_snp_data, snp_lineage):
    scores = defaultdict(int)
    for node_id, data in node_snp_data.items():
        if node_id not in snp_lineage:
            continue
        lineage = snp_lineage[node_id]['lineage']
        for ancestor in lineage:
            scores[ancestor] += data['count']
    return {k: v for k, v in scores.items() if v > 0}


def trace_paths(scores, node_dmp, node_snp_data):
    path_data = []
    children_map = defaultdict(list)
    for node_id, info in node_dmp.items():
        if info['parent'] is not None:
            children_map[info['parent']].append(node_id)
    def trace_recursive(current_node, path):
        children = children_map.get(current_node, [])
        valid_children = [c for c in children if scores.get(c, 0) > 0]
        if not valid_children:
            path_str = "->".join(map(str, path))
            path_scores = [scores.get(n, 0) for n in path]
            snp_counts = [node_snp_data.get(n, {}).get('count', 0) for n in path]
            path_data.append((path_str, ",".join(map(str, path_scores)), ",".join(map(str, snp_counts))))
            return
        for child in valid_children:
            trace_recursive(child, path + [child])
    if 1 in scores or any(scores.get(c, 0) > 0 for c in children_map.get(1, [])):
        trace_recursive(1, [1])
    return path_data


def cluster_and_rank_paths(path_data, node_dmp, branch_lengths, jaccard_percent=0.8):
    if not path_data:
        return []
    paths_info = []
    for path_str, scores_str, snp_str in path_data:
        nodes = path_str.split('->')
        scores = [int(s) for s in scores_str.split(',')]
        snp_counts = [int(s) for s in snp_str.split(',')]
        total_snp = sum(snp_counts)
        total_length = sum(branch_lengths.get(node_dmp.get(int(nodes[i]), {}).get('name', nodes[i]), 1e-6) for i in range(1, len(nodes)))
        density = total_snp / total_length if total_length > 0 else 0
        snp_branch = snp_counts[1:] if len(snp_counts) > 1 else []
        zero_ratio = sum(1 for s in snp_branch if s == 0) / len(snp_branch) if snp_branch else 0
        paths_info.append({
            'path': path_str, 'scores_str': scores_str, 'snp_str': snp_str, 'nodes': nodes, 'scores': scores,
            'total_snp': total_snp, 'density': density, 'total_length': total_length, 'zero_ratio': zero_ratio
        })
    n = len(paths_info)
    if n == 1:
        paths_info[0].update({'jaccard_cluster': 0, 'importance_rank': 1, 'snp_rank': 1, 'snp_proportion': 1.0})
        return paths_info
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        nodes_i = paths_info[i]['nodes']
        prefix_len_i = int(len(nodes_i) * jaccard_percent) + 1
        prefix_i = set(nodes_i[:prefix_len_i])
        sim_matrix[i, i] = 1.0
        for j in range(i + 1, n):
            nodes_j = paths_info[j]['nodes']
            prefix_len_j = int(len(nodes_j) * jaccard_percent) + 1
            prefix_j = set(nodes_j[:prefix_len_j])
            inter = len(prefix_i & prefix_j)
            union = len(prefix_i | prefix_j)
            sim = inter / union if union > 0 else 0
            sim_matrix[i, j] = sim
            sim_matrix[j, i] = sim
    dist_matrix = 1 - sim_matrix
    if np.all(dist_matrix == 0):
        jaccard_clusters = np.zeros(n, dtype=int)
    else:
        condensed_dist = squareform(dist_matrix)
        link = linkage(condensed_dist, method='average')
        jaccard_clusters = fcluster(link, t=0.5, criterion='distance') - 1
    for idx, p in enumerate(paths_info):
        p['jaccard_cluster'] = int(jaccard_clusters[idx])
    sort_key = lambda p: tuple(-int(s) for s in p['scores_str'].split(','))
    all_sorted = sorted(paths_info, key=sort_key)
    for i, p in enumerate(all_sorted):
        p['importance_rank'] = i + 1
    sorted_by_snp = sorted(paths_info, key=lambda p: p['total_snp'], reverse=True)
    max_snp = sorted_by_snp[0]['total_snp'] if sorted_by_snp else 0
    for i, p in enumerate(sorted_by_snp):
        p['snp_rank'] = i + 1
    for p in paths_info:
        p['snp_proportion'] = p['total_snp'] / max_snp if max_snp > 0 else 0.0
    return paths_info


def find_tied_best_paths(path_results):
    """
    Find paths with identical scores (tied for best).
    
    Returns list of paths that have the same score as the best path.
    """
    if not path_results:
        return []
    
    # Get best path
    best_path = next((p for p in path_results if p['importance_rank'] == 1), None)
    if not best_path:
        return []
    
    best_scores = best_path['scores_str']
    best_total_snp = best_path['total_snp']
    
    # Find all paths with identical scores
    tied_paths = []
    for p in path_results:
        if p['scores_str'] == best_scores and p['total_snp'] == best_total_snp:
            tied_paths.append(p)
    
    return tied_paths


def find_lca_of_genomes(genomes, original_genome_to_lineage):
    """
    Find the LCA (Lowest Common Ancestor) of multiple genomes.
    
    Returns:
        (lca_node, lca_lineage, lca_depth)
    """
    if not genomes:
        return None, [], 0
    
    # Get lineages for all genomes
    lineages = []
    for g in genomes:
        if g in original_genome_to_lineage:
            lineages.append(original_genome_to_lineage[g])
    
    if not lineages:
        return None, [], 0
    
    if len(lineages) == 1:
        return lineages[0][-1], lineages[0], len(lineages[0])
    
    # Find common prefix (LCA)
    common = []
    min_len = min(len(l) for l in lineages)
    for i in range(min_len):
        nodes_at_i = set(l[i] for l in lineages)
        if len(nodes_at_i) == 1:
            common.append(lineages[0][i])
        else:
            break
    
    if common:
        return common[-1], common, len(common)
    return None, [], 0


def check_trace_result_with_ties(path_results, lineage_stat, original_genome_to_lineage):
    """
    Check trace result, handling tied paths by finding their LCA.
    
    Returns:
        (is_terminal, final_node, genomes, trace_level, lca_info)
        
    trace_level: 'Genome' if traced to single genome, 'Node' if traced to LCA of multiple
    lca_info: dict with LCA details if applicable
    """
    if not path_results:
        return True, None, None, 'Unknown', None
    
    # Find tied best paths
    tied_paths = find_tied_best_paths(path_results)
    
    if len(tied_paths) == 1:
        # Single best path - normal case
        best_path = tied_paths[0]
        nodes = best_path['path'].split('->')
        final_node = nodes[-1]
        
        if final_node in lineage_stat:
            genome_count = lineage_stat[final_node]['count']
            genomes = lineage_stat[final_node]['genomes']
            
            if genome_count == 1:
                return True, final_node, genomes, 'Genome', None
            elif genome_count > 1 and final_node != '1':
                return False, final_node, genomes, 'Node', {
                    'reason': 'multiple_genomes_in_subtree',
                    'genome_count': genome_count
                }
        return True, final_node, None, 'Unknown', None
    
    else:
        # Multiple tied paths - find LCA of all endpoint genomes
        all_endpoint_genomes = []
        for p in tied_paths:
            nodes = p['path'].split('->')
            final_node = nodes[-1]
            if final_node in lineage_stat:
                all_endpoint_genomes.extend(lineage_stat[final_node]['genomes'])
        
        # Remove duplicates
        all_endpoint_genomes = list(set(all_endpoint_genomes))
        
        if len(all_endpoint_genomes) == 0:
            best_path = tied_paths[0]
            nodes = best_path['path'].split('->')
            return True, nodes[-1], None, 'Unknown', None
        
        if len(all_endpoint_genomes) == 1:
            # All paths lead to same genome
            genome = all_endpoint_genomes[0]
            if genome in original_genome_to_lineage:
                lineage = original_genome_to_lineage[genome]
                return True, lineage[-1], [genome], 'Genome', None
            return True, None, [genome], 'Genome', None
        
        # Find LCA of all endpoint genomes
        lca_node, lca_lineage, lca_depth = find_lca_of_genomes(all_endpoint_genomes, original_genome_to_lineage)
        
        lca_info = {
            'reason': 'tied_paths',
            'tied_path_count': len(tied_paths),
            'all_endpoint_genomes': all_endpoint_genomes,
            'lca_node': lca_node,
            'lca_depth': lca_depth
        }
        
        return True, lca_node, all_endpoint_genomes, 'Node', lca_info


def get_genome_global_path(genome_name, original_genome_to_lineage, node_iteration_map=None):
    if not genome_name or genome_name not in original_genome_to_lineage:
        return "N/A", "N/A"
    lineage = original_genome_to_lineage[genome_name]
    if node_iteration_map:
        path_parts = [f"{n}[R{node_iteration_map[n]}]" if n in node_iteration_map else str(n) for n in lineage]
        global_path = "->".join(path_parts)
    else:
        global_path = "->".join(str(n) for n in lineage)
    global_path_plain = "->".join(str(n) for n in lineage)
    return global_path, global_path_plain


def get_node_global_path(node_id, original_genome_to_lineage, node_iteration_map=None):
    """Get GlobalPath for a node (not genome)."""
    # Try to find lineage by node ID
    node_key = str(node_id)
    if node_key in original_genome_to_lineage:
        lineage = original_genome_to_lineage[node_key]
        if node_iteration_map:
            path_parts = [f"{n}[R{node_iteration_map[n]}]" if n in node_iteration_map else str(n) for n in lineage]
            global_path = "->".join(path_parts)
        else:
            global_path = "->".join(str(n) for n in lineage)
        global_path_plain = "->".join(str(n) for n in lineage)
        return global_path, global_path_plain
    return "N/A", "N/A"


def find_common_ancestor_in_original_db(genomes, original_genome_to_lineage):
    if not genomes:
        return None, []
    lineages = [original_genome_to_lineage[g] for g in genomes if g in original_genome_to_lineage]
    if not lineages:
        return None, []
    common = []
    min_len = min(len(l) for l in lineages)
    for i in range(min_len):
        nodes_at_i = set(l[i] for l in lineages)
        if len(nodes_at_i) == 1:
            common.append(lineages[0][i])
        else:
            break
    return (common[-1], common) if common else (None, [])


def compute_global_path_scores(sample_snps, original_db_dir, global_path_plain,
                                 cached_global_scores=None, cached_node_snp_data=None):
    """
    Compute cumulative SNP scores and per-node SNP counts along the global path
    (original database path), using the sample's SNP set.

    Args:
        sample_snps: set of sample SNP names
        original_db_dir: path to original (non-sub) database directory
        global_path_plain: '->'-joined node ID string (e.g. '1->5->12->47')
        cached_global_scores: precomputed {node_id: cumulative_score} dict (optional)
        cached_node_snp_data: precomputed {node_id: {'count': n}} dict (optional)

    Returns:
        (GlobalPathScores: str, GlobalNodeSNPCounts: str)
        Both are comma-separated strings aligned to global_path_plain nodes.
        Returns ('N/A', 'N/A') on any error.
    """
    try:
        global_path_nodes = [int(n) for n in global_path_plain.split('->')
                             if n.strip().isdigit()]
        if not global_path_nodes:
            return 'N/A', 'N/A'

        if cached_global_scores is not None and cached_node_snp_data is not None:
            global_scores  = cached_global_scores
            node_snp_data_ = cached_node_snp_data
        else:
            snp_to_node, combined_snps = load_lca_output(original_db_dir)
            snp_lineage = load_snp_lineage(original_db_dir)
            node_snp_data_, _, _ = snp_to_node_mapping(sample_snps, snp_to_node, combined_snps)
            global_scores = calculate_cumulative_scores(node_snp_data_, snp_lineage)

        score_list   = [global_scores.get(n, 0) for n in global_path_nodes]
        snp_cnt_list = [node_snp_data_.get(n, {}).get('count', 0) for n in global_path_nodes]

        return (','.join(str(s) for s in score_list),
                ','.join(str(c) for c in snp_cnt_list))
    except Exception:
        return 'N/A', 'N/A'



def build_iteration_markers(iteration_history, original_genome_to_lineage):
    node_iteration_map = {}
    for iter_info in iteration_history:
        iter_num = iter_info['iteration']
        iter_genomes = iter_info['genomes']
        if iter_genomes:
            lca_node, _ = find_common_ancestor_in_original_db(iter_genomes, original_genome_to_lineage)
            if lca_node is not None and lca_node not in node_iteration_map:
                node_iteration_map[lca_node] = iter_num
    return node_iteration_map


# =============================================================================
# SNP LCA Detail Output
# =============================================================================

def save_snp_lca_detail(snp_detail, snp_lineage, lineage_stat, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['SNP', 'NodeID', 'SNPType', 'NodeGenome', 'NodeLineage', 'TreeDepth', 'SubtreeGenomes', 'SubtreeGenomeCount'])
        for item in snp_detail:
            snp = item['snp']
            node_id = item['node_id']
            snp_type = item['type']
            node_info = snp_lineage.get(node_id, {})
            node_genome = node_info.get('genome', str(node_id))
            lineage = node_info.get('lineage', [])
            lineage_str = "->".join(str(n) for n in lineage)
            tree_depth = node_info.get('depth', len(lineage))
            stat_info = lineage_stat.get(str(node_id), {})
            subtree_genomes = stat_info.get('genomes', [])
            subtree_count = stat_info.get('count', 0)
            subtree_str = ";".join(subtree_genomes[:10])
            if len(subtree_genomes) > 10:
                subtree_str += f"...(+{len(subtree_genomes)-10})"
            writer.writerow([snp, node_id, snp_type, node_genome, lineage_str, tree_depth, subtree_str, subtree_count])


def save_all_paths_detail(path_results, lineage_stat, original_genome_to_lineage, output_file, iteration=1):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Iteration', 'Rank', 'Path', 'EndpointNode', 'EndpointGenome', 'EndpointGlobalPath',
                        'Scores', 'NodeSNPCounts', 'TotalSNP', 'JaccardCluster', 'SubtreeGenomeCount', 'SubtreeGenomes'])
        for p in sorted(path_results, key=lambda x: x['importance_rank']):
            nodes = p['path'].split('->')
            endpoint_node = nodes[-1]
            stat_info = lineage_stat.get(endpoint_node, {})
            endpoint_genomes = stat_info.get('genomes', [])
            endpoint_genome = endpoint_genomes[0] if endpoint_genomes else endpoint_node
            genome_count = stat_info.get('count', 0)
            global_path, _ = get_genome_global_path(endpoint_genome, original_genome_to_lineage)
            subtree_str = ";".join(endpoint_genomes[:5])
            if len(endpoint_genomes) > 5:
                subtree_str += f"...(+{len(endpoint_genomes)-5})"
            writer.writerow([iteration, p['importance_rank'], p['path'], endpoint_node, endpoint_genome, global_path,
                           p['scores_str'], p['snp_str'], p['total_snp'], p['jaccard_cluster'], genome_count, subtree_str])


# =============================================================================
# Visualization Functions (simplified)
# =============================================================================


def create_iteration_progress_figure(iteration_results, output_dir, sample_name):
    """Create iteration progress figure with only A and B panels."""
    n_iterations = len(iteration_results)
    if n_iterations == 0:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'Tracing Progress: {sample_name}', fontsize=14, fontweight='bold')
    
    iterations = list(range(1, n_iterations + 1))
    genome_counts = [r['genome_count'] for r in iteration_results]
    mapped_snps = [r['mapped_snps'] for r in iteration_results]
    
    # Panel A: Candidate Narrowing
    ax1 = axes[0]
    max_genomes = max(genome_counts) if genome_counts else 1
    for i, (iter_num, gc) in enumerate(zip(iterations, genome_counts)):
        width = gc / max_genomes * 0.8
        left = 0.5 - width/2
        color = '#27AE60' if gc == 1 else '#3498DB' if gc < 10 else '#E74C3C'
        ax1.add_patch(FancyBboxPatch((left, n_iterations - i - 0.8), width, 0.6, 
                                      boxstyle="round,pad=0.02", facecolor=color, 
                                      edgecolor='white', linewidth=2, alpha=0.8))
        ax1.text(0.5, n_iterations - i - 0.5, f'{gc}', fontsize=12, ha='center', 
                va='center', color='white', fontweight='bold')
        ax1.text(-0.1, n_iterations - i - 0.5, f'R{iter_num}', fontsize=9, ha='right', va='center')
    ax1.set_xlim(-0.3, 1.3)
    ax1.set_ylim(-0.5, n_iterations + 0.5)
    ax1.axis('off')
    ax1.set_title('A. Candidate Narrowing', fontsize=11, fontweight='bold', pad=10)
    
    # Panel B: Statistics
    ax2 = axes[1]
    ax2_twin = ax2.twinx()
    x = np.arange(len(iterations))
    ax2.bar(x - 0.175, genome_counts, 0.35, color='#3498DB', alpha=0.8, label='Genomes')
    ax2_twin.bar(x + 0.175, mapped_snps, 0.35, color='#E74C3C', alpha=0.8, label='SNPs')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Genomes', color='#3498DB')
    ax2_twin.set_ylabel('SNPs', color='#E74C3C')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'R{i}' for i in iterations])
    ax2.set_title('B. Statistics', fontsize=11, fontweight='bold', pad=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'iteration_progress.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, 'iteration_progress.png'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def create_all_visualizations(iteration_results, output_dir, sample_name):
    """Create only iteration_progress figure."""
    try:
        create_iteration_progress_figure(iteration_results, output_dir, sample_name)
    except Exception as e:
        print(f"  [WARNING] create_all_visualizations ({sample_name}): {type(e).__name__}: {e}")


# =============================================================================
# Output Functions
# =============================================================================

def save_path_scores(path_results, work_dir, iteration=1):
    filepath = os.path.join(work_dir, 'path_scores.tsv')
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Path', 'Scores', 'NodeSNPCounts', 'JaccardClusterID', 'TotalSNP', 'Density', 'ZeroRatio', 'TotalLength', 'ImportanceRank', 'SnpRank', 'SnpProportion'])
        for p in sorted(path_results, key=lambda x: x['importance_rank']):
            writer.writerow([p['path'], p['scores_str'], p['snp_str'], p['jaccard_cluster'], p['total_snp'], f"{p['density']:.4f}", f"{p['zero_ratio']:.4f}", f"{p['total_length']:.4f}", p['importance_rank'], p['snp_rank'], f"{p['snp_proportion']:.4f}"])
    return filepath


def save_node_list(final_node, genomes, work_dir):
    filepath = os.path.join(work_dir, f"{final_node}.Node.list")
    with open(filepath, 'w') as f:
        for genome in genomes:
            f.write(f"{genome}\n")
    return filepath


def save_trace_result(results, output_file):
    """
    Save Trace.result.tsv with new column structure (Issue 13).

    Columns:
    Sample, Genomes, GenomesNumber, FinalResult, TraceLevel, ConfidenceScore,
    ConfidenceLevel, ContainmentRankScore, Reads, Coverage, Depth,
    FractionMappedReads, TotalSNP, TracingSnpUsed, SnpRank, GlobalPath,
    GlobalPathScores, GlobalNodeSNPCounts, LCA_Round, SubPath, SubPathScores,
    NodeSNPCounts, ZeroRatio, JaccardClusterID
    """
    if not results:
        return

    # Sort by SnpRank (primary), FinalResult='Yes' first (secondary)
    sorted_results = sorted(results, key=lambda x: (
        0 if x.get('FinalResult') == 'Yes' else 1,
        x.get('SnpRank', 999)
    ))

    fieldnames = [
        'Sample', 'Genomes', 'GenomesNumber', 'FinalResult', 'TraceLevel',
        'ConfidenceScore', 'ConfidenceLevel',
        'ContainmentRankScore', 'Reads', 'Coverage', 'Depth',
        'FractionMappedReads', 'TotalSNP', 'TracingSnpUsed',
        'SnpRank', 'GlobalPath', 'GlobalPathScores', 'GlobalNodeSNPCounts',
        'LCA_Round', 'SubPath', 'SubPathScores', 'NodeSNPCounts',
        'ZeroRatio', 'JaccardClusterID',
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in sorted_results:
            # 创建安全副本，不修改调用方的原始 dict
            safe_row = {fn: row.get(fn, 'N/A') for fn in fieldnames}
            writer.writerow(safe_row)



def save_iteration_stats(iteration_results, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Iteration', 'DatabaseGenomes', 'DatabaseNodes', 'DatabaseSNPs', 'MaxDepth',
                        'MappedSNPs', 'CandidateGenomes', 'FinalNode', 'TraceLevel', 'BestGlobalPath', 'ConfidenceScore', 'ElapsedTime'])
        for i, r in enumerate(iteration_results):
            writer.writerow([i + 1, r.get('db_stats', {}).get('total_genomes', 'N/A'),
                           r.get('db_stats', {}).get('total_nodes', 'N/A'), r.get('db_stats', {}).get('total_snps', 'N/A'),
                           r.get('db_stats', {}).get('max_depth', 'N/A'), r.get('mapped_snps', 'N/A'),
                           r.get('genome_count', 'N/A'), r.get('final_node', 'N/A'),
                           r.get('trace_level', 'N/A'), r.get('best_global_path', 'N/A'),
                           f"{r.get('confidence_score', 0):.4f}", r.get('elapsed_time', 'N/A')])


def generate_summary_report(iteration_results, output_file, sample_name, total_time, containment_info=None):
    if not iteration_results:
        return
    final_result = iteration_results[-1]
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("                    mTracePatho v1.0 Source Tracing Report\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Sample Name: {sample_name}\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Runtime: {total_time}\n\n")
        
        f.write("-" * 80 + "\nCONFIDENCE SCORE\n" + "-" * 80 + "\n")
        f.write(f"   Total Score: {final_result.get('confidence_score', 0):.4f}\n")
        f.write(f"   Level: {final_result.get('confidence_level', 'Unknown')}\n")
        components = final_result.get('confidence_components', {})
        f.write(f"   Components (accuracy-focused):\n")
        f.write(f"      - SNP Mapping Rate (30%): {components.get('snp_mapping_rate', 0):.4f}\n")
        f.write(f"        (proportion of sample SNPs mapped to database)\n")
        f.write(f"      - Path Dominance (35%): {components.get('path_dominance', 0):.4f}\n")
        f.write(f"        (best path advantage over alternatives)\n")
        f.write(f"      - SNP Consistency (20%): {components.get('snp_consistency', 0):.4f}\n")
        f.write(f"        (SNP distribution along path)\n")
        f.write(f"      - Path Uniqueness (15%): {components.get('path_uniqueness', 0):.4f}\n")
        f.write(f"        (fewer competing paths = less ambiguity)\n\n")
        f.write("   Threshold Recommendations:\n")
        f.write("      - High confidence: >= 0.8 (reliable result)\n")
        f.write("      - Medium confidence: 0.5-0.8 (review recommended)\n")
        f.write("      - Low confidence: < 0.5 (use with caution)\n\n")
        
        if containment_info:
            f.write("-" * 80 + "\nCONTAINMENT PRE-FILTERING\n" + "-" * 80 + "\n")
            f.write(f"   Method: {containment_info.get('method', 'N/A')}\n")
            f.write(f"   Selected: {containment_info.get('n_selected', 'N/A')}\n\n")
        
        f.write("-" * 80 + "\nINPUT SUMMARY\n" + "-" * 80 + "\n")
        f.write(f"   Total SNPs: {iteration_results[0]['total_sample_snps']}\n")
        f.write(f"   Mapped SNPs: {iteration_results[0]['mapped_snps']}\n")
        f.write(f"   Initial genomes: {iteration_results[0]['genome_count']}\n\n")
        
        f.write("-" * 80 + "\nITERATIVE TRACING\n" + "-" * 80 + "\n\n")
        for i, r in enumerate(iteration_results):
            f.write(f"   Iteration {i + 1}: {r['mapped_snps']} SNPs -> {r['genome_count']} genomes\n")
            f.write(f"      Level: {r.get('trace_level', 'N/A')}, Node: {r['final_node']}\n")
        
        f.write("\n" + "-" * 80 + "\nFINAL RESULT\n" + "-" * 80 + "\n")
        f.write(f"   Iterations: {len(iteration_results)}\n")
        f.write(f"   Trace Level: {final_result.get('trace_level', 'Unknown')}\n")
        f.write(f"   Final genomes: {final_result['genome_count']}\n")
        f.write(f"   GlobalPath: {final_result.get('best_global_path', 'N/A')}\n\n")
        
        final_genomes = final_result.get('final_genomes', [])
        lca_info = final_result.get('lca_info')
        
        if lca_info:
            f.write(f"   LCA Information:\n")
            f.write(f"      Reason: {lca_info.get('reason', 'N/A')}\n")
            if lca_info.get('tied_path_count'):
                f.write(f"      Tied paths: {lca_info['tied_path_count']}\n")
            f.write(f"      LCA Node: {lca_info.get('lca_node', 'N/A')}\n\n")
        
        if final_genomes:
            f.write(f"   Target Genome(s): ({len(final_genomes)} total)\n")
            for g in final_genomes[:10]:
                f.write(f"      - {g}\n")
            if len(final_genomes) > 10:
                f.write(f"      ... and {len(final_genomes) - 10} more\n")
        
        f.write("\n" + "=" * 80 + "\n")


# =============================================================================
# Main Pipeline
# =============================================================================

def run_tracing_iteration(sample_snps, db_dir, work_dir, jaccard_percent, original_genome_to_lineage, global_max_depth):
    snp_to_node, combined_snps = load_lca_output(db_dir)
    snp_lineage = load_snp_lineage(db_dir)
    node_dmp = load_node_dmp(db_dir)
    lineage_stat = load_lineage_stat(db_dir)
    branch_lengths = load_tree_branch_lengths(db_dir, node_dmp)
    db_stats = get_database_stats(db_dir)
    
    node_snp_data, mapped_snps_set, snp_detail = snp_to_node_mapping(sample_snps, snp_to_node, combined_snps)
    mapped_snps = len(mapped_snps_set)
    
    if snp_detail:
        save_snp_lca_detail(snp_detail, snp_lineage, lineage_stat, os.path.join(work_dir, 'snp_lca_detail.tsv'))
    
    if not node_snp_data:
        return None, True, None, None, lineage_stat, db_stats, mapped_snps, "", "", mapped_snps_set, [], 'Unknown', None
    
    scores = calculate_cumulative_scores(node_snp_data, snp_lineage)
    path_data = trace_paths(scores, node_dmp, node_snp_data)
    
    if not path_data:
        return None, True, None, None, lineage_stat, db_stats, mapped_snps, "", "", mapped_snps_set, [], 'Unknown', None
    
    path_results = cluster_and_rank_paths(path_data, node_dmp, branch_lengths, jaccard_percent)
    
    save_all_paths_detail(path_results, lineage_stat, original_genome_to_lineage, os.path.join(work_dir, 'all_paths_detail.tsv'))
    
    # Check for tied paths and handle LCA
    is_terminal, final_node, genomes, trace_level, lca_info = check_trace_result_with_ties(
        path_results, lineage_stat, original_genome_to_lineage
    )
    
    best_path = next((p for p in path_results if p['importance_rank'] == 1), None)
    best_path_snp_str = best_path['snp_str'] if best_path else ""
    best_path_scores_str = best_path['scores_str'] if best_path else ""
    
    return (path_results, is_terminal, final_node, genomes, lineage_stat, db_stats, 
            mapped_snps, best_path_snp_str, best_path_scores_str, mapped_snps_set, snp_detail,
            trace_level, lca_info)



def _zero_ratio_from_gp_scores(gp_scores_str: str) -> str:
    """Compute fraction of zero values in GlobalPathScores string."""
    try:
        if not gp_scores_str or gp_scores_str == 'N/A':
            return '0.0000'
        vals = [int(x.strip()) for x in gp_scores_str.split(',') if x.strip()]
        if not vals:
            return '0.0000'
        zero_count = sum(1 for v in vals if v == 0)
        return f"{zero_count / len(vals):.4f}"
    except Exception:
        return '0.0000'


def merge_genome_info_with_globalpath(path_results, lineage_stat, sample_name, 
                                       original_genome_to_lineage, node_iteration_map,
                                       is_final_round, final_trace_level, final_confidence_score, 
                                       final_confidence_level, total_sample_snps=0,
                                       containment_score_map=None,
                                       containment_rank_map=None,
                                       global_scores=None,
                                       global_node_snp_data=None):
    """
    Merge path results with GlobalPath information.

    Issue 1 fix: When tied paths exist and containment scores are also equal,
    count distinct genomes — if >1 → TraceLevel=Node (not Genome).

    Column structure (Issue 13):
    Sample, Genomes, GenomesNumber, FinalResult, TraceLevel, ConfidenceScore, ConfidenceLevel,
    ContainmentRankScore, Reads, Coverage, Depth, FractionMappedReads, TotalSNP, TracingSnpUsed,
    SnpRank, GlobalPath, GlobalPathScores, GlobalNodeSNPCounts, LCA_Round, SubPath, SubPathScores,
    NodeSNPCounts, ZeroRatio, JaccardClusterID
    """
    if containment_score_map is None:
        containment_score_map = {}
    if containment_rank_map is None:
        containment_rank_map = {}
    if global_scores is None:
        global_scores = {}
    if global_node_snp_data is None:
        global_node_snp_data = {}

    merged = []

    # ── Determine which genome should be marked FinalResult=Yes ──
    final_genome_id = None
    final_trace_level_override = None

    if is_final_round and path_results:
        best = next((p for p in path_results if p['importance_rank'] == 1), None)
        if best:
            best_scores = best['scores_str']
            tied_paths = [p for p in path_results if p['scores_str'] == best_scores]

            if len(tied_paths) <= 1:
                final_genome_id = None
            else:
                tied_genomes_scores = []
                for tp in tied_paths:
                    nodes = tp['path'].split('->')
                    ep = nodes[-1]
                    ep_info = lineage_stat.get(ep, {})
                    ep_genomes = ep_info.get('genomes', [])
                    if ep_genomes:
                        for g in ep_genomes:
                            c_score = containment_score_map.get(g, 0)
                            tied_genomes_scores.append((g, c_score, tp))

                if tied_genomes_scores:
                    tied_genomes_scores.sort(key=lambda x: x[1], reverse=True)
                    best_g, best_cs, _ = tied_genomes_scores[0]
                    second_cs = tied_genomes_scores[1][1] if len(tied_genomes_scores) > 1 else best_cs

                    final_genome_id = best_g
                    if best_cs > second_cs:
                        final_trace_level_override = 'Node'
                    else:
                        # Issue 1: check distinct genome count to decide Node vs Genome
                        all_tied_genome_names = list(set(g for g, _, _ in tied_genomes_scores))
                        final_trace_level_override = 'Node' if len(all_tied_genome_names) > 1 else 'Genome'

    for p in path_results:
        nodes = p['path'].split('->')
        endpoint_node = nodes[-1]
        stat_info = lineage_stat.get(endpoint_node, {})
        endpoint_genomes = stat_info.get('genomes', [])
        endpoint_genome = endpoint_genomes[0] if endpoint_genomes else None
        genome_count = len(endpoint_genomes)
        
        if endpoint_genome and endpoint_genome in original_genome_to_lineage:
            global_path, global_path_plain = get_genome_global_path(endpoint_genome, original_genome_to_lineage, node_iteration_map)
        else:
            global_path, global_path_plain = get_node_global_path(endpoint_node, original_genome_to_lineage, node_iteration_map)
        
        if is_final_round:
            if final_genome_id is not None:
                is_final = 'Yes' if (endpoint_genome == final_genome_id or
                                     final_genome_id in endpoint_genomes) else 'No'
            else:
                is_final = 'Yes' if p['importance_rank'] == 1 else 'No'
        else:
            is_final = 'No'
        
        path_trace_level = 'Genome' if genome_count == 1 else 'Node'
        if is_final == 'Yes':
            row_trace_level = final_trace_level_override if final_trace_level_override else final_trace_level
        else:
            row_trace_level = path_trace_level

        # ContainmentRankScore
        cont_parts = []
        for g in (endpoint_genomes or []):
            if g in containment_rank_map:
                rank, score = containment_rank_map[g]
                cont_parts.append(f"{rank}:{score:.6f}")
        cont_rank_score = ','.join(cont_parts) if cont_parts else 'N/A'

        # GlobalPathScores & GlobalNodeSNPCounts
        if global_path_plain and global_path_plain != 'N/A':
            try:
                gp_nodes  = [int(n) for n in global_path_plain.split('->') if n.strip().isdigit()]
                gp_scores = [global_scores.get(n, 0) for n in gp_nodes]
                gp_snps   = [global_node_snp_data.get(n, {}).get('count', 0) for n in gp_nodes]
                gp_scores_str = ','.join(str(s) for s in gp_scores)
                gp_snps_str   = ','.join(str(c) for c in gp_snps)
            except Exception:
                gp_scores_str, gp_snps_str = 'N/A', 'N/A'
        else:
            gp_scores_str, gp_snps_str = 'N/A', 'N/A'

        merged.append({
            'Sample': sample_name,
            'Genomes': ';'.join(endpoint_genomes) if endpoint_genomes else endpoint_node,
            'GenomesNumber': genome_count,
            'FinalResult': is_final,
            'TraceLevel': row_trace_level,
            'ConfidenceScore': f"{final_confidence_score:.4f}" if is_final == 'Yes' else 'N/A',
            'ConfidenceLevel': final_confidence_level if is_final == 'Yes' else 'N/A',
            'ContainmentRankScore': cont_rank_score,
            'Reads': 'N/A',
            'Coverage': 'N/A',
            'Depth': 'N/A',
            'FractionMappedReads': 'N/A',
            'TotalSNP': total_sample_snps,
            'TracingSnpUsed': 'N/A',
            'SnpRank': p['snp_rank'],
            'GlobalPath': global_path,
            'GlobalPathScores': gp_scores_str,
            'GlobalNodeSNPCounts': gp_snps_str,
            'LCA_Round': 0,  # set later
            'SubPath': p['path'],
            'SubPathScores': p['scores_str'],
            'NodeSNPCounts': p['snp_str'],
            'ZeroRatio': _zero_ratio_from_gp_scores(gp_scores_str),
            'JaccardClusterID': p['jaccard_cluster'],
        })
    return merged




def process_single_sample(sample_name, sample_snps, db_dir, output_dir, jaccard_percent, 
                          max_iterations, keep_temp, no_visualization, total_var, jaccard_nwk,
                          original_genome_to_lineage, global_max_depth, use_containment=False, 
                          min_containment=0.0, top_n_cutoff=99999, min_candidates=3, genome_snps_map=None,
                          phylo_gap_ratio=999999999999, phylo_share_thresh=0.7,
                          scoring_mode='containment_weighted', cont_alpha=0.5, cont_ref_size=None,
                          score_cutoff_ratio=0.5, min_snp=5,
                          min_unique_snp_n_cutoff=1):
    """
    Process a single sample through the full tracing pipeline.

    New in v1.0:
        min_snp: minimum non-root mapped SNPs required for a reliable trace.
                 If mapped SNPs < min_snp, TraceLevel is marked 'SNP<=N' instead.
    """
    sample_timer = Timer()
    sample_dir = os.path.join(output_dir, 'result', sample_name)
    figures_dir = os.path.join(sample_dir, 'figures')
    os.makedirs(sample_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    # BUG#11: In multiprocessing workers genome_snps_map is None (not pickled);
    # fall back to the process-local copy loaded by _worker_init().
    if genome_snps_map is None and _worker_genome_snps_map is not None:
        genome_snps_map = _worker_genome_snps_map

    # ── Precompute global (original db) SNP scores for GlobalPathScores column ──
    try:
        _snp_to_node_orig, _combined_orig = load_lca_output(db_dir)
        _snp_lineage_orig = load_snp_lineage(db_dir)
        _node_snp_data_orig, _, _ = snp_to_node_mapping(sample_snps, _snp_to_node_orig, _combined_orig)
        _global_scores = calculate_cumulative_scores(_node_snp_data_orig, _snp_lineage_orig)
    except Exception:
        _snp_to_node_orig = {}
        _global_scores = {}
        _node_snp_data_orig = {}
    # ── Precompute genome→terminal-node mapping for unique-SNP filter ────────
    try:
        _, _, _genome_to_node_orig, _ = load_snp_lineage_with_genome_map(db_dir)
    except Exception:
        _genome_to_node_orig = {}

    containment_info = None
    current_db_dir = db_dir
    single_genome_result = None
    selected_genomes = []
    containment_results = []
    containment_rank_map = {}   # genome → (rank, score) from containment_results



    if use_containment and genome_snps_map:
        # Step 1: compute containment for all genomes (keep full list for logging)
        containment_results_all = calculate_containment(
            sample_snps, genome_snps_map,
            scoring_mode=scoring_mode,
            cont_alpha=cont_alpha,
            cont_ref_size=cont_ref_size,
        )
        if scoring_mode != 'containment':
            ref_used = containment_results_all[0].get('ref_size', 'N/A') if containment_results_all else 'N/A'
            print(f"  [Scoring] {sample_name}: mode={scoring_mode}"
                  + (f", alpha={cont_alpha}, ref_size={ref_used}" if scoring_mode in ('containment_weighted', 'cosine_weighted', 'jaccard_weighted') else ""))
        containment_results = containment_results_all  # will be replaced if filter activates

        # Step 2: [NEW] phylogenetic pre-filter — remove distant-branch genomes
        #         before top_n_cutoff selection to prevent PTC score hijacking
        containment_results, phylo_info = phylo_prefilter_candidates(
            containment_results_all, original_genome_to_lineage,
            phylo_gap_ratio=phylo_gap_ratio,
            phylo_share_thresh=phylo_share_thresh
        )
        if phylo_info.get('activated'):
            print(f"  [PhyloFilter] {sample_name}: activated "
                  f"(gap_ratio={phylo_info['gap_ratio']}, "
                  f"removed={phylo_info['n_removed']}/{phylo_info['n_before']}, "
                  f"kept={phylo_info['n_after']})")
        else:
            print(f"  [PhyloFilter] {sample_name}: not activated "
                  f"({phylo_info.get('reason', '')})")

        # Step 3: top_n_cutoff selection on the (possibly filtered) list
        selected_genomes, containment_info = select_candidates_by_containment(
            containment_results, min_containment=min_containment, 
            top_n_cutoff=top_n_cutoff, min_candidates=min_candidates,
            score_cutoff_ratio=score_cutoff_ratio
        )
        # Attach phylo_info into containment_info for downstream reporting
        if containment_info is None:
            containment_info = {}
        containment_info['phylo_prefilter'] = phylo_info
        containment_info['scoring_mode']    = scoring_mode
        containment_info['cont_alpha']      = cont_alpha

        # Log hard filter result
        hfc = containment_info.get('hard_filtered_count', 0)
        if hfc > 0:
            print(f"  [HardFilter] {sample_name}: score_cutoff_ratio={score_cutoff_ratio}, "
                  f"removed={hfc}, kept={len(selected_genomes)}")
        elif score_cutoff_ratio > 0:
            print(f"  [HardFilter] {sample_name}: score_cutoff_ratio={score_cutoff_ratio}, "
                  f"no genomes removed (all above threshold)")

        # ── Unique-SNP presence filter (Step 4: after score_cutoff_ratio) ────
        if min_unique_snp_n_cutoff > 0 and _snp_to_node_orig and _genome_to_node_orig:
            selected_genomes, usf_info = filter_by_unique_snp_presence(
                selected_genomes,
                sample_snps,
                _snp_to_node_orig,
                _genome_to_node_orig,
                min_unique_snp_n_cutoff=min_unique_snp_n_cutoff,
            )
            n_removed = usf_info['n_before'] - usf_info['n_after']
            if n_removed > 0 and not usf_info.get('skipped'):
                removed_names = [g for g, _, _ in usf_info['removed_genomes']]
                print(f"  [UniqueSNPFilter] {sample_name}: removed {n_removed} genome(s) "
                      f"with 0 unique-node SNPs in sample: {removed_names}")
            elif usf_info.get('skipped'):
                pass  # warning already printed inside the function
            else:
                print(f"  [UniqueSNPFilter] {sample_name}: no genomes removed "
                      f"(all have ≥{min_unique_snp_n_cutoff} unique-node SNP(s) in sample)")
            containment_info['unique_snp_filter'] = usf_info

        # Build containment_rank_map for ContainmentRankScore column
        for i, r in enumerate(containment_results):
            _score_val = r.get('score', r.get('containment', 0))
            containment_rank_map[r['genome']] = (i + 1, round(_score_val, 6))

        # Save post-filter TSV (what enters LCA) and full filter decision log
        save_containment_results(containment_results, containment_info,
                                 os.path.join(sample_dir, 'containment_results.tsv'), sample_name)
        save_containment_filter_log(containment_results_all, phylo_info, containment_info,
                                    os.path.join(sample_dir, 'containment_results.filter.log'), sample_name)
        plot_containment_curve(containment_results, containment_info, os.path.join(figures_dir, 'containment_curve.png'), sample_name)
        
        # 特殊情况：只选中1个基因组，直接返回该基因组作为结果
        if len(selected_genomes) == 1:
            single_genome = selected_genomes[0]
            single_genome_result = {
                'genome': single_genome,
                'reason': 'single_containment_match',
                'containment': containment_results[0]['containment'] if containment_results else 1.0,
                'shared_snps': containment_results[0]['shared_snps'] if containment_results else len(sample_snps)
            }
            # 不构建子数据库，后面直接生成结果
        elif len(selected_genomes) >= 2 and len(selected_genomes) < len(genome_snps_map):
            containment_db_dir = os.path.join(sample_dir, 'containment_filtered_db')
            genome_list_file = os.path.join(sample_dir, 'containment_candidates.list')
            with open(genome_list_file, 'w') as f:
                for g in selected_genomes:
                    f.write(f"{g}\n")
            try:
                build_sub_database(genome_list_file, total_var, jaccard_nwk, containment_db_dir, sample_snps=sample_snps)
                current_db_dir = containment_db_dir
            except Exception as e:
                print(f"  [WARNING] {sample_name}: containment sub-database build failed "
                      f"({type(e).__name__}: {e}), falling back to original database")
                # current_db_dir remains db_dir (original)
    
    trace_results = []
    iteration_results = []
    iteration_history = []
    temp_dirs_to_clean = []
    
    # 处理单基因组Containment直接命中的情况
    if single_genome_result:
        single_genome = single_genome_result['genome']
        
        # 获取该基因组的GlobalPath
        if single_genome in original_genome_to_lineage:
            lineage = original_genome_to_lineage[single_genome]
            global_path_plain = "->".join(str(n) for n in lineage)
            global_path = global_path_plain + "[R0]"  # R0表示直接命中
            path_depth = len(lineage)
        else:
            global_path_plain = "N/A"
            global_path = "N/A"
            path_depth = 0
        
        # 计算confidence score - 单基因组直接命中
        mapped_snps = single_genome_result['shared_snps']
        snp_mapping_rate = mapped_snps / len(sample_snps) if sample_snps else 0

        # path_dominance: 基于 containment score 竞争差异（非硬编码 1.0）
        if containment_results and len(containment_results) >= 2:
            top1_sc = containment_results[0].get('score', containment_results[0]['containment'])
            top2_sc = containment_results[1].get('score', containment_results[1]['containment'])
            path_dominance = (top1_sc - top2_sc) / top1_sc if top1_sc > 0 else 0.0
        else:
            path_dominance = 1.0

        # path_uniqueness: 基于 cutoff 前的候选数
        n_above_cutoff = containment_info.get('n_before_hard_filter', 1) if containment_info else 1
        if n_above_cutoff <= 1:
            path_uniqueness = 1.0
        elif n_above_cutoff <= 5:
            path_uniqueness = 0.7
        elif n_above_cutoff <= 10:
            path_uniqueness = 0.5
        else:
            path_uniqueness = 0.3

        # snp_consistency: 单基因组无路径信息，保守估计
        snp_consistency = 0.5

        confidence_score = (0.30 * snp_mapping_rate +
                            0.35 * path_dominance +
                            0.20 * snp_consistency +
                            0.15 * path_uniqueness)
        confidence_level = get_confidence_level(confidence_score)
        
        confidence_components = {
            'snp_mapping_rate': snp_mapping_rate,
            'path_dominance': path_dominance,
            'snp_consistency': snp_consistency,
            'path_uniqueness': path_uniqueness
        }
        
        # 生成结果记录 (新列结构)
        # ContainmentRankScore for single genome
        _sg_cont_rs = 'N/A'
        if single_genome in containment_rank_map:
            _r, _s = containment_rank_map[single_genome]
            _sg_cont_rs = f"{_r}:{_s:.6f}"

        # GlobalPathScores for single genome
        if single_genome in original_genome_to_lineage:
            _sg_lineage = original_genome_to_lineage[single_genome]
            _sg_gp_plain = '->'.join(str(n) for n in _sg_lineage)
            _sg_gp_scores = ','.join(str(_global_scores.get(n, 0)) for n in _sg_lineage)
            _sg_gp_snpcnts = ','.join(str(_node_snp_data_orig.get(n, {}).get('count', 0)) for n in _sg_lineage)
        else:
            _sg_gp_plain = 'N/A'
            _sg_gp_scores = 'N/A'
            _sg_gp_snpcnts = 'N/A'

        _sg_trace_level = 'Genome'
        if min_snp > 0 and mapped_snps < min_snp:
            _sg_trace_level = f'SNP<={min_snp}'

        trace_results.append({
            'Sample': sample_name,
            'Genomes': single_genome,
            'GenomesNumber': 1,
            'FinalResult': 'Yes',
            'TraceLevel': _sg_trace_level,
            'ConfidenceScore': f"{confidence_score:.4f}",
            'ConfidenceLevel': confidence_level,
            'ContainmentRankScore': _sg_cont_rs,
            'Reads': 'N/A',
            'Coverage': 'N/A',
            'Depth': 'N/A',
            'FractionMappedReads': 'N/A',
            'TotalSNP': len(sample_snps),
            'TracingSnpUsed': 'N/A',
            'SnpRank': 1,
            'GlobalPath': global_path,
            'GlobalPathScores': _sg_gp_scores,
            'GlobalNodeSNPCounts': _sg_gp_snpcnts,
            'LCA_Round': 0,
            'SubPath': 'N/A',
            'SubPathScores': str(mapped_snps),
            'NodeSNPCounts': str(mapped_snps),
            'ZeroRatio': _zero_ratio_from_gp_scores(_sg_gp_scores),
            'JaccardClusterID': 0,
        })
        
        iteration_results.append({
            'level': 0, 'db_stats': {}, 'mapped_snps': mapped_snps,
            'genome_count': 1, 'final_node': single_genome,
            'best_global_path': global_path, 'best_global_path_plain': global_path_plain,
            'elapsed_time': '0.0s', 'best_path_snp_str': str(mapped_snps),
            'best_path_scores_str': str(mapped_snps), 'total_sample_snps': len(sample_snps),
            'final_genomes': [single_genome], 'trace_level': 'Genome',
            'lca_info': None, 'confidence_score': confidence_score,
            'confidence_level': confidence_level, 'confidence_components': confidence_components,
            'local_path_nodes': []
        })
        
        total_sample_time = sample_timer.format_time(sample_timer.elapsed())
        
        # 保存结果
        save_trace_result(trace_results, os.path.join(sample_dir, 'Trace.result.tsv'))
        save_iteration_stats(iteration_results, os.path.join(sample_dir, 'iteration_statistics.tsv'))
        generate_summary_report(iteration_results, os.path.join(sample_dir, 'tracing_report.txt'), 
                               sample_name, total_sample_time, containment_info)
        
        return {
            'sample_name': sample_name, 'iterations': 0,
            'final_genomes': 1, 'time': total_sample_time,
            'target': single_genome, 'trace_level': 'Genome',
            'confidence_score': confidence_score, 'confidence_level': confidence_level,
            'containment_candidates': 1
        }
    
    # 正常LCA溯源流程
    level = 0
    prev_genome_count = float('inf')
    current_sample_snps = sample_snps.copy()
    
    while level < max_iterations:
        iter_start_time = time.time()
        iter_dir = os.path.join(sample_dir, f'LCA_Iteration_{level + 1}')
        os.makedirs(iter_dir, exist_ok=True)
        
        result = run_tracing_iteration(current_sample_snps, current_db_dir, iter_dir, 
                                       jaccard_percent, original_genome_to_lineage, global_max_depth)
        
        if result[0] is None:
            path_results, is_terminal, final_node, genomes = None, True, None, None
            lineage_stat, db_stats, mapped_snps = result[4], result[5], result[6]
            best_path_snp_str, best_path_scores_str = result[7], result[8]
            mapped_snps_set = result[9] if len(result) > 9 else set()
            trace_level = result[11] if len(result) > 11 else 'Unknown'
            lca_info = result[12] if len(result) > 12 else None
        else:
            (path_results, is_terminal, final_node, genomes, lineage_stat, db_stats, 
             mapped_snps, best_path_snp_str, best_path_scores_str, mapped_snps_set, _,
             trace_level, lca_info) = result
        
        iter_elapsed = sample_timer.format_time(time.time() - iter_start_time)
        
        if path_results is None:
            if not iteration_results:
                # ── ContainmentFallback: all SNPs map to root node ──────────────
                # When SNPs are too conserved / all map to the tree root, LCA
                # returns no usable paths.  If containment pre-selection was
                # active, the selected_genomes are the best candidates we can
                # provide.  We do an LCA of those candidates in the *original*
                # full-database tree and report that node as the result.
                fallback_used = False
                if selected_genomes and use_containment:
                    lca_node, lca_lineage = find_common_ancestor_in_original_db(
                        selected_genomes, original_genome_to_lineage)
                    if lca_node is not None and lca_lineage:
                        fallback_genome_count = len(selected_genomes)
                        if fallback_genome_count == 1:
                            fb_trace_level = 'Genome'
                            fb_path = "->".join(str(n) for n in
                                                original_genome_to_lineage.get(selected_genomes[0],
                                                                               lca_lineage))
                        else:
                            fb_trace_level = 'Node'
                            fb_path = "->".join(str(n) for n in lca_lineage)

                        best_cont   = containment_results[0].get('score',
                                          containment_results[0].get('containment', 0)
                                      ) if containment_results else 0
                        shared_snps = (containment_results[0]['shared_snps']
                                       if containment_results else 0)

                        fb_snp_rate = shared_snps / len(sample_snps) if sample_snps else 0
                        fb_dom      = 1.0 / fallback_genome_count
                        fb_cons     = best_cont
                        fb_uniq     = 1.0 / max(1, fallback_genome_count)
                        conf_score  = (0.30 * fb_snp_rate + 0.35 * fb_dom
                                       + 0.20 * fb_cons + 0.15 * fb_uniq)
                        conf_level  = get_confidence_level(conf_score)

                        fb_global_path = fb_path + '[ContainmentFallback]'

                        iteration_results.append({
                            'level': level + 1, 'db_stats': db_stats,
                            'mapped_snps': shared_snps,
                            'genome_count': fallback_genome_count,
                            'final_node': lca_node,
                            'best_global_path': fb_global_path,
                            'best_global_path_plain': fb_path,
                            'elapsed_time': iter_elapsed,
                            'best_path_snp_str': str(shared_snps),
                            'best_path_scores_str': str(shared_snps),
                            'total_sample_snps': len(sample_snps),
                            'final_genomes': selected_genomes,
                            'trace_level': fb_trace_level, 'lca_info': None,
                            'confidence_score': conf_score,
                            'confidence_level': conf_level,
                            'confidence_components': {
                                'snp_mapping_rate': fb_snp_rate,
                                'path_dominance':   fb_dom,
                                'snp_consistency':  fb_cons,
                                'path_uniqueness':  fb_uniq
                            }
                        })
                        # ContainmentRankScore for fallback
                        _fb_cont_parts = []
                        for _fg in selected_genomes:
                            if _fg in containment_rank_map:
                                _fr, _fs = containment_rank_map[_fg]
                                _fb_cont_parts.append(f"{_fr}:{_fs:.6f}")
                        _fb_cont_rs = ','.join(_fb_cont_parts) if _fb_cont_parts else 'N/A'

                        # GlobalPathScores for fallback path
                        _fb_gp_scores = ','.join(
                            str(_global_scores.get(n, 0))
                            for n in [int(x) for x in fb_path.split('->') if x.strip().isdigit()]
                        ) if fb_path != 'N/A' else 'N/A'
                        _fb_gp_snpcnts = ','.join(
                            str(_node_snp_data_orig.get(n, {}).get('count', 0))
                            for n in [int(x) for x in fb_path.split('->') if x.strip().isdigit()]
                        ) if fb_path != 'N/A' else 'N/A'

                        _fb_trace_level = fb_trace_level
                        if min_snp > 0 and shared_snps < min_snp:
                            _fb_trace_level = f'SNP<={min_snp}'

                        trace_results.append({
                            'Sample': sample_name,
                            'Genomes': ';'.join(selected_genomes),
                            'GenomesNumber': fallback_genome_count,
                            'FinalResult': 'Yes',
                            'TraceLevel': _fb_trace_level,
                            'ConfidenceScore': f"{conf_score:.4f}",
                            'ConfidenceLevel': conf_level,
                            'ContainmentRankScore': _fb_cont_rs,
                            'Reads': 'N/A',
                            'Coverage': 'N/A',
                            'Depth': 'N/A',
                            'FractionMappedReads': 'N/A',
                            'TotalSNP': len(sample_snps),
                            'TracingSnpUsed': 'N/A',
                            'SnpRank': 1,
                            'GlobalPath': fb_global_path,
                            'GlobalPathScores': _fb_gp_scores,
                            'GlobalNodeSNPCounts': _fb_gp_snpcnts,
                            'LCA_Round': level + 1,
                            'SubPath': 'N/A',
                            'SubPathScores': str(shared_snps),
                            'NodeSNPCounts': str(shared_snps),
                            'ZeroRatio': _zero_ratio_from_gp_scores(_fb_gp_scores),
                            'JaccardClusterID': 0,
                        })
                        fallback_used = True

                if not fallback_used:
                    iteration_results.append({
                        'level': level + 1, 'db_stats': db_stats,
                        'mapped_snps': mapped_snps,
                        'genome_count': 0, 'final_node': 'N/A',
                        'best_global_path': 'N/A', 'best_global_path_plain': 'N/A',
                        'elapsed_time': iter_elapsed, 'best_path_snp_str': '',
                        'best_path_scores_str': '',
                        'total_sample_snps': len(sample_snps), 'final_genomes': [],
                        'trace_level': 'Unknown', 'lca_info': None,
                        'confidence_score': 0, 'confidence_level': 'Low',
                        'confidence_components': {}
                    })
            break
        
        genome_count = len(genomes) if genomes else 0
        iteration_history.append({'iteration': level + 1, 'genomes': genomes if genomes else [], 'final_node': final_node})
        
        best_path = next((p for p in path_results if p['importance_rank'] == 1), None)
        local_path_nodes = [int(n) for n in best_path['path'].split('->')] if best_path else []
        
        # Get GlobalPath based on trace level
        if trace_level == 'Genome' and genomes:
            best_endpoint_genome = genomes[0]
            best_global_path, best_global_path_plain = get_genome_global_path(best_endpoint_genome, original_genome_to_lineage)
        elif lca_info and lca_info.get('lca_node'):
            best_global_path, best_global_path_plain = get_node_global_path(lca_info['lca_node'], original_genome_to_lineage)
        else:
            best_global_path, best_global_path_plain = "N/A", "N/A"
        
        path_depth = len(best_global_path_plain.split('->')) if best_global_path_plain != 'N/A' else 0
        
        # Calculate confidence score
        confidence_score, confidence_components = calculate_confidence_score(
            mapped_snps, len(sample_snps), best_path, path_results,
            genome_count, global_max_depth, path_depth, trace_level
        )
        confidence_level = get_confidence_level(confidence_score)
        
        iteration_results.append({
            'level': level + 1, 'db_stats': db_stats, 'mapped_snps': mapped_snps, 'genome_count': genome_count,
            'final_node': final_node, 'elapsed_time': iter_elapsed, 'best_path_snp_str': best_path_snp_str,
            'best_path_scores_str': best_path_scores_str, 'total_sample_snps': len(sample_snps),
            'final_genomes': genomes, 'path_results': path_results, 'lineage_stat': lineage_stat,
            'local_path_nodes': local_path_nodes,
            'best_global_path': best_global_path, 'best_global_path_plain': best_global_path_plain,
            'trace_level': trace_level, 'lca_info': lca_info,
            'confidence_score': confidence_score, 'confidence_level': confidence_level,
            'confidence_components': confidence_components
        })
        
        save_path_scores(path_results, iter_dir, level + 1)
        
        if is_terminal or genome_count >= prev_genome_count:
            break
        
        prev_genome_count = genome_count
        node_list_file = save_node_list(final_node, genomes, iter_dir)
        sub_db_dir = os.path.join(iter_dir, f'{final_node}_Node_DB')
        temp_dirs_to_clean.append(sub_db_dir)
        
        try:
            build_sub_database(node_list_file, total_var, jaccard_nwk, sub_db_dir, sample_snps=current_sample_snps)
        except Exception as e:
            print(f"  [WARNING] {sample_name}: iteration {level+1} sub-database build "
                  f"failed ({type(e).__name__}: {e}), stopping iteration")
            break
        
        current_db_dir = sub_db_dir
        if mapped_snps_set:
            # BUG#5: strip combo-SNP strings (contain '+') to prevent O(N²) pollution next round
            current_sample_snps = {s for s in mapped_snps_set if '+' not in s}
        level += 1
    
    total_sample_time = sample_timer.format_time(sample_timer.elapsed())
    
    node_iteration_map = build_iteration_markers(iteration_history, original_genome_to_lineage)
    
    # Get final result info
    final_result = iteration_results[-1] if iteration_results else {}
    final_trace_level = final_result.get('trace_level', 'Unknown')
    final_confidence_score = final_result.get('confidence_score', 0)
    final_confidence_level = final_result.get('confidence_level', 'Low')
    
    n_iters = len(iteration_results)
    # Build containment score map for tie-breaking in final result assignment
    _cont_score_map = {}
    if containment_results:
        for cr in containment_results:
            _cont_score_map[cr['genome']] = cr.get('score', cr.get('containment', 0))

    for i, r in enumerate(iteration_results):
        if r.get('path_results'):
            merged = merge_genome_info_with_globalpath(
                r['path_results'], r['lineage_stat'], sample_name,
                original_genome_to_lineage, node_iteration_map, i == n_iters - 1,
                final_trace_level, final_confidence_score, final_confidence_level,
                total_sample_snps=len(sample_snps),
                containment_score_map=_cont_score_map,
                containment_rank_map=containment_rank_map,
                global_scores=_global_scores,
                global_node_snp_data=_node_snp_data_orig,
            )
            for row in merged:
                row['LCA_Round'] = i + 1
                trace_results.append(row)

    # Issue 12: apply min_snp flag to FinalResult='Yes' rows
    if min_snp > 0:
        # Use mapped_snps from the first iteration as proxy for non-root SNP count
        _first_mapped = iteration_results[0].get('mapped_snps', 0) if iteration_results else 0
        if _first_mapped < min_snp:
            for row in trace_results:
                if row.get('FinalResult') == 'Yes':
                    _orig_level = row.get('TraceLevel', '')
                    if _orig_level not in ('Error', 'Unknown') and not _orig_level.startswith('SNP<='):
                        row['TraceLevel'] = f'SNP<={min_snp}'
    
    if iteration_results:
        save_trace_result(trace_results, os.path.join(sample_dir, 'Trace.result.tsv'))
        save_iteration_stats(iteration_results, os.path.join(sample_dir, 'iteration_statistics.tsv'))
        generate_summary_report(iteration_results, os.path.join(sample_dir, 'tracing_report.txt'), sample_name, total_sample_time, containment_info)
        if not no_visualization:
            create_all_visualizations(iteration_results, figures_dir, sample_name)
    
    if not keep_temp:
        for temp_dir in temp_dirs_to_clean:
            try:
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
            except Exception as e:
                print(f"  [WARNING] Failed to remove temp dir {temp_dir}: {e}")
    
    return {
        'sample_name': sample_name, 'iterations': len(iteration_results),
        'final_genomes': final_result.get('genome_count', 0), 'time': total_sample_time,
        'target': final_result.get('final_genomes', ['N/A'])[0] if final_result.get('final_genomes') else 'N/A',
        'trace_level': final_trace_level,
        'confidence_score': final_confidence_score,
        'confidence_level': final_confidence_level,
        'containment_candidates': containment_info.get('n_selected', 'N/A') if containment_info else 'N/A'
    }


def process_sample_wrapper(args):
    try:
        return process_single_sample(*args)
    except Exception as e:
        return {'sample_name': args[0], 'iterations': 0, 'final_genomes': 0, 'time': 'Error', 
                'target': f'Error: {str(e)}', 'trace_level': 'Error', 'confidence_score': 0, 
                'confidence_level': 'Low', 'containment_candidates': 'N/A'}


def main():
    parser = argparse.ArgumentParser(description='mTracePatho v1.0 - Pathogen source tracing via LCA algorithm')
    parser.add_argument('--sample_snp', required=True, help='Input TSV: col1=sample_name, col2=comma-separated SNPs')
    parser.add_argument('--db_dir', required=True, help='Database directory')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--jaccard_percent', type=float, default=0.8, help='Jaccard threshold (default: 0.8)')
    parser.add_argument('--max_iterations', type=int, default=10, help='Max LCA iterations (default: 10)')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Parallel workers (default: 1)')
    parser.add_argument('--keep_temp', action='store_true', help='Keep intermediate files')
    parser.add_argument('--no_visualization', action='store_true', help='Skip visualization output')
    parser.add_argument('--use_containment', action='store_true', default=True, help='Enable containment pre-selection (default: on)')
    parser.add_argument('--min_containment', type=float, default=0.0, help='Minimum containment threshold (default: 0.0)')
    parser.add_argument('--top_n_cutoff', type=int, default=99999, help='Top N distinct containment score levels to select (default: 99999 = all)')
    parser.add_argument('--min_candidates', type=int, default=3, help='Minimum candidate genomes to guarantee (default: 3)')
    parser.add_argument('--score_cutoff_ratio', type=float, default=0.5,
                        help='Hard cutoff: score >= top1 × ratio; 1 survivor → direct result (default: 0.5)')
    parser.add_argument('--phylo_gap_ratio', type=float, default=999999999999,
                        help='PhyloFilter activation: top1/level2 ratio threshold (default: disabled)')
    parser.add_argument('--phylo_share_thresh', type=float, default=0.7,
                        help='PhyloFilter: min shared lineage fraction to keep candidate (default: 0.7)')
    parser.add_argument('--scoring_mode', type=str, default='containment_weighted',
                        choices=['containment', 'containment_weighted', 'jaccard', 'cosine',
                                 'cosine_weighted', 'jaccard_weighted'],
                        help='Containment scoring mode (default: containment_weighted)')
    parser.add_argument('--cont_alpha', type=float, default=0.5,
                        help='Weighted scoring exponent: score = metric × (ref/size)^alpha (default: 0.5)')
    parser.add_argument('--cont_ref_size', type=float, default=None,
                        help='Reference genome size for weighted scoring (default: median GenomeSNPs)')
    parser.add_argument('--min_unique_snp_n_cutoff', type=int, default=1,
                        help='After score_cutoff_ratio filter, remove any genome whose '
                             'terminal-node unique SNPs have fewer than N hits in the sample '
                             '(default: 1). Set 0 to disable.')
    parser.add_argument('--min_snp', type=int, default=5,
                        help='Minimum non-root mapped SNPs for reliable tracing; '
                             'if below threshold TraceLevel is marked SNP<=N (default: 5)')
    
    args = parser.parse_args()
    
    timer = Timer()
    result_dir = os.path.join(args.output_dir, 'result')
    os.makedirs(result_dir, exist_ok=True)
    
    # BUG#6: use context manager so file handle is closed even on exception
    with TimingLogger(os.path.join(result_dir, 'timing_log.txt')) as timing_logger:

        print("\n" + "=" * 60)
        print("     mTracePatho v1.0 - Pathogen Source Tracing")
        print("=" * 60)
        print("")

        if not os.path.exists(args.sample_snp):
            sys.exit(f"Error: Input file not found: {args.sample_snp}")

        required_files = ['snp_lingage.tsv', 'node.tsv', 'LCA.output.tsv', 'lingage.stat.tsv',
                          'numbered_tree.tree', 'total.var.tsv', 'Jaccard.fastme.input.nwk']
        for f in required_files:
            if not os.path.exists(os.path.join(args.db_dir, f)):
                sys.exit(f"Error: Missing database file: {f}")

        timing_logger.log("Loading database...")
        _, original_genome_to_lineage, _, global_max_depth = load_snp_lineage_with_genome_map(args.db_dir)
        timing_logger.log(f"Loaded {len(original_genome_to_lineage)} genome lineages, max depth: {global_max_depth}")

        genome_snps_map = None
        if args.use_containment:
            print("Loading SNP matrix for containment...")
            timing_logger.log("Loading SNP matrix for containment...")
            genome_snps_map, _, _ = load_snp_matrix_for_containment(args.db_dir)
            timing_logger.log(f"Loaded SNP data for {len(genome_snps_map)} genomes")
            print(f"  Loaded SNP data for {len(genome_snps_map)} genomes")

        samples = load_samples_from_tsv(args.sample_snp)
        if not samples:
            sys.exit("Error: No valid samples found")

        print(f"Processing {len(samples)} sample(s) with {args.threads} worker(s)...")
        timing_logger.log(f"Processing {len(samples)} samples")

        total_var = os.path.join(args.db_dir, 'total.var.tsv')
        jaccard_nwk = os.path.join(args.db_dir, 'Jaccard.fastme.input.nwk')

        batch_results = []

        if args.threads == 1:
            for sample_name, sample_snps in samples:
                result = process_single_sample(
                    sample_name, sample_snps, args.db_dir, args.output_dir,
                    args.jaccard_percent, args.max_iterations, args.keep_temp,
                    args.no_visualization, total_var, jaccard_nwk, original_genome_to_lineage,
                    global_max_depth,
                    args.use_containment, args.min_containment, args.top_n_cutoff,
                    args.min_candidates,
                    genome_snps_map, args.phylo_gap_ratio, args.phylo_share_thresh,
                    args.scoring_mode, args.cont_alpha, args.cont_ref_size,
                    args.score_cutoff_ratio, args.min_snp,
                    args.min_unique_snp_n_cutoff,        # BUG#1 fix
                )
                batch_results.append(result)
                conf_str = f"{result['confidence_score']:.2f}({result['confidence_level']})"
                timing_logger.log_sample(result['sample_name'], result['iterations'], result['final_genomes'], result['target'], result['time'], conf_str)
                print(f"  {result['sample_name']}: {result['trace_level']} -> {result['final_genomes']}G, Conf: {conf_str}")
        else:
            # BUG#11: pass None for genome_snps_map; workers load it via _worker_init()
            sample_args = [
                (sample_name, sample_snps, args.db_dir, args.output_dir,
                 args.jaccard_percent, args.max_iterations, args.keep_temp,
                 args.no_visualization, total_var, jaccard_nwk, original_genome_to_lineage,
                 global_max_depth,
                 args.use_containment, args.min_containment, args.top_n_cutoff,
                 args.min_candidates,
                 None,                                   # genome_snps_map loaded by _worker_init
                 args.phylo_gap_ratio, args.phylo_share_thresh,
                 args.scoring_mode, args.cont_alpha, args.cont_ref_size,
                 args.score_cutoff_ratio, args.min_snp,
                 args.min_unique_snp_n_cutoff)            # BUG#1 fix
                for sample_name, sample_snps in samples
            ]
            _init_args = (args.db_dir,) if args.use_containment else None
            with ProcessPoolExecutor(
                    max_workers=args.threads,
                    initializer=_worker_init if args.use_containment else None,
                    initargs=_init_args) as executor:
                futures = {executor.submit(process_sample_wrapper, arg): arg[0] for arg in sample_args}
                for future in as_completed(futures):
                    result = future.result()
                    batch_results.append(result)
                    conf_str = f"{result['confidence_score']:.2f}({result['confidence_level']})"
                    timing_logger.log_sample(result['sample_name'], result['iterations'], result['final_genomes'], result['target'], result['time'], conf_str)
                    print(f"  {result['sample_name']}: {result['trace_level']} -> {result['final_genomes']}G, Conf: {conf_str}")

        sample_order = {name: i for i, (name, _) in enumerate(samples)}
        batch_results.sort(key=lambda x: sample_order.get(x['sample_name'], len(samples)))

        with open(os.path.join(result_dir, 'batch_summary.tsv'), 'w', newline='') as f:
            fieldnames = ['sample_name', 'iterations', 'final_genomes', 'target',
                          'trace_level', 'confidence_score', 'confidence_level', 'time']
            if args.use_containment:
                fieldnames.append('containment_candidates')
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            for r in batch_results:
                row = {k: r.get(k, 'N/A') for k in fieldnames}
                writer.writerow(row)

        total_time = timer.format_time(timer.elapsed())
        timing_logger.log(f"\nTotal runtime: {total_time}")
        # TimingLogger.close() called automatically by __exit__ (BUG#6)

        high_conf = sum(1 for r in batch_results if r.get('confidence_level') == 'High')
        med_conf = sum(1 for r in batch_results if r.get('confidence_level') == 'Medium')
        low_conf = sum(1 for r in batch_results if r.get('confidence_level') == 'Low')
        genome_level = sum(1 for r in batch_results if r.get('trace_level') == 'Genome')
        node_level = sum(1 for r in batch_results if r.get('trace_level') == 'Node')

        print("\n" + "=" * 60)
        print(f"COMPLETE - {len(samples)} samples in {total_time}")
        print(f"\nConfidence: High={high_conf} Medium={med_conf} Low={low_conf}")
        print(f"TraceLevel: Genome={genome_level} Node={node_level}")
        print(f"\nResults: {result_dir}")
        print("=" * 60 + "\n")


if __name__ == '__main__':
    main()
