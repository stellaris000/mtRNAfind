#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
annotation_free_adaptive_depth.py

Independent annotation-free mtRNA endpoint-pair discovery display pipeline with
optional high-confidence display filters and depth-adaptive parameter presets.

mtfind version: 0.311

Runtime input:
- raw SSC endpoint file only.

This version is intentionally annotation-free:
- It does not read an annotation GTF.
- It does not read any old result folder.
- It regenerates raw endpoint-pairs, sparse 2D connected components, post-cluster
  superclusters, high-confidence display GTF, audit tables, and 2D plots from SSC.
- Optional annotation/reference files are used only for plot overlays or
  post-hoc benchmark, never for candidate generation.

Biological boundary:
- GTF output is a coordinate-candidate display track, not validated transcript annotation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle

from gtf_coordinate_benchmark import add_benchmark_cli, read_coordinate_file, run_gtf_coordinate_benchmark


MTFIND_VERSION = "0.311"
NOISE_BOUNDARY_SEMANTICS = (
    "raw_support < min_signal_support is noise; "
    "raw_support >= min_signal_support is signal"
)


@dataclass
class Params:
    mtfind_version: str = MTFIND_VERSION
    sample_id: str = "BC04"
    chrom_name: str = "MT"
    mt_length: int = 16569
    density_min_raw_pair_support: int = 3
    density_seed_min_raw_pair_support: int = 5
    density_backfill_max_raw_pair_support: int = 4
    density_connect_start_bp: int = 5
    density_connect_end_bp: int = 5
    density_min_cluster_support: int = 20
    display_min_cluster_support: int = 200
    display_min_top_raw_fraction: float = 0.0
    display_min_representative_raw_support: int = 0
    display_max_start_span_bp: int = 0
    display_max_end_span_bp: int = 0
    compact_delta_start_bp: int = 50
    compact_delta_end_bp: int = 50
    shared_endpoint_delta_bp: int = 20
    shadow_other_endpoint_delta_bp: int = 300
    containment_min: float = 0.80
    plot_min_raw_pair_support: int = 5
    plot_max_raw_points: int = 40000
    plot_annotation_reference: str = ""
    plot_annotation_reference_format: str = "auto"
    run_mode: str = "full"
    adaptive_depth_mode: str = "fixed"
    adaptive_depth_total_support: int = 0
    adaptive_depth_profile: str = "NA"
    adaptive_raw_sqrt_divisor: float = 10.0
    adaptive_cluster_fraction: float = 0.015
    adaptive_display_fraction: float = 0.035
    adaptive_representative_fraction: float = 0.006
    adaptive_low_cluster_fraction: float = 0.0
    adaptive_low_display_fraction: float = 0.0
    adaptive_low_representative_fraction: float = 0.010
    adaptive_low_total_support_max: int = 10000
    adaptive_medium_total_support_max: int = 100000
    # ── noise-robust depth estimation ──
    noise_detection_method: str = "auto"           # "kneedle" | "mad" | "quantile" | "auto"
    min_signal_support: int = 0                    # runtime: minimum support classified as signal
    noise_floor: int = 0                           # runtime compatibility alias for min_signal_support
    noise_boundary_semantics: str = NOISE_BOUNDARY_SEMANTICS
    noise_knee_curvature: float = 0.0              # runtime: knee significance (kneedle mode)
    noise_mad_k: float = 3.0                       # k multiplier for MAD (default 3 ≈ 3σ)
    noise_min_knee_curvature: float = 0.03         # min curvature to accept knee; below → fallback to MAD
    signal_total_support: int = 0                  # runtime: signal-only total support
    noise_total_support: int = 0                   # runtime: noise total support
    noise_fraction: float = 0.0                    # runtime: noise / total support
    signal_pair_count: int = 0                     # runtime: signal pair count
    noise_pair_count: int = 0                      # runtime: noise pair count
    density_effective_seed_min_raw_pair_support: int = 0
    density_effective_backfill_max_raw_pair_support: int = 0
    density_seed_pair_count: int = 0
    density_backfill_candidate_pair_count: int = 0
    density_backfill_assigned_pair_count: int = 0
    density_backfill_unassigned_pair_count: int = 0
    density_backfill_ambiguous_pair_count: int = 0
    density_backfill_assigned_support: int = 0
    # ── multi-factor display scoring (replaces binary absolute threshold) ──
    display_score_mode: str = "scored"             # "scored" | "binary" (legacy)
    display_score_threshold: float = 0.40          # minimum composite score to include in GTF
    display_score_weight_support: float = 0.35     # weight of absolute support sufficiency
    display_score_weight_purity: float = 0.15      # weight of signal purity (top_raw_fraction)
    display_score_weight_peak: float = 0.30        # weight of peak signal strength (vs noise floor)
    display_score_weight_rank: float = 0.20        # weight of relative rank (vs median cluster)
    display_score_relative_factor: float = 0.60    # effective_threshold = max(base, median_score × this)


class DSU:
    def __init__(self, n:int):
        self.parent=list(range(n))
        self.rank=[0]*n
    def find(self,x:int)->int:
        while self.parent[x]!=x:
            self.parent[x]=self.parent[self.parent[x]]
            x=self.parent[x]
        return x
    def union(self,a:int,b:int)->None:
        ra,rb=self.find(a),self.find(b)
        if ra==rb: return
        if self.rank[ra]<self.rank[rb]: self.parent[ra]=rb
        elif self.rank[ra]>self.rank[rb]: self.parent[rb]=ra
        else:
            self.parent[rb]=ra
            self.rank[ra]+=1


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def ensure_dirs(root: Path) -> Dict[str, Path]:
    dirs={
        'root': root,
        'inputs': root/'inputs',
        'code': root/'code',
        'tables': root/'results'/'tables',
        'gtf': root/'results'/'gtf',
        'plots': root/'results'/'plots',
        'logs': root/'results'/'logs',
    }
    for d in dirs.values(): d.mkdir(parents=True, exist_ok=True)
    return dirs


def normalize_chrom_name(chrom: str) -> str:
    return {'chrM': 'MT', 'M': 'MT', 'chrMT': 'MT'}.get(str(chrom), str(chrom))


def load_raw_endpoint_pairs(ssc_path: Path, p: Params) -> pd.DataFrame:
    names=['read_id','chr','strand','start','end','unused6','unused7','unused8']
    df=pd.read_csv(ssc_path, sep='	', header=None, names=names, usecols=list(range(8)), dtype=str)
    df['chr_norm']=df['chr'].map(normalize_chrom_name)
    target=normalize_chrom_name(p.chrom_name)
    if target.lower() == 'auto':
        vc=df.loc[df['strand'].isin(['+','-']), 'chr_norm'].value_counts()
        if vc.empty:
            raise ValueError('No + or - strand rows found in SSC input.')
        target=str(vc.index[0])
    p.chrom_name=target
    df=df[(df['chr_norm']==target) & (df['strand'].isin(['+','-']))].copy()
    df['start']=pd.to_numeric(df['start'], errors='coerce')
    df['end']=pd.to_numeric(df['end'], errors='coerce')
    df=df.dropna(subset=['start','end'])
    df['start']=df['start'].astype(int)
    df['end']=df['end'].astype(int)
    df=df[(df['start']>=1) & (df['end']>=df['start']) & (df['end']<=p.mt_length)].copy()
    df['chr']=target
    pairs=df.groupby(['strand','start','end'], as_index=False).agg(raw_support=('read_id','size'))
    pairs=pairs.sort_values(['strand','start','end']).reset_index(drop=True)
    pairs['raw_pair_id']=[f'RP_{i:08d}' for i in range(1,len(pairs)+1)]
    pairs['length']=pairs['end']-pairs['start']+1
    return pairs


def clamp_int(value: float, low: int, high: int) -> int:
    return max(low, min(high, int(round(value))))


def quantile_or_zero(values: pd.Series, q: float) -> int:
    if values is None or len(values) == 0:
        return 0
    return int(round(float(values.quantile(q))))


def detect_noise_floor(raw_pairs: pd.DataFrame, p: Params) -> dict:
    """Detect the minimum signal support using a three-tier adaptive strategy.

    Tier 1 ─ Kneedle: knee detection on the cumulative support curve.
             The returned boundary is the minimum signal support: pairs below
             it are noise and pairs at or above it are signal.
             This is purely data-driven ─ no fixed thresholds.
    Tier 2 ─ MAD: Median Absolute Deviation on log10(support).
             Robust to extreme values; naturally adapts to depth.
             Fallback when the knee is too shallow (unimodal distribution).
    Tier 3 ─ Quantile: Q80 of support as last resort.
             Only used when MAD ≈ 0 (degenerate distribution).

    Returns a dict with keys:
        min_signal_support, method, knee_curvature (if kneedle),
        mad_log_median, mad_log_mad (if MAD), all_profile_data
    """
    if len(raw_pairs) == 0:
        return {
            "min_signal_support": 0, "noise_floor": 0, "method": "empty",
            "knee_curvature": 0.0, "mad_log_median": 0.0, "mad_log_mad": 0.0,
            "all_profile_data": {"noise_boundary_semantics": NOISE_BOUNDARY_SEMANTICS}
        }

    support = raw_pairs['raw_support'].values.astype(np.int64)
    total_all = int(support.sum())
    n_all = len(support)

    # Always compute quantiles on the full distribution for audit
    q50_full = int(round(float(np.percentile(support, 50))))
    q75_full = int(round(float(np.percentile(support, 75))))
    q90_full = int(round(float(np.percentile(support, 90))))
    q95_full = int(round(float(np.percentile(support, 95))))
    max_support = int(support.max())

    method_requested = str(p.noise_detection_method or "auto").lower()
    min_signal_support: int = 0
    method_used: str = ""
    knee_curvature: float = 0.0
    mad_log_median: float = 0.0
    mad_log_mad: float = 0.0

    # ── Internal: Kneedle knee detection ──
    def _kneedle_detect(support_arr: np.ndarray, total: int) -> tuple:
        """Find the knee point on the cumulative support curve.

        Builds normalised (x, y) where
          x[i] = i / n          (cumulative fraction of pairs)
          y[i] = cumsum(sorted_support)[i] / total  (cumulative fraction of support)

        The knee is the point farthest *below* the diagonal y = x
        (i.e. where cumulative contribution is most delayed).
        Returns (min_signal_support, curvature) where curvature = max(|y - x|).
        """
        n = len(support_arr)
        if n < 5:
            return 0, 0.0
        sorted_s = np.sort(support_arr)
        cumsum = np.cumsum(sorted_s).astype(np.float64)
        x = np.arange(1, n + 1, dtype=np.float64) / n
        y = cumsum / max(total, 1)
        diff = y - x                         # negative in noise zone, catches up in signal zone
        knee_idx = int(np.argmin(diff))      # point farthest below diagonal
        curvature = float(np.abs(diff[knee_idx]))
        # The knee support is the minimum signal support.  Values strictly
        # below it are noise; the boundary value itself remains signal.
        min_signal = max(int(sorted_s[knee_idx]), 1)
        return min_signal, curvature

    # ── Internal: MAD-based threshold ──
    def _mad_detect(support_arr: np.ndarray, k: float) -> tuple:
        """Robust threshold via Median Absolute Deviation on log10(support).

        signal_threshold = 10 ** (median_log + k * mad_log)

        Returns (min_signal_support, median_log, mad_log).
        """
        log_s = np.log10(support_arr.astype(np.float64) + 1.0)
        median_log = float(np.median(log_s))
        mad_log = float(np.median(np.abs(log_s - median_log)))
        if mad_log < 1e-9:
            return 0, median_log, mad_log
        threshold_log = median_log + k * mad_log
        min_signal = max(int(round(10.0 ** threshold_log)), 1)
        return min_signal, median_log, mad_log

    # ── Tier selection ──
    if method_requested in ("auto", "kneedle"):
        nf, curv = _kneedle_detect(support, total_all)
        knee_curvature = curv
        if method_requested == "kneedle" or curv >= p.noise_min_knee_curvature:
            min_signal_support = nf
            method_used = "kneedle"
            # Also compute MAD for audit
            _, mad_log_median, mad_log_mad = _mad_detect(support, p.noise_mad_k)

    if method_used == "" and method_requested in ("auto", "mad"):
        nf, mad_log_median, mad_log_mad = _mad_detect(support, p.noise_mad_k)
        if nf > 0:
            min_signal_support = nf
            method_used = "mad"
        elif method_requested == "mad":
            min_signal_support = 1
            method_used = "mad_zero"

    if method_used == "":
        # Tier 3: quantile fallback
        min_signal_support = max(int(round(float(np.percentile(support, 80)))), 1)
        method_used = "quantile_fallback"

    # ── Split signal / noise ──
    signal_mask = support >= min_signal_support
    signal_total = int(support[signal_mask].sum())
    noise_total = total_all - signal_total
    noise_fraction = noise_total / max(total_all, 1)
    signal_n = int(signal_mask.sum())
    noise_n = n_all - signal_n

    # ── Save to Params ──
    p.min_signal_support = min_signal_support
    # Retain the historical field name for downstream compatibility.  Its
    # value is explicitly a minimum signal support, not the maximum noise.
    p.noise_floor = min_signal_support
    p.noise_boundary_semantics = NOISE_BOUNDARY_SEMANTICS
    p.noise_detection_method = method_used
    p.noise_knee_curvature = knee_curvature
    p.signal_total_support = signal_total
    p.noise_total_support = noise_total
    p.noise_fraction = noise_fraction
    p.signal_pair_count = signal_n
    p.noise_pair_count = noise_n

    profile = {
        "raw_pair_count": n_all,
        "total_support": total_all,
        "signal_total_support": signal_total,
        "noise_total_support": noise_total,
        "noise_fraction": round(noise_fraction, 6),
        "signal_pair_count": signal_n,
        "noise_pair_count": noise_n,
        "min_signal_support": min_signal_support,
        "noise_floor": min_signal_support,
        "noise_floor_alias_for": "min_signal_support",
        "noise_boundary_semantics": NOISE_BOUNDARY_SEMANTICS,
        "noise_detection_method": method_used,
        "noise_knee_curvature": round(knee_curvature, 6),
        "noise_min_knee_curvature": p.noise_min_knee_curvature,
        "noise_mad_k": p.noise_mad_k,
        "noise_mad_log_median": round(mad_log_median, 6),
        "noise_mad_log_mad": round(mad_log_mad, 6),
        "raw_support_q50": q50_full,
        "raw_support_q75": q75_full,
        "raw_support_q90": q90_full,
        "raw_support_q95": q95_full,
        "raw_support_max": max_support,
    }
    return {
        "min_signal_support": min_signal_support,
        "noise_floor": min_signal_support,
        "method": method_used,
        "knee_curvature": knee_curvature,
        "mad_log_median": mad_log_median, "mad_log_mad": mad_log_mad,
        "all_profile_data": profile
    }


def apply_adaptive_depth_params(raw_pairs: pd.DataFrame, p: Params) -> None:
    """Mutate Params after raw-pair loading using noise-robust depth estimation.

    Noise detection runs first to separate signal from background.
    All adaptive thresholds are then derived from signal_total_support only,
    avoiding inflation by random fragmentation / PCR duplicate noise.
    """
    mode = str(p.adaptive_depth_mode or "fixed").lower()

    # ── Noise detection (always runs; even in "fixed" mode for audit) ──
    noise_result = detect_noise_floor(raw_pairs, p)
    p.adaptive_depth_profile = json.dumps(noise_result["all_profile_data"], sort_keys=True)

    # Use signal-only total for all threshold derivations
    total = p.signal_total_support if p.signal_total_support > 0 else (
        int(raw_pairs['raw_support'].sum()) if len(raw_pairs) else 0
    )
    p.adaptive_depth_total_support = total
    # The adaptive floor remains separate from the 0.311 seed cutoff.  The
    # latter is enforced in connected_components_2d after this calibration.
    min_signal_support = p.min_signal_support

    if mode == "fixed":
        return
    if mode == "auto":
        if total < p.adaptive_low_total_support_max:
            mode = "low"
        elif total < p.adaptive_medium_total_support_max:
            mode = "medium"
        else:
            mode = "high"
    if mode not in {"low", "medium", "high"}:
        raise ValueError(f"Unsupported adaptive depth mode: {p.adaptive_depth_mode}")

    # ── Noise-filtered parameter calibration ──
    # With noise removed, signal points are cleaner → tighten spatial radii.
    # Fractions slightly increased to compensate for signal-only total.
    # Low-mode display clamp floor lowered (50→20) so weak-but-real clusters survive.
    if mode == "low":
        cluster_fraction = p.adaptive_low_cluster_fraction if p.adaptive_low_cluster_fraction > 0 else p.adaptive_cluster_fraction
        display_fraction = p.adaptive_low_display_fraction if p.adaptive_low_display_fraction > 0 else p.adaptive_display_fraction
        rep_fraction = p.adaptive_low_representative_fraction if p.adaptive_low_representative_fraction > 0 else p.adaptive_representative_fraction
        p.density_min_raw_pair_support = clamp_int(math.sqrt(max(total, 1)) / p.adaptive_raw_sqrt_divisor, max(3, min_signal_support // 2), 25)
        p.density_min_cluster_support = clamp_int(total * cluster_fraction, 20, 200)
        p.display_min_cluster_support = clamp_int(total * display_fraction, 20, 1000)
        p.display_min_representative_raw_support = clamp_int(total * rep_fraction, 10, 1000)
        p.display_min_top_raw_fraction = max(float(p.display_min_top_raw_fraction), 0.10)
        p.density_connect_start_bp = max(p.density_connect_start_bp, 25)
        p.density_connect_end_bp = max(p.density_connect_end_bp, 25)
        p.compact_delta_start_bp = max(p.compact_delta_start_bp, 120)
        p.compact_delta_end_bp = max(p.compact_delta_end_bp, 120)
        p.shared_endpoint_delta_bp = max(p.shared_endpoint_delta_bp, 50)
        p.shadow_other_endpoint_delta_bp = min(max(p.shadow_other_endpoint_delta_bp, 400), 1000)
        p.containment_min = min(float(p.containment_min), 0.75)
        p.plot_min_raw_pair_support = min(p.plot_min_raw_pair_support, 1)
    elif mode == "medium":
        p.density_min_raw_pair_support = clamp_int(math.sqrt(max(total, 1)) / p.adaptive_raw_sqrt_divisor, max(3, min_signal_support // 2), 25)
        p.density_min_cluster_support = clamp_int(total * p.adaptive_cluster_fraction, 20, 200)
        p.display_min_cluster_support = clamp_int(total * p.adaptive_display_fraction, 30, 1000)
        p.display_min_representative_raw_support = clamp_int(total * p.adaptive_representative_fraction, 10, 1000)
        p.display_min_top_raw_fraction = max(float(p.display_min_top_raw_fraction), 0.10)
        p.density_connect_start_bp = max(p.density_connect_start_bp, 20)
        p.density_connect_end_bp = max(p.density_connect_end_bp, 20)
        p.compact_delta_start_bp = max(p.compact_delta_start_bp, 100)
        p.compact_delta_end_bp = max(p.compact_delta_end_bp, 100)
        p.shared_endpoint_delta_bp = max(p.shared_endpoint_delta_bp, 40)
        p.shadow_other_endpoint_delta_bp = min(max(p.shadow_other_endpoint_delta_bp, 400), 1000)
        p.containment_min = min(float(p.containment_min), 0.78)
        p.plot_min_raw_pair_support = min(p.plot_min_raw_pair_support, 3)
    else:
        p.density_min_raw_pair_support = clamp_int(math.sqrt(max(total, 1)) / p.adaptive_raw_sqrt_divisor, max(3, min_signal_support // 2), 25)
        p.density_min_cluster_support = clamp_int(total * p.adaptive_cluster_fraction, 20, 200)
        p.display_min_cluster_support = clamp_int(total * p.adaptive_display_fraction, 50, 1000)
        p.display_min_representative_raw_support = clamp_int(total * p.adaptive_representative_fraction, 10, 1000)
        p.display_min_top_raw_fraction = max(float(p.display_min_top_raw_fraction), 0.10)
        p.density_connect_start_bp = max(p.density_connect_start_bp, 15)
        p.density_connect_end_bp = max(p.density_connect_end_bp, 15)
        p.compact_delta_start_bp = max(p.compact_delta_start_bp, 80)
        p.compact_delta_end_bp = max(p.compact_delta_end_bp, 80)
        p.shared_endpoint_delta_bp = max(p.shared_endpoint_delta_bp, 30)
        p.shadow_other_endpoint_delta_bp = min(max(p.shadow_other_endpoint_delta_bp, 400), 1000)
        p.containment_min = min(float(p.containment_min), 0.82)
    p.adaptive_depth_mode = mode


def classify_shape(start_span:int, end_span:int) -> str:
    if start_span <=20 and end_span <=20: return 'compact_endpoint_pair_cluster'
    if start_span >60 and end_span <=20: return 'same_end_like_region'
    if end_span >60 and start_span <=20: return 'same_start_like_region'
    return 'complex_region'


def connected_components_2d(raw_pairs: pd.DataFrame, p: Params) -> pd.DataFrame:
    """Cluster high-support seeds, then backfill nearby low-support signal pairs.

    Only seed pairs establish 2D connectivity.  A 1–4-support pair can add
    support to exactly one adjacent seed component, but never joins components
    or propagates membership through other low-support pairs.
    """
    p.density_effective_seed_min_raw_pair_support = 0
    p.density_effective_backfill_max_raw_pair_support = 0
    p.density_seed_pair_count = 0
    p.density_backfill_candidate_pair_count = 0
    p.density_backfill_assigned_pair_count = 0
    p.density_backfill_unassigned_pair_count = 0
    p.density_backfill_ambiguous_pair_count = 0
    p.density_backfill_assigned_support = 0
    if raw_pairs is None or len(raw_pairs) == 0:
        return pd.DataFrame()

    min_signal_support = max(int(p.min_signal_support or p.noise_floor or 1), 1)
    effective_seed_min = max(
        int(p.density_min_raw_pair_support),
        int(p.density_seed_min_raw_pair_support),
        min_signal_support,
    )
    backfill_max = min(max(int(p.density_backfill_max_raw_pair_support), 0), effective_seed_min - 1)
    p.density_effective_seed_min_raw_pair_support = effective_seed_min
    p.density_effective_backfill_max_raw_pair_support = backfill_max

    signal_pairs = raw_pairs[raw_pairs['raw_support'] >= min_signal_support].copy()
    seed_pairs = signal_pairs[signal_pairs['raw_support'] >= effective_seed_min].copy()
    if backfill_max >= min_signal_support:
        backfill_pairs = signal_pairs[
            (signal_pairs['raw_support'] < effective_seed_min)
            & (signal_pairs['raw_support'] <= backfill_max)
        ].copy()
    else:
        backfill_pairs = signal_pairs.iloc[0:0].copy()

    p.density_seed_pair_count = int(len(seed_pairs))
    p.density_backfill_candidate_pair_count = int(len(backfill_pairs))
    if len(seed_pairs) == 0:
        p.density_backfill_unassigned_pair_count = int(len(backfill_pairs))
        return pd.DataFrame()

    rows = []
    cluster_no = 1
    for strand, seed_sub0 in seed_pairs.groupby('strand', sort=True):
        seed_sub = seed_sub0.sort_values(['start', 'end', 'raw_support'], ascending=[True, True, False]).reset_index(drop=True)
        n = len(seed_sub)
        dsu = DSU(n)
        grid = defaultdict(list)
        starts = seed_sub['start'].to_numpy(dtype=int)
        ends = seed_sub['end'].to_numpy(dtype=int)
        bx_size = p.density_connect_start_bp + 1
        by_size = p.density_connect_end_bp + 1
        for i, (start, end) in enumerate(zip(starts, ends)):
            bx = start // bx_size
            by = end // by_size
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for j in grid.get((bx + dx, by + dy), []):
                        if (
                            abs(start - starts[j]) <= p.density_connect_start_bp
                            and abs(end - ends[j]) <= p.density_connect_end_bp
                        ):
                            dsu.union(i, j)
            grid[(bx, by)].append(i)

        roots = defaultdict(list)
        for i in range(n):
            roots[dsu.find(i)].append(i)
        root_items = list(roots.items())
        root_to_component_index = {root: index for index, (root, _) in enumerate(root_items)}
        assigned_backfill = [[] for _ in root_items]

        strand_backfill = backfill_pairs[backfill_pairs['strand'].eq(strand)].sort_values(
            ['start', 'end', 'raw_support']
        )
        for low_pair in strand_backfill.to_dict('records'):
            low_start = int(low_pair['start'])
            low_end = int(low_pair['end'])
            bx = low_start // bx_size
            by = low_end // by_size
            adjacent_roots = set()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for j in grid.get((bx + dx, by + dy), []):
                        if (
                            abs(low_start - starts[j]) <= p.density_connect_start_bp
                            and abs(low_end - ends[j]) <= p.density_connect_end_bp
                        ):
                            adjacent_roots.add(dsu.find(j))
            if len(adjacent_roots) == 1:
                component_index = root_to_component_index[next(iter(adjacent_roots))]
                assigned_backfill[component_index].append(low_pair)
                p.density_backfill_assigned_pair_count += 1
                p.density_backfill_assigned_support += int(low_pair['raw_support'])
            elif len(adjacent_roots) > 1:
                p.density_backfill_ambiguous_pair_count += 1

        for component_index, (_, members) in enumerate(root_items):
            seed_comp = seed_sub.iloc[members].copy()
            backfill_records = assigned_backfill[component_index]
            if backfill_records:
                backfill_comp = pd.DataFrame(backfill_records)
                comp = pd.concat([seed_comp, backfill_comp], ignore_index=True)
            else:
                backfill_comp = seed_comp.iloc[0:0].copy()
                comp = seed_comp.copy()

            seed_support = int(seed_comp['raw_support'].sum())
            backfill_support = int(backfill_comp['raw_support'].sum())
            cluster_support = seed_support + backfill_support
            if cluster_support < p.density_min_cluster_support:
                continue
            rep = seed_comp.sort_values(
                ['raw_support', 'length', 'start', 'end'], ascending=[False, False, True, True]
            ).iloc[0]
            start_min = int(comp['start'].min())
            start_max = int(comp['start'].max())
            end_min = int(comp['end'].min())
            end_max = int(comp['end'].max())
            rep_start = int(rep['start'])
            rep_end = int(rep['end'])
            backfill_rule = 'unique_seed_component_within_connect_radius'
            backfill_support_rule = (
                f'backfill_signal_support_{min_signal_support}_to_{backfill_max}'
                if backfill_max >= min_signal_support
                else 'backfill_none_below_minimum_signal'
            )
            rows.append({
                'mtfind_version': MTFIND_VERSION,
                'cluster_id': f'DF{cluster_no:06d}',
                'cluster_origin': 'density_cluster',
                'candidate_source': 'annotation_free_independent_density_cluster',
                'strand': strand,
                'representative_method': 'highest_support_seed_endpoint_pair',
                'segmentation_method': 'seeded_2d_components_with_low_support_backfill_no_annotation',
                'density_component_id': cluster_no,
                'coordinate_basis': 'raw_representative_from_seed_component',
                'filter_reason': (
                    f'cluster_support_ge_{p.density_min_cluster_support};'
                    f'seed_pair_support_ge_{effective_seed_min};'
                    f'{backfill_support_rule};'
                    f'backfill_{backfill_rule}'
                ),
                'representative_raw_pair_id': rep['raw_pair_id'],
                'representative_start': rep_start,
                'representative_end': rep_end,
                'representative_length': rep_end - rep_start + 1,
                'representative_raw_support': int(rep['raw_support']),
                'seed_support': seed_support,
                'backfill_support': backfill_support,
                'cluster_support': cluster_support,
                'top_raw_fraction': round(float(rep['raw_support']) / float(cluster_support), 6) if cluster_support else 0,
                'seed_member_raw_pair_count': int(len(seed_comp)),
                'backfill_member_raw_pair_count': int(len(backfill_comp)),
                'member_raw_pair_count': int(len(comp)),
                'minimum_signal_support': min_signal_support,
                'effective_seed_min_raw_pair_support': effective_seed_min,
                'backfill_max_raw_pair_support': backfill_max,
                'backfill_assignment_rule': backfill_rule,
                'member_raw_pair_ids': ';'.join(comp['raw_pair_id'].astype(str).tolist()),
                'start_min': start_min,
                'start_max': start_max,
                'end_min': end_min,
                'end_max': end_max,
                'start_span': start_max - start_min,
                'end_span': end_max - end_min,
                'start_shift_min_from_rep': start_min - rep_start,
                'start_shift_max_from_rep': start_max - rep_start,
                'end_shift_min_from_rep': end_min - rep_end,
                'end_shift_max_from_rep': end_max - rep_end,
                'shape_class': classify_shape(start_max - start_min, end_max - end_min),
                'exact_raw_pair_exists': True,
                'candidate_id': f'DF2D_{cluster_no:06d}',
                'biological_status': 'coordinate_candidate_not_validated',
            })
            cluster_no += 1

    p.density_backfill_unassigned_pair_count = (
        p.density_backfill_candidate_pair_count - p.density_backfill_assigned_pair_count
    )
    return pd.DataFrame(rows)


def interval_overlap(a_start:int,a_end:int,b_start:int,b_end:int)->int:
    return max(0,min(a_end,b_end)-max(a_start,b_start)+1)


def merge_clusters(df: pd.DataFrame, p: Params, prefix: str='AFSC') -> Tuple[pd.DataFrame,pd.DataFrame]:
    if df is None or len(df)==0: return pd.DataFrame(), pd.DataFrame()
    d=df.copy().reset_index(drop=True)
    for c in [
        'representative_start', 'representative_end', 'representative_length',
        'representative_raw_support', 'seed_support', 'backfill_support',
        'cluster_support', 'seed_member_raw_pair_count',
        'backfill_member_raw_pair_count', 'member_raw_pair_count', 'start_min',
        'start_max', 'end_min', 'end_max', 'start_span', 'end_span',
    ]:
        if c in d.columns: d[c]=pd.to_numeric(d[c], errors='coerce')
    merge_edges=[]; group_rows=[]
    for strand in sorted(d['strand'].dropna().unique()):
        idx=np.asarray(d.index[d['strand'].eq(strand)].tolist(), dtype=int)
        n=len(idx)
        if n==0: continue
        dsu=DSU(n)
        rep_s=d.loc[idx,'representative_start'].to_numpy(dtype=int)
        rep_e=d.loc[idx,'representative_end'].to_numpy(dtype=int)
        rep_len=np.maximum(1,rep_e-rep_s+1)
        cl_ids=d.loc[idx,'cluster_id'].astype(str).to_numpy()
        cand_ids=d.loc[idx,'candidate_id'].astype(str).to_numpy()
        # ── spatial grid for O(n·k) instead of O(n²) all-pairs ──
        # Use the max of all merge-distance parameters as grid cell size.
        # Any pair within merge range will be in same or adjacent cells.
        max_delta = max(p.compact_delta_start_bp, p.compact_delta_end_bp,
                        p.shared_endpoint_delta_bp, p.shadow_other_endpoint_delta_bp)
        bin_size = max_delta + 1
        # Small-n fallback: O(n²) is faster when n < 20 (avoid grid overhead)
        if n < 20:
            for ai in range(n):
                for bi in range(ai+1,n):
                    ds=abs(rep_s[ai]-rep_s[bi]); de=abs(rep_e[ai]-rep_e[bi])
                    ov=max(0,min(rep_e[ai],rep_e[bi])-max(rep_s[ai],rep_s[bi])+1)
                    cont=ov/max(1,min(rep_len[ai],rep_len[bi]))
                    ok=False; reason='no_merge'
                    if ds <= p.compact_delta_start_bp and de <= p.compact_delta_end_bp:
                        ok=True; reason='compact_representative_close'
                    elif ds <= p.shared_endpoint_delta_bp and de <= p.shadow_other_endpoint_delta_bp and cont >= p.containment_min:
                        ok=True; reason='same_start_shadow_or_series_merge'
                    elif de <= p.shared_endpoint_delta_bp and ds <= p.shadow_other_endpoint_delta_bp and cont >= p.containment_min:
                        ok=True; reason='same_end_shadow_or_series_merge'
                    if ok:
                        dsu.union(ai,bi)
                        merge_edges.append({'strand':strand,'cluster_id_a':cl_ids[ai],'candidate_id_a':cand_ids[ai],'cluster_id_b':cl_ids[bi],'candidate_id_b':cand_ids[bi],'merge_reason':reason,'delta_start':int(ds),'delta_end':int(de),'containment':round(float(cont),6)})
        else:
            grid=defaultdict(list)
            for i in range(n):
                bx = int(rep_s[i]) // bin_size
                by = int(rep_e[i]) // bin_size
                grid[(bx, by)].append(i)
            for ai in range(n):
                bx = int(rep_s[ai]) // bin_size
                by = int(rep_e[ai]) // bin_size
                seen = set()
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for bi in grid.get((bx + dx, by + dy), []):
                            if bi <= ai or bi in seen:
                                continue
                            seen.add(bi)
                            ds = abs(rep_s[ai] - rep_s[bi]); de = abs(rep_e[ai] - rep_e[bi])
                            ov = max(0, min(rep_e[ai], rep_e[bi]) - max(rep_s[ai], rep_s[bi]) + 1)
                            cont = ov / max(1, min(rep_len[ai], rep_len[bi]))
                            ok = False; reason = 'no_merge'
                            if ds <= p.compact_delta_start_bp and de <= p.compact_delta_end_bp:
                                ok = True; reason = 'compact_representative_close'
                            elif ds <= p.shared_endpoint_delta_bp and de <= p.shadow_other_endpoint_delta_bp and cont >= p.containment_min:
                                ok = True; reason = 'same_start_shadow_or_series_merge'
                            elif de <= p.shared_endpoint_delta_bp and ds <= p.shadow_other_endpoint_delta_bp and cont >= p.containment_min:
                                ok = True; reason = 'same_end_shadow_or_series_merge'
                            if ok:
                                dsu.union(ai, bi)
                                merge_edges.append({'strand': strand, 'cluster_id_a': cl_ids[ai], 'candidate_id_a': cand_ids[ai], 'cluster_id_b': cl_ids[bi], 'candidate_id_b': cand_ids[bi], 'merge_reason': reason, 'delta_start': int(ds), 'delta_end': int(de), 'containment': round(float(cont), 6)})
        roots=defaultdict(list)
        for local,global_i in enumerate(idx): roots[dsu.find(local)].append(int(global_i))
        group_rows.extend(roots.values())
    rows=[]
    for n,members in enumerate(group_rows, start=1):
        sub=d.loc[members].copy()
        rep=sub.sort_values(['representative_raw_support','cluster_support','representative_length','representative_start','representative_end'], ascending=[False,False,False,True,True]).iloc[0]
        start_min=int(sub['start_min'].min()); start_max=int(sub['start_max'].max())
        end_min=int(sub['end_min'].min()); end_max=int(sub['end_max'].max())
        rep_start=int(rep['representative_start']); rep_end=int(rep['representative_end'])
        support=int(sub['cluster_support'].sum()); rep_support=int(rep['representative_raw_support'])
        seed_support = int(sub['seed_support'].sum()) if 'seed_support' in sub.columns else support
        backfill_support = int(sub['backfill_support'].sum()) if 'backfill_support' in sub.columns else 0
        seed_member_count = int(sub['seed_member_raw_pair_count'].sum()) if 'seed_member_raw_pair_count' in sub.columns else int(sub['member_raw_pair_count'].sum())
        backfill_member_count = int(sub['backfill_member_raw_pair_count'].sum()) if 'backfill_member_raw_pair_count' in sub.columns else 0
        rows.append({
            'mtfind_version':MTFIND_VERSION,
            'supercluster_id':f'{prefix}{n:06d}',
            'cluster_origin':'annotation_free_supercluster',
            'candidate_source':'annotation_free_independent_main',
            'strand':rep['strand'],
            'representative_raw_pair_id':rep.get('representative_raw_pair_id','NA'),
            'representative_start':rep_start,
            'representative_end':rep_end,
            'representative_length':rep_end-rep_start+1,
            'representative_raw_support':rep_support,
            'seed_support':seed_support,
            'backfill_support':backfill_support,
            'supercluster_support':support,
            'cluster_support':support,
            'top_raw_fraction':round(rep_support/support,6) if support else 0,
            'member_primary_cluster_count':int(len(sub)),
            'seed_member_raw_pair_count':seed_member_count,
            'backfill_member_raw_pair_count':backfill_member_count,
            'member_raw_pair_count':int(sub['member_raw_pair_count'].sum()),
            'minimum_signal_support':int(p.min_signal_support),
            'effective_seed_min_raw_pair_support':int(p.density_effective_seed_min_raw_pair_support),
            'backfill_max_raw_pair_support':int(p.density_effective_backfill_max_raw_pair_support),
            'backfill_assignment_rule':'unique_seed_component_within_connect_radius',
            'member_cluster_ids':';'.join(sub['cluster_id'].astype(str).tolist()),
            'member_candidate_ids':';'.join(sub['candidate_id'].astype(str).tolist()),
            'start_min':start_min,
            'start_max':start_max,
            'end_min':end_min,
            'end_max':end_max,
            'start_span':start_max-start_min,
            'end_span':end_max-end_min,
            'start_shift_min_from_rep':start_min-rep_start,
            'start_shift_max_from_rep':start_max-rep_start,
            'end_shift_min_from_rep':end_min-rep_end,
            'end_shift_max_from_rep':end_max-rep_end,
            'shape_class':'supercluster_from_primary_density_clusters',
            'coordinate_basis':'raw_representative_from_supercluster',
            'representative_method':'highest_support_raw_endpoint_pair_across_member_clusters',
            'segmentation_method':'annotation_free_post_cluster_merge',
            'biological_status':'coordinate_candidate_not_validated',
            'annotation_used':'false',
        })
    sc=pd.DataFrame(rows)
    if len(sc):
        sc['candidate_id']=[f'AF_CAND_{i:06d}' for i in range(1,len(sc)+1)]
        # Build boolean conditions once (vectorized)
        if p.display_score_mode == "scored":
            # ── Multi-factor scoring with relative thresholds ──
            support_arr = sc['cluster_support'].astype(float)
            rep_arr = sc['representative_raw_support'].astype(float)
            top_frac_arr = sc['top_raw_fraction'].astype(float)

            median_support = float(support_arr.median()) if len(support_arr) > 0 else 1.0
            abs_threshold = max(float(p.display_min_cluster_support), 1.0)
            min_signal_support = max(p.min_signal_support or p.noise_floor, 1)

            # ── Combined confidence: two independent indicators of whether
            #     the noise/signal separation is meaningful on this chromosome.
            #     ① Knee curvature: high → clear bimodal → reliable minimum-signal boundary
            #     ② Signal fraction: low → most pairs are noise → good separation
            #        (If >50% of pairs pass the boundary → separation is weak.)
            raw_pair_count = max(p.signal_pair_count + p.noise_pair_count, 1)
            signal_fraction = p.signal_pair_count / raw_pair_count
            curvature_conf = min(1.0, max(0.20, p.noise_knee_curvature / 0.30))
            separation_conf = min(1.0, 0.50 / max(signal_fraction, 0.01))
            combined_confidence = curvature_conf * separation_conf

            # Factor 1: Support sufficiency (relative to adaptive absolute threshold)
            f_support = (support_arr / abs_threshold).clip(0.0, 1.0)

            # Factor 2: Signal purity (rep dominates cluster = clean peak)
            f_purity = (top_frac_arr / 0.30).clip(0.0, 1.0)

            # Factor 3: Peak signal strength
            f_peak = (rep_arr / (min_signal_support * 3.0)).clip(0.0, 1.0)

            # Factor 4: Relative rank among sibling clusters on this chromosome
            f_rank = (support_arr / max(median_support * 0.30, 1.0)).clip(0.0, 1.0)

            # Raw composite score
            score_raw = (
                f_support * p.display_score_weight_support +
                f_purity  * p.display_score_weight_purity +
                f_peak    * p.display_score_weight_peak +
                f_rank    * p.display_score_weight_rank
            )
            # Global confidence multiplier: weak noise separation → entire score discounted.
            # This ensures single-candidate chromosomes (f_purity=f_rank=1.0 always)
            # don't pass on structural advantage alone.
            score = score_raw * combined_confidence

            # ── Relative threshold ──
            median_score = float(score.median()) if len(score) > 0 else 0.0
            effective_threshold = max(
                p.display_score_threshold,
                median_score * p.display_score_relative_factor
            )
            include = score >= effective_threshold

            sc['display_score'] = score.round(4)

            # Diagnostic reason strings
            reasons = []
            for i in range(len(sc)):
                tag = (
                    f"sc_{score.iloc[i]:.3f}"
                    f"=raw{score_raw.iloc[i]:.3f}×cf{combined_confidence:.2f}"
                    f"_s{f_support.iloc[i]:.2f}"
                    f"_p{f_purity.iloc[i]:.2f}"
                    f"_k{f_peak.iloc[i]:.2f}"
                    f"_r{f_rank.iloc[i]:.2f}"
                )
                if not include.iloc[i]:
                    fail_parts = []
                    if f_support.iloc[i] < 0.5:
                        fail_parts.append(f"sup{support_arr.iloc[i]:.0f}vs{abs_threshold:.0f}")
                    if f_purity.iloc[i] < 0.5:
                        fail_parts.append(f"pur{top_frac_arr.iloc[i]:.3f}")
                    if f_peak.iloc[i] < 0.5:
                        fail_parts.append(f"pk{rep_arr.iloc[i]:.0f}ms{min_signal_support}")
                    if f_rank.iloc[i] < 0.5:
                        fail_parts.append(f"rk{support_arr.iloc[i]:.0f}md{median_support:.0f}")
                    if fail_parts:
                        tag += "_" + "_".join(fail_parts)
                    tag += f"_below_{effective_threshold:.3f}"
                reasons.append(tag)
            sc['display_filter_reason'] = reasons
            sc['include_high_confidence_gtf'] = include
        else:
            # ── Legacy binary mode (backward compatible) ──
            cond_cluster = sc['cluster_support'].astype(int) >= p.display_min_cluster_support
            cond_top = pd.Series(True, index=sc.index)
            if p.display_min_top_raw_fraction > 0:
                cond_top = sc['top_raw_fraction'].astype(float) >= p.display_min_top_raw_fraction
            cond_rep = pd.Series(True, index=sc.index)
            if p.display_min_representative_raw_support > 0:
                cond_rep = sc['representative_raw_support'].astype(int) >= p.display_min_representative_raw_support
            cond_start_span = pd.Series(True, index=sc.index)
            if p.display_max_start_span_bp > 0:
                cond_start_span = sc['start_span'].astype(int) <= p.display_max_start_span_bp
            cond_end_span = pd.Series(True, index=sc.index)
            if p.display_max_end_span_bp > 0:
                cond_end_span = sc['end_span'].astype(int) <= p.display_max_end_span_bp
            include = cond_cluster & cond_top & cond_rep & cond_start_span & cond_end_span
            reasons = []
            n_sc = len(sc)
            for i in range(n_sc):
                fail = []
                if not cond_cluster.iloc[i]:
                    fail.append(f"cluster_support_lt_{p.display_min_cluster_support}")
                if p.display_min_top_raw_fraction > 0 and not cond_top.iloc[i]:
                    fail.append(f"top_raw_fraction_lt_{p.display_min_top_raw_fraction}")
                if p.display_min_representative_raw_support > 0 and not cond_rep.iloc[i]:
                    fail.append(f"representative_raw_support_lt_{p.display_min_representative_raw_support}")
                if p.display_max_start_span_bp > 0 and not cond_start_span.iloc[i]:
                    fail.append(f"start_span_gt_{p.display_max_start_span_bp}")
                if p.display_max_end_span_bp > 0 and not cond_end_span.iloc[i]:
                    fail.append(f"end_span_gt_{p.display_max_end_span_bp}")
                reasons.append('pass' if not fail else ';'.join(fail))
            sc['display_filter_reason'] = reasons
            sc['include_high_confidence_gtf'] = include
    return sc, pd.DataFrame(merge_edges)


def clean_val(v: object) -> str:
    if pd.isna(v): return 'NA'
    if isinstance(v,float) and v.is_integer(): return str(int(v))
    if isinstance(v,(np.integer,)): return str(int(v))
    if isinstance(v,(np.floating,)): return f'{float(v):.6g}'
    return str(v).replace('"',"'")


def make_attrs(row: pd.Series, layer: str) -> str:
    keys=['mtfind_version','candidate_id','supercluster_id','cluster_origin','candidate_source','coordinate_basis','representative_method','segmentation_method','representative_start','representative_end','representative_length','seed_support','backfill_support','cluster_support','supercluster_support','representative_raw_support','top_raw_fraction','member_primary_cluster_count','seed_member_raw_pair_count','backfill_member_raw_pair_count','member_raw_pair_count','minimum_signal_support','effective_seed_min_raw_pair_support','backfill_max_raw_pair_support','backfill_assignment_rule','start_min','start_max','end_min','end_max','start_span','end_span','start_shift_min_from_rep','start_shift_max_from_rep','end_shift_min_from_rep','end_shift_max_from_rep','shape_class','display_filter_reason','biological_status','annotation_used']
    fields={'gene_id':clean_val(row.get('supercluster_id','NA')),'transcript_id':f"{clean_val(row.get('candidate_id','NA'))}.t1",'gtf_layer':layer}
    for k in keys:
        if k in row.index: fields[k]=clean_val(row[k])
    return ' '.join([f'{k} "{v}";' for k,v in fields.items()])


def write_gtf(df: pd.DataFrame, path: Path, layer: str, source: str, chrom: str) -> None:
    if df is None or len(df)==0:
        path.write_text('', encoding='utf-8'); return
    if 'include_high_confidence_gtf' in df.columns:
        keep = df[df['include_high_confidence_gtf'].astype(bool)]
    else:
        keep = df
    if len(keep)==0:
        path.write_text('', encoding='utf-8'); return
    keep=keep.sort_values(['strand','representative_start','representative_end']).reset_index(drop=True)
    with path.open('w', encoding='utf-8') as f:
        for _,r in keep.iterrows():
            f.write('\t'.join(map(str,[chrom,source,'exon',int(r['representative_start']),int(r['representative_end']),'.',r['strand'],'.',make_attrs(r,layer)]))+'\n')


def load_plot_annotation_reference(path_text: str, fmt: str, chrom_name: str) -> pd.DataFrame:
    if not path_text:
        return pd.DataFrame(columns=['record_id','chrom','strand','start','end'])
    path=Path(path_text)
    if not path.exists():
        raise FileNotFoundError(f'Plot annotation reference not found: {path}')
    ref=read_coordinate_file(path, fmt=fmt, id_prefix='plot_ref')
    if ref.empty:
        return ref
    target=normalize_chrom_name(chrom_name)
    ref=ref[ref['chrom'].map(normalize_chrom_name).eq(target)].copy()
    return ref


def plot_overview(raw_pairs: pd.DataFrame, main: pd.DataFrame, outdir: Path, p: Params, annotation_ref: pd.DataFrame | None = None) -> None:
    raw=raw_pairs[raw_pairs['raw_support']>=p.plot_min_raw_pair_support].copy()
    if len(raw)>p.plot_max_raw_points:
        raw=raw.sort_values(['raw_support','start','end'], ascending=[False,True,True]).head(p.plot_max_raw_points)
    main_keep=main[main['include_high_confidence_gtf'].astype(bool)].copy() if len(main) else pd.DataFrame(columns=['strand','representative_start','representative_end'])
    ann=annotation_ref.copy() if annotation_ref is not None and len(annotation_ref) else pd.DataFrame()
    for strand in ['+','-']:
        fig,ax=plt.subplots(figsize=(8,8),dpi=180)
        rw=raw[raw['strand']==strand]
        if len(rw):
            sizes=np.clip(np.log10(rw['raw_support']+1)*3.5,1.5,18)
            ax.scatter(rw['start'],rw['end'],s=sizes,alpha=0.18,linewidths=0,label='raw endpoint-pairs')
        aw=ann[ann['strand']==strand] if len(ann) else pd.DataFrame()
        if len(aw):
            ax.scatter(aw['start'],aw['end'],s=80,marker='o',facecolors='none',edgecolors='tab:red',linewidths=1.2,label='annotation/reference')
        mw=main_keep[main_keep['strand']==strand]
        if len(mw): ax.scatter(mw['representative_start'],mw['representative_end'],s=48,marker='x',label='annotation-free high-confidence')
        ax.set_xlim(1,p.mt_length); ax.set_ylim(1,p.mt_length)
        ax.set_xlabel('start coordinate'); ax.set_ylabel('end coordinate')
        ax.set_title(f'{p.sample_id} annotation-free adaptive-depth endpoint-pair overview, strand {strand}')
        handles, labels = ax.get_legend_handles_labels()
        if handles: ax.legend(loc='best',fontsize=7)
        ax.grid(alpha=0.15); fig.tight_layout()
        suffix='plus' if strand=='+' else 'minus'
        fig.savefig(outdir/f'overview_annotation_free_adaptive_depth_{suffix}.png')
        plt.close(fig)


def plot_slabs(raw_pairs: pd.DataFrame, main: pd.DataFrame, out_pdf: Path, p: Params, annotation_ref: pd.DataFrame | None = None) -> None:
    windows=[]; start=1; window=4000; overlap=500; step=window-overlap
    while start<=p.mt_length:
        end=min(p.mt_length,start+window-1); windows.append((start,end))
        if end>=p.mt_length: break
        start+=step
    raw=raw_pairs[raw_pairs['raw_support']>=p.plot_min_raw_pair_support].copy()
    main_keep=main[main['include_high_confidence_gtf'].astype(bool)].copy() if len(main) else pd.DataFrame(columns=['strand','representative_start','representative_end'])
    ann=annotation_ref.copy() if annotation_ref is not None and len(annotation_ref) else pd.DataFrame()
    with PdfPages(out_pdf) as pdf:
        for strand in ['+','-']:
            for w0,w1 in windows:
                rw=raw[(raw['strand']==strand)&raw['start'].between(w0,w1)]
                mw=main_keep[(main_keep['strand']==strand)&main_keep['representative_start'].between(w0,w1)]
                aw=ann[(ann['strand']==strand)&ann['start'].between(w0,w1)] if len(ann) else pd.DataFrame()
                vals=[]
                if len(rw): vals += rw['end'].astype(int).tolist()
                if len(mw): vals += mw['representative_end'].astype(int).tolist()
                if len(aw): vals += aw['end'].astype(int).tolist()
                ymin,ymax=(max(1,min(vals)-200),min(p.mt_length,max(vals)+200)) if vals else (1,p.mt_length)
                fig,ax=plt.subplots(figsize=(10,6),dpi=140)
                if len(rw):
                    sizes=np.clip(np.log10(rw['raw_support']+1)*4.0,1.5,20)
                    ax.scatter(rw['start'],rw['end'],s=sizes,alpha=0.18,linewidths=0,label='raw')
                if len(aw):
                    ax.scatter(aw['start'],aw['end'],s=80,marker='o',facecolors='none',edgecolors='tab:red',linewidths=1.2,label='annotation/reference')
                if len(mw): ax.scatter(mw['representative_start'],mw['representative_end'],s=48,marker='x',label='main')
                ax.set_xlim(w0,w1); ax.set_ylim(ymin,ymax)
                ax.set_xlabel('start coordinate'); ax.set_ylabel('end coordinate')
                ax.set_title(f'{p.sample_id} annotation-free adaptive-depth strand {strand}, start slab {w0}-{w1}')
                handles, labels = ax.get_legend_handles_labels()
                if handles: ax.legend(fontsize=7,loc='best')
                ax.grid(alpha=0.15); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def audit_summary(rows: List[Tuple[str,object]], path: Path) -> pd.DataFrame:
    df=pd.DataFrame(rows, columns=['metric','value']); df.to_csv(path,sep='\t',index=False); return df


def finalize_manifest(root: Path, params: Params) -> None:
    rows=[('params_json',json.dumps(asdict(params),sort_keys=True,ensure_ascii=False))]
    for path in sorted(root.rglob('*')):
        if path.is_file() and path.relative_to(root).as_posix()!='package_manifest.tsv':
            rows.append((path.relative_to(root).as_posix(), sha256_file(path)))
    pd.DataFrame(rows, columns=['item','sha256_or_value']).to_csv(root/'package_manifest.tsv',sep='\t',index=False)



def write_run_parameters(params: Params, path_tsv: Path, path_json: Path | None = None) -> pd.DataFrame:
    """Write the actual effective runtime parameters used by this run."""
    rows = []
    for k, v in asdict(params).items():
        rows.append({"parameter": k, "value": v})
    df = pd.DataFrame(rows)
    path_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path_tsv, sep='\t', index=False)
    if path_json is not None:
        path_json.parent.mkdir(parents=True, exist_ok=True)
        path_json.write_text(json.dumps(asdict(params), indent=2, sort_keys=True, ensure_ascii=False), encoding='utf-8')
    return df


def add_param_cli(ap: argparse.ArgumentParser) -> None:
    """Expose all tunable Params fields to the terminal while preserving defaults."""
    ap.add_argument('--density-min-raw-pair-support', type=int, default=Params.density_min_raw_pair_support,
                    help=f'Default: {Params.density_min_raw_pair_support}. Additional raw-support floor for seed endpoint pairs.')
    ap.add_argument('--density-seed-min-raw-pair-support', type=int, default=Params.density_seed_min_raw_pair_support,
                    help=f'Default: {Params.density_seed_min_raw_pair_support}. Minimum raw support for a pair to establish 2D seed connectivity.')
    ap.add_argument('--density-backfill-max-raw-pair-support', type=int, default=Params.density_backfill_max_raw_pair_support,
                    help=f'Default: {Params.density_backfill_max_raw_pair_support}. Maximum raw support eligible for non-connecting 1–4-style backfill.')
    ap.add_argument('--density-connect-start-bp', type=int, default=Params.density_connect_start_bp,
                    help=f'Default: {Params.density_connect_start_bp}. Start-axis connection radius for 2D connected components.')
    ap.add_argument('--density-connect-end-bp', type=int, default=Params.density_connect_end_bp,
                    help=f'Default: {Params.density_connect_end_bp}. End-axis connection radius for 2D connected components.')
    ap.add_argument('--density-min-cluster-support', type=int, default=Params.density_min_cluster_support,
                    help=f'Default: {Params.density_min_cluster_support}. Minimum total support for retaining a primary density cluster.')
    ap.add_argument('--display-min-cluster-support', type=int, default=Params.display_min_cluster_support,
                    help=f'Default: {Params.display_min_cluster_support}. Minimum supercluster support written to high-confidence GTF.')
    ap.add_argument('--display-min-top-raw-fraction', type=float, default=Params.display_min_top_raw_fraction,
                    help=f'Default: {Params.display_min_top_raw_fraction}. Minimum representative_raw_support/supercluster_support for writing high-confidence GTF. 0 disables.')
    ap.add_argument('--display-min-representative-raw-support', type=int, default=Params.display_min_representative_raw_support,
                    help=f'Default: {Params.display_min_representative_raw_support}. Minimum raw support of the representative endpoint pair for writing high-confidence GTF. 0 disables.')
    ap.add_argument('--display-max-start-span-bp', type=int, default=Params.display_max_start_span_bp,
                    help=f'Default: {Params.display_max_start_span_bp}. Maximum supercluster start-coordinate span for writing high-confidence GTF. 0 disables.')
    ap.add_argument('--display-max-end-span-bp', type=int, default=Params.display_max_end_span_bp,
                    help=f'Default: {Params.display_max_end_span_bp}. Maximum supercluster end-coordinate span for writing high-confidence GTF. 0 disables.')
    ap.add_argument('--compact-delta-start-bp', type=int, default=Params.compact_delta_start_bp,
                    help=f'Default: {Params.compact_delta_start_bp}. Representative start distance for compact supercluster merge.')
    ap.add_argument('--compact-delta-end-bp', type=int, default=Params.compact_delta_end_bp,
                    help=f'Default: {Params.compact_delta_end_bp}. Representative end distance for compact supercluster merge.')
    ap.add_argument('--shared-endpoint-delta-bp', type=int, default=Params.shared_endpoint_delta_bp,
                    help=f'Default: {Params.shared_endpoint_delta_bp}. Tolerance for declaring shared start or shared end between clusters.')
    ap.add_argument('--shadow-other-endpoint-delta-bp', type=int, default=Params.shadow_other_endpoint_delta_bp,
                    help=f'Default: {Params.shadow_other_endpoint_delta_bp}. Other-endpoint tolerance for same-start/same-end shadow merge.')
    ap.add_argument('--containment-min', type=float, default=Params.containment_min,
                    help=f'Default: {Params.containment_min}. Minimum interval containment/overlap fraction for shadow/series merge.')
    ap.add_argument('--plot-min-raw-pair-support', type=int, default=Params.plot_min_raw_pair_support,
                    help=f'Default: {Params.plot_min_raw_pair_support}. Raw pair support cutoff for plots; ignored in simple mode.')
    ap.add_argument('--plot-max-raw-points', type=int, default=Params.plot_max_raw_points,
                    help=f'Default: {Params.plot_max_raw_points}. Maximum raw endpoint-pair points in plots; ignored in simple mode.')
    ap.add_argument('--plot-annotation-reference', type=Path, default=None,
                    help='Optional GTF/SSC reference to overlay on full-mode 2D plots. This is display-only and never used for discovery. If omitted, --benchmark-reference is used when benchmark is enabled.')
    ap.add_argument('--plot-annotation-reference-format', choices=['auto','gtf','ssc'], default=Params.plot_annotation_reference_format,
                    help=f'Default: {Params.plot_annotation_reference_format}. Parser for --plot-annotation-reference.')
    ap.add_argument('--adaptive-depth-mode', choices=['fixed','auto','low','medium','high'], default=Params.adaptive_depth_mode,
                    help='Default: fixed. fixed preserves explicit CLI/default parameters. auto profiles each chrom and selects low/medium/high depth-adaptive parameters.')
    ap.add_argument('--adaptive-depth-enable', action='store_true',
                    help='Shortcut for --adaptive-depth-mode auto unless --adaptive-depth-mode is explicitly set.')
    ap.add_argument('--adaptive-raw-sqrt-divisor', type=float, default=Params.adaptive_raw_sqrt_divisor,
                    help=f'Default: {Params.adaptive_raw_sqrt_divisor}. raw threshold = clamp(round(sqrt(total_support)/divisor), 3, 25).')
    ap.add_argument('--adaptive-cluster-fraction', type=float, default=Params.adaptive_cluster_fraction,
                    help=f'Default: {Params.adaptive_cluster_fraction}. cluster threshold = clamp(round(total_support*fraction), 20, 200).')
    ap.add_argument('--adaptive-display-fraction', type=float, default=Params.adaptive_display_fraction,
                    help=f'Default: {Params.adaptive_display_fraction}. display threshold = clamp(round(total_support*fraction), 50, 1000).')
    ap.add_argument('--adaptive-representative-fraction', type=float, default=Params.adaptive_representative_fraction,
                    help=f'Default: {Params.adaptive_representative_fraction}. representative raw threshold = clamp(round(total_support*fraction), 10, 1000).')
    ap.add_argument('--adaptive-low-cluster-fraction', type=float, default=Params.adaptive_low_cluster_fraction,
                    help=f'Default: {Params.adaptive_low_cluster_fraction}. Optional low-depth-only cluster fraction override; 0 disables.')
    ap.add_argument('--adaptive-low-display-fraction', type=float, default=Params.adaptive_low_display_fraction,
                    help=f'Default: {Params.adaptive_low_display_fraction}. Optional low-depth-only display fraction override; 0 disables.')
    ap.add_argument('--adaptive-low-representative-fraction', type=float, default=Params.adaptive_low_representative_fraction,
                    help=f'Default: {Params.adaptive_low_representative_fraction}. Optional low-depth-only representative fraction override; 0 disables.')
    ap.add_argument('--adaptive-low-total-support-max', type=int, default=Params.adaptive_low_total_support_max,
                    help=f'Default: {Params.adaptive_low_total_support_max}. auto mode selects low depth below this total support.')
    ap.add_argument('--adaptive-medium-total-support-max', type=int, default=Params.adaptive_medium_total_support_max,
                    help=f'Default: {Params.adaptive_medium_total_support_max}. auto mode selects medium depth below this total support.')
    # ── noise-robust depth estimation ──
    ap.add_argument('--noise-detection-method', choices=['auto','kneedle','mad','quantile'], default=Params.noise_detection_method,
                    help='Default: auto. Noise/signal boundary detection. auto tries kneedle → MAD → quantile fallback.')
    ap.add_argument('--noise-mad-k', type=float, default=Params.noise_mad_k,
                    help=f'Default: {Params.noise_mad_k}. k multiplier for MAD threshold (≈ k-σ equivalent).')
    ap.add_argument('--noise-min-knee-curvature', type=float, default=Params.noise_min_knee_curvature,
                    help=f'Default: {Params.noise_min_knee_curvature}. Min curvature to accept knee; below → fallback to MAD.')
    # ── multi-factor display scoring ──
    ap.add_argument('--display-score-mode', choices=['scored','binary'], default=Params.display_score_mode,
                    help='Default: scored. scored = multi-factor composite; binary = legacy absolute thresholds.')
    ap.add_argument('--display-score-threshold', type=float, default=Params.display_score_threshold,
                    help=f'Default: {Params.display_score_threshold}. Min composite score [0-1] to include in GTF.')
    ap.add_argument('--display-score-weight-support', type=float, default=Params.display_score_weight_support,
                    help=f'Default: {Params.display_score_weight_support}. Weight of absolute support sufficiency.')
    ap.add_argument('--display-score-weight-purity', type=float, default=Params.display_score_weight_purity,
                    help=f'Default: {Params.display_score_weight_purity}. Weight of signal purity (top_raw_fraction).')
    ap.add_argument('--display-score-weight-peak', type=float, default=Params.display_score_weight_peak,
                    help=f'Default: {Params.display_score_weight_peak}. Weight of peak signal strength (vs noise floor).')
    ap.add_argument('--display-score-weight-rank', type=float, default=Params.display_score_weight_rank,
                    help=f'Default: {Params.display_score_weight_rank}. Weight of relative rank (vs median cluster).')
    ap.add_argument('--display-score-relative-factor', type=float, default=Params.display_score_relative_factor,
                    help=f'Default: {Params.display_score_relative_factor}. effective_threshold = max(base, median_score × this). 0 disables.')


def build_params_from_args(args: argparse.Namespace) -> Params:
    return Params(
        sample_id=args.sample_id,
        chrom_name=args.chrom,
        mt_length=args.sequence_length,
        density_min_raw_pair_support=args.density_min_raw_pair_support,
        density_seed_min_raw_pair_support=args.density_seed_min_raw_pair_support,
        density_backfill_max_raw_pair_support=args.density_backfill_max_raw_pair_support,
        density_connect_start_bp=args.density_connect_start_bp,
        density_connect_end_bp=args.density_connect_end_bp,
        density_min_cluster_support=args.density_min_cluster_support,
        display_min_cluster_support=args.display_min_cluster_support,
        display_min_top_raw_fraction=args.display_min_top_raw_fraction,
        display_min_representative_raw_support=args.display_min_representative_raw_support,
        display_max_start_span_bp=args.display_max_start_span_bp,
        display_max_end_span_bp=args.display_max_end_span_bp,
        compact_delta_start_bp=args.compact_delta_start_bp,
        compact_delta_end_bp=args.compact_delta_end_bp,
        shared_endpoint_delta_bp=args.shared_endpoint_delta_bp,
        shadow_other_endpoint_delta_bp=args.shadow_other_endpoint_delta_bp,
        containment_min=args.containment_min,
        plot_min_raw_pair_support=args.plot_min_raw_pair_support,
        plot_max_raw_points=args.plot_max_raw_points,
        plot_annotation_reference=str(args.plot_annotation_reference) if args.plot_annotation_reference is not None else '',
        plot_annotation_reference_format=args.plot_annotation_reference_format,
        run_mode=args.run_mode,
        adaptive_depth_mode='auto' if args.adaptive_depth_enable and args.adaptive_depth_mode == Params.adaptive_depth_mode else args.adaptive_depth_mode,
        adaptive_raw_sqrt_divisor=args.adaptive_raw_sqrt_divisor,
        adaptive_cluster_fraction=args.adaptive_cluster_fraction,
        adaptive_display_fraction=args.adaptive_display_fraction,
        adaptive_representative_fraction=args.adaptive_representative_fraction,
        adaptive_low_cluster_fraction=args.adaptive_low_cluster_fraction,
        adaptive_low_display_fraction=args.adaptive_low_display_fraction,
        adaptive_low_representative_fraction=args.adaptive_low_representative_fraction,
        adaptive_low_total_support_max=args.adaptive_low_total_support_max,
        adaptive_medium_total_support_max=args.adaptive_medium_total_support_max,
        noise_detection_method=args.noise_detection_method,
        noise_mad_k=args.noise_mad_k,
        noise_min_knee_curvature=args.noise_min_knee_curvature,
        display_score_mode=args.display_score_mode,
        display_score_threshold=args.display_score_threshold,
        display_score_weight_support=args.display_score_weight_support,
        display_score_weight_purity=args.display_score_weight_purity,
        display_score_weight_peak=args.display_score_weight_peak,
        display_score_weight_rank=args.display_score_weight_rank,
        display_score_relative_factor=args.display_score_relative_factor,
    )

def main() -> None:
    ap=argparse.ArgumentParser(
        description=f'mtfind {MTFIND_VERSION}: annotation-free endpoint-pair discovery display pipeline with optional high-confidence display filters and depth-adaptive parameters. '
                    'All tunable parameters are exposed through CLI; omitted parameters use defaults.'
    )
    ap.add_argument('--version', action='version', version=f'mtfind {MTFIND_VERSION}')
    ap.add_argument('--ssc', type=Path, required=True, help='Raw SSC/SSC2 endpoint file')
    ap.add_argument('--outdir', type=Path, required=True, help='Output directory')
    ap.add_argument('--sample-id', default=Params.sample_id, help=f'Default: {Params.sample_id}')
    ap.add_argument('--chrom', default=Params.chrom_name, help='Target chromosome/contig in SSC. Use auto to select the most frequent + or - strand contig.')
    ap.add_argument('--sequence-length', '--mt-length', dest='sequence_length', type=int, default=Params.mt_length,
                    help=f'Default: {Params.mt_length}. Target contig length used for filtering and plot axes.')
    ap.add_argument('--run-mode', choices=['simple','data','full'], default='full',
                    help='simple: generate GTF + actual runtime parameter files only. data: generate all non-plot outputs. full: generate all tables, plots, logs, manifest, and parameter files.')
    add_param_cli(ap)
    add_benchmark_cli(ap)
    args=ap.parse_args()
    p=build_params_from_args(args)
    if not p.plot_annotation_reference and args.benchmark_enable and args.benchmark_reference is not None:
        p.plot_annotation_reference=str(args.benchmark_reference)
        p.plot_annotation_reference_format=args.benchmark_reference_format
    dirs=ensure_dirs(args.outdir)

    # Copy inputs/code only in data/full modes. Simple mode intentionally keeps output minimal.
    if p.run_mode in {'data', 'full'}:
        ssc_dst=dirs['inputs']/args.ssc.name
        code_dst=dirs['code']/'annotation_free_adaptive_depth.py'
        benchmark_module_src=Path(__file__).with_name('gtf_coordinate_benchmark.py')
        benchmark_module_dst=dirs['code']/'gtf_coordinate_benchmark.py'
        if args.ssc.resolve()!=ssc_dst.resolve(): shutil.copy2(args.ssc,ssc_dst)
        if Path(__file__).resolve()!=code_dst.resolve(): shutil.copy2(Path(__file__).resolve(),code_dst)
        if benchmark_module_src.exists() and benchmark_module_src.resolve()!=benchmark_module_dst.resolve(): shutil.copy2(benchmark_module_src,benchmark_module_dst)

    raw_pairs=load_raw_endpoint_pairs(args.ssc,p)
    apply_adaptive_depth_params(raw_pairs,p)
    primary_clusters=connected_components_2d(raw_pairs,p)
    superclusters, merge_edges=merge_clusters(primary_clusters,p,'AFSC')
    high=superclusters[superclusters['include_high_confidence_gtf'].astype(bool)].copy() if len(superclusters) else pd.DataFrame()

    # Always write GTF and actual runtime parameters.
    default_prediction_gtf = dirs['gtf']/f'{p.sample_id}.annotation_free.adaptive_depth.gtf'
    write_gtf(superclusters, default_prediction_gtf,'annotation_free_adaptive_depth','annotation_free', p.chrom_name)
    write_run_parameters(p, dirs['tables']/ 'run_parameters.tsv', dirs['tables']/ 'run_parameters.json')

    summary_rows=[
        ('mtfind_version',MTFIND_VERSION),
        ('run_mode',p.run_mode),
        ('runtime_dependency_on_previous_results','false'),
        ('runtime_inputs','raw_ssc_only'),
        ('target_chromosome',p.chrom_name),
        ('target_sequence_length',p.mt_length),
        ('annotation_used','false'),
        ('adaptive_depth_mode',p.adaptive_depth_mode),
        ('adaptive_depth_total_support',p.adaptive_depth_total_support),
        ('adaptive_depth_profile',p.adaptive_depth_profile),
        ('adaptive_raw_sqrt_divisor',p.adaptive_raw_sqrt_divisor),
        ('adaptive_cluster_fraction',p.adaptive_cluster_fraction),
        ('adaptive_display_fraction',p.adaptive_display_fraction),
        ('adaptive_representative_fraction',p.adaptive_representative_fraction),
        ('adaptive_low_cluster_fraction',p.adaptive_low_cluster_fraction),
        ('adaptive_low_display_fraction',p.adaptive_low_display_fraction),
        ('adaptive_low_representative_fraction',p.adaptive_low_representative_fraction),
        ('adaptive_low_total_support_max',p.adaptive_low_total_support_max),
        ('adaptive_medium_total_support_max',p.adaptive_medium_total_support_max),
        ('noise_detection_method',p.noise_detection_method),
        ('min_signal_support',p.min_signal_support),
        ('noise_floor',p.noise_floor),
        ('noise_floor_alias_for','min_signal_support'),
        ('noise_boundary_semantics',p.noise_boundary_semantics),
        ('noise_knee_curvature',p.noise_knee_curvature),
        ('noise_mad_k',p.noise_mad_k),
        ('signal_total_support',p.signal_total_support),
        ('noise_total_support',p.noise_total_support),
        ('noise_fraction',p.noise_fraction),
        ('signal_pair_count',p.signal_pair_count),
        ('noise_pair_count',p.noise_pair_count),
        ('raw_endpoint_pair_count',len(raw_pairs)),
        ('raw_endpoint_pair_total_support',int(raw_pairs['raw_support'].sum()) if len(raw_pairs) else 0),
        ('density_min_raw_pair_support',p.density_min_raw_pair_support),
        ('density_seed_min_raw_pair_support',p.density_seed_min_raw_pair_support),
        ('density_effective_seed_min_raw_pair_support',p.density_effective_seed_min_raw_pair_support),
        ('density_backfill_max_raw_pair_support',p.density_backfill_max_raw_pair_support),
        ('density_effective_backfill_max_raw_pair_support',p.density_effective_backfill_max_raw_pair_support),
        ('density_seed_pair_count',p.density_seed_pair_count),
        ('density_backfill_candidate_pair_count',p.density_backfill_candidate_pair_count),
        ('density_backfill_assigned_pair_count',p.density_backfill_assigned_pair_count),
        ('density_backfill_unassigned_pair_count',p.density_backfill_unassigned_pair_count),
        ('density_backfill_ambiguous_pair_count',p.density_backfill_ambiguous_pair_count),
        ('density_backfill_assigned_support',p.density_backfill_assigned_support),
        ('density_connect_start_bp',p.density_connect_start_bp),
        ('density_connect_end_bp',p.density_connect_end_bp),
        ('density_min_cluster_support',p.density_min_cluster_support),
        ('primary_density_cluster_count',len(primary_clusters)),
        ('supercluster_count',len(superclusters)),
        ('supercluster_merge_edge_count',len(merge_edges)),
        ('display_min_cluster_support',p.display_min_cluster_support),
        ('display_min_top_raw_fraction',p.display_min_top_raw_fraction),
        ('display_min_representative_raw_support',p.display_min_representative_raw_support),
        ('display_max_start_span_bp',p.display_max_start_span_bp),
        ('display_max_end_span_bp',p.display_max_end_span_bp),
        ('compact_delta_start_bp',p.compact_delta_start_bp),
        ('compact_delta_end_bp',p.compact_delta_end_bp),
        ('shared_endpoint_delta_bp',p.shared_endpoint_delta_bp),
        ('shadow_other_endpoint_delta_bp',p.shadow_other_endpoint_delta_bp),
        ('containment_min',p.containment_min),
        ('plot_min_raw_pair_support',p.plot_min_raw_pair_support),
        ('plot_max_raw_points',p.plot_max_raw_points),
        ('high_confidence_gtf_count',int(superclusters['include_high_confidence_gtf'].sum()) if len(superclusters) else 0),
        ('plot_annotation_reference',p.plot_annotation_reference if p.plot_annotation_reference else 'NA'),
        ('plot_annotation_reference_format',p.plot_annotation_reference_format),
        ('biological_interpretation_from_coordinates','false'),
    ]
    if args.benchmark_enable:
        if args.benchmark_reference is None:
            raise ValueError('--benchmark-enable requires --benchmark-reference')
        prediction_gtf = args.benchmark_prediction_gtf if args.benchmark_prediction_gtf is not None else default_prediction_gtf
        benchmark_outdir = args.benchmark_output_dir if args.benchmark_output_dir is not None else (args.outdir/'results'/'benchmark')
        benchmark_result = run_gtf_coordinate_benchmark(
            prediction_gtf=prediction_gtf,
            reference_file=args.benchmark_reference,
            outdir=benchmark_outdir,
            reference_format=args.benchmark_reference_format,
            tolerances=args.benchmark_tolerances,
            primary_tolerance=args.benchmark_primary_tolerance,
        )
        for key in ['tolerance_bp','predicted_gtf_count','reference_count','TP','FP','FN','precision','recall','F1','off_reference_candidate_count']:
            summary_rows.append((f'benchmark_{key}', benchmark_result.get(key, 'NA')))
        summary_rows.append(('benchmark_enabled','true'))
        summary_rows.append(('benchmark_reference', str(args.benchmark_reference)))
        summary_rows.append(('benchmark_prediction_gtf', str(prediction_gtf)))
    else:
        summary_rows.append(('benchmark_enabled','false'))

    audit_summary(summary_rows, dirs['tables']/ 'run_summary.tsv')

    if p.run_mode in {'data', 'full'}:
        raw_pairs.to_csv(dirs['tables']/ 'raw_endpoint_pair_count_table.tsv',sep='\t',index=False)
        primary_clusters.to_csv(dirs['tables']/ 'annotation_free_primary_density_clusters.tsv',sep='\t',index=False)
        superclusters.to_csv(dirs['tables']/ 'annotation_free_supercluster_table.tsv',sep='\t',index=False)
        merge_edges.to_csv(dirs['tables']/ 'annotation_free_supercluster_merge_edges.tsv',sep='\t',index=False)
        high.to_csv(dirs['tables']/ 'annotation_free_high_confidence_display_table.tsv',sep='\t',index=False)
        pd.DataFrame(summary_rows, columns=['metric','value']).to_csv(dirs['logs']/ 'run_summary.txt',sep='\t',index=False)
        if p.run_mode == 'full':
            annotation_ref=load_plot_annotation_reference(p.plot_annotation_reference, p.plot_annotation_reference_format, p.chrom_name)
            plot_overview(raw_pairs,superclusters,dirs['plots'],p,annotation_ref)
            plot_slabs(raw_pairs,superclusters,dirs['plots']/f'{p.sample_id}.annotation_free.adaptive_depth.start_slab.pdf',p,annotation_ref)
        cmd_parts = [
            "python3 code/annotation_free_adaptive_depth.py",
            f"  --ssc inputs/{args.ssc.name}",
            "  --outdir <output_dir>",
            f"  --sample-id {p.sample_id}",
            f"  --chrom {p.chrom_name}",
            f"  --sequence-length {p.mt_length}",
            f"  --run-mode {p.run_mode}",
            f"  --adaptive-depth-mode {p.adaptive_depth_mode}",
            f"  --adaptive-raw-sqrt-divisor {p.adaptive_raw_sqrt_divisor}",
            f"  --adaptive-cluster-fraction {p.adaptive_cluster_fraction}",
            f"  --adaptive-display-fraction {p.adaptive_display_fraction}",
            f"  --adaptive-representative-fraction {p.adaptive_representative_fraction}",
            f"  --adaptive-low-cluster-fraction {p.adaptive_low_cluster_fraction}",
            f"  --adaptive-low-display-fraction {p.adaptive_low_display_fraction}",
            f"  --adaptive-low-representative-fraction {p.adaptive_low_representative_fraction}",
            f"  --adaptive-low-total-support-max {p.adaptive_low_total_support_max}",
            f"  --adaptive-medium-total-support-max {p.adaptive_medium_total_support_max}",
            f"  --density-min-raw-pair-support {p.density_min_raw_pair_support}",
            f"  --density-seed-min-raw-pair-support {p.density_seed_min_raw_pair_support}",
            f"  --density-backfill-max-raw-pair-support {p.density_backfill_max_raw_pair_support}",
            f"  --density-connect-start-bp {p.density_connect_start_bp}",
            f"  --density-connect-end-bp {p.density_connect_end_bp}",
            f"  --density-min-cluster-support {p.density_min_cluster_support}",
            f"  --display-min-cluster-support {p.display_min_cluster_support}",
            f"  --display-min-top-raw-fraction {p.display_min_top_raw_fraction}",
            f"  --display-min-representative-raw-support {p.display_min_representative_raw_support}",
            f"  --display-max-start-span-bp {p.display_max_start_span_bp}",
            f"  --display-max-end-span-bp {p.display_max_end_span_bp}",
            f"  --compact-delta-start-bp {p.compact_delta_start_bp}",
            f"  --compact-delta-end-bp {p.compact_delta_end_bp}",
            f"  --shared-endpoint-delta-bp {p.shared_endpoint_delta_bp}",
            f"  --shadow-other-endpoint-delta-bp {p.shadow_other_endpoint_delta_bp}",
            f"  --containment-min {p.containment_min}",
            f"  --plot-annotation-reference {p.plot_annotation_reference}" if p.plot_annotation_reference else "",
            f"  --plot-annotation-reference-format {p.plot_annotation_reference_format}",
        ]
        cmd_parts=[x for x in cmd_parts if x]
        if args.benchmark_enable:
            cmd_parts.extend([
                "  --benchmark-enable",
                f"  --benchmark-reference {args.benchmark_reference}",
                f"  --benchmark-reference-format {args.benchmark_reference_format}",
                f"  --benchmark-tolerances {args.benchmark_tolerances}",
                f"  --benchmark-primary-tolerance {args.benchmark_primary_tolerance}",
            ])
        run_command = "# From package root annotation-free adaptive-depth:\n" + " \\\n".join(cmd_parts) + "\n"
        (dirs['code']/ 'run_command.txt').write_text(run_command,encoding='utf-8')
        readme=(
            f"# mtfind {MTFIND_VERSION}: annotation-free adaptive-depth for {p.sample_id}\n\n"
            "This version runs discovery from raw SSC/SSC2 only. Annotation/reference inputs are optional display overlays or post-hoc benchmarks and never participate in candidate generation. "
            "All tunable parameters, including high-confidence display filters, are exposed by terminal CLI; omitted values use defaults. "
            "`--adaptive-depth-mode auto` profiles each chrom and derives conservative depth-aware support thresholds.\n\n"
            "In 0.311, Kneedle produces `min_signal_support`: raw endpoint pairs below this boundary are noise, and pairs at or above it are signal. "
            "Only high-support pairs establish 2D connectivity; nearby 1–4-support signal pairs can backfill exactly one seed component and never bridge components.\n\n"
            "Run modes:\n"
            "- `--run-mode simple`: generate only GTF + actual runtime parameter files + minimal run summary.\n"
            "- `--run-mode data`: generate GTF, tables, logs, manifest, copied inputs/code, and parameter files; no plots are rendered.\n"
            "- `--run-mode full`: generate GTF, tables, plots, logs, manifest, copied inputs/code, parameter files, and optional annotation/reference overlays in 2D plots.\n\n"
            "GTF is a coordinate-candidate display track, not validated transcript annotation.\n\nOptional benchmark module:\n- Use `--benchmark-enable --benchmark-reference <file>` to compare the generated GTF against a GTF/SSC coordinate reference after discovery has finished.\n- Benchmark and plot annotations are post-hoc only and never participate in candidate generation.\n\n"
            "See `code/run_command.txt`.\n"
        )
        (args.outdir/'README_annotation_free_adaptive_depth.md').write_text(readme,encoding='utf-8')
        finalize_manifest(args.outdir,p)

    print(pd.DataFrame(summary_rows, columns=['metric','value']).to_string(index=False))

if __name__=='__main__':
    main()
