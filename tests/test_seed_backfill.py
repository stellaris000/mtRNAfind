#!/usr/bin/env python3
"""Focused regression tests for mtfind 0.311 seed/backfill clustering."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd


CODE_DIR = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE_DIR))
import annotation_free_adaptive_depth as mtfind  # noqa: E402


def raw_pairs(rows: list[tuple[str, int, int, int, str]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows, columns=["strand", "start", "end", "raw_support", "raw_pair_id"])
    frame["length"] = frame["end"] - frame["start"] + 1
    return frame


class SeedBackfillTests(unittest.TestCase):
    def test_kneedle_boundary_is_minimum_signal_support(self) -> None:
        frame = pd.DataFrame({"raw_support": [1, 1, 2, 3, 8, 13, 21]})
        params = mtfind.Params(noise_detection_method="kneedle")

        result = mtfind.detect_noise_floor(frame, params)

        boundary = result["min_signal_support"]
        self.assertEqual(params.min_signal_support, boundary)
        self.assertEqual(params.noise_floor, boundary)
        self.assertEqual(
            params.signal_pair_count,
            int((frame["raw_support"] >= boundary).sum()),
        )
        self.assertEqual(
            params.noise_pair_count,
            int((frame["raw_support"] < boundary).sum()),
        )
        self.assertEqual(
            result["all_profile_data"]["noise_boundary_semantics"],
            mtfind.NOISE_BOUNDARY_SEMANTICS,
        )

    def test_low_support_pairs_backfill_one_seed_without_bridging(self) -> None:
        frame = raw_pairs([
            ("+", 100, 200, 10, "RP_A_SEED"),
            ("+", 110, 210, 10, "RP_B_SEED"),
            ("+", 102, 202, 4, "RP_A_4"),
            ("+", 103, 203, 1, "RP_A_1"),
            ("+", 108, 208, 2, "RP_B_2"),
            ("+", 105, 205, 4, "RP_AMBIGUOUS"),
        ])
        params = mtfind.Params(
            min_signal_support=1,
            noise_floor=1,
            density_min_raw_pair_support=5,
            density_seed_min_raw_pair_support=5,
            density_backfill_max_raw_pair_support=4,
            density_connect_start_bp=5,
            density_connect_end_bp=5,
            density_min_cluster_support=5,
        )

        clusters = mtfind.connected_components_2d(frame, params).sort_values("representative_start")

        self.assertEqual(len(clusters), 2)
        self.assertEqual(clusters["cluster_support"].tolist(), [15, 12])
        self.assertEqual(clusters["seed_support"].tolist(), [10, 10])
        self.assertEqual(clusters["backfill_support"].tolist(), [5, 2])
        self.assertEqual(params.density_backfill_candidate_pair_count, 4)
        self.assertEqual(params.density_backfill_assigned_pair_count, 3)
        self.assertEqual(params.density_backfill_unassigned_pair_count, 1)
        self.assertEqual(params.density_backfill_ambiguous_pair_count, 1)

    def test_pairs_below_minimum_signal_are_not_backfilled(self) -> None:
        frame = raw_pairs([
            ("+", 100, 200, 10, "RP_SEED"),
            ("+", 102, 202, 4, "RP_SIGNAL_4"),
            ("+", 103, 203, 2, "RP_NOISE_2"),
        ])
        params = mtfind.Params(
            min_signal_support=3,
            noise_floor=3,
            density_min_raw_pair_support=5,
            density_seed_min_raw_pair_support=5,
            density_backfill_max_raw_pair_support=4,
            density_connect_start_bp=5,
            density_connect_end_bp=5,
            density_min_cluster_support=5,
        )

        clusters = mtfind.connected_components_2d(frame, params)

        self.assertEqual(len(clusters), 1)
        self.assertEqual(int(clusters.iloc[0]["cluster_support"]), 14)
        self.assertEqual(params.density_backfill_candidate_pair_count, 1)
        self.assertEqual(params.density_backfill_assigned_pair_count, 1)


if __name__ == "__main__":
    unittest.main()
