import sys
import unittest
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import analysis_v3


class WilsonIntervalTests(unittest.TestCase):
    def test_9_of_397(self):
        low, high = analysis_v3.wilson_interval_pct(9, 397)
        self.assertAlmostEqual(low, 1.1972, places=3)
        self.assertAlmostEqual(high, 4.2517, places=3)

    def test_20_of_397(self):
        low, high = analysis_v3.wilson_interval_pct(20, 397)
        self.assertAlmostEqual(low, 3.2846, places=3)
        self.assertAlmostEqual(high, 7.6528, places=3)

    def test_interval_contains_observed_proportion(self):
        low, high = analysis_v3.wilson_interval_pct(9, 397)
        p = 100 * 9 / 397
        self.assertLessEqual(low, p)
        self.assertGreaterEqual(high, p)


class CompetingRiskTests(unittest.TestCase):
    def test_aalen_johansen_with_competing_event(self):
        records = [(2.0, "restart"), (3.0, "icu_exit"), (5.0, "censor")]
        self.assertAlmostEqual(analysis_v3.cif_aj(records, 2.0), 1/3)
        self.assertAlmostEqual(analysis_v3.cif_aj(records, 3.0), 1/3)

    def test_restricted_mean_without_recorded_restart(self):
        records = [(2.0, "restart"), (24.0, "censor")]
        self.assertAlmostEqual(analysis_v3.rmst_without_restart(records, 24.0), 13.0)


class DeterministicClusterBootstrapTests(unittest.TestCase):
    def setUp(self):
        self.records = [
            {"subject_id": 101, "time": 2.0, "kind": "restart"},
            {"subject_id": 102, "time": 3.0, "kind": "icu_exit"},
            {"subject_id": 103, "time": 5.0, "kind": "censor"},
            {"subject_id": 104, "time": 4.0, "kind": "restart"},
            {"subject_id": 104, "time": 8.0, "kind": "censor"},
        ]

    def test_seeded_bootstrap_is_repeatable(self):
        a = analysis_v3.bootstrap_cluster(self.records, 24, reps=500, seed=77, metric="cif")
        b = analysis_v3.bootstrap_cluster(self.records, 24, reps=500, seed=77, metric="cif")
        self.assertEqual(a, b)

    def test_seeded_bootstrap_is_invariant_to_row_order(self):
        a = analysis_v3.bootstrap_cluster(self.records, 24, reps=500, seed=77, metric="cif")
        b = analysis_v3.bootstrap_cluster(list(reversed(self.records)), 24,
                                          reps=500, seed=77, metric="cif")
        self.assertEqual(a, b)

    def test_cluster_bootstrap_preserves_all_records_per_subject(self):
        # Subject 104 contributes two records; the bootstrap must resample that
        # subject as one cluster rather than treating the two stays independently.
        a = analysis_v3.bootstrap_cluster(self.records, 24, reps=500, seed=91, metric="cif")
        reordered = [self.records[4], self.records[2], self.records[0],
                     self.records[3], self.records[1]]
        b = analysis_v3.bootstrap_cluster(reordered, 24, reps=500, seed=91, metric="cif")
        self.assertEqual(a, b)


class VolumeIntegrationTests(unittest.TestCase):
    def test_independent_orders_are_summed(self):
        t0 = datetime(2020, 1, 1, 0, 0, 0)
        t1 = t0 + timedelta(hours=1)
        intervals = [(t0, t1, 10.0, "order-A"), (t0, t1, 5.0, "order-B")]
        self.assertAlmostEqual(
            analysis_v3.integrate_volume_by_order(intervals, t0, t1), 15.0
        )

    def test_overlapping_rows_within_same_order_are_deduplicated_by_max_rate(self):
        t0 = datetime(2020, 1, 1, 0, 0, 0)
        t1 = t0 + timedelta(hours=1)
        intervals = [(t0, t1, 10.0, "same-link"), (t0, t1, 20.0, "same-link")]
        self.assertAlmostEqual(
            analysis_v3.integrate_volume_by_order(intervals, t0, t1), 20.0
        )


if __name__ == "__main__":
    unittest.main()
