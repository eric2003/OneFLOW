#!/usr/bin/env python3
"""Regression tests for the MPI rank-utilization timing gate."""

import pathlib
import subprocess
import sys
import tempfile
import unittest


VERIFIER = pathlib.Path(__file__).with_name("verify-rank-utilization.py")


class RankUtilizationVerifierTest(unittest.TestCase):
    def run_verifier(self, rows_by_rank, expected_ranks=None):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            for rank, rows in rows_by_rank.items():
                path = root / f"breakdown-{rank}.tsv"
                lines = ["category\telapsed_ms"]
                lines.extend(f"{category}\t{elapsed}" for category, elapsed in rows)
                path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            ranks = expected_ranks
            if ranks is None:
                ranks = len(rows_by_rank)
            return subprocess.run(
                [
                    sys.executable,
                    str(VERIFIER),
                    "--expected-ranks",
                    str(ranks),
                    "--label",
                    "unit-test",
                    "--pattern",
                    str(root / "breakdown-*.tsv"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )

    def test_all_ranks_with_solver_work_pass(self):
        result = self.run_verifier(
            {
                0: [("initialization", 20.0), ("gradient", 100.0)],
                1: [("initialization", 20.0), ("calc_time_step", 80.0)],
            }
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("RANK_UTILIZATION_PASS", result.stdout)
        self.assertIn("active_ranks=2", result.stdout)

    def test_missing_rank_file_fails(self):
        result = self.run_verifier(
            {0: [("gradient", 100.0)]},
            expected_ranks=2,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("missing=[1]", result.stderr)

    def test_initialization_and_mpi_only_rank_fails(self):
        result = self.run_verifier(
            {
                0: [("gradient", 100.0)],
                1: [("initialization", 40.0), ("mpi_interface", 10.0)],
            }
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("active_ranks=1", result.stderr)
        self.assertIn("idle_ranks=[1]", result.stderr)

    def test_tiny_scheduler_overhead_does_not_count_as_active(self):
        result = self.run_verifier(
            {
                0: [("gradient", 100.0)],
                1: [("calc_time_step", 0.009)],
            }
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("idle_ranks=[1]", result.stderr)


if __name__ == "__main__":
    unittest.main()
