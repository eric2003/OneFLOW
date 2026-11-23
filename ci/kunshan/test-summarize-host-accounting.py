import importlib.util
import pathlib
import tempfile
import unittest


spec = importlib.util.spec_from_file_location(
    "accounting", pathlib.Path(__file__).with_name("summarize-host-accounting.py"))
accounting = importlib.util.module_from_spec(spec)
spec.loader.exec_module(accounting)


class HostAccountingTest(unittest.TestCase):
    def fixture(self, root, rank=0, steps=1000, total=1000.0):
        path = root / ("rank-%d.tsv.host.tsv" % rank)
        path.write_text("scope\tcount\tinclusive_ms\texclusive_ms\n"
                        "mg_run\t1\t%s\t0\nmg_inner_iteration\t%d\t%s\t%s\n" %
                        (total, steps, total, total))
        return path

    def test_missing_rank_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            self.fixture(root)
            with self.assertRaisesRegex(ValueError, "rank files"):
                accounting.read_profiles(str(root / "*.tsv"), 2, 1000)

    def test_iteration_count_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            self.fixture(root, steps=1)
            with self.assertRaisesRegex(ValueError, "iteration count"):
                accounting.read_profiles(str(root / "*.tsv"), 1, 1000)

    def test_double_counting_and_nonfinite_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            path = self.fixture(root)
            for value in ("50", "nan"):
                path.write_text(path.read_text().split("extra")[0] +
                                "extra\t1\t%s\t%s\n" % (value, value))
                with self.assertRaises(ValueError):
                    accounting.read_profiles(str(root / "*.tsv"), 1, 1000)

    def test_startup_subtraction_and_critical_rank(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            self.fixture(root, rank=0, total=900)
            self.fixture(root, rank=1, total=1000)
            long = accounting.read_profiles(str(root / "*.tsv"), 2, 1000)
            self.fixture(root, rank=0, steps=1, total=0.9)
            self.fixture(root, rank=1, steps=1, total=1.0)
            short = accounting.read_profiles(str(root / "*.tsv"), 2, 1)
            report = accounting.summarize(long, short, 1000, 1, 2500, 1501)
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["critical_rank"], 1)
            self.assertAlmostEqual(report["wall_ms_per_step"], 1)
            self.assertAlmostEqual(report["unexplained_fraction"], 0)
            report = accounting.summarize(long, short, 1000, 1, 3500, 1501)
            self.assertEqual(report["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
