from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from .vencat_reference import (
    vencat,
    vencat_c,
    vencat_eric,
    write_vencat_reference,
)


class VencatReferenceTest(unittest.TestCase):
    def test_retained_formula_variants(self) -> None:
        self.assertAlmostEqual(vencat(0.0), 1.0)
        self.assertAlmostEqual(vencat(1.0), 0.75)
        self.assertAlmostEqual(vencat_c(1.0, 2.0), 0.75)
        self.assertAlmostEqual(vencat_eric(1.0), 1.0)

    def test_writer_preserves_ten_zone_reference_layout(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "nested" / "vencat.dat"
            written = write_vencat_reference(output)

            self.assertEqual(written, output)
            lines = output.read_text(encoding="ascii").splitlines()
            self.assertEqual(
                lines[:5],
                [
                    'title="THE FLOW FIELD OF ONEFLOW"',
                    "variables=",
                    '"x"',
                    '"vencat"',
                    '"gvencat"',
                ],
            )
            zone_indices = [
                index for index, line in enumerate(lines)
                if line.startswith(" zone  i = ")
            ]
            self.assertEqual(len(zone_indices), 10)

            for zone_index, start in enumerate(zone_indices):
                end = (
                    zone_indices[zone_index + 1]
                    if zone_index + 1 < len(zone_indices)
                    else len(lines)
                )
                self.assertEqual(lines[start], " zone  i = 101")
                rows = [line.split() for line in lines[start + 1:end]]
                self.assertEqual(len(rows), 101)
                self.assertTrue(all(len(row) == 3 for row in rows))

                first = [float(value) for value in rows[0]]
                last = [float(value) for value in rows[-1]]
                coefficient = 2.0 + zone_index
                self.assertAlmostEqual(first[0], 0.0)
                self.assertAlmostEqual(first[1], 1.0)
                self.assertAlmostEqual(first[2], 0.0)
                self.assertAlmostEqual(last[0], 2.0)
                expected = vencat_c(2.0, coefficient)
                self.assertAlmostEqual(last[1], expected)
                self.assertAlmostEqual(last[2], 2.0 * expected)

    def test_writer_rejects_invalid_point_count_and_coefficients(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "vencat.dat"
            with self.assertRaises(ValueError):
                write_vencat_reference(output, point_count=1)
            with self.assertRaises(ValueError):
                write_vencat_reference(output, c_values=())
            with self.assertRaises(ValueError):
                write_vencat_reference(output, c_values=(0.0,))


if __name__ == "__main__":
    unittest.main()
