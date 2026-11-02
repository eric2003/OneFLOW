"""Tests for the evidence-backed parameter consumer index."""
from __future__ import annotations

import unittest

from tools.oneflow_gui.parameter_consumer_catalog import (
    ConsumerReadKind,
    ParameterConsumerCatalog,
)


class ParameterConsumerCatalogTest(unittest.TestCase):
    def setUp(self):
        self.catalog = ParameterConsumerCatalog()

    def test_known_parameter_returns_source_backed_consumer(self):
        consumers = self.catalog.find("MAXSTEPS")

        self.assertEqual(len(consumers), 1)
        self.assertEqual(consumers[0].source_path, "codes/global/src/Iteration.cpp")
        self.assertEqual(consumers[0].symbol, "Iteration::Init")
        self.assertEqual(consumers[0].read_kind, ConsumerReadKind.REQUIRED)
        self.assertIn("GetDataValue", consumers[0].read_expression)
        self.assertIn("termination", consumers[0].behavior)

    def test_parameter_can_have_multiple_runtime_read_sites(self):
        consumers = self.catalog.find("time_integral")

        self.assertEqual(
            {(item.source_path, item.symbol) for item in consumers},
            {
                ("codes/global/src/Ctrl.cpp", "Ctrl::Init"),
                ("codes/ns/src/NsCtrl.cpp", "NsCtrl::Init"),
            },
        )

    def test_lazy_consumer_is_not_misreported_as_eager(self):
        consumers = self.catalog.find("ptTol")

        self.assertEqual(len(consumers), 1)
        self.assertEqual(consumers[0].read_kind, ConsumerReadKind.LAZY)
        self.assertIn("first Tolerance::GetTol", consumers[0].behavior)

    def test_unknown_parameter_has_no_invented_consumer(self):
        self.assertEqual(self.catalog.find("vendorPrivateParameter"), ())
        self.assertFalse(self.catalog.has_verified_consumer("vendorPrivateParameter"))

    def test_multigrid_consumer_records_normalization_behavior(self):
        consumers = self.catalog.find("nmg")

        self.assertEqual(len(consumers), 1)
        self.assertIn("nmg <= 1", consumers[0].behavior)
        self.assertIn("smoothing", consumers[0].behavior)


if __name__ == "__main__":
    unittest.main()
