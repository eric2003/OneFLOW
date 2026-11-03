"""Tests for the evidence-backed parameter consumer index."""
from __future__ import annotations

import unittest
from pathlib import Path

from tools.verification.curve_catalog import list_theory_entries

from tools.oneflow_gui.parameter_consumer_catalog import (
    ConsumerReadKind,
    ParameterConsumerCatalog,
    validate_parameter_consumer_sources,
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
                ("codes/ins/src/INsCtrl.cpp", "INsCtrl::Init"),
                ("codes/turb/src/TurbCtrl.cpp", "TurbCtrl::Init"),
            },
        )

    def test_rk_stage_tracks_all_solver_control_read_sites(self):
        consumers = self.catalog.find("rk_stage")

        self.assertEqual(len(consumers), 4)
        self.assertTrue(all(item.read_kind is ConsumerReadKind.REQUIRED for item in consumers))
        self.assertTrue(all("1 / (stage + 1)" in item.behavior for item in consumers))

    def test_time_integration_behavior_links_to_dispatch_and_capability_sources(self):
        consumer = self.catalog.find("time_integral")[0]

        self.assertIn("MULTI_STAGE", consumer.behavior)
        self.assertIn("TimeIntegral.cpp", " ".join(consumer.behavior_source_paths))
        self.assertIn("EvaluateEulerRungeKuttaCapability", consumer.behavior_symbols)
        self.assertEqual(consumer.knowledge_ids, ("time_integral_rk_stage_weight",))
        self.assertEqual(
            self.catalog.find("rk_stage")[0].knowledge_ids,
            ("time_integral_rk_stage_weight",),
        )


    def test_cfl_ramp_parameters_link_to_iteration_index_and_schedule(self):
        for key, phrase in (
            ("cflst", "starting value"),
            ("cfled", "ending value"),
            ("ncfl", "ramp length"),
        ):
            with self.subTest(parameter=key):
                consumers = self.catalog.find(key)
                self.assertEqual(len(consumers), 1)
                consumer = consumers[0]
                self.assertEqual(consumer.source_path, "codes/global/src/Iteration.cpp")
                self.assertEqual(consumer.symbol, "Iteration::Init")
                self.assertIn(phrase, consumer.behavior)
                self.assertIn("CflAtIteration", consumer.behavior_symbols)
                self.assertIn("CflIterationIndex", consumer.behavior_symbols)

    def test_time_step_strategy_and_global_dt_are_traced_separately(self):
        strategy = self.catalog.find("timestepModel")
        self.assertEqual(len(strategy), 1)
        self.assertIn("0 local, 1 global", strategy[0].behavior)
        self.assertIn("UTimeStep::CalcLgTimeStep", strategy[0].behavior_symbols)

        global_steps = self.catalog.find("global_dt")
        self.assertEqual(len(global_steps), 4)
        self.assertEqual(
            {item.source_path for item in global_steps},
            {
                "codes/global/src/Ctrl.cpp",
                "codes/ns/src/NsCtrl.cpp",
                "codes/ins/src/INsCtrl.cpp",
                "codes/turb/src/TurbCtrl.cpp",
            },
        )
        self.assertTrue(all("Ctrl::pdt" in item.behavior for item in global_steps))
        self.assertTrue(all("UTimeStep::CalcGlobalTimeStep" in item.behavior_symbols for item in global_steps))


    def test_output_frequency_parameters_link_to_scheduler_consumers(self):
        for key, symbol in (
            ("nFieldSave", "kDumpRestartTaskName"),
            ("nVisualSave", "kVisualizationTaskName"),
            ("nResSave", "Iteration::ResOk"),
            ("nForceSave", "Iteration::ForceOk"),
        ):
            with self.subTest(parameter=key):
                consumers = self.catalog.find(key)
                self.assertEqual(len(consumers), 1)
                self.assertEqual(consumers[0].source_path, "codes/global/src/Iteration.cpp")
                self.assertIn(symbol, consumers[0].behavior_symbols)

    def test_catalog_read_and_behavior_evidence_resolves_to_current_sources(self):
        repository_root = Path(__file__).resolve().parents[2]

        issues = validate_parameter_consumer_sources(repository_root)

        self.assertEqual(issues, [], "\n".join(issues))

    def test_parameter_knowledge_references_resolve_in_formula_catalog(self):
        formula_ids = {entry.id for entry in list_theory_entries()}
        referenced_ids = {
            knowledge_id
            for consumer in self.catalog.consumers
            for knowledge_id in consumer.knowledge_ids
        }

        self.assertTrue(referenced_ids)
        self.assertEqual(referenced_ids - formula_ids, set())

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
