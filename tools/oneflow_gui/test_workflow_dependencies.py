"""Tests for explicit workflow case dependency semantics."""

from __future__ import annotations

import unittest

from tools.oneflow_gui.workflow_dependencies import validate_case_dependencies


class WorkflowDependenciesTest(unittest.TestCase):
    def test_classifies_independent_producer_consumer_and_chain_cases(self) -> None:
        graph = validate_case_dependencies(
            [
                {"id": "standalone", "directory": "grid/standalone"},
                {"id": "mesh", "directory": "grid/generated"},
                {"id": "solve", "directory": "test/cfd/solver", "depends_on": ["mesh"]},
                {"id": "post", "directory": "tools/post", "depends_on": ["solve"]},
            ]
        )
        self.assertEqual(graph.role("standalone"), "independent")
        self.assertEqual(graph.role("mesh"), "producer")
        self.assertEqual(graph.role("solve"), "consumer_and_producer")
        self.assertEqual(graph.role("post"), "consumer")
        self.assertEqual(graph.ordered_ids, ("standalone", "mesh", "solve", "post"))

    def test_dependencies_override_declaration_order(self) -> None:
        graph = validate_case_dependencies(
            [
                {"id": "solve", "directory": "solver", "depends_on": ["mesh"]},
                {"id": "mesh", "directory": "generator"},
                {"id": "independent", "directory": "independent"},
            ]
        )
        self.assertLess(graph.ordered_ids.index("mesh"), graph.ordered_ids.index("solve"))
        self.assertEqual(set(graph.ordered_ids), {"mesh", "solve", "independent"})

    def test_rejects_unknown_dependency(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown case"):
            validate_case_dependencies(
                [{"id": "solve", "directory": "solver", "depends_on": ["missing"]}]
            )

    def test_rejects_dependency_cycle(self) -> None:
        with self.assertRaisesRegex(ValueError, "contain a cycle"):
            validate_case_dependencies(
                [
                    {"id": "mesh", "directory": "grid", "depends_on": ["solve"]},
                    {"id": "solve", "directory": "solver", "depends_on": ["mesh"]},
                ]
            )

    def test_rejects_duplicate_ids_self_dependency_and_duplicate_edges(self) -> None:
        with self.assertRaisesRegex(ValueError, "Duplicate case id"):
            validate_case_dependencies(
                [
                    {"id": "same", "directory": "a"},
                    {"id": "same", "directory": "b"},
                ]
            )
        with self.assertRaisesRegex(ValueError, "cannot depend on itself"):
            validate_case_dependencies(
                [{"id": "self", "directory": "a", "depends_on": ["self"]}]
            )
        with self.assertRaisesRegex(ValueError, "duplicate dependency"):
            validate_case_dependencies(
                [
                    {"id": "mesh", "directory": "a"},
                    {"id": "solve", "directory": "b", "depends_on": ["mesh", "mesh"]},
                ]
            )

    def test_rejects_empty_graph_and_missing_directory(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one case"):
            validate_case_dependencies([])
        with self.assertRaisesRegex(ValueError, "non-empty directory"):
            validate_case_dependencies([{"id": "case"}])


if __name__ == "__main__":
    unittest.main()
