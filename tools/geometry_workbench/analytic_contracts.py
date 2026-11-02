"""Analytic intermediate contracts (no C++).

These evaluate published closed-form expressions for workbench tables.
They are not a reimplementation of production solvers.
"""

from __future__ import annotations

from typing import Any


def sa_fv1(chi: float, cv1: float = 7.1) -> float:
    chi3 = chi * chi * chi
    cv13 = cv1 * cv1 * cv1
    return chi3 / (chi3 + cv13)


def run_analytic_sa_fv1(scenario: dict[str, Any]) -> dict[str, Any]:
    cv1 = float(scenario.get("input", {}).get("cv1", 7.1))
    chis = scenario.get("input", {}).get("chi") or []
    fv1 = [sa_fv1(float(c), cv1) for c in chis]
    return {"fv1": fv1}


def run_analytic(scenario: dict[str, Any]) -> dict[str, Any] | None:
    runner = scenario.get("runner")
    if runner == "analytic_sa_fv1":
        return run_analytic_sa_fv1(scenario)
    return None
