"""Evidence-backed index of selected legacy parameter runtime consumers.

This registry records source-level reads that have been checked in the current
repository. It is intentionally a small, curated first slice, not an automated
claim that every parameter consumer has been discovered.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import re


class ConsumerReadKind(str, Enum):
    REQUIRED = "required"
    OPTIONAL = "optional"
    LAZY = "lazy"


@dataclass(frozen=True, slots=True)
class ParameterConsumer:
    """One source-backed runtime read of a parameter declaration."""

    parameter_name: str
    source_path: str
    symbol: str
    read_expression: str
    read_kind: ConsumerReadKind
    behavior: str
    applicability: str
    evidence_note: str
    behavior_source_paths: tuple[str, ...] = ()
    behavior_symbols: tuple[str, ...] = ()
    knowledge_ids: tuple[str, ...] = ()


# Keep source paths repository-relative and symbols precise enough to navigate.
# Entries describe reads, not every downstream use of the resulting member.
_CONSUMERS = (
    ParameterConsumer(
        "maxSteps", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<int>("maxSteps")', ConsumerReadKind.REQUIRED,
        "Initializes Iteration::maxSteps and maxIterSteps; SimuIterState::Running "
        "uses maxSteps as the outer-iteration termination limit when iexitflag == 0.",
        "Solver run initialization; step-count termination mode.",
        "Direct typed read and termination predicate are present in Iteration.cpp.",
    ),
    ParameterConsumer(
        "time_integral", "codes/global/src/Ctrl.cpp", "Ctrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the selected time-integration method in Ctrl::time_integral.",
        "Global solver control initialization.",
        "Direct typed read is present in Ctrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "time_integral", "codes/ns/src/NsCtrl.cpp", "NsCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the time-integration method in the compressible Navier-Stokes control.",
        "Compressible Navier-Stokes control initialization.",
        "Direct typed read is present in NsCtrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "time_integral", "codes/ins/src/INsCtrl.cpp", "INsCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the time-integration selector for incompressible Navier-Stokes control.",
        "Incompressible Navier-Stokes control initialization.",
        "Direct typed read is present in INsCtrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "time_integral", "codes/turb/src/TurbCtrl.cpp", "TurbCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the time-integration selector for turbulence control.",
        "Turbulence control initialization.",
        "Direct typed read is present in TurbCtrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "rk_stage", "codes/global/src/Ctrl.cpp", "Ctrl::Init",
        'GetDataValue<int>("rk_stage")', ConsumerReadKind.REQUIRED,
        "Builds rk_coef as 1 / (stage + 1); TimeIntegral::RungeKutta consumes the coefficients.",
        "Global solver control initialization; multi-stage integration path.",
        "Direct read and coefficient construction are present in Ctrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/global/src/Ctrl.cpp", "codes/ns/src/NsCtrl.cpp", "codes/ins/src/INsCtrl.cpp", "codes/turb/src/TurbCtrl.cpp"),
        ("Ctrl::Init", "NsCtrl::Init", "INsCtrl::Init", "TurbCtrl::Init", "TimeIntegral::RungeKutta"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "rk_stage", "codes/ns/src/NsCtrl.cpp", "NsCtrl::Init",
        'GetDataValue<int>("rk_stage")', ConsumerReadKind.REQUIRED,
        "Builds rk_coef as 1 / (stage + 1).",
        "Compressible Navier-Stokes control initialization.",
        "Direct read and coefficient construction are present in NsCtrl.cpp.",
        (), (), ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "rk_stage", "codes/ins/src/INsCtrl.cpp", "INsCtrl::Init",
        'GetDataValue<int>("rk_stage")', ConsumerReadKind.REQUIRED,
        "Builds rk_coef as 1 / (stage + 1).",
        "Incompressible Navier-Stokes control initialization.",
        "Direct read and coefficient construction are present in INsCtrl.cpp.",
        (), (), ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "rk_stage", "codes/turb/src/TurbCtrl.cpp", "TurbCtrl::Init",
        'GetDataValue<int>("rk_stage")', ConsumerReadKind.REQUIRED,
        "Builds rk_coef as 1 / (stage + 1).",
        "Turbulence control initialization.",
        "Direct read and coefficient construction are present in TurbCtrl.cpp.",
        (), (), ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "mach_ref", "codes/ns/src/NsCom.cpp", "NsCom::Init",
        'GetDataValue<Real>("mach_ref")', ConsumerReadKind.REQUIRED,
        "Loads the reference Mach number; NsCom::Init uses it in "
        "const_cp = 1 / ((gama_ref - 1) * mach_ref^2).",
        "Compressible Navier-Stokes common-state initialization.",
        "Direct typed read and derived const_cp expression are present in NsCom.cpp.",
    ),
    ParameterConsumer(
        "vismodel", "codes/global/src/Ctrl.cpp", "VisModel::Init",
        'GetDataValue<int>("vismodel")', ConsumerReadKind.REQUIRED,
        "Selects the viscous-model identifier. Other code branches on the resulting "
        "vis_model.vismodel, including NS model validation and wall-distance policy.",
        "Solver setup; downstream behavior depends on the selected model.",
        "Direct typed read is present in Ctrl.cpp; downstream examples are in "
        "NsCom.cpp and MultiBlock.cpp.",
    ),
    ParameterConsumer(
        "ptTol", "codes/global/src/Tolerance.cpp", "Tolerance::Init",
        'GetDataValue<Real>("ptTol")', ConsumerReadKind.LAZY,
        "Initializes the process-level geometric tolerance on the first Tolerance::GetTol call.",
        "Only code paths that request the shared geometric tolerance.",
        "Lazy initialization is explicit in Tolerance::GetTol; Reset clears its initialized flag.",
    ),
    ParameterConsumer(
        "nmg", "codes/multigrid/src/Multigrid.cpp", "MG::Init",
        'GetDataValue<int>("nmg")', ConsumerReadKind.REQUIRED,
        "Sets the multigrid level count. When nmg <= 1, MG::Init forces one pre-smoothing "
        "step and zero post-smoothing steps.",
        "Multigrid initialization.",
        "Direct typed read and the nmg <= 1 normalization are present in Multigrid.cpp.",
    ),
)


class ParameterConsumerCatalog:
    """Query a curated index of verified reads without guessing unknown consumers."""

    def __init__(self, consumers: tuple[ParameterConsumer, ...] = _CONSUMERS) -> None:
        self._consumers = consumers

    @property
    def consumers(self) -> tuple[ParameterConsumer, ...]:
        return self._consumers

    def find(self, parameter_name: str) -> tuple[ParameterConsumer, ...]:
        """Return verified reads for a key using case-insensitive matching."""
        key = parameter_name.casefold()
        return tuple(
            consumer for consumer in self._consumers
            if consumer.parameter_name.casefold() == key
        )

    def has_verified_consumer(self, parameter_name: str) -> bool:
        """Whether this curated index contains at least one verified read."""
        return bool(self.find(parameter_name))


def validate_parameter_consumer_sources(
    repository_root: Path,
    catalog: ParameterConsumerCatalog | None = None,
) -> list[str]:
    """Check that direct reads and named downstream behavior symbols still exist."""
    root = Path(repository_root).resolve()
    consumer_catalog = catalog or ParameterConsumerCatalog()
    issues: list[str] = []
    for consumer in consumer_catalog.consumers:
        source_path = root / consumer.source_path
        if not source_path.is_file():
            issues.append(f"{consumer.parameter_name}: missing consumer source: {consumer.source_path}")
            continue
        source = re.sub(r"\s+", "", source_path.read_text(encoding="utf-8", errors="replace"))
        if consumer.read_expression.replace(" ", "") not in source:
            issues.append(f"{consumer.parameter_name}: read expression not found in {consumer.source_path}")
        if consumer.symbol.replace(" ", "") not in source:
            issues.append(f"{consumer.parameter_name}: consumer symbol not found in {consumer.source_path}")
        behavior_sources = []
        for relative_path in consumer.behavior_source_paths:
            behavior_path = root / relative_path
            if not behavior_path.is_file():
                issues.append(f"{consumer.parameter_name}: missing behavior source: {relative_path}")
                continue
            behavior_sources.append(re.sub(r"\s+", "", behavior_path.read_text(encoding="utf-8", errors="replace")))
        combined_behavior = "\n".join(behavior_sources)
        for symbol in consumer.behavior_symbols:
            if symbol.replace(" ", "") not in combined_behavior:
                issues.append(f"{consumer.parameter_name}: behavior symbol not found: {symbol}")
    return issues
