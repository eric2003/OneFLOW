"""Evidence-backed index of selected legacy parameter runtime consumers.

This registry records source-level reads that have been checked in the current
repository. It is intentionally a small, curated first slice, not an automated
claim that every parameter consumer has been discovered.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


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
    ),
    ParameterConsumer(
        "time_integral", "codes/ns/src/NsCtrl.cpp", "NsCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the time-integration method in the compressible Navier-Stokes control.",
        "Compressible Navier-Stokes control initialization.",
        "Direct typed read is present in NsCtrl.cpp.",
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
