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
        "maxTime", "codes/global/src/Ctrl.cpp", "Ctrl::Init",
        'GetDataValue<Real>("maxTime")', ConsumerReadKind.REQUIRED,
        "Loads the simulation-time stop threshold into Ctrl::maxTime and initializes "
        "currTime to zero. SimuIterState::Running compares currTime >= maxTime when "
        "iexitflag != 0; otherwise the outer loop uses maxSteps.",
        "Solver run initialization; time-based outer-loop termination mode.",
        "Direct typed read and time-threshold predicate are present in Ctrl.cpp and Iteration.cpp.",
        ("codes/global/src/Iteration.cpp",),
        ("SimuIterState::Running",),
    ),
    ParameterConsumer(
        "time_integral", "codes/global/src/Ctrl.cpp", "Ctrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the selected method in Ctrl::time_integral. TimeIntegral::Init dispatches MULTI_STAGE to RungeKutta and SIMPLE to Simple; all other values fall back to Lusgs. The context-aware Relaxation overload applies the same dispatch.",
        "Global solver control initialization.",
        "Direct typed read is present in Ctrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "time_integral", "codes/ns/src/NsCtrl.cpp", "NsCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the selected method in Ctrl::time_integral. TimeIntegral::Init dispatches MULTI_STAGE to RungeKutta and SIMPLE to Simple; all other values fall back to Lusgs. The context-aware Relaxation overload applies the same dispatch.",
        "Compressible Navier-Stokes control initialization.",
        "Direct typed read is present in NsCtrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "time_integral", "codes/ins/src/INsCtrl.cpp", "INsCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the selected method in Ctrl::time_integral. TimeIntegral::Init dispatches MULTI_STAGE to RungeKutta and SIMPLE to Simple; all other values fall back to Lusgs. The context-aware Relaxation overload applies the same dispatch.",
        "Incompressible Navier-Stokes control initialization.",
        "Direct typed read is present in INsCtrl.cpp.",
        ("codes/solver/src/TimeIntegral.cpp", "codes/solver/src/EulerRungeKuttaCapability.cpp", "codes/solver/include/EulerRungeKuttaCapability.h"),
        ("TimeIntegral::Init", "TimeIntegral::Relaxation", "TimeIntegral::RungeKutta", "EvaluateEulerRungeKuttaCapability"),
        ("time_integral_rk_stage_weight",),
    ),
    ParameterConsumer(
        "time_integral", "codes/turb/src/TurbCtrl.cpp", "TurbCtrl::Init",
        'GetDataValue<int>("time_integral")', ConsumerReadKind.REQUIRED,
        "Stores the selected method in Ctrl::time_integral. TimeIntegral::Init dispatches MULTI_STAGE to RungeKutta and SIMPLE to Simple; all other values fall back to Lusgs. The context-aware Relaxation overload applies the same dispatch.",
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

    ParameterConsumer(
        "idualtime", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'FindDataEntry("idualtime")', ConsumerReadKind.OPTIONAL,
        "Optionally loads the dual-time selector into Iteration::dualtime. Missing entries retain the default value 0. TimeStep::CalcCfl uses the outer-step index for physical-time mode and the inner-step index when dualtime == 1.",
        "CFL ramp iteration indexing; optional compatibility input.",
        "The guarded FindDataEntry read and default-zero behavior are present in Iteration.cpp.",
        ("codes/ns/src/TimeStep.cpp", "codes/ns/src/TimeStepPolicy.cpp", "codes/ns/include/TimeStepPolicy.h"),
        ("TimeStep::CalcCfl", "CflIterationIndex"),
    ),
    ParameterConsumer(
        "cflst", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<Real>("cflst")', ConsumerReadKind.REQUIRED,
        "Loads the starting value of the CFL ramp. TimeStep::CalcCfl selects the ramp index and CflAtIteration evaluates the schedule.",
        "Compressible NS time-step calculation; ramp index depends on physical-time or dual-time mode.",
        "Direct typed read is present in Iteration.cpp; the schedule is implemented by the pure TimeStepPolicy functions.",
        ("codes/ns/src/TimeStep.cpp", "codes/ns/src/TimeStepPolicy.cpp", "codes/ns/include/TimeStepPolicy.h"),
        ("TimeStep::CalcCfl", "CflIterationIndex", "CflAtIteration"),
    ),
    ParameterConsumer(
        "cfled", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<Real>("cfled")', ConsumerReadKind.REQUIRED,
        "Loads the ending value of the CFL ramp. CflAtIteration returns this value once the configured ramp length is reached.",
        "Compressible NS time-step calculation; ramp index depends on physical-time or dual-time mode.",
        "Direct typed read is present in Iteration.cpp; the ramp endpoint behavior is implemented in TimeStepPolicy.cpp.",
        ("codes/ns/src/TimeStep.cpp", "codes/ns/src/TimeStepPolicy.cpp", "codes/ns/include/TimeStepPolicy.h"),
        ("TimeStep::CalcCfl", "CflIterationIndex", "CflAtIteration"),
    ),
    ParameterConsumer(
        "ncfl", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<int>("ncfl")', ConsumerReadKind.REQUIRED,
        "Loads the CFL ramp length in iterations. CflAtIteration interpolates while iteration < ncfl and then returns cfled.",
        "Compressible NS time-step calculation; physical-time uses outerSteps and dual-time uses innerSteps.",
        "Direct typed read is present in Iteration.cpp; indexing and interpolation are implemented in TimeStepPolicy.cpp.",
        ("codes/ns/src/TimeStep.cpp", "codes/ns/src/TimeStepPolicy.cpp", "codes/ns/include/TimeStepPolicy.h"),
        ("TimeStep::CalcCfl", "CflIterationIndex", "CflAtIteration"),
    ),
    ParameterConsumer(
        "timestepModel", "codes/ns/src/NsCom.cpp", "NsCom::Init",
        'GetDataValue<int>("timestepModel")', ConsumerReadKind.REQUIRED,
        "Selects the unstructured solver time-step strategy: 0 local, 1 global, and the remaining branch uses the local/global hybrid path (the case-file comment documents 2 as lg).",
        "Compressible unstructured Navier-Stokes time-step calculation.",
        "The setting is read in NsCom.cpp and dispatched in UTimeStep::CalcTimeStep.",
        ("codes/uns/src/UTimeStep.cpp",),
        ("UTimeStep::CalcTimeStep", "UTimeStep::CalcLocalTimeStep", "UTimeStep::CalcGlobalTimeStep", "UTimeStep::CalcLgTimeStep"),
    ),
    ParameterConsumer(
        "global_dt", "codes/global/src/Ctrl.cpp", "Ctrl::Init",
        'GetDataValue<Real>("global_dt")', ConsumerReadKind.REQUIRED,
        "Loads the configured global time step into Ctrl::pdt. UTimeStep::CalcGlobalTimeStep applies pdt; the local/global hybrid path replaces pdt with the minimum local step.",
        "Solver control initialization; used by the global strategy and as mutable state in the local/global hybrid strategy.",
        "Direct typed read is present in Ctrl.cpp; strategy-specific use is present in UTimeStep.cpp.",
        ("codes/uns/src/UTimeStep.cpp",),
        ("UTimeStep::CalcGlobalTimeStep", "UTimeStep::CalcLgTimeStep"),
    ),
    ParameterConsumer(
        "global_dt", "codes/ns/src/NsCtrl.cpp", "NsCtrl::Init",
        'GetDataValue<Real>("global_dt")', ConsumerReadKind.REQUIRED,
        "Loads the configured global time step into Ctrl::pdt. UTimeStep::CalcGlobalTimeStep applies pdt; the local/global hybrid path replaces pdt with the minimum local step.",
        "Compressible Navier-Stokes control initialization.",
        "Direct typed read is present in NsCtrl.cpp; strategy-specific use is present in UTimeStep.cpp.",
        ("codes/uns/src/UTimeStep.cpp",),
        ("UTimeStep::CalcGlobalTimeStep", "UTimeStep::CalcLgTimeStep"),
    ),
    ParameterConsumer(
        "global_dt", "codes/ins/src/INsCtrl.cpp", "INsCtrl::Init",
        'GetDataValue<Real>("global_dt")', ConsumerReadKind.REQUIRED,
        "Loads the configured global time step into Ctrl::pdt. UTimeStep::CalcGlobalTimeStep applies pdt; the local/global hybrid path replaces pdt with the minimum local step.",
        "Incompressible Navier-Stokes control initialization.",
        "Direct typed read is present in INsCtrl.cpp; strategy-specific use is present in UTimeStep.cpp.",
        ("codes/uns/src/UTimeStep.cpp",),
        ("UTimeStep::CalcGlobalTimeStep", "UTimeStep::CalcLgTimeStep"),
    ),
    ParameterConsumer(
        "global_dt", "codes/turb/src/TurbCtrl.cpp", "TurbCtrl::Init",
        'GetDataValue<Real>("global_dt")', ConsumerReadKind.REQUIRED,
        "Loads the configured global time step into Ctrl::pdt. UTimeStep::CalcGlobalTimeStep applies pdt; the local/global hybrid path replaces pdt with the minimum local step.",
        "Turbulence control initialization.",
        "Direct typed read is present in TurbCtrl.cpp; strategy-specific use is present in UTimeStep.cpp.",
        ("codes/uns/src/UTimeStep.cpp",),
        ("UTimeStep::CalcGlobalTimeStep", "UTimeStep::CalcLgTimeStep"),
    ),


    ParameterConsumer(
        "inflow_intensity", "codes/turb/src/TurbCom.cpp", "TurbCom::Init",
        'GetDataValue<Real>("inflow_intensity")', ConsumerReadKind.REQUIRED,
        "Loads the turbulence inflow intensity. In the two-equation inflow initializer, a positive value sets k_e = 1.5*(0.01*intensity)^2; non-positive values use the empirical intensity fallback.",
        "Turbulence inflow initialization for the two-equation model.",
        "The value is read during TurbCom initialization and consumed by TurbCom::InitInflow.",
        ("codes/turb/src/TurbCom.cpp",),
        ("TurbCom::InitInflow",),
    ),
    ParameterConsumer(
        "inflow_viscosity", "codes/turb/src/TurbCom.cpp", "TurbCom::Init",
        'GetDataValue<Real>("inflow_viscosity")', ConsumerReadKind.REQUIRED,
        "Loads the turbulence inflow eddy viscosity. The two-equation inflow initializer uses it as mutoo to derive the reference k and omega values together with Reynolds number and reference density.",
        "Turbulence inflow initialization for the two-equation model.",
        "The value is read during TurbCom initialization and consumed by TurbCom::InitInflow.",
        ("codes/turb/src/TurbCom.cpp",),
        ("TurbCom::InitInflow",),
    ),
    ParameterConsumer(
        "turb_cfl_ratio", "codes/turb/src/TurbCom.cpp", "TurbCom::Init",
        'GetDataValue<Real>("turb_cfl_ratio")', ConsumerReadKind.REQUIRED,
        "Scales the turbulence-equation local time-step used in the unsteady spectral contribution. UTurbSpectrum::CalcUnsteadySpectrum multiplies the selected flow/turbulence time-step by turb_cfl_ratio before forming the implicit unsteady term; the dual-time branch also includes the physical-time contribution through ctrl.pdt1.",
        "Unstructured turbulence spectrum assembly; single-step and dual-time branches.",
        "Direct typed read is present in TurbCom.cpp; both single-step and dual-time uses are present in UTurbSpectrum.cpp.",
        ("codes/uturb/src/UTurbSpectrum.cpp",),
        ("UTurbSpectrum::CalcUnsteadySpectrum",),
    ),
    ParameterConsumer(
        "max_vis_ratio", "codes/turb/src/TurbCom.cpp", "TurbCom::Init",
        'GetDataValue<Real>("max_vis_ratio")', ConsumerReadKind.REQUIRED,
        "Loads max_vis_ratio into the turbulence common-state object. The downstream numerical use and intended physical meaning have not yet been established by this audit.",
        "Turbulence initialization; downstream effect not yet audited.",
        "The typed read is present in TurbCom::Init; no downstream consumer is claimed here.",
    ),
    ParameterConsumer(
        "rprod", "codes/turb/src/TurbCom.cpp", "TurbCom::Init",
        'GetDataValue<Real>("rprod")', ConsumerReadKind.REQUIRED,
        "Supplies the production multiplier in the Spalart-Allmaras limited-diffusion policy: SaLimitedDiffusion returns min(rprod * production, diffusion). The audited call is in the SA source-term path.",
        "Spalart-Allmaras turbulence source-term assembly.",
        "Direct typed read, call site, and policy implementation are present in TurbCom.cpp and TurbSaPolicy.cpp.",
        ("codes/turb/src/TurbCom.cpp", "codes/turb/src/TurbSaPolicy.cpp"),
        ("SaLimitedDiffusion",),
    ),
    ParameterConsumer(
        "nFieldSave", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<int>("nFieldSave")', ConsumerReadKind.REQUIRED,
        "Loads the field/restart output interval. Representative NS, incompressible NS, turbulence, and multigrid paths schedule restart output using outer-step or solver-specific iteration checks.",
        "Solver output scheduling; exact behavior depends on solver path and start strategy.",
        "The value is loaded in Iteration.cpp and used by solver implementation task scheduling.",
        ("codes/ns/src/NsSolverImp.cpp", "codes/ins/src/INsSolverImp.cpp", "codes/turb/src/TurbSolverImp.cpp", "codes/multigrid/src/Multigrid.cpp"),
        ("kDumpRestartTaskName",),
    ),
    ParameterConsumer(
        "nVisualSave", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<int>("nVisualSave")', ConsumerReadKind.REQUIRED,
        "Loads the visualization output interval. Representative compressible and incompressible NS paths schedule visualization tasks using this interval.",
        "Solver output scheduling; exact behavior depends on solver path and start strategy.",
        "The value is loaded in Iteration.cpp and used by NS solver implementation task scheduling.",
        ("codes/ns/src/NsSolverImp.cpp", "codes/ins/src/INsSolverImp.cpp"),
        ("kVisualizationTaskName", "kDumpPressureCoeffTaskName"),
    ),
    ParameterConsumer(
        "nResSave", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<int>("nResSave")', ConsumerReadKind.REQUIRED,
        "Loads the residual output interval. Iteration::ResOk uses innerSteps when startStrategy == 2 and outerSteps otherwise.",
        "Global iteration output scheduling.",
        "Direct read and interval predicate are present in Iteration.cpp.",
        ("codes/global/src/Iteration.cpp",),
        ("Iteration::ResOk",),
    ),
    ParameterConsumer(
        "nForceSave", "codes/global/src/Iteration.cpp", "Iteration::Init",
        'GetDataValue<int>("nForceSave")', ConsumerReadKind.REQUIRED,
        "Loads the force output interval. Iteration::ForceOk uses innerSteps when startStrategy == 2 and outerSteps otherwise.",
        "Global iteration output scheduling.",
        "Direct read and interval predicate are present in Iteration.cpp.",
        ("codes/global/src/Iteration.cpp",),
        ("Iteration::ForceOk",),
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
