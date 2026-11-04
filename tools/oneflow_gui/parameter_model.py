"""Logical CFD/grid parameter model composed from one or more legacy files.

Physical file boundaries are treated as storage details. The model groups typed
declarations by meaning while retaining each declaration's source location.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import re
from typing import Any, Iterable

from .parameter_document import ParameterDocument


class ParameterDomain(str, Enum):
    """Logical ownership of a parameter, independent of its physical file."""

    SHARED = "shared"
    CFD = "cfd"
    GRID = "grid"
    OTHER = "other"


class ParameterApplicability(str, Enum):
    """Task families in which a parameter is known to apply."""

    BOTH = "both"
    SOLVE = "solve"
    GRID = "grid"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class ParameterClassification:
    """Stable, conservative classification used by Workbench presentation."""

    domain: ParameterDomain
    group: str
    applicability: ParameterApplicability


@dataclass(frozen=True, slots=True)
class ParameterEntry:
    """A typed parameter value together with its original declaration source."""

    name: str
    declared_type: str
    raw_value: str
    value: Any
    source_path: Path
    file_order: int
    declaration_order: int
    classification: ParameterClassification

    @property
    def source_file(self) -> str:
        return self.source_path.name


@dataclass(frozen=True, slots=True)
class UnsupportedParameterDeclaration:
    """A declaration-like statement outside the supported scalar syntax."""

    name: str
    declared_type: str
    source_path: Path
    line_number: int
    reason: str
    raw_text: str

@dataclass(frozen=True, slots=True)
class ParameterGroup:
    """A display group containing effective values in declaration order."""

    domain: ParameterDomain
    name: str
    entries: tuple[ParameterEntry, ...]


_SHARED_KEYS = frozenset({"simutask", "dimension", "pttol"})

_GRID_GROUPS = {
    "Grid definition": frozenset({
        "gridobj", "topotype", "multiblock", "axis_dir", "igene",
        "codeofoversetgrid", "ignorenobc",
    }),
    "Grid formats and files": frozenset({
        "fileformat", "gridscale", "gridtrans",
        "sourcegridtype", "targetgridtype", "iplot3d",
        "sourcegridfilename", "sourcegridbcname", "targetgridfilename",
        "patchfilename", "gridlayoutfilename",
        "boundaryregionnamemapfilename", "volumeboundarynamemapfilename",
        "oversetgridfilename",
    }),
    "Partitioning": frozenset({
        "partition_type", "partition_c2n", "npartproc", "us_type",
        "ori_uns_file", "part_uns_file", "ori_layout_file",
        "ori_c2n_file", "part_c2n_file", "parallelstrategy",
    }),
}

_CFD_GROUPS = {
    "Run control": frozenset({
        "startstrategy", "iexitflag", "maxsteps", "maxtime",
        "numberofminimumsubiterationsteps",
        "numberofmaximumsubiterationsteps",
    }),
    "Solver system": frozenset({"nequ", "ntequ", "ntmodel"}),
    "Time stepping": frozenset({
        "timestepmodel", "global_dt", "idualtime", "ncfl", "cflst",
        "cfled", "max_time_ratio", "vissrmodel", "vistimestepmodel",
    }),
    "Output control": frozenset({
        "nfieldsave", "nforcesave", "nressave", "nvisualsave",
        "addvisualizationsteps", "showfield", "idump",
    }),
    "Output files": frozenset({
        "flowfile", "aerofile", "resfile", "wallaerofile", "heatfluxfile",
        "visualfile", "turbresfile", "turbflowfile",
    }),
    "Input and output files": frozenset({
        "flowinputroot", "flowoutputroot", "gridfilename",
    }),
    "Numerical methods": frozenset({
        "time_integral", "rk_stage", "lineartwostepmethods", "ischeme",
        "ivischeme", "ieigenfix", "centropy1", "centropy2", "ilim",
        "vencat_coef", "lusgsfacepressuremode",
        "turbulentresidualsmoothornot", "modifyturbulentresidualstrategy",
    }),
    "Flow and reference conditions": frozenset({
        "mach_ref", "aoa_degree", "sideslip_degree", "reynolds", "lref",
        "aref", "xref", "yref", "zref",
    }),
    "Thermodynamic reference state": frozenset({
        "gasinfostrategy", "machstrategy", "elevation", "dref_dim",
        "tref_dim", "pref_dim", "vref_dim", "reylref_dim", "gama_ref",
    }),
    "Transport properties": frozenset({"prl", "prt", "schmidtl", "schmidtt"}),
    "Boundary and inflow conditions": frozenset({
        "inflowtype", "manufacturedinflowtype", "inflowvelocity",
        "isowallbc", "twall_dim",
    }),
    "Physical models": frozenset({
        "vismodel", "visname", "chemmodel", "gasmodelfile",
        "transition_model", "des_model", "easm_model",
    }),
    "Turbulence": frozenset({
        "turbulentrelaxationcoefficient", "turbulentsourcemodel",
        "inflow_intensity", "inflow_viscosity",
        "mgturbulentviscosity", "mgfactor", "mgturbfactor",
        "coefficientofsades", "coefficientofsstdeske",
        "coefficientofsstdeskw", "numberofanisotropicstress",
        "coefficientofminimumturbulentquantity", "ref_sa", "ref_sst",
        "iturb_visflux", "iprod_sa", "iprod_sst", "max_vis_ratio",
        "rprod", "monitormaximumturbulentviscositycoefficient",
        "numberofabnormalturbulentcellslimit", "turb_cfl_ratio",
        "turb_ilim", "tns_ilim", "ft2_flag",
    }),
    "Multigrid": frozenset({
        "nmg", "nmgpresmoothingsteps", "nmgpostsmoothingsteps",
        "nmgsmoothingsteps", "mgprolongationmethod", "mgcycletype",
        "mgiterMode".casefold(), "mglaminarviscosity",
    }),
    "Wall distance": frozenset({"walldist_file", "ireadwdst"}),
}

_GRID_KEY_GROUP = {
    key: group for group, keys in _GRID_GROUPS.items() for key in keys
}
_CFD_KEY_GROUP = {
    key: group for group, keys in _CFD_GROUPS.items() for key in keys
}


_DECLARATION_HEAD = re.compile(
    r"^[ \t]*(?P<type>[A-Za-z_][A-Za-z0-9_]*)[ \t]+"
    r"(?P<name>[A-Za-z_][A-Za-z0-9_]*)(?P<shape>[ \t]*\[[^\]\r\n]*\])?[ \t]*=",
    re.MULTILINE,
)


def _normalize_key(name: str) -> str:
    return name.casefold()




def classify_parameter(name: str) -> ParameterClassification:
    """Classify by known parameter semantics, never by the active task alone.

    Task selection describes the current run; it is not evidence of parameter
    ownership. Unknown keys therefore remain explicitly unclassified.
    """
    key = _normalize_key(name)
    if key in _SHARED_KEYS:
        return ParameterClassification(
            ParameterDomain.SHARED, "General", ParameterApplicability.BOTH
        )
    if key in _GRID_KEY_GROUP:
        return ParameterClassification(
            ParameterDomain.GRID, _GRID_KEY_GROUP[key], ParameterApplicability.GRID
        )
    if key in _CFD_KEY_GROUP:
        return ParameterClassification(
            ParameterDomain.CFD, _CFD_KEY_GROUP[key], ParameterApplicability.SOLVE
        )
    return ParameterClassification(
        ParameterDomain.OTHER, "Unclassified", ParameterApplicability.UNKNOWN
    )


def _display_group_key(entry: ParameterEntry) -> tuple[ParameterDomain, str]:
    return entry.classification.domain, entry.classification.group


class ParameterModel:
    """Combine declarations from multiple files into one logical parameter view.

    The source files remain authoritative and are never rewritten by this model.
    Repeated names follow the runtime database's last-declaration-wins behavior;
    shadowed declarations remain available for inspection and diagnostics.
    """

    def __init__(
        self,
        task: str | None,
        declarations: tuple[ParameterEntry, ...],
        entries: tuple[ParameterEntry, ...],
        shadowed_entries: tuple[ParameterEntry, ...],
        unsupported_declarations: tuple[UnsupportedParameterDeclaration, ...] = (),
    ) -> None:
        self.task = task
        self.declarations = declarations
        self.entries = entries
        self.shadowed_entries = shadowed_entries
        self.unsupported_declarations = unsupported_declarations

    @classmethod
    def from_documents(
        cls,
        documents: Iterable[ParameterDocument],
    ) -> "ParameterModel":
        """Build one logical model from documents in control-manifest order."""
        document_list = tuple(documents)
        parsed: list[tuple[int, int, ParameterDocument, Any]] = []
        unsupported: list[UnsupportedParameterDeclaration] = []
        task = None
        for file_order, document in enumerate(document_list):
            parsed_lines = {
                document.text.count("\n", 0, declaration.value_start) + 1
                for declaration in document.declarations
            }
            for match in _DECLARATION_HEAD.finditer(document.text):
                line_number = document.text.count("\n", 0, match.start()) + 1
                shape = match.group("shape")
                if shape is not None:
                    reason = "Array declarations are not supported by the scalar parameter view."
                elif line_number not in parsed_lines:
                    reason = "Declaration type or syntax is not supported by the scalar parser."
                else:
                    continue
                line_start = document.text.rfind("\n", 0, match.start()) + 1
                line_end = document.text.find("\n", match.start())
                if line_end < 0:
                    line_end = len(document.text)
                unsupported.append(UnsupportedParameterDeclaration(
                    name=match.group("name"),
                    declared_type=match.group("type"),
                    source_path=document.path,
                    line_number=line_number,
                    reason=reason,
                    raw_text=document.text[line_start:line_end].strip(),
                ))
            for declaration_order, declaration in enumerate(document.declarations):
                parsed.append((file_order, declaration_order, document, declaration))
                if declaration.name.casefold() == "simutask":
                    task = str(declaration.value)

        declarations = tuple(
            ParameterEntry(
                name=declaration.name,
                declared_type=declaration.declared_type,
                raw_value=declaration.raw_value,
                value=declaration.value,
                source_path=document.path,
                file_order=file_order,
                declaration_order=declaration_order,
                classification=classify_parameter(declaration.name),
            )
            for file_order, declaration_order, document, declaration in parsed
        )

        effective_positions: dict[str, int] = {}
        effective: list[ParameterEntry] = []
        shadowed: list[ParameterEntry] = []
        for entry in declarations:
            key = entry.name.casefold()
            previous_position = effective_positions.get(key)
            if previous_position is None:
                effective_positions[key] = len(effective)
                effective.append(entry)
            else:
                shadowed.append(effective[previous_position])
                effective[previous_position] = entry

        return cls(task, declarations, tuple(effective), tuple(shadowed), tuple(unsupported))

    @property
    def groups(self) -> tuple[ParameterGroup, ...]:
        """Return grouped effective values while preserving first-seen group order."""
        grouped: dict[tuple[ParameterDomain, str], list[ParameterEntry]] = {}
        for entry in self.entries:
            grouped.setdefault(_display_group_key(entry), []).append(entry)
        return tuple(
            ParameterGroup(domain, name, tuple(entries))
            for (domain, name), entries in grouped.items()
        )

    def get(self, name: str) -> ParameterEntry | None:
        """Return the effective value using case-insensitive key matching."""
        key = name.casefold()
        return next((entry for entry in self.entries if entry.name.casefold() == key), None)
