#!/usr/bin/env python3
"""Generate a 3-D tetrahedral volume mesh directly from a closed STL shell."""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from grid.stl_to_gmsh.converter import Stl2DError, convert, inspect_msh4  # noqa: E402


DEFAULT_SPHERE_STL = (
    ROOT / "grid" / "stl_to_gmsh" / "3d" / "sphere"
    / "automatic_tetrahedra" / "python" / "input" / "model.stl"
)
DEFAULT_MESH_SIZE = 0.20


def default_output_path() -> Path:
    """Return the platform-independent OneFLOW user-data output location."""
    configured_home = os.environ.get("ONEFLOW_HOME", "").strip()
    user_root = Path(configured_home).expanduser() if configured_home else Path.home() / ".oneflow"
    return user_root / "mesh-generation" / "sphere" / "sphere.msh"


def discover_gmsh() -> str | None:
    """Find Gmsh without requiring beginners to enter its executable path."""
    configured = os.environ.get("GMSH_EXECUTABLE", "").strip()
    if configured and Path(configured).is_file():
        return str(Path(configured).expanduser().resolve())

    on_path = shutil.which("gmsh") or shutil.which("gmsh.exe")
    if on_path:
        return on_path

    candidates: list[Path] = []
    if os.name == "nt":
        for variable in ("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA"):
            base = os.environ.get(variable)
            if base:
                candidates.extend(sorted(Path(base).glob("Gmsh*/gmsh.exe")))
        software = Path("C:/software")
        if software.is_dir():
            candidates.extend(sorted(software.glob("gmsh*/gmsh.exe")))
            candidates.extend(sorted(software.glob("gmsh*/gmsh.EXE")))
    else:
        candidates.extend((Path("/usr/local/bin/gmsh"), Path("/usr/bin/gmsh")))

    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    return None


def generate(
    stl: Path, output: Path, gmsh: str | None, mesh_size_min: float, mesh_size_max: float,
    boundary_name: str = "wall", volume_name: str = "fluid_domain",
    surface_angle_degrees: float = 40.0, curve_angle_degrees: float = 180.0,
) -> dict:
    from grid.stl_to_gmsh.run_report import build_run_report, format_run_summary, stage_result, write_run_report

    stl = stl.expanduser().resolve()
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    geo = output.with_suffix(".geo")
    report_path = output.with_suffix(".report.json")
    # A previous successful report must never describe this new attempt.
    report_path.unlink(missing_ok=True)
    stages = {
        "input": stage_result("NOT_RUN"), "geometry": stage_result("NOT_RUN"),
        "gmsh": stage_result("NOT_RUN"), "mesh_structure": stage_result("NOT_RUN"),
        "physical_groups": stage_result("NOT_RUN"),
        "mesh_quality": stage_result("NOT_CHECKED",
            "Tetrahedron quality acceptance criteria are not implemented yet."),
        "cfd_validation": stage_result("NOT_RUN",
            "Run Mesh-to-CFD Validation separately; it is outside this mesh-generation command."),
    }
    geometry_report = None
    mesh = None
    active_stage = "input"

    def make_report(exit_code: int, error: Exception | None = None) -> dict:
        metrics = {}
        if mesh:
            metrics = {"nodes":mesh["nodes"], "elements":mesh["elements"],
                       "surface_triangles":mesh["triangle_elements"],
                       "tetrahedra":mesh["tetrahedron_elements"], "msh_format":mesh["format"]}
            quality = mesh.get("triangle_quality", {})
            if "minimum_quality" in quality:
                metrics["minimum_surface_triangle_quality"] = quality["minimum_quality"]
        message = str(error) if error is not None else None
        if isinstance(error, subprocess.CalledProcessError):
            detail = error.stderr or error.stdout or ""
            detail_lines = detail.strip().splitlines()
            if detail_lines:
                message += "\n" + "\n".join(detail_lines[-12:])
        report = build_run_report(
            "STL -> Gmsh", stages, exit_code, metrics=metrics,
            artifacts={"geometry":str(geo),"mesh":str(output),"report":str(report_path)},
            error=({"type":type(error).__name__,"message":message} if error is not None else None),
        )
        report.update({
            "input_stl":str(stl), "geometry":str(geo), "mesh":str(output),
            "gmsh_executable":gmsh, "mesh_size_min":mesh_size_min,
            "mesh_size_max":mesh_size_max, "geometry_validation":geometry_report,
            "mesh_summary":mesh, "validated":exit_code == 0,
        })
        return report

    try:
        if not stl.is_file():
            raise FileNotFoundError(f"STL input does not exist: {stl}")
        if output.suffix.lower() != ".msh":
            raise ValueError("output path must use the .msh extension")
        if output == stl:
            raise ValueError("output mesh must not overwrite the input STL")
        if (not math.isfinite(mesh_size_min) or mesh_size_min <= 0
                or not math.isfinite(mesh_size_max) or mesh_size_max <= 0):
            raise ValueError("mesh sizes must be positive and finite")
        if mesh_size_min > mesh_size_max:
            raise ValueError("--mesh-size-min must not exceed --mesh-size-max")
        if not gmsh:
            raise ValueError("Gmsh was not found; install it, add it to PATH, or set GMSH_EXECUTABLE")
        stages["input"] = stage_result("PASS", "Input file and mesh-size parameters are valid.")

        active_stage = "geometry"
        manifest_path = stl.with_suffix(".geometry.json")
        geometry_kind = None
        if manifest_path.is_file():
            try:
                geometry_kind = json.loads(manifest_path.read_text(encoding="utf-8")).get("geometry")
            except (OSError, json.JSONDecodeError) as error:
                raise ValueError(f"Cannot read STL geometry manifest: {manifest_path}: {error}") from error
        region_manifest = stl.parent / "region_manifest.json"
        if region_manifest.is_file():
            volume_mesh_kind = "region_aware_closed_surface_automatic_tetrahedra"
        else:
            volume_mesh_kind = ("cylinder_internal_flow" if geometry_kind == "cylinder_internal_flow"
                                else "closed_surface_automatic_tetrahedra")
        config = {
            "volume_mesh":volume_mesh_kind, "mesh_dimension":3,
            "region_manifest":str(region_manifest) if region_manifest.is_file() else None,
            "mesh_size_min":mesh_size_min, "mesh_size_max":mesh_size_max,
            "surface_name":volume_name, "boundary_name":boundary_name,
            "surface_angle_degrees":surface_angle_degrees,
            "curve_angle_degrees":curve_angle_degrees,
        }
        with tempfile.TemporaryDirectory(prefix="oneflow-stl-to-gmsh-") as temporary:
            sidecar = Path(temporary) / "layout.json"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            geometry_report = convert(stl, sidecar, geo, tolerance=1.0e-10)
        stages["geometry"] = stage_result("PASS", "STL topology and geometry checks passed.",
            {key:geometry_report[key] for key in ("triangles","euler_characteristic","signed_volume")
             if key in geometry_report})

        active_stage = "gmsh"
        subprocess.run([gmsh,str(geo),"-3","-format","msh4","-o",str(output)],
                       check=True,cwd=output.parent,capture_output=True,text=True)
        stages["gmsh"] = stage_result("PASS", "Gmsh completed with exit code 0.")

        active_stage = "mesh_structure"
        mesh = inspect_msh4(output)
        if mesh["tetrahedron_elements"] <= 0:
            raise RuntimeError("Gmsh output contains no TETRA_4 volume elements")
        if mesh["triangle_elements"] <= 0:
            raise RuntimeError("Gmsh output contains no TRI_3 boundary elements")
        stages["mesh_structure"] = stage_result("PASS",
            "MSH 4.1 is readable; volume and boundary elements are present.",
            {"nodes":mesh["nodes"],"elements":mesh["elements"],
             "tetrahedra":mesh["tetrahedron_elements"],
             "surface_triangles":mesh["triangle_elements"]})

        active_stage = "physical_groups"
        groups = {(item["dimension"],item["name"]) for item in mesh["physical_groups"]}
        expected = ({(2,"wall"),(2,"inflow"),(2,"outflow"),(3,"fluid_domain")}
                    if volume_mesh_kind in {"cylinder_internal_flow", "region_aware_closed_surface_automatic_tetrahedra"}
                    else {(2,boundary_name),(3,volume_name)})
        if groups != expected:
            raise RuntimeError(f"Physical Groups {sorted(groups)} do not match {sorted(expected)}")
        group_counts = {
            (item["dimension"], item["name"]): item.get("element_count", 0)
            for item in mesh["physical_groups"]
        }
        empty_groups = sorted(
            f"{dimension}:{name}" for dimension, name in expected
            if group_counts.get((dimension, name), 0) <= 0
        )
        if empty_groups:
            raise RuntimeError(f"Physical Groups contain no mesh elements: {empty_groups}")
        stages["physical_groups"] = stage_result("PASS",
            "All expected Physical Groups contain mesh elements.",
            {"groups": {
                f"{dimension}:{name}": group_counts[(dimension, name)]
                for dimension, name in sorted(expected)
            }})

        report = make_report(0)
        write_run_report(report_path, report)
        print(format_run_summary(report))
        return report
    except Exception as error:
        stages[active_stage] = stage_result("FAIL", str(error))
        report = make_report(2, error)
        write_run_report(report_path, report)
        print(format_run_summary(report))
        raise
def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a tetrahedral volume mesh from a closed, genus-zero STL shell. "
            "Defaults run the bundled unit-sphere example without extra arguments."
        )
    )
    parser.add_argument(
        "--stl", type=Path, default=DEFAULT_SPHERE_STL,
        help=f"input closed STL surface (default: {DEFAULT_SPHERE_STL})",
    )
    parser.add_argument(
        "--output", type=Path, default=default_output_path(),
        help=f"output .msh file (default: {default_output_path()})",
    )
    parser.add_argument(
        "--mesh-size", type=float,
        help=f"set minimum and maximum size to one value (default: {DEFAULT_MESH_SIZE})",
    )
    parser.add_argument("--mesh-size-min", type=float, help="minimum target characteristic length")
    parser.add_argument("--mesh-size-max", type=float, help="maximum target characteristic length")
    parser.add_argument(
        "--gmsh", default=discover_gmsh(),
        help="Gmsh executable (auto-detected from GMSH_EXECUTABLE, PATH and common install folders)",
    )
    parser.add_argument("--boundary-name", default="wall")
    parser.add_argument("--volume-name", default="fluid_domain")
    parser.add_argument("--surface-angle", type=float, default=40.0)
    parser.add_argument("--curve-angle", type=float, default=180.0)
    args = parser.parse_args()

    if args.mesh_size_min is not None or args.mesh_size_max is not None:
        if args.mesh_size is not None:
            parser.error("use --mesh-size OR the --mesh-size-min/--mesh-size-max pair")
        if args.mesh_size_min is None or args.mesh_size_max is None:
            parser.error("provide both --mesh-size-min and --mesh-size-max")
        mesh_size_min, mesh_size_max = args.mesh_size_min, args.mesh_size_max
    elif args.mesh_size is not None:
        mesh_size_min = mesh_size_max = args.mesh_size
    else:
        mesh_size_min = mesh_size_max = DEFAULT_MESH_SIZE

    try:
        generate(
            args.stl, args.output, args.gmsh,
            mesh_size_min, mesh_size_max,
            args.boundary_name, args.volume_name,
            args.surface_angle, args.curve_angle,
        )
    except (
        OSError, ValueError, Stl2DError,
        subprocess.CalledProcessError, RuntimeError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
