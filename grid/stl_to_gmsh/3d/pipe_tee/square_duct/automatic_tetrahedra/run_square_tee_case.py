#!/usr/bin/env python3
"""Generate a square T-junction STL and automatic TETRA_4 volume mesh."""
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[6]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from grid.stl_to_gmsh.generate_volume_mesh import generate as generate_volume_mesh
GENERATOR=ROOT/"grid/stl_to_gmsh/3d/pipe_tee/square_duct/geometry/python/generate_square_tee.py"
def run_case(case_directory:Path,gmsh:str,mesh_size:float=0.75,operating_mode:str="main_inlets")->Path:
 case_directory=case_directory.expanduser().resolve(); inp=case_directory/"input"; out=case_directory/"results"
 inp.mkdir(parents=True,exist_ok=True); out.mkdir(parents=True,exist_ok=True)
 stl=inp/"square_tee.stl"; mesh=out/"square_tee.msh"
 subprocess.run([sys.executable,str(GENERATOR),"--output",str(stl),"--operating-mode",operating_mode],check=True,cwd=ROOT)
 report=generate_volume_mesh(stl,mesh,gmsh,mesh_size,mesh_size)
 if report.get("validated") is not True: raise RuntimeError("Square T-junction STL-derived volume meshing failed.")
 return mesh
def main():
 parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--case-directory",type=Path,required=True)
 parser.add_argument("--gmsh",required=True); parser.add_argument("--mesh-size",type=float,default=0.75)
 parser.add_argument("--operating-mode",choices=("benchmark","main_inlets","branch_inlet"),default="main_inlets")
 args=parser.parse_args()
 try: mesh=run_case(args.case_directory,args.gmsh,args.mesh_size,args.operating_mode)
 except (OSError,ValueError,RuntimeError,subprocess.CalledProcessError) as error: parser.exit(2,f"Square T-junction Case failed: {error}\n")
 print(f"PASS: generated {mesh}"); return 0
if __name__=="__main__": raise SystemExit(main())
