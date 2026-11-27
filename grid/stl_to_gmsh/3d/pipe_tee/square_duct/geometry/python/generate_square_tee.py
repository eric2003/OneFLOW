#!/usr/bin/env python3
"""Generate an exact, manually triangulated square-duct T-junction STL."""
from __future__ import annotations
import argparse, hashlib, json, math, struct
from pathlib import Path

ROOT = Path(__file__).resolve().parents[7]
DEFAULTS = {"main_length": 4.0, "branch_length": 2.0, "duct_width": 1.0}

def cross(a,b,c):
    u=tuple(b[i]-a[i] for i in range(3)); v=tuple(c[i]-a[i] for i in range(3))
    return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])

def signed_volume(faces):
    return sum(sum(f[0][k]*cross(*f)[k] for k in range(3))/6.0 for f in faces)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def write_stl(path,faces,label):
    with path.open("wb") as stream:
        stream.write(("OneFLOW square T-junction "+label).encode("ascii")[:80].ljust(80,b"\0"))
        stream.write(struct.pack("<I",len(faces)))
        for face in faces:
            n=cross(*face); mag=math.sqrt(sum(v*v for v in n))
            n=tuple(v/mag for v in n) if mag else (0.0,0.0,0.0)
            stream.write(struct.pack("<12fH",*n,*face[0],*face[1],*face[2],0))

def generate(output: Path, operating_mode="main_inlets"):
    if operating_mode not in {"benchmark","main_inlets","branch_inlet"}:
        raise ValueError(f"Unsupported operating mode: {operating_mode}")
    output=output.expanduser().resolve()
    if output==ROOT or ROOT in output.parents:
        raise ValueError("Generated STL assets must stay outside the OneFLOW source tree.")
    output.parent.mkdir(parents=True,exist_ok=True)
    w=DEFAULTS["duct_width"]; h=w/2; x0=-2.0; x1=2.0; z1=h+2.0; y0=-h; y1=h
    # Exact T-shaped cross-section in X-Z; the branch starts at the main duct wall.
    polygon=[(x0,-h),(x1,-h),(x1,h),(h,h),(h,z1),(-h,z1),(-h,h),(x0,h)]
    groups={"wall":[],"main_left":[],"main_right":[],"branch":[]}
    def p(x,z,y): return (float(x),float(y),float(z))
    def quad(target,a,b,c,d): target.extend(((a,b,c),(a,c,d)))
    port_edges={1:"main_right",4:"branch",7:"main_left"}
    for i,start in enumerate(polygon):
        end=polygon[(i+1)%len(polygon)]
        target=groups[port_edges[i]] if i in port_edges else groups["wall"]
        quad(target,p(*start,y0),p(*start,y1),p(*end,y1),p(*end,y0))
    # Two adjoining rectangles tile the T-shaped y-normal walls without gaps.
    for y in (y0,y1):
        for xa,za,xb,zb in ((x0,-h,x1,h),(-h,h,h,z1)):
            a,b,c,d=p(xa,za,y),p(xb,za,y),p(xb,zb,y),p(xa,zb,y)
            if y==y0: quad(groups["wall"],a,b,c,d)       # outward normal -Y
            else: quad(groups["wall"],a,d,c,b)           # outward normal +Y
    if signed_volume(sum(groups.values(),[])) < 0:
        for name in groups: groups[name]=[(a,c,b) for a,b,c in groups[name]]
    ports={"main_left":groups["main_left"],"main_right":groups["main_right"],"branch":groups["branch"]}
    if operating_mode=="main_inlets":
        regions={"wall":groups["wall"],"inflow":ports["main_left"]+ports["main_right"],"outflow":ports["branch"]}
        named={"main_left_inflow":("inflow",ports["main_left"]),"main_right_inflow":("inflow",ports["main_right"]),"branch_outflow":("outflow",ports["branch"])}
    elif operating_mode=="branch_inlet":
        regions={"wall":groups["wall"],"inflow":ports["branch"],"outflow":ports["main_left"]+ports["main_right"]}
        named={"branch_inflow":("inflow",ports["branch"]),"main_left_outflow":("outflow",ports["main_left"]),"main_right_outflow":("outflow",ports["main_right"])}
    else:
        regions={"wall":groups["wall"],"inflow":ports["main_left"]+ports["branch"],"outflow":ports["main_right"]}
        named={"main_inflow":("inflow",ports["main_left"]),"branch_inflow":("inflow",ports["branch"]),"main_outflow":("outflow",ports["main_right"])}
    faces=regions["wall"]+regions["inflow"]+regions["outflow"]
    edges={}
    for f in faces:
        for a,b in ((f[0],f[1]),(f[1],f[2]),(f[2],f[0])):
            key=tuple(sorted((a,b))); edges[key]=edges.get(key,0)+1
    bad=sum(n!=2 for n in edges.values())
    volume=signed_volume(faces)
    if bad: raise ValueError(f"Square T-junction shell has {bad} non-manifold/open edges.")
    if volume<=0: raise ValueError(f"Square T-junction shell has non-positive signed volume: {volume}")
    write_stl(output,faces,"closed-fluid-domain")
    records={}
    for name,patch_faces in regions.items():
        path=output.parent/f"{output.stem}.{name}.stl"; write_stl(path,patch_faces,name)
        records[name]={"file":path.name,"sha256":sha(path),"triangle_count":len(patch_faces)}
    subs={}
    for name,(parent,patch_faces) in named.items():
        path=output.parent/f"{output.stem}.{name}.stl"; write_stl(path,patch_faces,name)
        subs[name]={"file":path.name,"sha256":sha(path),"triangle_count":len(patch_faces),"parent_region":parent}
    region_path=output.parent/"region_manifest.json"
    region_path.write_text(json.dumps({"schema":"oneflow.stl-surface-regions","version":1,"closed_surface_sha256":sha(output),"regions":records,"subregions":subs},indent=2,sort_keys=True)+"\n",encoding="utf-8")
    surface={"triangle_count":len(faces),"region_triangle_counts":{k:len(v) for k,v in regions.items()},"port_subregion_triangle_counts":{k:len(v) for k,(_p,v) in named.items()},"invalid_edge_count":bad,"signed_volume_m3":volume}
    geometry={"schema":"oneflow.generated-geometry","version":1,"geometry":"square_duct_pipe_tee","operating_mode":operating_mode,"units":"m","parameters":{"main_length":4.0,"branch_length":2.0,"cross_section":[w,w],"main_axis":"+X","branch_axis":"+Z","junction":"sharp exact planar intersection"},"surface":surface,"region_manifest":region_path.name,"stl_file":output.name,"stl_sha256":sha(output),"generator":{"kind":"python","method":"explicit_planar_facets_no_boolean_intersection"}}
    sidecar=output.with_suffix(".geometry.json"); sidecar.write_text(json.dumps(geometry,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return {"stl":output,"geometry_manifest":sidecar,"region_manifest":region_path,"surface":surface}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--operating-mode",choices=("benchmark","main_inlets","branch_inlet"),default="main_inlets")
    args=parser.parse_args()
    result=generate(args.output,args.operating_mode)
    print(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in result.items()},indent=2))
    return 0
if __name__=="__main__": raise SystemExit(main())
