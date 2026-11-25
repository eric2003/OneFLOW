// C++17 reference generator for the OneFLOW 90-degree circular pipe bend.
// The closed STL is a fluid-domain shell; inlet/outlet caps are not solid walls.
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double pi = 3.14159265358979323846;
using V = std::array<double,3>;
struct Ring { V center,tangent,normal; };
struct Face { std::array<std::size_t,3> v; std::string patch; };
V add(V a,V b){return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
V sub(V a,V b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
V mul(V a,double s){return {a[0]*s,a[1]*s,a[2]*s};}
V cross(V a,V b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double dot(V a,V b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm(V a){return std::sqrt(dot(a,a));}
struct Config { double radius=.5,bend_radius=1.,inlet=1.5,outlet=1.5,angle=90.; int circum=48,bend=24,inlet_div=12,outlet_div=12; std::string units="m"; };
void validate(const Config& c){
  if(c.units!="m"&&c.units!="mm"&&c.units!="cm"&&c.units!="in") throw std::invalid_argument("units must be one of m, mm, cm, in");
  if(!(std::isfinite(c.radius)&&std::isfinite(c.bend_radius)&&std::isfinite(c.inlet)&&std::isfinite(c.outlet)&&std::isfinite(c.angle))) throw std::invalid_argument("parameters must be finite");
  if(c.radius<=0||c.bend_radius<=c.radius||c.inlet<=0||c.outlet<=0||c.angle<=0||c.angle>180) throw std::invalid_argument("require radius>0, bend_radius>radius, positive lengths and angle in (0,180]");
  if(c.circum<8||c.bend<2||c.inlet_div<1||c.outlet_div<1) throw std::invalid_argument("insufficient surface divisions");
}
std::pair<std::vector<V>,std::vector<Face>> build(const Config& c){
  validate(c); const double theta=c.angle*pi/180.; std::vector<Ring> rings;
  for(int i=0;i<=c.inlet_div;++i){double s=-c.inlet+c.inlet*i/double(c.inlet_div);rings.push_back({{s,0,0},{1,0,0},{0,1,0}});}
  for(int i=1;i<=c.bend;++i){double p=theta*i/double(c.bend);rings.push_back({{c.bend_radius*std::sin(p),c.bend_radius*(1-std::cos(p)),0},{std::cos(p),std::sin(p),0},{-std::sin(p),std::cos(p),0}});}
  V end{c.bend_radius*std::sin(theta),c.bend_radius*(1-std::cos(theta)),0},t{std::cos(theta),std::sin(theta),0},n{-std::sin(theta),std::cos(theta),0};
  for(int i=1;i<=c.outlet_div;++i) rings.push_back({add(end,mul(t,c.outlet*i/double(c.outlet_div))),t,n});
  std::vector<V> points; for(const auto& r:rings) for(int j=0;j<c.circum;++j){double a=2*pi*j/c.circum;points.push_back(add(r.center,add(mul(r.normal,c.radius*std::cos(a)),V{0,0,c.radius*std::sin(a)})));}
  std::vector<Face> faces;
  for(std::size_t i=0;i+1<rings.size();++i) for(int j=0;j<c.circum;++j){std::size_t a=i*c.circum+j,b=i*c.circum+(j+1)%c.circum,d=(i+1)*c.circum+j,e=(i+1)*c.circum+(j+1)%c.circum;faces.push_back({{a,b,e},"wall"});faces.push_back({{a,e,d},"wall"});}
  const std::size_t inlet=points.size();points.push_back(rings.front().center);const std::size_t outlet=points.size();points.push_back(rings.back().center);const std::size_t last=(rings.size()-1)*c.circum;
  for(int j=0;j<c.circum;++j){std::size_t k=(j+1)%c.circum;faces.push_back({{inlet,k,std::size_t(j)},"inflow"});faces.push_back({{outlet,last+std::size_t(j),last+k},"outflow"});}
  double volume=0;for(const auto& f:faces){V a=points[f.v[0]],b=points[f.v[1]],d=points[f.v[2]];volume+=dot(a,cross(b,d))/6.;}
  if(volume<0)for(auto& f:faces)std::swap(f.v[1],f.v[2]);
  return {points,faces};
}
void write(const std::filesystem::path& path,const Config& c,const std::vector<V>& p,const std::vector<Face>& f){
  std::filesystem::create_directories(path.parent_path());std::ofstream stl(path);if(!stl)throw std::runtime_error("cannot open STL output");
  stl<<"solid oneflow_pipe_bend_90_degree\n"<<std::setprecision(9);
  for(const auto& q:f){V a=p[q.v[0]],b=p[q.v[1]],d=p[q.v[2]],n=cross(sub(b,a),sub(d,a));double m=norm(n);if(m)n=mul(n,1/m);stl<<" facet normal "<<n[0]<<' '<<n[1]<<' '<<n[2]<<"\n  outer loop\n";for(V x:{a,b,d})stl<<"   vertex "<<x[0]<<' '<<x[1]<<' '<<x[2]<<"\n";stl<<"  endloop\n endfacet\n";}
  stl<<"endsolid oneflow_pipe_bend_90_degree\n";stl.close();
  std::size_t wall=0,inflow=0,outflow=0;for(const auto& x:f){if(x.patch=="wall")++wall;else if(x.patch=="inflow")++inflow;else ++outflow;}
  auto sidecar=path;sidecar.replace_extension(".geometry-report.json");std::ofstream meta(sidecar);if(!meta)throw std::runtime_error("cannot open geometry report");
  meta<<std::setprecision(17)<<"{\n  \"schema\": \"oneflow.pipe-bend-geometry-report\",\n  \"version\": 1,\n  \"asset_id\": \"geometry.pipe_bend.90_degree\",\n  \"geometry\": \"pipe_bend_90_degree\",\n  \"units\": \""<<c.units<<"\",\n  \"parameters\": {\"radius\": "<<c.radius<<", \"bend_radius\": "<<c.bend_radius<<", \"inlet_length\": "<<c.inlet<<", \"outlet_length\": "<<c.outlet<<", \"angle_degrees\": "<<c.angle<<", \"circumferential_divisions\": "<<c.circum<<", \"bend_divisions\": "<<c.bend<<", \"inlet_divisions\": "<<c.inlet_div<<", \"outlet_divisions\": "<<c.outlet_div<<"},\n  \"surface\": {\"vertex_count\": "<<p.size()<<", \"triangle_count\": "<<f.size()<<"},\n  \"patch_triangle_counts\": {\"wall\": "<<wall<<", \"inflow\": "<<inflow<<", \"outflow\": "<<outflow<<"},\n  \"intended_gmsh_physical_groups\": {\"surface\": [\"wall\", \"inflow\", \"outflow\"], \"volume\": [\"fluid_domain\"]},\n  \"notes\": \"Inlet/outlet caps close the fluid-domain shell; they are boundary patches, not solid walls.\"\n}\n";
}
}
std::filesystem::path find_repository_root(){
  auto current=std::filesystem::weakly_canonical(std::filesystem::current_path());
  for(auto candidate=current;!candidate.empty();){
    if(std::filesystem::exists(candidate/".git")&&std::filesystem::exists(candidate/"grid"))return candidate;
    auto parent=candidate.parent_path();
    if(parent==candidate)break;
    candidate=parent;
  }
  return current;
}
int main(int argc,char** argv){
 try{
  Config c;const char* home=std::getenv("USERPROFILE");if(!home)home=std::getenv("HOME");
  std::filesystem::path output;
  if(home) output=std::filesystem::path(home)/".oneflow/mesh-generation/generated_stl/pipe_bend_90_degree/pipe_bend_90_degree.stl";
  else output="pipe_bend_90_degree.stl";
  for(int i=1;i<argc;++i){std::string a=argv[i];auto val=[&](){if(i+1>=argc)throw std::invalid_argument("missing value for "+a);return std::string(argv[++i]);};
   if(a=="--output")output=val();else if(a=="--radius")c.radius=std::stod(val());else if(a=="--bend-radius")c.bend_radius=std::stod(val());else if(a=="--inlet-length")c.inlet=std::stod(val());else if(a=="--outlet-length")c.outlet=std::stod(val());else if(a=="--angle-degrees")c.angle=std::stod(val());else if(a=="--circumferential-divisions")c.circum=std::stoi(val());else if(a=="--bend-divisions")c.bend=std::stoi(val());else if(a=="--inlet-divisions")c.inlet_div=std::stoi(val());else if(a=="--outlet-divisions")c.outlet_div=std::stoi(val());else if(a=="--units")c.units=val();else throw std::invalid_argument("unknown option: "+a);
  }
  output=std::filesystem::absolute(output);
  auto root=find_repository_root();auto target=std::filesystem::weakly_canonical(output.parent_path());
  auto rel=target.lexically_relative(root);if(!rel.empty()&&*rel.begin()!=std::filesystem::path(".."))throw std::invalid_argument("output must be outside the current source working tree");
  auto [points,faces]=build(c);write(output,c,points,faces);std::cout<<"PASS: "<<output<<" ("<<points.size()<<" vertices, "<<faces.size()<<" triangles)\n";
 }catch(const std::exception& e){std::cerr<<"pipe-bend generation failed: "<<e.what()<<'\n';return 2;}return 0;
}
