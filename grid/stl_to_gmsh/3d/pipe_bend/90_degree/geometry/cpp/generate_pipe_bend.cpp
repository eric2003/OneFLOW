// C++17 reference generator for the OneFLOW 90-degree circular pipe bend.
// The closed STL is a fluid-domain shell; inlet/outlet caps are not solid walls.
#include <array>
#include <cstdint>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <iterator>
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

std::string sha256(const std::string& input) {
  static constexpr std::uint32_t k[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
  };
  auto rotr=[](std::uint32_t x,int n){return (x>>n)|(x<<(32-n));};
  std::vector<unsigned char> data(input.begin(),input.end());
  const std::uint64_t bits=static_cast<std::uint64_t>(data.size())*8;
  data.push_back(0x80); while(data.size()%64!=56)data.push_back(0);
  for(int i=7;i>=0;--i)data.push_back(static_cast<unsigned char>(bits>>(i*8)));
  std::uint32_t h[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  for(std::size_t offset=0;offset<data.size();offset+=64){
    std::uint32_t w[64]{};
    for(int i=0;i<16;++i){std::size_t p=offset+4*i;w[i]=(std::uint32_t(data[p])<<24)|(std::uint32_t(data[p+1])<<16)|(std::uint32_t(data[p+2])<<8)|data[p+3];}
    for(int i=16;i<64;++i){auto s0=rotr(w[i-15],7)^rotr(w[i-15],18)^(w[i-15]>>3);auto s1=rotr(w[i-2],17)^rotr(w[i-2],19)^(w[i-2]>>10);w[i]=w[i-16]+s0+w[i-7]+s1;}
    auto a=h[0],b=h[1],cc=h[2],d=h[3],e=h[4],ff=h[5],g=h[6],hh=h[7];
    for(int i=0;i<64;++i){auto s1=rotr(e,6)^rotr(e,11)^rotr(e,25);auto ch=(e&ff)^(~e&g);auto t1=hh+s1+ch+k[i]+w[i];auto s0=rotr(a,2)^rotr(a,13)^rotr(a,22);auto maj=(a&b)^(a&cc)^(b&cc);auto t2=s0+maj;hh=g;g=ff;ff=e;e=d+t1;d=cc;cc=b;b=a;a=t1+t2;}
    h[0]+=a;h[1]+=b;h[2]+=cc;h[3]+=d;h[4]+=e;h[5]+=ff;h[6]+=g;h[7]+=hh;
  }
  std::ostringstream out;out<<std::hex<<std::setfill('0');for(auto value:h)out<<std::setw(8)<<value;return out.str();
}
std::string read_bytes(const std::filesystem::path& path) {
  std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("cannot read generated asset");
  return std::string(std::istreambuf_iterator<char>(in),{});
}
void write_patch(const std::filesystem::path& path,const std::vector<V>& p,
                 const std::vector<Face>& f,const std::string& patch) {
  std::vector<Face> selected;for(const auto& face:f)if(face.patch==patch)selected.push_back(face);
  std::ofstream out(path);if(!out)throw std::runtime_error("cannot open region STL");
  out<<"solid oneflow_pipe_bend_"<<patch<<"\n"<<std::setprecision(9);
  for(const auto& face:selected){V a=p[face.v[0]],b=p[face.v[1]],d=p[face.v[2]],n=cross(sub(b,a),sub(d,a));double m=norm(n);if(m)n=mul(n,1/m);
    out<<" facet normal "<<n[0]<<' '<<n[1]<<' '<<n[2]<<"\n  outer loop\n";
    for(V x:{a,b,d})out<<"   vertex "<<x[0]<<' '<<x[1]<<' '<<x[2]<<"\n";
    out<<"  endloop\n endfacet\n";
  }
  out<<"endsolid oneflow_pipe_bend_"<<patch<<"\n";
}

void write(const std::filesystem::path& path,const Config& c,const std::vector<V>& p,const std::vector<Face>& f){
  std::filesystem::create_directories(path.parent_path());std::ofstream stl(path);if(!stl)throw std::runtime_error("cannot open STL output");
  stl<<"solid oneflow_pipe_bend_90_degree\n"<<std::setprecision(9);
  for(const auto& q:f){V a=p[q.v[0]],b=p[q.v[1]],d=p[q.v[2]],n=cross(sub(b,a),sub(d,a));double m=norm(n);if(m)n=mul(n,1/m);stl<<" facet normal "<<n[0]<<' '<<n[1]<<' '<<n[2]<<"\n  outer loop\n";for(V x:{a,b,d})stl<<"   vertex "<<x[0]<<' '<<x[1]<<' '<<x[2]<<"\n";stl<<"  endloop\n endfacet\n";}
  stl<<"endsolid oneflow_pipe_bend_90_degree\n";stl.close();
  std::size_t wall=0,inflow=0,outflow=0;for(const auto& x:f){if(x.patch=="wall")++wall;else if(x.patch=="inflow")++inflow;else ++outflow;}
  std::map<std::string,std::filesystem::path> region_files;
  std::map<std::string,std::size_t> region_counts{{"wall",wall},{"inflow",inflow},{"outflow",outflow}};
  for(const std::string name:{"wall","inflow","outflow"}){
    auto region=path;region.replace_filename(path.stem().string()+"."+name+".stl");
    write_patch(region,p,f,name);region_files[name]=region;
  }
  auto region_manifest=path.parent_path()/"region_manifest.json";
  std::ofstream regions(region_manifest);if(!regions)throw std::runtime_error("cannot open region manifest");
  regions<<"{\n  \"schema\": \"oneflow.stl-surface-regions\",\n  \"version\": 1,\n  \"closed_surface_sha256\": \""<<sha256(read_bytes(path))<<"\",\n  \"regions\": {\n";
  int index=0;for(const std::string name:{"wall","inflow","outflow"}){
    if(index++)regions<<",\n";
    regions<<"    \""<<name<<"\": { \"file\": \""<<region_files[name].filename().string()
      <<"\", \"sha256\": \""<<sha256(read_bytes(region_files[name]))
      <<"\", \"triangle_count\": "<<region_counts[name]<<" }";
  }
  regions<<"\n  }\n}\n";regions.close();
  std::filesystem::path source_path(__FILE__);
  std::string source_hash;
  if(std::filesystem::is_regular_file(source_path))source_hash=sha256(read_bytes(source_path));
  auto sidecar=path;sidecar.replace_extension(".geometry.json");std::ofstream meta(sidecar);
  if(!meta)throw std::runtime_error("cannot open geometry manifest");
  meta<<std::setprecision(17)<<"{\n  \"schema\": \"oneflow.stl-geometry-manifest\",\n  \"version\": 1,\n"
      <<"  \"generator\": { \"name\": \"oneflow_pipe_bend_geometry\", \"version\": \"1.0\", \"source\": \"generate_pipe_bend.cpp\", \"source_sha256\": \""<<source_hash<<"\" },\n"
      <<"  \"asset_id\": \"geometry.pipe_bend.90_degree\",\n  \"geometry\": \"pipe_bend_90_degree\",\n"
      <<"  \"units\": \""<<c.units<<"\",\n  \"parameters\": {"
      <<"\"radius\": "<<c.radius<<", \"bend_radius\": "<<c.bend_radius
      <<", \"inlet_length\": "<<c.inlet<<", \"outlet_length\": "<<c.outlet
      <<", \"angle_degrees\": "<<c.angle<<", \"circumferential_divisions\": "<<c.circum
      <<", \"bend_divisions\": "<<c.bend<<", \"inlet_divisions\": "<<c.inlet_div
      <<", \"outlet_divisions\": "<<c.outlet_div<<"},\n"
      <<"  \"surface\": { \"vertex_count\": "<<p.size()<<", \"triangle_count\": "<<f.size()<<" },\n"
      <<"  \"region_manifest\": \""<<region_manifest.filename().string()<<"\",\n"
      <<"  \"region_triangle_counts\": { \"wall\": "<<wall<<", \"inflow\": "<<inflow<<", \"outflow\": "<<outflow<<" },\n"
      <<"  \"stl_file\": \""<<path.filename().string()<<"\", \"stl_sha256\": \""<<sha256(read_bytes(path))<<"\",\n"
      <<"  \"intended_gmsh_physical_groups\": { \"surface\": [\"wall\",\"inflow\",\"outflow\"], \"volume\": [\"fluid_domain\"] },\n"
      <<"  \"notes\": \"Inlet/outlet caps close the fluid-domain shell; they are boundary patches, not solid walls.\"\n}\n";

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
