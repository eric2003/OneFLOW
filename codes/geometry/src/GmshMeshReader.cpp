#include "GmshMeshReader.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>
namespace ONEFLOW { namespace {
long long integer(const std::string& s,const std::string& field) {
 std::size_t n=0; long long v=0;
 try { v=std::stoll(s,&n); } catch(...) { throw std::runtime_error("Invalid Gmsh "+field+": "+s); }
 if(n!=s.size()) throw std::runtime_error("Invalid Gmsh "+field+": "+s); return v;
}
double real(const std::string& s,const std::string& field) {
 std::size_t n=0; double v=0;
 try { v=std::stod(s,&n); } catch(...) { throw std::runtime_error("Invalid Gmsh "+field+": "+s); }
 if(n!=s.size()||!std::isfinite(v)) throw std::runtime_error("Invalid Gmsh "+field+": "+s); return v;
}
int nodeCount(int t) {
 switch(t) {
 case 1:return 2;case 2:return 3;case 3:return 4;case 4:return 4;case 5:return 8;
 case 6:return 6;case 7:return 5;case 8:return 3;case 9:return 6;case 10:return 9;
 case 11:return 10;case 12:return 27;case 13:return 18;case 14:return 14;case 15:return 1;
 case 16:return 8;case 17:return 20;case 18:return 15;case 19:return 13;case 20:return 9;
 case 21:return 10;case 22:return 12;case 23:return 15;case 24:return 15;case 25:return 21;
 case 26:return 4;case 27:return 5;case 28:return 6;case 29:return 20;case 30:return 35;case 31:return 56;
 default:throw std::runtime_error("Unsupported Gmsh element type "+std::to_string(t));
 }
}
class Reader {
 std::vector<std::string> t;
public:
 explicit Reader(const std::string& s){std::istringstream in(s);std::string x;while(in>>x)t.push_back(x);}
 std::size_t section(const std::string& s)const{for(std::size_t i=0;i<t.size();++i)if(t[i]==s)return i+1;throw std::runtime_error("Missing Gmsh section "+s);}
 std::string next(std::size_t& p,const std::string& f)const{if(p>=t.size())throw std::runtime_error("Truncated Gmsh "+f);return t[p++];}
 long long i(std::size_t& p,const std::string& f)const{return integer(next(p,f),f);}
 double r(std::size_t& p,const std::string& f)const{return real(next(p,f),f);}
 void end(std::size_t p,const std::string& s)const{if(p>=t.size()||t[p]!=s)throw std::runtime_error("Malformed Gmsh terminator "+s);}
};
}
GmshMesh ReadGmsh41AsciiMesh(const std::string& content) {
 Reader r(content); GmshMesh m; std::size_t p=r.section("$MeshFormat");
 double version=r.r(p,"version"); long long binary=r.i(p,"file type"),size=r.i(p,"data size");
 if(version<4.1||version>=4.2||binary!=0||(size!=4&&size!=8))throw std::runtime_error("Expected MSH 4.1 ASCII");
 r.end(p,"$EndMeshFormat");m.version="4.1";
 std::unordered_set<long long> nodeIds;p=r.section("$Nodes");
 long long nb=r.i(p,"node blocks"),nn=r.i(p,"node count"),nmin=r.i(p,"minimum node tag"),nmax=r.i(p,"maximum node tag");
 if(nb<0||nn<1||nmin<1||nmax<nmin)throw std::runtime_error("Invalid Gmsh node header");
 m.nodes.reserve(static_cast<std::size_t>(nn));long long omin=std::numeric_limits<long long>::max(),omax=0;
 for(long long b=0;b<nb;++b) {
  long long dim=r.i(p,"node entity dimension"),ent=r.i(p,"node entity tag"),param=r.i(p,"parametric flag"),count=r.i(p,"node block size");
  if(dim<0||dim>3||ent<1||(param!=0&&param!=1)||count<0)throw std::runtime_error("Invalid Gmsh node block");
  std::vector<long long> ids;ids.reserve(static_cast<std::size_t>(count));
  for(long long k=0;k<count;++k){long long id=r.i(p,"node tag");if(id<1||!nodeIds.insert(id).second)throw std::runtime_error("Duplicate or invalid Gmsh node tag");ids.push_back(id);omin=std::min(omin,id);omax=std::max(omax,id);}
  for(long long k=0;k<count;++k){GmshNode n;n.tag=ids[static_cast<std::size_t>(k)];n.entityDimension=static_cast<int>(dim);n.entityTag=ent;for(double& x:n.coordinates)x=r.r(p,"coordinate");for(long long a=0;a<(param?dim:0);++a)n.parametricCoordinates.push_back(r.r(p,"parametric coordinate"));m.nodeTagToIndex.emplace(n.tag, m.nodes.size());
   m.nodes.push_back(std::move(n));}
 }
 if(static_cast<long long>(m.nodes.size())!=nn||omin!=nmin||omax!=nmax)throw std::runtime_error("Gmsh node totals or tag range mismatch");r.end(p,"$EndNodes");
 std::unordered_set<long long> elemIds;p=r.section("$Elements");
 long long eb=r.i(p,"element blocks"),ne=r.i(p,"element count"),emin=r.i(p,"minimum element tag"),emax=r.i(p,"maximum element tag");
 if(eb<0||ne<1||emin<1||emax<emin)throw std::runtime_error("Invalid Gmsh element header");
 m.elements.reserve(static_cast<std::size_t>(ne));long long oeMin=std::numeric_limits<long long>::max(),oeMax=0;
 for(long long b=0;b<eb;++b) {
  long long dim=r.i(p,"element entity dimension"),ent=r.i(p,"element entity tag"),tv=r.i(p,"element type"),count=r.i(p,"element block size");
  if(dim<0||dim>3||ent<1||tv<1||tv>std::numeric_limits<int>::max()||count<0)throw std::runtime_error("Invalid Gmsh element block");
  int type=static_cast<int>(tv),arity=nodeCount(type);
  for(long long k=0;k<count;++k){GmshElement e;e.tag=r.i(p,"element tag");if(e.tag<1||!elemIds.insert(e.tag).second)throw std::runtime_error("Duplicate or invalid Gmsh element tag");oeMin=std::min(oeMin,e.tag);oeMax=std::max(oeMax,e.tag);e.type=type;e.entityDimension=static_cast<int>(dim);e.entityTag=ent;for(int a=0;a<arity;++a){long long id=r.i(p,"element node tag");if(!nodeIds.count(id))throw std::runtime_error("Gmsh element references unknown node");e.nodeTags.push_back(id);}m.elements.push_back(std::move(e));}
 }
 if(static_cast<long long>(m.elements.size())!=ne||oeMin!=emin||oeMax!=emax)throw std::runtime_error("Gmsh element totals or tag range mismatch");r.end(p,"$EndElements");return m;
}
}


namespace ONEFLOW {

GmshMeshReadReport BuildGmshMeshReadReport(const GmshMesh& mesh) {
    GmshMeshReadReport report;
    report.version = mesh.version;
    report.nodeCount = mesh.nodes.size();
    report.elementCount = mesh.elements.size();
    if (!mesh.nodes.empty()) {
        report.minimumNodeTag = mesh.nodes.front().tag;
        report.maximumNodeTag = mesh.nodes.front().tag;
        report.minimumCoordinates = mesh.nodes.front().coordinates;
        report.maximumCoordinates = mesh.nodes.front().coordinates;
        for (const GmshNode& node : mesh.nodes) {
            report.minimumNodeTag = std::min(report.minimumNodeTag, node.tag);
            report.maximumNodeTag = std::max(report.maximumNodeTag, node.tag);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                report.minimumCoordinates[axis] =
                    std::min(report.minimumCoordinates[axis], node.coordinates[axis]);
                report.maximumCoordinates[axis] =
                    std::max(report.maximumCoordinates[axis], node.coordinates[axis]);
            }
        }
    }
    if (!mesh.elements.empty()) {
        report.minimumElementTag = mesh.elements.front().tag;
        report.maximumElementTag = mesh.elements.front().tag;
        for (const GmshElement& element : mesh.elements) {
            report.minimumElementTag = std::min(report.minimumElementTag, element.tag);
            report.maximumElementTag = std::max(report.maximumElementTag, element.tag);
            ++report.elementTypeCounts[element.type];
            if (element.type != 2 || element.nodeTags.size() != 3) continue;
            ++report.triangleCount;
            const auto a = mesh.nodeTagToIndex.find(element.nodeTags[0]);
            const auto b = mesh.nodeTagToIndex.find(element.nodeTags[1]);
            const auto c = mesh.nodeTagToIndex.find(element.nodeTags[2]);
            if (a == mesh.nodeTagToIndex.end() ||
                b == mesh.nodeTagToIndex.end() ||
                c == mesh.nodeTagToIndex.end()) {
                throw std::runtime_error("Cannot report TRI_3 with unresolved node tag");
            }
            const auto& p0 = mesh.nodes[a->second].coordinates;
            const auto& p1 = mesh.nodes[b->second].coordinates;
            const auto& p2 = mesh.nodes[c->second].coordinates;
            const double ux = p1[0] - p0[0];
            const double uy = p1[1] - p0[1];
            const double uz = p1[2] - p0[2];
            const double vx = p2[0] - p0[0];
            const double vy = p2[1] - p0[1];
            const double vz = p2[2] - p0[2];
            const double cx = uy * vz - uz * vy;
            const double cy = uz * vx - ux * vz;
            const double cz = ux * vy - uy * vx;
            const double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
            const double signedXYArea = 0.5 * (ux * vy - uy * vx);
            const double scale = std::max({
                report.maximumCoordinates[0] - report.minimumCoordinates[0],
                report.maximumCoordinates[1] - report.minimumCoordinates[1],
                report.maximumCoordinates[2] - report.minimumCoordinates[2]
            });
            const double areaTolerance =
                64.0 * std::numeric_limits<double>::epsilon() * scale * scale;
            if (area <= areaTolerance) ++report.degenerateTriangleCount;
            if (signedXYArea > areaTolerance) ++report.positiveXYTriangleCount;
            else if (signedXYArea < -areaTolerance) ++report.negativeXYTriangleCount;
            if (report.triangleCount == 1) {
                report.minimumTriangleArea = area;
                report.maximumTriangleArea = area;
            } else {
                report.minimumTriangleArea = std::min(report.minimumTriangleArea, area);
                report.maximumTriangleArea = std::max(report.maximumTriangleArea, area);
            }
        }
    }
    return report;
}

std::string GmshMeshReadReportToJson(const GmshMeshReadReport& report) {
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n"
        << "  \"schema\": \"oneflow.gmsh-read-report\",\n"
        << "  \"schema_version\": 1,\n"
        << "  \"status\": \"READ_SUCCESS\",\n"
        << "  \"format_version\": \"" << report.version << "\",\n"
        << "  \"node_count\": " << report.nodeCount << ",\n"
        << "  \"element_count\": " << report.elementCount << ",\n"
        << "  \"node_tag_range\": [" << report.minimumNodeTag << ", "
        << report.maximumNodeTag << "],\n"
        << "  \"element_tag_range\": [" << report.minimumElementTag << ", "
        << report.maximumElementTag << "],\n"
        << "  \"coordinate_bounds\": {\n"
        << "    \"minimum\": [" << report.minimumCoordinates[0] << ", "
        << report.minimumCoordinates[1] << ", " << report.minimumCoordinates[2] << "],\n"
        << "    \"maximum\": [" << report.maximumCoordinates[0] << ", "
        << report.maximumCoordinates[1] << ", " << report.maximumCoordinates[2] << "]\n"
        << "  },\n"
        << "  \"element_type_counts\": {";
    bool first = true;
    for (const auto& entry : report.elementTypeCounts) {
        if (!first) out << ",";
        out << "\n    \"" << entry.first << "\": " << entry.second;
        first = false;
    }
    if (!first) out << "\n  ";
    out << "},\n"
        << "  \"tri_3_geometry\": {\n"
        << "    \"count\": " << report.triangleCount << ",\n"
        << "    \"degenerate_count\": " << report.degenerateTriangleCount << ",\n"
        << "    \"positive_xy_orientation_count\": " << report.positiveXYTriangleCount << ",\n"
        << "    \"negative_xy_orientation_count\": " << report.negativeXYTriangleCount << ",\n"
        << "    \"minimum_area\": " << report.minimumTriangleArea << ",\n"
        << "    \"maximum_area\": " << report.maximumTriangleArea << "\n"
        << "  },\n"
        << "  \"validation\": {\n"
        << "    \"unique_source_tags\": \"PASS\",\n"
        << "    \"node_references\": \"PASS\",\n"
        << "    \"finite_coordinates\": \"PASS\",\n"
        << "    \"connectivity_arity\": \"PASS\",\n"
        << "    \"physical_groups\": \"NOT_YET_PARSED\",\n"
        << "    \"cgns_round_trip\": \"NOT_YET_AVAILABLE\"\n"
        << "  }\n"
        << "}\n";
    return out.str();
}


} // namespace ONEFLOW
