#include "GmshMeshReader.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <limits>
#include <iomanip>
#include <map>
#include <set>
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
 explicit Reader(const std::string& s) {
  for (std::size_t p=0; p<s.size();) {
   while (p<s.size() && std::isspace(static_cast<unsigned char>(s[p]))) ++p;
   if (p==s.size()) break;
   std::string token;
   if (s[p]=='"') {
    ++p;
    bool closed=false;
    while (p<s.size()) {
     char c=s[p++];
     if (c=='"') { closed=true; break; }
     if (c=='\\' && p<s.size()) c=s[p++];
     token.push_back(c);
    }
    if (!closed) throw std::runtime_error("Unterminated quoted Gmsh string");
   } else {
    while (p<s.size() && !std::isspace(static_cast<unsigned char>(s[p]))) token.push_back(s[p++]);
   }
   t.push_back(std::move(token));
  }
 }
 bool hasSection(const std::string& s) const { for(const auto& x:t) if(x==s) return true; return false; }
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
 if (r.hasSection("$PhysicalNames")) {
  p=r.section("$PhysicalNames");
  const long long count=r.i(p,"physical name count");
  if(count<0) throw std::runtime_error("Invalid Gmsh physical name count");
  for(long long i=0;i<count;++i) {
   const long long dim=r.i(p,"physical name dimension");
   const long long tag=r.i(p,"physical name tag");
   const std::string name=r.next(p,"physical name");
   if(dim<0||dim>3||tag<1||name.empty()) throw std::runtime_error("Invalid Gmsh physical name");
   if(!m.physicalNames.emplace(std::make_pair(static_cast<int>(dim),tag),name).second)
    throw std::runtime_error("Duplicate Gmsh physical name");
  }
  r.end(p,"$EndPhysicalNames");
 }
 if (r.hasSection("$Entities")) {
  p=r.section("$Entities");
  long long counts[4];
  for(int dim=0;dim<4;++dim) { counts[dim]=r.i(p,"entity count"); if(counts[dim]<0) throw std::runtime_error("Invalid Gmsh entity count"); }
  for(int dim=0;dim<4;++dim) for(long long i=0;i<counts[dim];++i) {
   const long long tag=r.i(p,"entity tag");
   if(tag<1) throw std::runtime_error("Invalid Gmsh entity tag");
   const int coordinateCount=dim==0?3:6;
   for(int axis=0;axis<coordinateCount;++axis) (void)r.r(p,"entity bounding coordinate");
   const long long physicalCount=r.i(p,"entity physical tag count");
   if(physicalCount<0) throw std::runtime_error("Invalid Gmsh entity physical tag count");
   std::vector<long long> tags;
   for(long long j=0;j<physicalCount;++j) {
    const long long physicalTag=r.i(p,"entity physical tag");
    if(physicalTag<1) throw std::runtime_error("Invalid Gmsh entity physical tag");
    tags.push_back(physicalTag);
   }
   if(dim>0) {
    const long long boundaryCount=r.i(p,"entity boundary count");
    if(boundaryCount<0) throw std::runtime_error("Invalid Gmsh entity boundary count");
    for(long long j=0;j<boundaryCount;++j) (void)r.i(p,"bounding entity tag");
   }
   if(!m.entityPhysicalTags.emplace(std::make_pair(dim,tag),std::move(tags)).second)
    throw std::runtime_error("Duplicate Gmsh entity");
  }
  r.end(p,"$EndEntities");
 }
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

GmshTetrahedralTopologyReport ValidateGmshTetrahedralTopology(const GmshMesh& mesh) {
    using FaceKey = std::array<long long, 3>;
    std::map<FaceKey, std::size_t> faceIncidence;
    std::set<FaceKey> surfaceFaces;
    GmshTetrahedralTopologyReport report;

    for (const GmshElement& element : mesh.elements) {
        if (element.entityDimension == 3) {
            if (element.type != 4 || element.nodeTags.size() != 4) {
                throw std::runtime_error(
                    "Tetrahedral topology validation requires first-order TETRA_4 volume cells");
            }
            std::set<long long> uniqueNodes(element.nodeTags.begin(), element.nodeTags.end());
            if (uniqueNodes.size() != 4) {
                throw std::runtime_error("TETRA_4 element contains repeated node tags");
            }
            ++report.tetrahedronCount;
            for (std::size_t omitted = 0; omitted < 4; ++omitted) {
                FaceKey face{};
                std::size_t target = 0;
                for (std::size_t node = 0; node < 4; ++node) {
                    if (node != omitted) face[target++] = element.nodeTags[node];
                }
                std::sort(face.begin(), face.end());
                const std::size_t count = ++faceIncidence[face];
                if (count > 2) {
                    throw std::runtime_error(
                        "Non-manifold tetrahedral face is shared by more than two cells");
                }
            }
        } else if (element.entityDimension == 2) {
            if (element.type != 2 || element.nodeTags.size() != 3) {
                throw std::runtime_error(
                    "Tetrahedral boundary validation requires first-order TRI_3 surface elements");
            }
            FaceKey face{{element.nodeTags[0], element.nodeTags[1], element.nodeTags[2]}};
            std::sort(face.begin(), face.end());
            if (!surfaceFaces.insert(face).second) {
                throw std::runtime_error("Duplicate TRI_3 boundary element detected");
            }
        }
    }

    if (report.tetrahedronCount == 0) {
        throw std::runtime_error("Tetrahedral topology validation found no TETRA_4 volume cells");
    }
    report.surfaceTriangleCount = surfaceFaces.size();
    for (const auto& entry : faceIncidence) {
        if (entry.second == 1) {
            ++report.boundaryFaceCount;
            if (surfaceFaces.find(entry.first) == surfaceFaces.end()) {
                std::ostringstream message;
                message << "TETRA_4 exterior face has no matching TRI_3 boundary element: "
                        << entry.first[0] << ", " << entry.first[1] << ", " << entry.first[2];
                throw std::runtime_error(message.str());
            }
        } else {
            ++report.internalFaceCount;
        }
    }
    for (const FaceKey& face : surfaceFaces) {
        if (faceIncidence.find(face) == faceIncidence.end()) {
            std::ostringstream message;
            message << "TRI_3 surface element does not match any TETRA_4 face: "
                    << face[0] << ", " << face[1] << ", " << face[2];
            throw std::runtime_error(message.str());
        }
    }
    return report;
}

GmshMeshReadReport BuildGmshMeshReadReport(const GmshMesh& mesh) {
    GmshMeshReadReport report;
    report.version = mesh.version;
    report.nodeCount = mesh.nodes.size();
    report.elementCount = mesh.elements.size();
    report.physicalGroupCount = mesh.physicalNames.size();
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
        << "  \"physical_group_count\": " << report.physicalGroupCount << ",\n"
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
        << "    \"physical_groups\": \"PARSED\",\n"
        << "    \"cgns_round_trip\": \"NOT_YET_AVAILABLE\"\n"
        << "  }\n"
        << "}\n";
    return out.str();
}


namespace {
std::string GmshElementTypeName(int type) {
    switch (type) {
        case 1: return "LINE_2"; case 2: return "TRI_3"; case 3: return "QUAD_4";
        case 4: return "TETRA_4"; case 5: return "HEXA_8"; case 6: return "PENTA_6";
        case 7: return "PYRA_5"; case 8: return "LINE_3"; case 9: return "TRI_6";
        case 10: return "QUAD_9"; case 11: return "TETRA_10"; case 12: return "HEXA_27";
        case 13: return "PENTA_18"; case 14: return "PYRA_14"; case 15: return "POINT";
        case 16: return "QUAD_8"; case 17: return "HEXA_20"; case 18: return "PENTA_15";
        case 19: return "PYRA_13"; case 20: return "TRI_9"; case 21: return "TRI_10";
        case 22: return "TRI_12"; case 23: return "TRI_15"; case 24: return "TRI_15I";
        case 25: return "TRI_21"; case 26: return "LINE_4"; case 27: return "LINE_5";
        case 28: return "LINE_6"; case 29: return "TETRA_20"; case 30: return "TETRA_35";
        case 31: return "TETRA_56"; default: return "UNKNOWN";
    }
}
std::string EscapeMarkdownCell(const std::string& value) {
    std::string result;
    for (char c : value) {
        if (c == '|' || c == '\\') result.push_back('\\');
        result.push_back((c == '\r' || c == '\n') ? ' ' : c);
    }
    return result;
}
}

std::string GmshMeshReportToMarkdown(
    const GmshMesh& mesh, const std::string& cgnsOutputFile) {
    const GmshMeshReadReport report = BuildGmshMeshReadReport(mesh);
    int dimension = 0;
    for (const auto& element : mesh.elements)
        dimension = std::max(dimension, element.entityDimension);
    const double scale = std::max({
        report.maximumCoordinates[0] - report.minimumCoordinates[0],
        report.maximumCoordinates[1] - report.minimumCoordinates[1],
        report.maximumCoordinates[2] - report.minimumCoordinates[2], 1.0});
    const bool planarXY = std::abs(report.maximumCoordinates[2] -
        report.minimumCoordinates[2]) <= 1.0e-12 * scale;
    std::map<int, std::size_t> entityCounts;
    for (const auto& entity : mesh.entityPhysicalTags) ++entityCounts[entity.first.first];
    std::map<std::pair<int, long long>, std::size_t> groupElementCounts;
    for (const auto& element : mesh.elements) {
        const auto entity = mesh.entityPhysicalTags.find(
            {element.entityDimension, element.entityTag});
        if (entity == mesh.entityPhysicalTags.end()) continue;
        for (long long tag : entity->second)
            ++groupElementCounts[{element.entityDimension, tag}];
    }

    std::ostringstream out;
    out << std::setprecision(10);
    out << "# Gmsh Mesh Inspection Report\n\n## Overall Summary\n\n";
    out << "- **Format version:** Gmsh MSH " << report.version << " (ASCII)\n";
    out << "- **Topological dimension:** " << dimension << "D\n";
    out << "- **Coordinate embedding:** " << (planarXY ? "planar in XY" : "3D coordinates") << "\n";
    out << "- **Physical groups:** " << mesh.physicalNames.size() << "\n";
    out << "- **Geometric entities:** " << mesh.entityPhysicalTags.size() << "\n";
    out << "- **Mesh size:** " << mesh.nodes.size() << " nodes, " << mesh.elements.size() << " elements\n";
    out << "- **Coordinate bounds:** X=[" << report.minimumCoordinates[0] << ", " << report.maximumCoordinates[0]
        << "], Y=[" << report.minimumCoordinates[1] << ", " << report.maximumCoordinates[1]
        << "], Z=[" << report.minimumCoordinates[2] << ", " << report.maximumCoordinates[2] << "]\n";
    out << "- **Source node-tag range:** [" << report.minimumNodeTag << ", " << report.maximumNodeTag << "]\n";
    out << "- **Source element-tag range:** [" << report.minimumElementTag << ", " << report.maximumElementTag << "]\n\n";

    out << "## Physical Groups\n\n| Dimension | Tag | Name | Associated elements |\n|---:|---:|---|---:|\n";
    for (const auto& group : mesh.physicalNames) {
        const auto count = groupElementCounts.find(group.first);
        out << "| " << group.first.first << " | " << group.first.second << " | "
            << EscapeMarkdownCell(group.second) << " | "
            << (count == groupElementCounts.end() ? 0 : count->second) << " |\n";
    }
    out << "\n## Geometric Entities ($Entities)\n\n| Dimension | Entity count |\n|---:|---:|\n";
    for (int dim = 0; dim <= 3; ++dim) {
        const auto count = entityCounts.find(dim);
        out << "| " << dim << " | " << (count == entityCounts.end() ? 0 : count->second) << " |\n";
    }
    out << "\n| Dimension | Entity tag | Physical tags |\n|---:|---:|---|\n";
    for (const auto& entity : mesh.entityPhysicalTags) {
        out << "| " << entity.first.first << " | " << entity.first.second << " | ";
        if (entity.second.empty()) out << "None";
        for (std::size_t i = 0; i < entity.second.size(); ++i) {
            if (i) out << ", ";
            out << entity.second[i];
        }
        out << " |\n";
    }
    out << "\n## Element Types and Counts\n\n| Gmsh type | Element family | Count |\n|---:|---|---:|\n";
    for (const auto& item : report.elementTypeCounts)
        out << "| " << item.first << " | " << GmshElementTypeName(item.first) << " | " << item.second << " |\n";

    auto cgnsElementTypeName = [](int type) -> const char* {
        switch (type) {
            case 1: return "BAR_2";
            case 2: return "TRI_3";
            case 3: return "QUAD_4";
            case 4: return "TETRA_4";
            case 5: return "HEXA_8";
            case 6: return "PENTA_6";
            case 7: return "PYRA_5";
            default: return "UNSUPPORTED";
        }
    };
    std::map<std::pair<int, std::string>, std::size_t> cgnsBoundaryCounts;
    for (const auto& element : mesh.elements) {
        if (element.entityDimension != dimension - 1) continue;
        const auto entity = mesh.entityPhysicalTags.find(
            {element.entityDimension, element.entityTag});
        if (entity != mesh.entityPhysicalTags.end() && !entity->second.empty())
            ++cgnsBoundaryCounts[{element.type, cgnsElementTypeName(element.type)}];
    }
    out << "\n## CGNS Boundary Element Sections\n\n"
        << "| Gmsh type | CGNS section type | Count |\n|---:|---|---:|\n";
    for (const auto& item : cgnsBoundaryCounts)
        out << "| " << item.first.first << " | " << item.first.second
            << " | " << item.second << " |\n";

    const std::size_t limit = 30;
    out << "\n## Nodes ($Nodes)\n\nSource node tags are identifiers, not dense array indices. Original tags are retained while an internal index map is built.\n\n";
    out << "| Source tag | X | Y | Z | Entity (dimension:tag) |\n|---:|---:|---:|---:|---|\n";
    for (std::size_t i = 0; i < std::min(mesh.nodes.size(), limit); ++i) {
        const auto& node = mesh.nodes[i];
        out << "| " << node.tag << " | " << node.coordinates[0] << " | "
            << node.coordinates[1] << " | " << node.coordinates[2] << " | "
            << node.entityDimension << ":" << node.entityTag << " |\n";
    }
    if (mesh.nodes.size() > limit)
        out << "\n_Only the first " << limit << " of " << mesh.nodes.size() << " nodes are shown._\n";

    out << "\n## Elements ($Elements)\n\nConnectivity uses original Gmsh node tags and preserves source ordering.\n\n";
    out << "| Element tag | Type | Entity (dimension:tag) | Ordered node tags |\n|---:|---|---|---|\n";
    for (std::size_t i = 0; i < std::min(mesh.elements.size(), limit); ++i) {
        const auto& element = mesh.elements[i];
        out << "| " << element.tag << " | " << GmshElementTypeName(element.type)
            << " (" << element.type << ") | " << element.entityDimension << ":"
            << element.entityTag << " | ";
        for (std::size_t j = 0; j < element.nodeTags.size(); ++j) {
            if (j) out << ", ";
            out << element.nodeTags[j];
        }
        out << " |\n";
    }
    if (mesh.elements.size() > limit)
        out << "\n_Only the first " << limit << " of " << mesh.elements.size() << " elements are shown._\n";

    out << "\n## Geometry and Mesh Checks\n\n";
    out << "- **Unique source tags:** PASS\n- **Element node references:** PASS\n- **Finite coordinates:** PASS\n";
    out << "- **TRI_3 count:** " << report.triangleCount << "\n";
    out << "- **Degenerate TRI_3 count:** " << report.degenerateTriangleCount << "\n";
    out << "- **Positive XY orientation:** " << report.positiveXYTriangleCount << "\n";
    out << "- **Negative XY orientation:** " << report.negativeXYTriangleCount << "\n";
    if (report.triangleCount)
        out << "- **TRI_3 area range:** [" << report.minimumTriangleArea << ", " << report.maximumTriangleArea << "]\n";

    out << "\n## CGNS Conversion\n\n";
    out << "- **Output file:** " << cgnsOutputFile << "\n";
    out << "- **Status:** CGNS writer completed successfully before this report was emitted.\n";
    out << "- **Boundary semantics:** Physical group names are retained as user-defined labels; solver-specific wall/inlet/outlet meanings are not inferred.\n";
    out << "- **Independent read-back:** Covered by the CGNS conversion test for zones, coordinates, cell and boundary sections, and BC point lists.\n";
    out << "\n## Current Scope and Limitations\n\n";
    out << "- Input: ASCII MSH 4.1.\n";
    out << "- Conversion: supported first-order cells and codimension-one boundary elements.\n";
    out << "- High-order conversion, periodic links, partition metadata, and ghost-cell metadata are not yet implemented.\n";
    out << "- Physical names are metadata, not solver boundary-condition assignments.\n";
    return out.str();
}


} // namespace ONEFLOW
