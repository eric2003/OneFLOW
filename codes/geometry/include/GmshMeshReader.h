#pragma once
#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>
namespace ONEFLOW {
struct GmshNode { long long tag{0}; std::array<double,3> coordinates{}; int entityDimension{0}; long long entityTag{0}; std::vector<double> parametricCoordinates; };
struct GmshElement { long long tag{0}; int type{0}; int entityDimension{0}; long long entityTag{0}; std::vector<long long> nodeTags; };
struct GmshMesh {
    std::string version;
    std::vector<GmshNode> nodes;
    std::vector<GmshElement> elements;
    // Maps sparse source tags to stable zero-based positions in nodes.
    std::map<long long, std::size_t> nodeTagToIndex;
    std::map<std::pair<int, long long>, std::string> physicalNames;
    std::map<std::pair<int, long long>, std::vector<long long>> entityPhysicalTags;
};
struct GmshMeshReadReport {
    std::string version;
    std::size_t nodeCount{0};
    std::size_t elementCount{0};
    std::size_t physicalGroupCount{0};
    long long minimumNodeTag{0};
    long long maximumNodeTag{0};
    long long minimumElementTag{0};
    long long maximumElementTag{0};
    std::map<int, std::size_t> elementTypeCounts;
    std::array<double, 3> minimumCoordinates{};
    std::array<double, 3> maximumCoordinates{};
    std::size_t triangleCount{0};
    std::size_t degenerateTriangleCount{0};
    std::size_t positiveXYTriangleCount{0};
    std::size_t negativeXYTriangleCount{0};
    double minimumTriangleArea{0.0};
    double maximumTriangleArea{0.0};
};
GmshMesh ReadGmsh41AsciiMesh(const std::string& content);
GmshMeshReadReport BuildGmshMeshReadReport(const GmshMesh& mesh);
std::string GmshMeshReadReportToJson(const GmshMeshReadReport& report);
std::string GmshMeshReportToMarkdown(const GmshMesh& mesh, const std::string& cgnsOutputFile);
}
