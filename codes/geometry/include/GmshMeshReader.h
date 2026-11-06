#pragma once
#include <array>
#include <string>
#include <vector>
namespace ONEFLOW {
struct GmshNode { long long tag{0}; std::array<double,3> coordinates{}; int entityDimension{0}; long long entityTag{0}; std::vector<double> parametricCoordinates; };
struct GmshElement { long long tag{0}; int type{0}; int entityDimension{0}; long long entityTag{0}; std::vector<long long> nodeTags; };
struct GmshMesh { std::string version; std::vector<GmshNode> nodes; std::vector<GmshElement> elements; };
GmshMesh ReadGmsh41AsciiMesh(const std::string& content);
}
