#pragma once

#include <array>
#include <string>
#include <vector>

namespace ONEFLOW
{

// Value-only per-zone metadata. It has no dependency on runtime grid ownership.
struct GridInfoZoneSnapshot
{
    std::string name;
    int dimension = 0;
    int nodeCount = 0;
    int faceCount = 0;
    int cellCount = 0;
    int boundaryFaceCount = 0;
};

// L0 mesh summary contract: aggregate counts and the bounding box of all stored
// coordinates. Snapshot creation rejects null zones and non-finite coordinates.
// Serialization is handled separately by GridInfoJsonWriter.
struct GridInfoSnapshot
{
    std::string caseRoot;
    std::string projectRoot;
    std::string inputRoot;
    std::string outputRoot;
    std::string meshInputRoot;
    std::string flowInputRoot;
    std::string meshOutputRoot;
    std::string flowOutputRoot;
    std::string sourceGridPath;
    std::string outputGridPath;
    int dimension = 0;
    int nodeCount = 0;
    int faceCount = 0;
    int cellCount = 0;
    int boundaryFaceCount = 0;
    bool hasBoundingBox = false;
    std::array< double, 3 > minimum{};
    std::array< double, 3 > maximum{};
    std::vector< GridInfoZoneSnapshot > zones;
};

} // namespace ONEFLOW
