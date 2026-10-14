#pragma once

#include <array>
#include <string>
#include <vector>

namespace ONEFLOW
{

struct GridInfoZoneSnapshot
{
    std::string name;
    int dimension = 0;
    int nodeCount = 0;
    int faceCount = 0;
    int cellCount = 0;
    int boundaryFaceCount = 0;
};

struct GridInfoSnapshot
{
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
