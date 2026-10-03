#pragma once

#include "HXDefine.h"

#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

enum class GridDistributionType
{
    Ratio,
    Distance,
    Tanh,
    Copy,
    Exponential
};

struct GridPointDefinition
{
    int id = 0;
    Real x = 0.0;
    Real y = 0.0;
    Real z = 0.0;
};

struct GridLineDefinition
{
    int id = 0;
    int p1 = 0;
    int p2 = 0;
};

struct GridCircleDefinition
{
    int id = 0;
    int p1 = 0;
    int pc = 0;
    int p2 = 0;
};

struct GridDimensionDefinition
{
    int id = 0;
    int pointCount = 0;
};

struct GridDistributionDefinition
{
    int lineId = 0;
    GridDistributionType type = GridDistributionType::Distance;
    Real value1 = 0.0;
    Real value2 = 0.0;
    std::vector< int > copyLineIds;
};

struct GridBoundaryDefinition
{
    int id = 0;
    int boundaryType = 0;
};

struct GridLineToFaceDefinition
{
    int faceId = 0;
    int position = 0;
    int lineId = 0;
};

struct GridFaceToBlockDefinition
{
    int blockId = 0;
    int position = 0;
    int faceId = 0;
};

class GridLayout
{
public:
    std::vector< GridPointDefinition > points;
    std::vector< GridLineDefinition > lines;
    std::vector< GridCircleDefinition > circles;
    std::vector< GridDimensionDefinition > dimensions;
    std::vector< GridDistributionDefinition > distributions;
    std::vector< GridBoundaryDefinition > boundaries;
    std::vector< GridLineToFaceDefinition > lineToFaces;
    std::vector< GridFaceToBlockDefinition > faceToBlocks;

    void Validate() const;
};

EndNameSpace
