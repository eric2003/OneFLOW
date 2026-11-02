/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "HXType.h"

#include <cstddef>
#include <vector>

namespace ONEFLOW
{

enum class BoundaryMemorySpace
{
    Host,
    Device
};

enum class BoundaryTransferMode
{
    HostStaged,
    DeviceDirect
};

// Endpoint ownership is explicit because a partition interface sends values
// from real cells and receives values into ghost cells.  Unspecified keeps
// the value-only contract usable by legacy callers that do not have a cell
// ownership split.
enum class BoundaryEndpointRole
{
    Unspecified,
    OwnedCell,
    GhostCell
};

struct BoundaryBufferView
{
    Real * values = nullptr;
    std::size_t valueCount = 0;
    BoundaryMemorySpace memorySpace = BoundaryMemorySpace::Host;
};

struct BoundaryEndpointMap
{
    int cellCount = 0;
    std::vector< int > cellIds;
    int ownedCellCount = -1;
    BoundaryEndpointRole role = BoundaryEndpointRole::Unspecified;

    void Validate() const;
};

// Describes one directed partition-boundary transfer. Values are packed in
// equation-major order: packed[equation * interfaceCount + interface].
struct BoundaryFieldExchangePlan
{
    int sourceZone = -1;
    int targetZone = -1;
    int sourceRank = -1;
    int targetRank = -1;
    int ghostLevel = 0;
    int nEquations = 0;
    int sourceCellCount = 0;
    int targetCellCount = 0;
    std::vector< int > sendCellIds;
    std::vector< int > receiveCellIds;

    void Validate() const;
    std::size_t InterfaceCount() const;
    std::size_t ValueCount() const;
};

// Backend-neutral description of one directed interface exchange.  The
// source endpoint must name owned cells and the target endpoint must name
// ghost cells; transfer code is responsible only for the selected buffers.
struct BoundaryExchangeContract
{
    int sourceZone = -1;
    int targetZone = -1;
    int sourceRank = -1;
    int targetRank = -1;
    int nEquations = 0;
    BoundaryEndpointMap source;
    BoundaryEndpointMap target;

    void Validate() const;
    std::size_t InterfaceCount() const;
    std::size_t ValueCount() const;
};

void ValidateBoundaryBuffer(
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    const BoundaryBufferView & buffer,
    BoundaryTransferMode mode );

void ValidateBoundaryExchangeBuffers(
    const BoundaryExchangeContract & contract,
    const BoundaryBufferView & send,
    const BoundaryBufferView & receive,
    BoundaryTransferMode mode );

void PackBoundaryHost(
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    const std::vector< const Real * > & sourceComponents,
    std::vector< Real > & packedValues );

void UnpackBoundaryHost(
    const BoundaryEndpointMap & endpoint,
    int nEquations,
    const std::vector< Real > & packedValues,
    const std::vector< Real * > & targetComponents );

void PackBoundaryHost(
    const BoundaryFieldExchangePlan & plan,
    const std::vector< const Real * > & sourceComponents,
    std::vector< Real > & packedValues );

void UnpackBoundaryHost(
    const BoundaryFieldExchangePlan & plan,
    const std::vector< Real > & packedValues,
    const std::vector< Real * > & targetComponents );

}
