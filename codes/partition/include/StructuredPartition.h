/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "HXDefine.h"
#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

// Axis-aligned, cell-count-balanced block partition.
//
// METIS minimizes the edge cut of a graph; it does not guarantee that every
// part ends up with the same interface count.  A weak-scaling measurement on
// top of a METIS partition therefore measures the partitioner as much as the
// solver: the per-part interface count can vary by more than 2x inside one
// run and grow with the part count.
//
// This partitioner cuts the bounding box into an nx*ny*nz array of blocks and
// places every axis cut at a cell-count quantile, so each part holds the same
// number of cells (within one cell) and the interface count is proportional to
// the block surface.  The cuts follow the cell distribution rather than the
// coordinates, so graded meshes stay balanced as well.
//
// The split is deterministic and does not depend on METIS.
struct StructuredSplit
{
    int nx = 1;
    int ny = 1;
    int nz = 1;
};

// Effective number of cells along one axis: (max - min) / median positive
// spacing.  This is the length the split heuristic wants, not the physical
// box size: a unit cube refined 4x more in z must be treated as four times
// longer, otherwise the blocks come out physically cubic and the interface
// count stops matching between configurations.  The median spacing is used
// so a few clustered cells do not dominate.  Returns 0 when the axis is
// degenerate (fewer than two distinct coordinates) or the input is invalid.
Real EffectiveAxisLength( const Real * coordinate, int nCells );

// Pick the factorization nx*ny*nz == nPart that minimizes the total internal
// interface area of a block array with effective lengths (lx, ly, lz), as
// returned by EffectiveAxisLength.  A zero or negative length means that the
// axis is not split.  nPart < 1 returns {1, 1, 1}.
StructuredSplit ChooseStructuredSplit( int nPart, Real lx, Real ly, Real lz );

// Fill zoneOfCell with the block id in [0, nx*ny*nz).  Every part receives
// either floor(nCells/nPart) or ceil(nCells/nPart) cells.  Returns false on
// invalid input (empty grid, null coordinates, inconsistent split).
bool AssignStructuredZones(
    int nCells,
    const Real * x,
    const Real * y,
    const Real * z,
    const StructuredSplit & split,
    IntField & zoneOfCell );

EndNameSpace
