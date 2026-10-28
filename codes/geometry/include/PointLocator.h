/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#pragma once
#include "HXDefine.h"
#include "HXAdtTree.h"
#include "GridHandles.h"

BeginNameSpace( ONEFLOW )

class Grid;

using AdtTree = HXAdtTree<int, Real>;
using AdtNode = HXAdtNode<int, Real>;

/**
 * Locate and deduplicate points using an ADT spatial index.
 *
 * The locator owns its spatial index and coordinate store. FindPoint is a
 * query that updates the last-result ID; AddPoint may also extend the store.
 * Matching tolerance is established by Initialize or InitializeSpecial.
 */
class PointLocator
{
public:
    PointLocator();
    // Virtual so IFaceLink tests can inject a lightweight double without an ADT.
    virtual ~PointLocator();

    // [Refactored] Disable copying to prevent double-free of unique_ptr.
    // unique_ptr is move-only, so we explicitly delete copy operations.
    PointLocator(const PointLocator&) = delete;
    PointLocator& operator=(const PointLocator&) = delete;
protected:
    int id;

    // [Refactored] Changed from raw pointer to std::unique_ptr.
    // Automatic cleanup, exception safety, no double-free possible.
    std::unique_ptr<AdtTree> coorTree;

    Real tolerance;
    RealField xCoor, yCoor, zCoor;
public:
    /**
     * Initialize the search bounds and ADT using an explicit tolerance.
     * The input bounds are read-only; the expanded bounds are stored locally.
     */
    void Initialize( const RealField & pmin, const RealField & pmax, Real toleranceIn );

    /// Initialize from one grid; tolerance is derived from its minimum spacing.
    void Initialize( Grid * grid );

    /// Initialize from one grid, then override the derived tolerance.
    void InitializeSpecial( Grid * grid, Real toleranceIn );

    /// Initialize from owned grid handles using the standard multi-grid policy.
    void Initialize( Grids & grids );

    /// Initialize from non-owning grid views using the standard multi-grid policy.
    void Initialize( const GridViews & grids );
public:
    // Virtual hooks used by mid-point / interface geometry paths. Production
    // builds the ADT; unit tests may override without grids or coorTree.
    /// Return the number of points currently stored.
    virtual int GetNPoint() const { return static_cast<int> (xCoor.size()); }

    /**
     * Find the closest stored point strictly inside the tolerance radius.
     * Returns INVALID_INDEX when no point matches.
     */
    virtual int FindPoint( Real xm, Real ym, Real zm );

    /**
     * Return an existing point ID within tolerance or insert a new point.
     * Insertion updates the coordinate store and ADT.
     */
    virtual int AddPoint( Real xm, Real ym, Real zm );

    /// Copy the coordinates of a stored point to the output arguments.
    virtual void GetPoint( int id, Real & xm, Real & ym, Real & zm ) const;

    /// Return the active point-matching tolerance.
    virtual Real GetTol() const { return tolerance; }
protected:
    int AddPoint( const RealField & coordinate );
    int FindPoint( const RealField & coordinate );
public:
    void GetFaceCoorList( const IntField & nodeId, RealField &xList, RealField &yList, RealField &zList ) const;
};

/**
 * Build an ADT from bounds without modifying the caller's bounds.
 * The supplied tolerance controls how far the copied bounds are expanded.
 */
void CreateStandardADT( const RealField & ptmin, const RealField & ptmax, std::unique_ptr<AdtTree>& adtTree, Real & tolerance );
void CreateStandardADT( Grid * grid, std::unique_ptr<AdtTree>& adtTree, Real & tolerance );
void CreateStandardADT( Grids & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance );
void CreateStandardADT( const GridViews & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance );
void CreateStandardADTByTolerance( Grids & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance );
void CreateStandardADTByTolerance( const GridViews & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance );

void ShiftMinMaxBox( RealField & pmin, RealField & pmax, Real tolerance );
void GetGridsMinMaxDistance( Grids & grids, Real & mindis, Real & maxdis );
void GetGridsMinMaxDistance( const GridViews & grids, Real & mindis, Real & maxdis );
Real CalcGridTolerance( Grids & grids );
Real CalcGridTolerance( const GridViews & grids );
void GetBoundingBoxOfMultiZoneGrids( Grids & grids, RealField & pmin, RealField & pmax );
void GetBoundingBoxOfMultiZoneGrids( const GridViews & grids, RealField & pmin, RealField & pmax );

EndNameSpace
