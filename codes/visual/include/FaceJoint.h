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
#include "Point.h"
#include <memory> // Added for std::unique_ptr

BeginNameSpace( ONEFLOW )

class FaceJointManager;
class PointLocator;
class FaceJoint;

class FaceJointManager
{
public:
    FaceJointManager();
    ~FaceJointManager();

    // FIX: Disable copy to prevent double-free of unique_ptrs
    FaceJointManager(const FaceJointManager&) = delete;
    FaceJointManager& operator=(const FaceJointManager&) = delete;

public:
    // FIX: Use std::unique_ptr for automatic memory management
    HXVector< std::unique_ptr< FaceJoint > > patch;
    std::unique_ptr< FaceJoint > global;

public:
    void ConstructPointIndex();
    void CalcNodeValue();
};

class WallVisual;

class FaceJoint
{
public:
    using PointType = Point< Real >;
    using PointField = HXVector< PointType >;
    using PointLink = HXVector< PointField >;

public:
    FaceJoint();
    ~FaceJoint();

    // FIX: Disable copy to prevent double-free of unique_ptrs
    FaceJoint(const FaceJoint&) = delete;
    FaceJoint& operator=(const FaceJoint&) = delete;

public:
    bool isValid;
    IntField l2g;

public:
    PointLink fvp; 
    LinkField fLink; 
    IntField  weightId;
    RealField fcv; 
    RealField fnv; 

public:
    RealField pmin, pmax;
    Real dismin, dismax;

    // FIX: Use std::unique_ptr for automatic memory management
    std::unique_ptr< PointLocator > ps;
    std::unique_ptr< WallVisual > wallVisual;

public:
    void CalcBoundBox();
    void ConstructPointIndex();
    void ConstructPointIndexMap( FaceJoint * globalBasicWall );
    void CalcNodeValue();
    void RemapNodeValue( FaceJoint * globalBasicWall );

public:
    int GetSize() { return fvp.size(); }

public:
    void AddFacePoint( int nSolidCells, FaceJoint::PointLink & ptLink );
    void AddFaceCenterValue( int nSolidCells, RealField & fcvIn );
    void Visual( std::fstream & file );
};

EndNameSpace
