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
#include "HXArray.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class BcRecord;

const int F_INNER = 0;
const int F_GHOST = 1;
const int F_TOTAL = 2;

class UnsGrid;
class FaceTopo;
class CellTopo;
class UGeom
{
public:
    UGeom() = default;
    ~UGeom() = default;
public:
    // Bind current zone geometry pointers into the global UGeom (ug).
    void Init();
    // Ensure BcRecord::bcInfo is built for the current zone.
    void EnsureBcTypeRegion();
    // Set face iteration range: inner / boundary / total.
    void SetStEd( int flag );
    // Debug dump of faces attached to a cell.
    void DumpCellFace( int cId );
public:
    UnsGrid * grid;
    int fId;
    int cId;
    int bcfId;
    int nFaces;
    int nBFaces;
    int nCells;
    int nTCell;
    int ist, ied;
    int lc;
    int rc;
    int ir, bcNameId, bctype, nRegion, nRBFace;
    int ireconface;
public:
    IntField * lcf;
    IntField * rcf;
    IntField * blankf;
    LinkField * c2f;

    RealField * xfn;
    RealField * yfn;
    RealField * zfn;
    RealField * vfn;
    RealField * farea;

    RealField * xfc;
    RealField * yfc;
    RealField * zfc;

    RealField * vfx;
    RealField * vfy;
    RealField * vfz;

    RealField * cvol, * cvol1, * cvol2;
    RealField * xcc;
    RealField * ycc;
    RealField * zcc;
    BcRecord * bcRecord;

private:
    // Reinitialization replaces these topology objects. When their identity
    // is unchanged, repeated stage setup only needs to refresh the range.
    const UnsGrid * boundGrid_ = nullptr;
    const FaceTopo * boundFaceTopo_ = nullptr;
    const CellTopo * boundCellTopo_ = nullptr;
    int boundFaces_ = -1;
    int boundCells_ = -1;
    int boundBoundaryFaces_ = -1;
};

extern UGeom ug;

void AddF2CField( MRField * cellField, MRField * faceField );
void AddF2CFieldDebug( MRField * cellField, MRField * faceField );

class HXDebug
{
public:
    HXDebug() = default;
    ~HXDebug() = default;
    static std::string fileName1, fileName2;
public:
    static void DumpResField( const std::string & fileName );
    static void DumpField( const std::string & fileName, const MRField * field );
    static std::string GetFullFileName( const std::string & fileName, int startStrategy );
    static void CompareFile( Real mindiff, int idump );
    static std::unique_ptr<MRField> ReadField( const std::string & fileName );
    static void DumpCellInfo( int iCell );
    static void CheckNANField( const MRField * field );
};

EndNameSpace
