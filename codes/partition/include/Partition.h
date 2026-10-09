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
#include <memory>
#include <optional>
#include "HXDefine.h"
#include "GridTypes.h"
#include "GridHandles.h"
#include "HXCgns.h"
#include <vector>
#include <string>

#ifdef ENABLE_METIS
#include "metis.h"
#endif




BeginNameSpace( ONEFLOW )

class Grid;
class UnsGrid;
class G2LMapping;

class L2GMapping
{
public:
    IntField l2g_node;
    IntField l2g_face;
    IntField l2g_cell;
public:
    void Alloc( UnsGrid & grid );
    void CalcL2G    ( UnsGrid & ggrid, int zid, UnsGrid & grid, G2LMapping & g2l );
    void CalcL2GNode( UnsGrid & ggrid, G2LMapping & g2l );
    void CalcL2GFace( UnsGrid & ggrid, G2LMapping & g2l );
    void CalcL2GCell( UnsGrid & ggrid, int zid, G2LMapping & g2l );
};

class G2LMapping
{
public:
    G2LMapping( UnsGrid & ggrid, int npartproc );
public:
    IntField g2l_node;
    IntField g2l_face;
    IntField g2l_cell;
    std::vector<idx_t> gc2lzone;
    const int npartproc;
public:
    void GenerateGC2Z( UnsGrid & ggrid );
#ifdef ENABLE_METIS
    void GetXadjAdjncy( UnsGrid & ggrid, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy );
    void PartByMetis( idx_t nCells, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy );
#endif
};

class Partition
{
private:
    std::string sourceFile;
    int partitionType;
    std::optional< G2LMapping > g2l;
    G2LMapping & GetG2LMapping();
public:
    explicit Partition( const GridConfig & config );
    ~Partition();
public:
    Grids grids;
public:
    int npartproc;
    L2GMapping l2g;
public:
    void Run();
    UnsGrid & ReadGrid( const std::string & sourceFile );
    void GenerateMultiZoneGrid( UnsGrid & ggrid );
    void CreatePart( UnsGrid & ggrid );
    void AllocPart();
    void BuildCalculationalGrid( UnsGrid & ggrid );
    void BuildCalculationalGrid( UnsGrid & ggrid, int zid );
public:
    void CalcG2lCell( UnsGrid & ggrid );
    void CalcG2lFace( UnsGrid & ggrid, int zid, UnsGrid & grid );
    void CalcG2lNode( UnsGrid & ggrid, UnsGrid & grid );
    int GetNCell( UnsGrid & ggrid, int zid );
    void CreateL2g( UnsGrid & ggrid, int zid, UnsGrid & grid );
    void SetCoor( UnsGrid & ggrid, UnsGrid & grid );
    void SetGeometricRelationship( UnsGrid & ggrid, int zid, UnsGrid & grid );
    void CalcF2N( UnsGrid & ggrid, UnsGrid & grid );
    void SetF2CAndBC( UnsGrid & ggrid, int zid, UnsGrid & grid );
    void SetInterface( UnsGrid & ggrid, int zid, UnsGrid & grid, int partitionType );
};

class FacePairBasic
{
public:
    FacePairBasic() {};
    ~FacePairBasic() {};
public:
    int zone_id;
    int face_id;
    int cell_id;
};

class FacePair
{
public:
    FacePair() {};
    ~FacePair() {};
public:
    FacePairBasic lf, rf;
};

bool FindMatch( UnsGrid & grid, FacePair & facePair );

EndNameSpace
