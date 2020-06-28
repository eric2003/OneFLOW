/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "GridDef.h"
#include "HXCgns.h"
#include <vector>
#include <string>
#include <fstream>

#ifdef ENABLE_METIS
#include "metis.h"
#endif


using namespace std;

BeginNameSpace( ONEFLOW )

class Grid;
class UnsGrid;
class G2LMapping;

class L2GMapping
{
public:
    L2GMapping();
    ~L2GMapping();
public:
    IntField l2g_node;
    IntField l2g_face;
    IntField l2g_cell;
public:
    G2LMapping * g2l;
public:
    void Alloc( UnsGrid * grid );
    void CalcL2G    ( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void CalcL2GNode( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void CalcL2GFace( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void CalcL2GCell( UnsGrid * ggrid, int zid, UnsGrid * grid );
};

class G2LMapping
{
public:
    G2LMapping( UnsGrid * ggrid );
    ~G2LMapping();
public:
    IntField g2l_node;
    IntField g2l_face;
    IntField g2l_cell;
    vector<idx_t> gc2lzone;
    UnsGrid * ggrid;
    int npartproc;
    LinkField c2c;
public:
    void GenerateGC2Z();
#ifdef ENABLE_METIS
    void GetXadjAdjncy( UnsGrid * ggrid, vector<idx_t>& xadj, vector<idx_t>& adjncy );
    void PartByMetis( idx_t nCell, vector<idx_t>& xadj, vector<idx_t>& adjncy );
#endif
    void DumpXadjAdjncy( UnsGrid * grid, IntField & xadj, IntField & adjncy );
    void DumpGC2Z( UnsGrid * grid );
    void ReadGC2Z( UnsGrid * grid );
};

class Partition
{
public:
    Partition();
    ~Partition();
public:
    Grids grids;
public:
    UnsGrid * uns_grid;
    int npartproc;
    int partition_type;
    int partition_c2n;
    G2LMapping * g2l;
    L2GMapping * l2g;
public:
    void Run();
    void ReadGrid();
    void GenerateMultiZoneGrid();
    void CreatePart();
    void AllocPart();
    void BuildCalcutationalGrid();
    void BuildCalcutationalGrid( int zid );
    void PreProcess();
    void PostProcess();
public:
    void CalcGC2N();
    void CalcG2lCell();
    void CalcG2lFace( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void CalcG2lNode( UnsGrid * ggrid, int zid, UnsGrid * grid );
    int GetNCell( UnsGrid * ggrid, int zid );
    void CreateL2g( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void SetCoor( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void SetGeometricRelationship( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void CalcF2N( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void SetF2CAndBC( UnsGrid * ggrid, int zid, UnsGrid * grid );
    void SetInterface( UnsGrid * ggrid, int zid, UnsGrid * grid );
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

bool FindMatch( UnsGrid * grid, FacePair * facePair );

EndNameSpace