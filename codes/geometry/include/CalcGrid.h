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

using namespace std;

BeginNameSpace( ONEFLOW )

class Grid;
class IFaceLink;
class NodeMesh;

class CalcGrid
{
public:
    CalcGrid();
    ~CalcGrid();
public:
    Grids grids;
    string gridFileName;
    IFaceLink * iFaceLink;
public:
    void BuildInterfaceLink();
    void Dump();
    void Post();
public:
    void Init( Grids & grids );
public:
    void GenerateOverset();
    void GenerateLink();
    void ResetGridScaleAndTranslate();
public:
    void ReconstructLink();
public:
    void ModifyBcType();
    void GenerateLgMapping();
    void ReGenerateLgMapping();
    void UpdateLgMapping();
    void UpdateOtherTopologyTerm();
    void ReconstructInterFace();
    void MatchInterfaceTopology();
    void ReconstructLink( int iZone );
};

string GetTargetGridFileName();
int GetIgnoreNoBc();

void GenerateMultiZoneCalcGrids( Grids & grids );
void ResetGridScaleAndTranslate( NodeMesh * nodeMesh );
void TurnZAxisToYAxis( NodeMesh * nodeMesh );

EndNameSpace