/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "Grid.h"

BeginNameSpace( ONEFLOW )

class DataBook;
class FaceTopo;
class FaceMesh;
class CellMesh;
class IFaceLink;
class VirtualFile;

class UnsGrid : public Grid
{
public:
    IMPLEMENT_GRID_CLONE( UnsGrid )
public:
    UnsGrid();
    ~UnsGrid() override;
public:
    void Init() override;
    FaceTopo * faceTopo;
    FaceMesh * faceMesh;
    CellMesh * cellMesh;
public:
    void Decode( DataBook * databook ) override;
    void Encode( DataBook * databook ) override;
public:
    void ReadGrid( DataBook * databook ) override;
    void WriteGrid( DataBook * databook ) override;
public:
    void ReadGridFaceTopology( DataBook * databook );
    void ReadBoundaryTopology( DataBook * databook );
    void WriteGridFaceTopology( DataBook * databook );
    void WriteBoundaryTopology( DataBook * databook );
    void WriteBoundaryTopology1D( DataBook * databook );
    void WriteGridFaceTopology1D( DataBook * databook );
public:
    void ModifyBcType( int bcType1, int bcType2 ) override;
    void GenerateLgMapping( IFaceLink * iFaceLink ) override;
    void ReGenerateLgMapping( IFaceLink * iFaceLink ) override;
    void UpdateOtherTopologyTerm( IFaceLink * iFaceLink ) override;
    void NormalizeBc();
public:
    void GetMinMaxDistance( Real & dismin, Real & dismax ) override;
    void WriteGrid( std::fstream & file );
public:
    void CalcMetrics() override;
    void CalcMetrics1D();
    void CalcMetrics2D();
    void CalcMetrics3D();
    void AllocMetrics();
private:
    void CalcFaceCenter1D();
    void CalcCellCenterVol1D();
    void CalcFaceNormal1D();
    void CalcGhostCellCenterVol1D();
private:
    void CalcFaceNormal2D();
    void CalcFaceCenter2D();
    void CalcCellCenterVol2D();
private:
    void CalcFaceNormal3D();
    void CalcFaceCenter3D();
    void CalcCellCenterVol3D();
};

UnsGrid * UnsGridCast( Grid * gridIn );

EndNameSpace
