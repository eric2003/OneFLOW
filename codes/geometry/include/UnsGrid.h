/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
	Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
    ~UnsGrid();
public:
    void Init();
    FaceTopo * faceTopo;
	FaceMesh * faceMesh;
	CellMesh * cellMesh;
public:
    void Decode( DataBook * databook );
    void Encode( DataBook * databook );
public:
    void ReadGrid( DataBook * databook );
    void WriteGrid( DataBook * databook );
public:
    void ReadGridFaceTopology( DataBook * databook );
    void ReadBoundaryTopology( DataBook * databook );
    void WriteGridFaceTopology( DataBook * databook );
    void WriteBoundaryTopology( DataBook * databook );
public:
    void ModifyBcType( int bcType1, int bcType2 );
    void GenerateLgMapping( IFaceLink * iFaceLink );
    void ReGenerateLgMapping( IFaceLink * iFaceLink );
    void UpdateOtherTopologyTerm( IFaceLink * iFaceLink );
    void NormalizeBc();
public:
    void GetMinMaxDistance( Real & dismin, Real & dismax );
    void WriteGrid( fstream & file );
private:
    void WriteGridFaceTopology( VirtualFile * vf );
    void WriteBoundaryTopology( VirtualFile * vf );
public:
    void CmpMetrics();
    void CmpMetrics1D();
    void CmpMetrics2D();
    void CmpMetrics3D();
    void AllocMetrics();
private:
	void ComputeFaceCenter1D();
	void ComputeCellCenterVol1D();
	void ComputeFaceNormal1D();
    void ComputeGhostCellCenterVol1D();
private:
	void ComputeFaceNormal2D();
	void ComputeFaceCenter2D();
	void ComputeCellCenterVol2D();
private:
	void ComputeFaceNormal3D();
	void ComputeFaceCenter3D();
	void ComputeCellCenterVol3D();
};

UnsGrid * UnsGridCast( Grid * gridIn );

EndNameSpace