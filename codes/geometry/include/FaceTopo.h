/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "Constant.h"
#include "HXDefine.h"
#include <vector>
#include <string>
#include <fstream>


BeginNameSpace( ONEFLOW )

class BcRecord;
class IFaceLink;
class InterFace;
class BcManager;
class Grid;

class FaceTopo
{
public:
    FaceTopo();
    ~FaceTopo();
public:
    int nCells;
    IntField fTypes;
    LinkField faces;
    LinkField c2f;

    IntField lCells, rCells;
    IntField lPosition, rPosition;

    HXSize_t nBFaces;
    BcManager * bcManager;
    Grid * grid;
public:
    LinkField facesNew;
    IntField lCellsNew, rCellsNew;
public:
    HXSize_t GetNFaces() { return fTypes.size();  }
    HXSize_t CalcTotalFaceNodes();
    HXSize_t GetNBFaces();
    void SetNBFaces( HXSize_t nBFaces );
public:
    void ModifyFaceNodeId( IFaceLink * iFaceLink );
    void SetNewFace2Node( IFaceLink * iFaceLink );
    void SetNewFace2Cell( IFaceLink * iFaceLink );
    void ModifyBoundaryInformation( IFaceLink * iFaceLink );
    void ResetNumberOfBoundaryCondition( IFaceLink * iFaceLink );
    void ConstructNewInterfaceMap( IFaceLink * iFaceLink );
    void UpdateOtherTopologyTerm();
    void GenerateI2B( InterFace * interFace );
public:
    bool GetSId( int iFace, int iPosition, int & sId );
    bool GetTId( int iFace, int iPosition, int & tId );
    void CalcC2C( LinkField & c2c );
};

EndNameSpace
