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
#include "BcSolver.h"
#include "UBcSolver.h"
#include "NsBcSolver.h"
#include <vector>

BeginNameSpace( ONEFLOW )

class UNsBcSolver : public NsBcSolver
{
public:
    UNsBcSolver();
    ~UNsBcSolver();
public:
    void Init();
    void CalcBcRegion();
    void CalcBc();
    void SetId( int bcfId );
    void PrepareData();
    void UpdateBc();
private:
    void CacheFieldPointers();

    std::vector< Real * > qData;
    std::vector< Real * > temprData;
    std::vector< Real * > bcQData;
    Real * gamaData = nullptr;

    const int * leftCellData = nullptr;
    const int * rightCellData = nullptr;
    const Real * xfnData = nullptr;
    const Real * yfnData = nullptr;
    const Real * zfnData = nullptr;
    const Real * vfxData = nullptr;
    const Real * vfyData = nullptr;
    const Real * vfzData = nullptr;
    const Real * vfnData = nullptr;
    const Real * faceAreaData = nullptr;
    const Real * faceCenterXData = nullptr;
    const Real * faceCenterYData = nullptr;
    const Real * faceCenterZData = nullptr;
    const Real * cellCenterXData = nullptr;
    const Real * cellCenterYData = nullptr;
    const Real * cellCenterZData = nullptr;
};

EndNameSpace
