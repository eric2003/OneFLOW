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
#include "TurbSrcFlux.h"

BeginNameSpace( ONEFLOW )

class UTurbSrcFlux : public TurbSrcFlux
{
public:
    UTurbSrcFlux();
    ~UTurbSrcFlux();
public:
    void Init();
    void InitVist();
    void CalcSrcFlux();
    void CalcSrcFlux1Equ();
    void CalcSrcFlux2Equ();

    void Update1Equ();
    void Update2Equ();
    void CalcGrad();
    void CalcGradDebug();
    void CalcVist();
    void CalcVist1Equ();
    void CalcVist2Equ();
    void CalcVistMax();
    void SetGhostCellVist();
    void ZeroSpectrum();
    void ReadTmp();
public:
    void PrepareCellValue();
    void PrepareCellValue1Equ();
    void CalcLengthScaleSa();
    void CalcLengthScaleSst();

    void CalcLengthScaleOfSaDes();
    void CalcLengthScaleOfSaDdes();
    void CalcLengthScaleOfSaIddes();

    void CalcLengthScaleOfSstDes();
    void CalcLengthScaleOfSstDdes();
    void CalcLengthScaleOfSstIddes();

    void CalcLengthScaleOfWallDist();
public:
    void CalcBlendingTerm();
    void CalcCrossingTerm();
    void CalcCdkwmax();
    void CalcCrossDiff();
    void CalcBlendField();
    void ModifyBlendingTerm();
};

void CalcLengthLesOfSa( RealField & lesLengthField );
void CalcLengthLesOfSst( RealField & lesLengthField );

RealField GetLengthScale();
RealField & GetLargestSpacing();
void CalcSubgridLengthScale( RealField & lenth_scale );
void CalcLowReynoldsCorrection( RealField & lowReynoldsCorr );

EndNameSpace