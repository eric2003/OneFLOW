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
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class TurbSrcFlux
{
public:
    TurbSrcFlux();
    ~TurbSrcFlux();
    typedef void ( TurbSrcFlux:: * SrcPointer )();
    SrcPointer srcFlux;
    SrcPointer cmpBeta;
    SrcPointer cmpProd;
public:
    void SetSrcFluxPointer();
    void CmpSrcSa();
    void CmpSrc2Equ();
    void CmpSrc2EquKwMenter();
    void CmpSrc2EquKwWilcox1998();
    void CmpSrc2EquKwWilcox2006();
    void CmpSrc2EquKwDefault();
    void CmpSrc2EquEasmKw2003();
    void CmpSrc2EquEasmKw2005();
    void CmpFbetaCoef();
    void CmpFbetaDefault();
    void CmpFbetaOfKwWilcox1998();
    void CmpFbetaOfKwWilcox2006();
    void CmpFbetaOfEasmKw2003();
public:
    void CmpVGrad();
    void CmpTransition();
    void CmpProdW();
    void CmpProdk();
    void CmpDissk();
    void LimitProdk();
    void CmpProdwKwMenter();
    void CmpProdwKwWilcox1998();
    void CmpProdwKwWilcox2006();
    void CmpProdwKwDefault();
    void CmpProdwEasmKw2003();
    void CmpProdwEasmKw2005();
    void ModifyPd();
    void CmpSrc();
public:
    void CmpCellVist1Equ();
    void CmpCellVist2Equ();
};


EndNameSpace