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
#include "TurbInvFlux.h"

BeginNameSpace( ONEFLOW )

class Limiter;
class LimField;


class UTurbInvFlux : public TurbInvFlux
{
public:
    UTurbInvFlux ();
    ~UTurbInvFlux();
public:
    Limiter * limiter, * nslimiter;
    LimField * limf;
    MRField * invflux;
public:
    void Alloc();
    void DeAlloc();
    void CmpFlux();
    void CmpInvFlux();
    void CmpInvFace();
    void CmpLimiter();
    void AddInvFlux();
    void PrepareFaceValue();
    void UpdateFaceInvFlux();
public:
    void GetQlQrField();
    void ReconstructFaceValueField();
    void BoundaryQlQrFixField();
};


EndNameSpace