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
#include "TurbUpdate.h"
BeginNameSpace( ONEFLOW )

class UTurbUpdate : public TurbUpdate
{
public:
    UTurbUpdate ();
    ~UTurbUpdate();
public:
    void UpdateFlowField( int sTid );
    void UpdateFlowField1Equ();
    void UpdateFlowField2Equ();
    void UpdateFlowField1EquStd();
    void UpdateFlowField2EquStd();

    void PrepareData1Equ();
    void PrepareData2Equ();
    void UpdateValue();
public:
    void DumpProbeInfo();
    void SmoothTurbulentPoint();
};

EndNameSpace