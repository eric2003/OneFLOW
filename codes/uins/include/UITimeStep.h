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
#include "ITimestep.h"

BeginNameSpace( ONEFLOW )
class UITimestep : public ITimestep
{
public:
    UITimestep();
    ~UITimestep();
public:
    void Init();
    void ReadTmp();
    void CmpTimestep();
    void CmpLocalTimestep();
    void CmpGlobalTimestep();
    void CmpLgTimestep();
    void CmpInvTimestep();
    void CmpVisTimestep();
    void CmpMinTimestep();
    void SetTimestep( Real timestep );
public:
    void CmpSpectrumField();
    void CmpInvSpectrumField();
    void CmpVisSpectrumField();
public:
    void SetId( int fId );
    void PrepareData();
    void PrepareVisData();
    void UpdateInvSpectrumField();
    void UpdateVisSpectrumField();
    void ModifyTimestep();
};


EndNameSpace