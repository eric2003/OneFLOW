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

#include "MultigridTaskReg.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterMultigridTask )

void RegisterMultigridTask()
{
    REGISTER_DATA_CLASS( RestrictAllQ );
    REGISTER_DATA_CLASS( RestrictDefect );
    REGISTER_DATA_CLASS( ModifyFineGrid );
    REGISTER_DATA_CLASS( ModifyCoarseGrid );
    REGISTER_DATA_CLASS( RecoverCoarseGrid );
    REGISTER_DATA_CLASS( RecoverResidual );
}

void RestrictAllQ( StringField & data )
{
    ;
}

void RestrictDefect( StringField & data )
{
    ;
}

void ModifyFineGrid( StringField & data )
{
    ;
}

void ModifyCoarseGrid( StringField & data )
{
    ;
}

void RecoverCoarseGrid( StringField & data )
{
    ;
}

void RecoverResidual( StringField & data )
{
    ;
}

EndNameSpace
