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
#include "HXDefine.h"
#include <string>

BeginNameSpace( ONEFLOW )

// Orchestrates field definition, interface registration, runtime allocation,
// and constant initialization for one solverType.
// Config root: system/<basicString>/alloc/
//   inner.txt / face.txt / bc.txt / unsteady.txt  -> FieldManager definitions
//   inter.txt / interDq.txt / interGrad.txt / interOverset.txt
//                                                 -> communication field names
//   init.txt                                      -> constant initial values
// Definition is solverType-scoped; allocation is per current Grid/DataStorage.
class FieldAllocator
{
public:
    static void Allocate(
        int solverType,
        const std::string & basicString );
};

EndNameSpace
