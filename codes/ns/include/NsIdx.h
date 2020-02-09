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
#include "Configure.h"

BeginNameSpace( ONEFLOW )

class IDX
{
public:
    IDX(){};
    ~IDX(){};
public:
    static const int IR = 0;
    static const int IU = 1;
    static const int IV = 2;
    static const int IW = 3;
    static const int IP = 4;

    static const int IRU = 1;
    static const int IRV = 2;
    static const int IRW = 3;
    static const int IRE = 4;

    static const int ITT = 0;
};

EndNameSpace
