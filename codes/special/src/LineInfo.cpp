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

#include "LineInfo.h"

BeginNameSpace( ONEFLOW )

LineInfo::LineInfo()
{
    this->p1 = 0;
    this->p2 = 0;
    this->id = 1;
    this->type = LINE;
}

LineInfo::LineInfo( int p1, int p2, int id )
{
    this->p1 = p1;
    this->p2 = p2;
    this->id = id;
    this->type = LINE;
}

LineInfo::~LineInfo()
{
    ;
}

EndNameSpace