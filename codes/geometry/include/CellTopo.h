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
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class FaceTopo;
class UnsGrid;

class CellTopo
{
public:
    CellTopo();
    ~CellTopo();
public:
    IntField cellType;
    IntField blank;
    LinkField cellToNode;
    LinkField c2f;
    LinkField c2c;
public:
    void PushElement( int p1, int p2, int p3, int elementType );
    void PushElement( int p1, int p2, int p3, int p4, int elementType );
public:
    void Alloc( int nCell );
    UInt GetNumberOfCells() { return cellType.size(); }
    void CalcC2f( FaceTopo * faceTopo );
    void CalcC2C( FaceTopo * faceTopo );
};

void CalcC2f( UnsGrid * grid );
void CalcC2C( UnsGrid * grid );

EndNameSpace