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
#include "HXCgns.h"
#include <vector>
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )
class FaceSolver;

class EState
{
public:
    EState();
    ~EState();
public:
    static int SPLIT   ;
    static int NOCHANGE;
    static int MERGE   ;
    static int DELETE  ;
    static int HIDDEN  ;
    static int ON      ;
    static int BC_FACE ;
    static int H_FACE  ;
    static int G_FACE  ;
};

class ElemFeature
{
public:
    ElemFeature();
    ~ElemFeature();
public:
    IntField * eType;       //��Ԫ����
    CgLinkField  eNodeId;     //��Ԫindex
public:
    FaceSolver * face_solver;
    void ScanElements();

};

EndNameSpace