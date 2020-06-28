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
#include "Point.h"

BeginNameSpace( ONEFLOW )

class CurveLine;

typedef Point< Real > PointType;

class DomainData
{
public:
    DomainData();
    ~DomainData();
public:
    RealField2D x;
    RealField2D y;
    RealField2D z;
    int ni, nj;
public:
    void Alloc();
    void Symmetry( DomainData * datain );
    void Join( DomainData * d1, DomainData * d2 );
};

class StrCurveLoop;

class Cylinder
{
public:
    Cylinder();
    ~Cylinder();
public:
    DomainData domain_data;
    DomainData symm_domain;
    DomainData final_domain;
    int nZone;

    StrCurveLoop * strCurveLoop;
public:
    Real beta;
public:
    void Run( int igene );
    void HalfCylinder();
    void QuarterCylinder();
    void GenePlate();
public:
    void SetBoundaryGrid();
    void GeneDomain();
public:
    void CalcCircleCenter( PointType & p1, PointType & p2, PointType & p0, PointType & pcenter );
public:
    void DumpGrid( const string & fileName, DomainData * domain );
    void DumpBcFile( const string & fileName, DomainData * domain, IntField & bcList );
    void ToTecplot( const string & fileName, DomainData * domain );
};

void ToTecplot( fstream & file, RealField2D & coor, int ni, int nj, int nk );
void DumpBc( fstream &file, int imin, int imax, int jmin, int jmax, int bcType );

EndNameSpace