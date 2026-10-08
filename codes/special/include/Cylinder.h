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
#include "Point.h"
#include "GridTypes.h"
#include "CurveLine.h"
#include <string>

BeginNameSpace( ONEFLOW )

using PointType = Point< Real >;

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
    void Symmetry( const DomainData & datain );
    void Join( const DomainData & d1, const DomainData & d2 );
};


class Cylinder
{
public:
    Cylinder();
    ~Cylinder();
    Cylinder( const Cylinder & ) = delete;
    Cylinder & operator = ( const Cylinder & ) = delete;
    Cylinder( Cylinder && ) = delete;
    Cylinder & operator = ( Cylinder && ) = delete;
public:
    DomainData domain_data;
    DomainData symm_domain;
    DomainData final_domain;
    int nZone;

    StrCurveLoop strCurveLoop;
public:
    Real beta;
public:
    void Run( const GridConfig & config, const std::string & caseDir );
    void HalfCylinder( const GridConfig & config, const std::string & caseDir );
    void QuarterCylinder( const GridConfig & config, const std::string & caseDir );
    void GenePlate();
public:
    void SetBoundaryGrid( const GridConfig & config, const std::string & caseDir );
    void GeneDomain();
public:
    void CalcCircleCenter( PointType & p1, PointType & p2, PointType & p0, PointType & pcenter );
public:
    void DumpGrid( const std::string & fileName, const std::string & caseDir, const DomainData & domain );
    void DumpBcFile( const std::string & fileName, const std::string & caseDir, const DomainData & domain, const IntField & bcList );
    void ToTecplot( const std::string & fileName, const std::string & caseDir, const DomainData & domain );
};

void ToTecplot( std::fstream & file, const RealField2D & coor, int ni, int nj, int nk );
void DumpBc( std::fstream &file, int imin, int imax, int jmin, int jmax, int bcType );

EndNameSpace
