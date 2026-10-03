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
#include "CalcCoor.h"
#include "SimpleDomain.h"
#include <map>
#include <memory>


BeginNameSpace( ONEFLOW )

class SDomain;
class Block2D;
class CurveMesh;

class SLine
{
public:
    SLine();
    ~SLine();
public:
    int line_id;
    int ni;
    RealField x1d, y1d, z1d;
    IntField ctrlpoints;
public:
    void SetDomainBcMesh( SDomain * sDomain );
    void SetBlkBcMesh( Block2D * blk2d );
    void ConstructCtrlPoints( const IntField & pointIdList );
    void Alloc();
    void CopyMesh( const CurveMesh & curveMesh );
    void ConstructPointToLineMap( const IntField & pointIdList, std::map< int, IntSet > & pointToLineMap ) const;

};

class MLine : public DomData
{
public:
    explicit MLine( CoorMap * coorMap );
    ~MLine();
public:
    int pos;
    IntField lineList;
    HXVector< std::unique_ptr< SLine > > slineList;
    CoorMap * coorMap;
public:
    std::map< int, IntSet > pointToLine;
public:
    void ConstructLineToDomainMap();
    void ConstructLineToDomainMap( int domain_id, std::map< int, IntSet > & lineToDomainMap );
    void ConstructPointToDomainMap();
    void ConstructPointToDomainMap( int domain_id, std::map< int, IntSet > & pointToDomainMap ) const;
    void ConstructPointToPointMap();
    void ConstructPointToPointMap( std::map< int, IntSet > & pointToPointMap ) const;
    void ConstructPointToLineMap( std::map< int, IntSet > & pointToLineMap ) const;
public:
    void AddSubLine( int line_id );
    void ConstructDomainTopo();
    void ConstructCtrlPoint();
    void ConstructSLineCtrlPoint();
    //void CalcCoor( CoorMap * localCoorMap );
    void CalcCoor();
    void SetDomainBcMesh( SDomain * sDomain );
    void CreateInpFaceList1D( HXVector< std::unique_ptr<Face2D> > &facelist );
    void SetBlkBcMesh( Block2D * blk2d );
};


EndNameSpace
