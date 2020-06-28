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
#include "Point.h"

BeginNameSpace( ONEFLOW )

class FaceJointManager;

class PointSearch;
class FaceJoint;

class FaceJointManager
{
public:
    FaceJointManager();
    ~FaceJointManager();
public:
    HXVector< FaceJoint * > patch;
    FaceJoint * global;
public:
    void ConstructPointIndex();
    void CalcNodeValue();
};

class WallVisual;

class FaceJoint
{
public:
    typedef Point< Real > PointType;
    typedef HXVector< PointType > PointField;
    typedef HXVector< PointField > PointLink;
public:
    FaceJoint();
    ~FaceJoint();
public:
    bool isValid;
    IntField l2g;
public:
    PointLink fvp; //face vertex point;
    LinkField fLink; //face link
    IntField  weightId;
    RealField fcv; //face center value
    RealField fnv; //face node value
public:
    RealField pmin, pmax;
    Real dismin, dismax;
    PointSearch * ps;
    WallVisual * wallVisual;
public:
    void CalcBoundBox();
    void ConstructPointIndex();
    void ConstructPointIndexMap( FaceJoint * globalBasicWall );
    void CalcNodeValue();
    void RemapNodeValue( FaceJoint * globalBasicWall );
public:
    int GetSize() { return fvp.size(); }
public:
    void AddFacePoint( int nSolidCell, FaceJoint::PointLink & ptLink );
    void AddFaceCenterValue( int nSolidCell, RealField & fcvIn );
    void Visual( fstream & file );
};

EndNameSpace