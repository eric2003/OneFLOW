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
#include "HXSort.h"

BeginNameSpace( ONEFLOW )

template < typename T >
void Reorder( T & field, IntField & order_map )
{
    T t_swap = field;
    int nElem = field.size();
    for ( int iElem = 0; iElem < nElem; ++ iElem )
    {
        int id = order_map[ iElem ];
        field[ iElem ] = t_swap[ id ];
    }
}

class PointDebug
{
public:
    PointDebug() {}
    ~PointDebug() {}
public:
    RealField xArray, yArray, zArray;
    Real tolerance;
public:
    void AddPoint( Real x, Real y, Real z );
    void GetPoint( int index, Real & x, Real & y, Real & z );
};

class WallVisual
{
public:
    WallVisual();
    ~WallVisual();
public:
    IntField lCell;
    IntField rCell;

    IntField lPos;
    IntField rPos;

    LinkField fLink;
    LinkField eLink;

    IntField faceType;
    IntField elementType;
    set< HXSort< IntField > > * faceSet;
public:
    RealField xN, yN, zN;
public:
    void CalcOrderMap( int & nBFace, IntField & orderMapping );
    void ConstructTopology();
    void ConstructTopology2D();
    void ConstructTopology3D();
    void BuildFaceTopo( IntField & faceNodeIndexArray, int loc_Face, int iCell, int face_type );
    void PushElement( int p1, int p2, int elementType );
    void PushElement( int p1, int p2, int p3, int elementType );
    void PushElement( int p1, int p2, int p3, int p4, int elementType );
public:
    void Visual( fstream & file, StringField & titleOfTecplot, RealField2D & qNodeField );
    void Visual3D( fstream & file, StringField & titleOfTecplot, RealField2D & qNodeField );
    void VisualLine( fstream & file, StringField & titleOfTecplot, RealField2D & qNodeField );
};

EndNameSpace