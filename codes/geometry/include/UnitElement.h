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
#include "Constant.h"
#include "HXVector.h"
#include "HXDefine.h"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class UnitElement
{
public:
    UnitElement();
    ~UnitElement();
public:
    int elementType;
    IntField nodeId;
    LinkField middlePointStruct; //Midpoint data structure
    LinkField childElementIndex;
    IntField childElementType;
    IntField faceTypeContainer;
    LinkField faceList;

    IntField compositeFaceElementType;
    LinkField compositeFaceList;
public:
    void Initialize( int elementType );
public:
    IntField & GetElementPhysicsFace(int iFace);
    int GetElementPhysicsFaceType(int iFace) const;
    int GetElementFaceNumber() const;
    int GetElementPhysicsFaceNumber() const;
    int GetElementType() const;
    int GetElementNodeNumbers( int elementType );

    int GetRefinedElementType( int elementType );
    int GetSimpleElementType ( int elementType );

    bool IsUnitElementType  ( int elementType );
    bool IsFaceElementType   ( int elementType );

    bool IsFaceElementTypeAtLeast( int elementType );
    bool IsBasicVolumeElementType( int elementType );

    IntField & GetElementFace(int iFace);
    int GetFaceType(int iFace) const;

    IntField & GetRelatedPointListForMiddlePointCalcutation(int iMiddlePoint);

    IntField & GetChildElementRelativeNodeIndex(int iChildElement);

    int GetChildElementNumbers() const;
    int GetChildElementType(int iChild);
public:
    void PushElementFace( int faceType, int p1 );
    void PushElementFace( int faceType, int p1, int p2 );
    void PushElementFace( int faceType, int p1, int p2, int p3 );
    void PushElementFace( int faceType, int p1, int p2, int p3, int p4 );
public:
    void PushMiddlePoint( int pm, int p1, int p2 );
    void PushMiddlePoint( int pm, int p1, int p2, int p3 );
    void PushMiddlePoint( int pm, int p1, int p2, int p3, int p4 );
    void PushMiddlePoint( int pm, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8 );
public:
    void PushChildElement( int elementType, int p1, int p2 );
    void PushChildElement( int elementType, int p1, int p2, int p3 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5, int p6 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8 );
public:
    void PushCompositeFace( int faceElementType, int p1, int p2, int p3 );
    void PushCompositeFace( int faceElementType, int p1, int p2, int p3, int p4, int p5, int p6 );
    void PushCompositeFace( int faceElementType, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8, int p9 );
};

EndNameSpace