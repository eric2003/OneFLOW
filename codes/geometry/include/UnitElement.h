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
#include "Constant.h"
#include "HXVector.h"
#include "HXDefine.h"
#include <vector>
#include <string>
#include <fstream>


BeginNameSpace( ONEFLOW )

class UnitElement
{
private:
    // Local topology storage is maintained through the member API so related
    // arrays cannot be edited independently by callers.
    int elementType = -1;
    IntField nodeId;
    LinkField middlePointNodeIndices;
    LinkField childElementIndex;
    IntField childElementType;
    IntField faceTypeContainer;
    LinkField faceList;

    IntField compositeFaceElementType;
    LinkField compositeFaceList;

public:
    UnitElement() = default;
    ~UnitElement() = default;

    void Initialize( int elementType );

private:
    void BuildTopology( int elementType );

public:
    const IntField & GetCompositeFace( int iFace ) const;
    int GetCompositeFaceType( int iFace ) const;
    int GetElementFaceCount() const;
    int GetCompositeFaceCount() const;

    // Compatibility wrappers; prefer the Count-suffixed APIs for new code.
    int GetElementFaceNumber() const;
    int GetCompositeFaceNumber() const;
    int GetElementType() const;
    static int GetElementNodeNumbers( int elementType );

    int GetRefinedElementType( int elementType ) const;
    int GetSimpleElementType ( int elementType ) const;

    bool IsUnitElementType  ( int elementType ) const;
    bool IsFaceElementType   ( int elementType ) const;

    bool IsFaceElementTypeAtLeast( int elementType ) const;
    bool IsBasicVolumeElementType( int elementType ) const;

    const IntField & GetElementFace( int iFace ) const;
    int GetFaceType( int iFace ) const;

    const IntField & GetRelatedPointListForMiddlePointCalculation( int iMiddlePoint ) const;

    // Compatibility wrappers; use the correctly spelled API for new code.
    const IntField & GetRelatedPointListForMiddlePointCalcutation( int iMiddlePoint ) const;

    const IntField & GetChildElementRelativeNodeIndex( int iChildElement ) const;

    int GetChildElementCount() const;
    // Compatibility wrapper; prefer GetChildElementCount() for new code.
    int GetChildElementNumbers() const;
    int GetChildElementType( int iChild ) const;

private:
    void PushElementFace( int faceType, int p1 );
    void PushElementFace( int faceType, int p1, int p2 );
    void PushElementFace( int faceType, int p1, int p2, int p3 );
    void PushElementFace( int faceType, int p1, int p2, int p3, int p4 );

    void PushMiddlePoint( int pm, int p1, int p2 );
    void PushMiddlePoint( int pm, int p1, int p2, int p3 );
    void PushMiddlePoint( int pm, int p1, int p2, int p3, int p4 );
    void PushMiddlePoint( int pm, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8 );

    void PushChildElement( int elementType, int p1, int p2 );
    void PushChildElement( int elementType, int p1, int p2, int p3 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5, int p6 );
    void PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8 );

    void PushCompositeFace( int faceElementType, int p1, int p2, int p3 );
    void PushCompositeFace( int faceElementType, int p1, int p2, int p3, int p4, int p5, int p6 );
    void PushCompositeFace( int faceElementType, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8, int p9 );
};

EndNameSpace
