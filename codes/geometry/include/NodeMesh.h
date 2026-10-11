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

BeginNameSpace( ONEFLOW )
class NodeMesh
{
public:
    struct Point3
    {
        Real x;
        Real y;
        Real z;
    };

public:
    NodeMesh();
public:
    RealField xN, yN, zN;
    RealField pmin, pmax;
    bool boxFlag;
public:
    HXSize_t GetNumberOfNodes() const { return xN.size(); }

    // Keep SoA as the owning storage while offering a value-type point interface.
    Point3 GetPoint( HXSize_t nodeIndex ) const;
    void SetPoint( HXSize_t nodeIndex, Real xp, Real yp, Real zp );

    void CreateNodes( int numberOfNodes );
    void SetCoordinates( RealField x, RealField y, RealField z );
    void CalcMinMaxBox();
    void AddPoint( Real xp, Real yp, Real zp );

    // Geometry mutations live with the coordinate owner so cached bounds stay valid.
    void ScaleAndTranslate( Real scale, Real tx, Real ty, Real tz );
    void TurnZAxisToYAxis();
};


EndNameSpace
