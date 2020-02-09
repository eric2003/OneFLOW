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

#include "LineMeshImp.h"
#include "LineMesh.h"
#include "CircleLineMesh.h"
#include "LineInfo.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

CurveMesh * CreateLineMesh( int lineType )
{
    CurveMesh * lineMesh = 0;

    if ( lineType == 0 )
    {
        lineMesh = new LineMesh();
    }
    else if ( lineType == 1 )
    {
        lineMesh = new CircleLineMesh();
    }

    return lineMesh;
}

CurveMesh * CreateLineMesh( CurveInfo * curveInfo )
{
    CurveMesh * lineMesh = 0;

    if ( curveInfo->type == 0 )
    {
        lineMesh = new LineMesh();
    }
    else if ( curveInfo->type == 1 )
    {
        lineMesh = new CircleLineMesh();
    }

    lineMesh->curveInfo = curveInfo;
    return lineMesh;

}

EndNameSpace