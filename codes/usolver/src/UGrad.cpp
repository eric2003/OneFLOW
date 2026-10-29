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

#include "UGrad.h"
#include "UGradPolicy.h"
#include "UCom.h"
#include "FieldManager.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "HXMath.h"
#include "Zone.h"
#include <iostream>
#include <iomanip>


BeginNameSpace( ONEFLOW )

void CalcGrad( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        const Real value = ArithmeticFaceAverage( q[ ug.lc ], q[ ug.rc ] );

        Real fnxa, fnya, fnza;
        AreaWeightedFaceNormal(
            ( * ug.xfn )[ ug.fId ], ( * ug.yfn )[ ug.fId ], ( * ug.zfn )[ ug.fId ],
            ( * ug.farea )[ ug.fId ], fnxa, fnya, fnza );

        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, 1.0,
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ] );

        if ( ug.fId < ug.nBFaces ) continue;
        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, -1.0,
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ScaleGradientByVolume(
            ( * ug.cvol )[ cId ], dqdx[ cId ], dqdy[ cId ], dqdz[ cId ] );
    }

    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        CopyGradient(
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ],
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }
}

void CalcGradGGCellWeight( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        const Real dxl = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
        const Real dyl = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
        const Real dzl = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

        const Real dxr = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
        const Real dyr = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
        const Real dzr = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

        Real cl, cr;
        InverseDistanceFaceWeights( dxl, dyl, dzl, dxr, dyr, dzr, cl, cr );
        const Real value = WeightedFaceAverage( q[ ug.lc ], q[ ug.rc ], cl, cr );

        Real fnxa, fnya, fnza;
        AreaWeightedFaceNormal(
            ( * ug.xfn )[ ug.fId ], ( * ug.yfn )[ ug.fId ], ( * ug.zfn )[ ug.fId ],
            ( * ug.farea )[ ug.fId ], fnxa, fnya, fnza );

        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, 1.0,
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ] );

        if ( ug.fId < ug.nBFaces ) continue;
        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, -1.0,
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ScaleGradientByVolume(
            ( * ug.cvol )[ cId ], dqdx[ cId ], dqdy[ cId ], dqdz[ cId ] );
    }

    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        CopyGradient(
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ],
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }
}


void CalcGradGGCellWeightDebug( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0.0;
    dqdy = 0.0;
    dqdz = 0.0;

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        const Real dxl = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
        const Real dyl = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
        const Real dzl = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

        const Real dxr = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
        const Real dyr = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
        const Real dzr = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

        Real cl, cr;
        InverseDistanceFaceWeights( dxl, dyl, dzl, dxr, dyr, dzr, cl, cr );
        const Real value = WeightedFaceAverage( q[ ug.lc ], q[ ug.rc ], cl, cr );

        Real fnxa, fnya, fnza;
        AreaWeightedFaceNormal(
            ( * ug.xfn )[ ug.fId ], ( * ug.yfn )[ ug.fId ], ( * ug.zfn )[ ug.fId ],
            ( * ug.farea )[ ug.fId ], fnxa, fnya, fnza );

        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, 1.0,
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ] );

        if ( ug.fId < ug.nBFaces ) continue;
        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, -1.0,
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ScaleGradientByVolume(
            ( * ug.cvol )[ cId ], dqdx[ cId ], dqdy[ cId ], dqdz[ cId ] );
    }

    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        CopyGradient(
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ],
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }

}

void CalcGradDebug( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        const Real value = ArithmeticFaceAverage( q[ ug.lc ], q[ ug.rc ] );

        Real fnxa, fnya, fnza;
        AreaWeightedFaceNormal(
            ( * ug.xfn )[ ug.fId ], ( * ug.yfn )[ ug.fId ], ( * ug.zfn )[ ug.fId ],
            ( * ug.farea )[ ug.fId ], fnxa, fnya, fnza );

        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, 1.0,
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ] );

        if ( ug.lc == 0 || ug.rc == 0 )
        {
            std::cout << " fId = " << fId << " value = " << value << " ql= " << q[ ug.lc ] << " qr = " << q[ ug.rc ] << std::endl;
        }

        if ( ug.fId < ug.nBFaces ) continue;
        AccumulateFaceGradient(
            fnxa, fnya, fnza, value, -1.0,
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ScaleGradientByVolume(
            ( * ug.cvol )[ cId ], dqdx[ cId ], dqdy[ cId ], dqdz[ cId ] );
    }

    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        CopyGradient(
            dqdx[ ug.lc ], dqdy[ ug.lc ], dqdz[ ug.lc ],
            dqdx[ ug.rc ], dqdy[ ug.rc ], dqdz[ ug.rc ] );
    }
}

EndNameSpace
