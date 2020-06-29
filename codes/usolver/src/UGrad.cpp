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

#include "UGrad.h"
#include "UCom.h"
#include "FieldImp.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "HXMath.h"
#include "Zone.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void CalcGrad( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        Real value = half * ( q[ ug.lc ] + q[ ug.rc ] );

        Real fnxa = ( * ug.xfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];
        Real fnya = ( * ug.yfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];
        Real fnza = ( * ug.zfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];

        dqdx[ ug.lc ] += fnxa * value;
        dqdy[ ug.lc ] += fnya * value;
        dqdz[ ug.lc ] += fnza * value;

        if ( ug.fId < ug.nBFace ) continue;
        dqdx[ ug.rc ] -= fnxa * value;
        dqdy[ ug.rc ] -= fnya * value;
        dqdz[ ug.rc ] -= fnza * value;
    }

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        Real ovol = one / ( * ug.cvol )[ cId ];
        dqdx[ cId ] *= ovol;
        dqdy[ cId ] *= ovol;
        dqdz[ cId ] *= ovol;
    }

    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        dqdx[ ug.rc ] = dqdx[ ug.lc ];
        dqdy[ ug.rc ] = dqdy[ ug.lc ];
        dqdz[ ug.rc ] = dqdz[ ug.lc ];
    }
}

void CalcGradGGCellWeight( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        if ( fId == 432 )
        {
            int kkk = 1;
        }
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        Real dxl = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
        Real dyl = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
        Real dzl = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

        Real dxr = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
        Real dyr = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
        Real dzr = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

        Real delt1  = DIST( dxl, dyl, dzl );
        Real delt2  = DIST( dxr, dyr, dzr );
        Real delta  = 1.0 / ( delt1 + delt2 + SMALL );

        Real cl = delt2 * delta;
        Real cr = delt1 * delta;

        Real value = cl * q[ ug.lc ] + cr * q[ ug.rc ];

        Real fnxa = ( * ug.xfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];
        Real fnya = ( * ug.yfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];
        Real fnza = ( * ug.zfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];

        dqdx[ ug.lc ] += fnxa * value;
        dqdy[ ug.lc ] += fnya * value;
        dqdz[ ug.lc ] += fnza * value;

        if ( ug.fId < ug.nBFace ) continue;
        dqdx[ ug.rc ] -= fnxa * value;
        dqdy[ ug.rc ] -= fnya * value;
        dqdz[ ug.rc ] -= fnza * value;
    }

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        Real ovol = one / ( * ug.cvol )[ cId ];
        dqdx[ cId ] *= ovol;
        dqdy[ cId ] *= ovol;
        dqdz[ cId ] *= ovol;
    }

    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        dqdx[ ug.rc ] = dqdx[ ug.lc ];
        dqdy[ ug.rc ] = dqdy[ ug.lc ];
        dqdz[ ug.rc ] = dqdz[ ug.lc ];
    }
}

void CalcGradDebug( RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz )
{
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( ug.lc == 0 || ug.rc == 0 )
        {
            int kkk = 1;
        }

        Real value = half * ( q[ ug.lc ] + q[ ug.rc ] );

        Real fnxa = ( * ug.xfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];
        Real fnya = ( * ug.yfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];
        Real fnza = ( * ug.zfn )[ ug.fId ] * ( * ug.farea )[ ug.fId ];

        dqdx[ ug.lc ] += fnxa * value;
        dqdy[ ug.lc ] += fnya * value;
        dqdz[ ug.lc ] += fnza * value;

        if ( ug.lc == 0 || ug.rc == 0 )
        {
            cout << " fId = " << fId << " value = " << value << " ql= " << q[ ug.lc ] << " qr = " << q[ ug.rc ] << endl;
        }

        if ( ug.fId < ug.nBFace ) continue;
        dqdx[ ug.rc ] -= fnxa * value;
        dqdy[ ug.rc ] -= fnya * value;
        dqdz[ ug.rc ] -= fnza * value;
    }

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        if ( cId == 0 )
        {
            int kkk = 1;
        }
        Real ovol = one / ( * ug.cvol )[ cId ];
        dqdx[ cId ] *= ovol;
        dqdy[ cId ] *= ovol;
        dqdz[ cId ] *= ovol;
    }

    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        dqdx[ ug.rc ] = dqdx[ ug.lc ];
        dqdy[ ug.rc ] = dqdy[ ug.lc ];
        dqdz[ ug.rc ] = dqdz[ ug.lc ];
    }
}

EndNameSpace