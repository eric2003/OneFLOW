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

#include "CgnsBc1to1.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsPeriod.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

int AbsoluteDiagonalId( int x, int y )
{
    if ( ABS( x ) == ABS( y ) )
    {
        return 1;
    }
    return 0;
}

CgnsBc1to1::CgnsBc1to1( CgnsZone * cgnsZone )
    : CgnsBcLink( cgnsZone )
{
}

CgnsBc1to1::~CgnsBc1to1()
{
}

void CgnsBc1to1::ReadCgnsBc1To1()
{
    int fileId = this->cgnsZone->cgnsBase->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zId = this->cgnsZone->zId;

    this->nConnPoints = 6;
    this->nConnDonorPoints = 6;

    this->connPoint.resize( nConnPoints );
    this->connDonorPoint.resize( nConnDonorPoints );

    this->donorPointSetType = CGNS_ENUMV( PointRange );
    this->donorDataType     = CGNS_ENUMV( Integer );

    CgnsTraits::char33 connName;
    CgnsTraits::char33 donorZoneName;


    cg_1to1_read( fileId, baseId, zId, this->bcId, connName, donorZoneName, & this->connPoint[ 0 ], & this->connDonorPoint[ 0 ], itranfrm );

    this->connName = connName;
    this->donorZoneName  = donorZoneName;

    cout << "\n";
    cout << "   connName      = " << connName << " donorZoneName = " << donorZoneName << "\n";
    cout << "   donorDataType = " << DataTypeName[ this->donorDataType ] << "\n";
    cout << "   donorPointSetType = " << PointSetTypeName[ this->donorPointSetType ] << "\n";
    cout << "   nConnPoints      = " << nConnPoints << "\n";
    cout << "   nConnDonorPoints = " << nConnDonorPoints << "\n";

    cout << "   range (this zone )= ";
    int width = 5;
    cout << setw( width ) << connPoint[ 0 ];
    cout << setw( width ) << connPoint[ 1 ];
    cout << setw( width ) << connPoint[ 2 ] << "\n";
    cout << "                       ";
    cout << setw( width ) << connPoint[ 3 ];
    cout << setw( width ) << connPoint[ 4 ];
    cout << setw( width ) << connPoint[ 5 ] << "\n";
    cout << "   range (donor zone)= ";
    cout << setw( width ) << connDonorPoint[ 0 ];
    cout << setw( width ) << connDonorPoint[ 1 ];
    cout << setw( width ) << connDonorPoint[ 2 ] << "\n";
    cout << "                       ";
    cout << setw( width ) << connDonorPoint[ 3 ];
    cout << setw( width ) << connDonorPoint[ 4 ];
    cout << setw( width ) << connDonorPoint[ 5 ] << "\n";
    cout << "   transform = " << itranfrm[ 0 ] << " " << itranfrm[ 1 ] << " " << itranfrm[ 2 ] << "\n";

    int transform[ 3 ][ 3 ];

    //For 3-D, the
    //transformation matrix T is constructed from Transform = [¡Àa, ¡Àb, ¡Àc ] as follows:
    //T =
    //[ sgn( a ) del( a ? 1 ) sgn( b ) del( b ? 1 ) sgn( c ) del( c ? 1 )]
    //[ sgn( a ) del( a ? 2 ) sgn( b ) del( b ? 2 ) sgn( c ) del( c ? 2 )]
    //[ sgn( a ) del( a ? 3 ) sgn( b ) del( b ? 3 ) sgn( c ) del( c ? 3 )]

    int celldim = this->cgnsZone->cgnsBase->celldim;

    for ( int i = 0; i < celldim; ++ i )
    {
        for ( int j = 0; j < celldim; ++ j )
        {
            transform[ i ][ j ] = SIGN( 1, itranfrm[ j ] ) * AbsoluteDiagonalId( itranfrm[ j ], i + 1 );
        }
    }
}


#endif

EndNameSpace