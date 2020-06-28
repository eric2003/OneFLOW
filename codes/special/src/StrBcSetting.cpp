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

#include "StrBcSetting.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "StrGrid.h"

BeginNameSpace( ONEFLOW )

StrBcSetting::StrBcSetting()
{
}

StrBcSetting::~StrBcSetting()
{
}

void StrBcSetting::ConstructBcMap()
{
    boundaryMap.insert( pair< string, int >( "interface", BC::INTERFACE ) );
    boundaryMap.insert( pair< string, int >( "extrapolation", BC::EXTRAPOLATION ) );
    boundaryMap.insert( pair< string, int >( "symmetry", BC::SYMMETRY ) );
    boundaryMap.insert( pair< string, int >( "inflow", BC::INFLOW ) );
    boundaryMap.insert( pair< string, int >( "outflow", BC::OUTFLOW ) );
    boundaryMap.insert( pair< string, int >( "solidsurface", BC::SOLID_SURFACE ) );
    boundaryMap.insert( pair< string, int >( "pole", BC::POLE ) );
    boundaryMap.insert( pair< string, int >( "no", BC::NO_BOUNDARY ) );
    boundaryMap.insert( pair< string, int >( "farfield", BC::FARFIELD ) );
    boundaryMap.insert( pair< string, int >( "overset", BC::OVERSET ) );
    boundaryMap.insert( pair< string, int >( "generic1", BC::GENERIC_1 ) );
    boundaryMap.insert( pair< string, int >( "generic2", BC::GENERIC_2 ) );
    boundaryMap.insert( pair< string, int >( "generic3", BC::GENERIC_3 ) );
}

int StrBcSetting::GetBcType( const string & bcTypeName )
{
    map< string, int >::iterator iter;
    iter = boundaryMap.find( bcTypeName );
    return iter->second;
}

void StrBcSetting::PushBc( int imin, int imax, int jmin, int jmax, int kmin, int kmax, int bcType )
{
    iminList.push_back( imin );
    imaxList.push_back( imax );
    jminList.push_back( jmin );
    jmaxList.push_back( jmax );
    kminList.push_back( kmin );
    kmaxList.push_back( kmax );
    bcTypeList.push_back( bcType );
}

int StrBcSetting::CalcNBcRegion()
{
    int nSize = this->bcTypeList.size();
    int nCount = 0;
    for ( int ir = 0; ir < nSize; ++ ir )
    {
        if ( this->bcTypeList[ ir ] > 0 )
        {
            ++ nCount;
        }
    }
    return nCount;
}

void StrBcSetting::SetBcRegion( StrGrid * grid )
{
    int nBcRegion = CalcNBcRegion();

    int iZone = grid->id;

    grid->bcRegionGroup->Create( nBcRegion );
    BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
    int ir_count = 0;
    int ir = 0;
    while ( ir < nBcRegion )
    {
        int imin = iminList[ ir_count ];
        int imax = imaxList[ ir_count ];
        int jmin = jminList[ ir_count ];
        int jmax = jmaxList[ ir_count ];
        int kmin = kminList[ ir_count ];
        int kmax = kmaxList[ ir_count ];
        int bcType = bcTypeList[ ir_count ];

        BcRegion * bcRegion = new BcRegion( iZone, ir );
        bcRegion->s->SetRegion( imin, imax, jmin, jmax, kmin, kmax );
        bcRegion->s->zid = iZone;
        bcRegion->bcType = bcType;
        bcRegionGroup->SetBcRegion( ir, bcRegion );

        if ( bcType < 0 )
        {
            ++ ir_count;
            imin = iminList[ ir_count ];
            imax = imaxList[ ir_count ];
            jmin = jminList[ ir_count ];
            jmax = jmaxList[ ir_count ];
            kmin = kminList[ ir_count ];
            kmax = kmaxList[ ir_count ];
            int zid = bcTypeList[ ir_count ];

            bcRegion->t->SetRegion( imin, imax, jmin, jmax, kmin, kmax );
            bcRegion->t->zid = zid;
        }
        ++ ir;
        ++ ir_count;
    }
}

EndNameSpace