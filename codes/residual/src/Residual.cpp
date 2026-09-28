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

#include "Residual.h"
#include "Parallel.h"
#include "HXMath.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

ResAver::ResAver()
{
    ;
}

ResAver::~ResAver()
{
    ;
}

void ResAver::Init( int nEqu )
{
    this->res.resize( nEqu );
}

void ResAver::Zero()
{
    this->res = 0;
    this->nCells = 0;
}

void ResAver::CalcAver( HXVector< ResData > & dataList )
{
    this->Zero();
    for ( int i = 0; i < dataList.size(); ++ i )
    {
        ResData & t = dataList[ i ];
        * this += t.resave;
    }

    int nEqu = this->res.size();
    RealField resSum( nEqu );
    resSum = this->res;
    int nCellSum = this->nCells;
    HXReduceInt( & this->nCells, & nCellSum, 1, PL_SUM );
    HXReduceReal( & this->res[ 0 ], & resSum[ 0 ], nEqu, PL_SUM );

    this->nCells = nCellSum;
    if ( nCellSum <= 0 )
    {
        this->res = 0;
        return;
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        this->res[ iEqu ] = sqrt( resSum[ iEqu ] / nCellSum );
    }
}

ResAver & ResAver::operator += ( const ResAver & rhs )
{
    int nEqu = this->res.size();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        this->res[ iEqu ] += rhs.res[ iEqu ];
    }
    this->nCells += rhs.nCells;
    return *this;
}

ResMax::ResMax()
{
    ;
}

ResMax::~ResMax()
{
    ;
}

void ResMax::Init( int nEqu )
{
    resmax.resize( nEqu );
    index .resize( nEqu );
    zid   .resize( nEqu );

    xcc.resize( nEqu );
    ycc.resize( nEqu );
    zcc.resize( nEqu );
    vol.resize( nEqu );
}

void ResMax::SwapMax( ResMax & rhs )
{
    int nEqu = resmax.size();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        Real a = this->resmax[ iEqu ];
        Real b = rhs.resmax[ iEqu ];
        if ( a < b )
        {
            this->zid[ iEqu ] = rhs.zid[ iEqu ];
            this->index[ iEqu ] = rhs.index[ iEqu ];
            this->resmax[ iEqu ] = rhs.resmax[ iEqu ];
            this->xcc[ iEqu ] = rhs.xcc[ iEqu ];
            this->ycc[ iEqu ] = rhs.ycc[ iEqu ];
            this->zcc[ iEqu ] = rhs.zcc[ iEqu ];
            this->vol[ iEqu ] = rhs.vol[ iEqu ];
        }
    }
}

void ResMax::CalcMax( HXVector< ResData > & dataList )
{
    const bool hasLocalData = ! dataList.empty();
    const int nEqu = resmax.size();

    if ( hasLocalData )
    {
        * this = dataList[ 0 ].resmax;
        for ( int i = 1; i < dataList.size(); ++ i )
        {
            this->SwapMax( dataList[ i ].resmax );
        }
    }
    else
    {
        resmax = 0;
        index = -1;
        zid = -1;
        xcc = 0;
        ycc = 0;
        zcc = 0;
        vol = 0;
    }

    RealField globalMax( nEqu );
    globalMax = resmax;
    HXReduceReal( & resmax[ 0 ], & globalMax[ 0 ], nEqu, PL_MAX );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        const Real localMax = resmax[ iEqu ];
        int ownerCandidate = Parallel::nProc;
        if ( hasLocalData && localMax == globalMax[ iEqu ] )
        {
            ownerCandidate = Parallel::pid;
        }

        int owner = ownerCandidate;
        HXReduceInt( & ownerCandidate, & owner, 1, PL_MIN );
        resmax[ iEqu ] = globalMax[ iEqu ];

        if ( owner >= Parallel::nProc )
        {
            index[ iEqu ] = -1;
            zid[ iEqu ] = -1;
            xcc[ iEqu ] = 0;
            ycc[ iEqu ] = 0;
            zcc[ iEqu ] = 0;
            vol[ iEqu ] = 0;
            continue;
        }

        HXBcast( & index[ iEqu ], 1, owner );
        HXBcast( & zid[ iEqu ], 1, owner );
        HXBcast( & xcc[ iEqu ], 1, owner );
        HXBcast( & ycc[ iEqu ], 1, owner );
        HXBcast( & zcc[ iEqu ], 1, owner );
        HXBcast( & vol[ iEqu ], 1, owner );
    }
}

int ResMax::CalcMaxId()
{
    int nEqu = resmax.size();
    int id = 0;
    Real s = resmax[ id ];
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        if ( s < resmax[ iEqu ] )
        {
            s = resmax[ iEqu ];
            id = iEqu;
        }
    }
    return id;
}

ResData::ResData()
{
    ;
}

ResData::~ResData()
{
    ;
}

void ResData::Init( int nEqu )
{
    resmax.Init( nEqu );
    resave.Init( nEqu );
}

Residual::Residual()
{
    ;
}

Residual::~Residual()
{
    ;
}

EndNameSpace
