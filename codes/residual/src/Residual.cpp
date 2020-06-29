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

#include "Residual.h"
#include "Parallel.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

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
    this->nCell = 0;
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
    int nCellSum = 0;
    HXReduceInt( & this->nCell, & nCellSum, 1, PL_SUM );
    HXReduceReal( & this->res[ 0 ], & resSum[ 0 ], nEqu, PL_SUM );

    this->nCell = nCellSum;
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
    this->nCell += rhs.nCell;
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
    ResData & t = dataList[ 0 ];
    * this = t.resmax;

    for ( int i = 0; i < dataList.size(); ++ i )
    {
        ResData & t = dataList[ i ];
        this->SwapMax( t.resmax );
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