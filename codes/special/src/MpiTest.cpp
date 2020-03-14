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

#include "MpiTest.h"
#include "mpi.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

MpiTest::MpiTest()
{
    this->HX_MPI_Init();
}

MpiTest::~MpiTest()
{
    this->HX_MPI_Final();
}

void MpiTest::Init()
{
}

void MpiTest::Run()
{
    this->Init();
    this->Test();
}

void MpiTest::HX_MPI_Init()
{
    MPI_Init( 0, 0 );
    MPI_Comm_size( MPI_COMM_WORLD, & nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, & myrank );
    cout << "myrank = " << myrank << " nprocs = " << nprocs << "\n";
}

void MpiTest::HX_MPI_Final()
{
    MPI_Finalize();
}

void MpiTest::Test()
{
    if ( this->nprocs == 1 ) return;

    int src = 0, dest = 1;

    int buf[ 20 ];

    if ( myrank == src || myrank == dest )
    {
        int partner = src;
        if ( myrank == partner )
        {
            partner = dest;
        }
        MPI_Sendrecv( MPI_BOTTOM, 0, MPI_INT, partner, 0, MPI_BOTTOM, 0, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }

    if ( myrank == src )
    {
        MPI_Send( buf, 20, MPI_INT, dest, 0, MPI_COMM_WORLD );
    }
    else if ( myrank == dest )
    {
        MPI_Recv( buf, 20, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }
}

EndNameSpace