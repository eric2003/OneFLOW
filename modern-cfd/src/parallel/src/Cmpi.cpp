/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "Cmpi.h"

int Cmpi::pid = 0;
int Cmpi::nproc = 1;
int Cmpi::server_id = 0;
int Cmpi::num_gpus = 0;
int Cmpi::num_cpus = 1;
ServerCout Cmpi::server_out;

Cmpi::Cmpi()
{
    ;
}

Cmpi::~Cmpi()
{
    ;
}

void Cmpi::Init(int argc, char **argv)
{
#ifdef PRJ_ENABLE_MPI
    MPI_Init( &argc, &argv ); 
    MPI_Comm_rank( MPI_COMM_WORLD, &Cmpi::pid ); 
    MPI_Comm_size( MPI_COMM_WORLD, &Cmpi::nproc );
#endif
}

void Cmpi::Finalize()
{
#ifdef PRJ_ENABLE_MPI
    MPI_Finalize();
#endif
}

bool Cmpi::IsServer()
{
    return Cmpi::pid == Cmpi::server_id;
}

ServerCout::ServerCout()
{
    ;
}

ServerCout::~ServerCout()
{
    ;
}