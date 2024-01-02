/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#ifdef PRJ_ENABLE_MPI
#include <mpi.h>
#endif

#include <iostream>

class ServerCout;

class Cmpi
{
public:
    Cmpi();
    ~Cmpi();
public:
    static void Init(int argc, char **argv);
    static void Finalize();
    static bool IsServer();
public:
    static int pid;
    static int nproc;
    static int server_id;
    static int num_gpus;
    static int num_cpus;
    static ServerCout server_out;
public:
    template <typename... Args>
    static void printf( char const * const format, Args&&... args )
    {
        if ( Cmpi::IsServer() )
        {
            std::printf( format, args... );
        }
    }
};

class ServerCout
{
public:
    ServerCout();
    ~ServerCout();
public:
    ServerCout&  operator<<(ServerCout&(* _Pfn)(ServerCout&) ) {
        return _Pfn(*this);
    }
};

template< typename T >
ServerCout & operator << ( ServerCout & sout, const T & value )
{
    if ( Cmpi::IsServer() )
    {
        std::cout << value;
    }
    return sout;
}




