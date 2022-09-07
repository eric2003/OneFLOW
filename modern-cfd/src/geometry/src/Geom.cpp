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
#include "Geom.h"
#include "Cmpi.h"
#include "Grid.h"
#include "Project.h"
#include "tools.h"
#ifdef PRJ_ENABLE_CGNS
#include "cgnstest.h"
#endif
#include <iostream>

void HXGenerateGrid( int ni, float xmin, float xmax, float * xcoor )
{
    float dx = ( xmax - xmin ) / ( ni - 1 );
    int ist = 0;
    int ied = ni + 1;
    for ( int i = 1; i <= ni; ++ i )
    {
        float xm = xmin + ( i - 1 ) * dx;

        xcoor[ i ] = xm;
    }
    xcoor[ ist ] = 2 * xcoor[ ist + 1 ] - xcoor[ ist + 2 ];
    xcoor[ ied ] = 2 * xcoor[ ied - 1 ] - xcoor[ ied - 2 ];
}

void CopyGrid( std::vector<float> xcoor, float * xcoor_global )
{
    int ni = xcoor.size();
    int ist = 0;
    int ied = ni + 1;
    //for ( int i = 1; i <= ni; ++ i )
    for ( int i = 0; i < ni; ++ i )
    {
        xcoor_global[ i + 1 ] = xcoor[ i ];
    }
    xcoor_global[ ist ] = 2 * xcoor_global[ ist + 1 ] - xcoor_global[ ist + 2 ];
    xcoor_global[ ied ] = 2 * xcoor_global[ ied - 1 ] - xcoor_global[ ied - 2 ];
}

int Geom_t::ni_ghost = 2;
int Geom_t::ni_global = 1;
int Geom_t::ni_global_total = -1;
float * Geom_t::xcoor_global = 0;
float Geom_t::dx = -1.0;
int Geom_t::gridobj = 1;
std::string Geom_t::gridName;

std::vector<int> Geom_t::zone_nis;
std::vector<int> Geom_t::proc_ids;
std::vector<int> Geom_t::zone_ids;

Geom_t::Geom_t()
{
}

Geom_t::~Geom_t()
{
}

void Geom_t::Init()
{
#ifdef PRJ_ENABLE_CGNS
    if ( Cmpi::pid == 0 )
    {
        //TestCgnsLink();
    }
#endif
    Geom_t::gridName = add_string( Project::prj_grid_dir, "/", Geom_t::gridName );
    //Geom_t::ni_global = 41;
    //Geom_t::ni_ghost = 2;
    //Geom_t::ni_global_total = Geom_t::ni_global + Geom_t::ni_ghost;

    //int nZones = Cmpi::nproc;
    //Geom_t::zone_nis.resize( nZones );
    //int grid_ni = ( Geom_t::ni_global + nZones - 1 ) / nZones;
    //int ni_last = Geom_t::ni_global - ( nZones - 1 ) * ( grid_ni - 1 );

    //for ( int i = 0; i < nZones - 1; ++ i )
    //{
    //    Geom_t::zone_nis[i] = grid_ni;
    //}
    //Geom_t::zone_nis[nZones - 1] = ni_last;
    //std::printf( "zone ni----------------------\n" );
    //for ( int i = 0; i < nZones; ++ i )
    //{
    //    std::printf( "%d ", Geom_t::zone_nis[i] );
    //}
    //std::printf( "\n" );
}

void Geom_t::Finalize()
{
    delete [] Geom_t::xcoor_global;
    Geom_t::xcoor_global = 0;
}

void Geom_t::DumpGrid( const std::string & fileName )
{
#ifdef PRJ_ENABLE_CGNS
    ::cgns_dump_grid( Geom_t::xcoor_global+1, Geom_t::ni_global, fileName );
#endif
}

void Geom_t::GenerateGrid()
{
    Geom_t::ni_global = 41;
    Geom_t::ni_ghost = 2;
    Geom_t::ni_global_total = Geom_t::ni_global + Geom_t::ni_ghost;

    int nZones = Cmpi::nproc;
    Geom_t::zone_nis.resize( nZones );
    int grid_ni = ( Geom_t::ni_global + nZones - 1 ) / nZones;
    int ni_last = Geom_t::ni_global - ( nZones - 1 ) * ( grid_ni - 1 );

    for ( int i = 0; i < nZones - 1; ++ i )
    {
        Geom_t::zone_nis[i] = grid_ni;
    }
    Geom_t::zone_nis[nZones - 1] = ni_last;
    std::printf( "zone ni----------------------\n" );
    for ( int i = 0; i < nZones; ++ i )
    {
        std::printf( "%d ", Geom_t::zone_nis[i] );
    }
    std::printf( "\n" );

    Geom_t::xcoor_global = new float[ Geom_t::ni_global_total ];
    float xmin = 0.0;
    float xmax = 2.0;
    ::HXGenerateGrid( Geom_t::ni_global, xmin, xmax, Geom_t::xcoor_global );
    if ( Cmpi::pid == 0 )
    {
        Geom_t::DumpGrid( Geom_t::gridName );
    }
}

void Geom_t::ReadGrid( const std::string &gridName )
{
    if ( Cmpi::pid == 0 )
    {
        std::vector<float> xcoor;
        ::cgns_read_grid( xcoor, gridName );
        int ni = xcoor.size();
        Geom_t::ni_global = ni;
        Geom_t::ni_ghost = 2;
        Geom_t::ni_global_total = Geom_t::ni_global + Geom_t::ni_ghost;
        Geom_t::xcoor_global = new float[ Geom_t::ni_global_total ];
        ::CopyGrid( xcoor, xcoor_global );

        int nZones = Cmpi::nproc;
        Geom_t::zone_nis.resize( nZones );
        int grid_ni = ( Geom_t::ni_global + nZones - 1 ) / nZones;
        int ni_last = Geom_t::ni_global - ( nZones - 1 ) * ( grid_ni - 1 );

        for ( int i = 0; i < nZones - 1; ++ i )
        {
            Geom_t::zone_nis[i] = grid_ni;
        }
        Geom_t::zone_nis[nZones - 1] = ni_last;
        std::printf( "zone ni----------------------\n" );
        for ( int i = 0; i < nZones; ++ i )
        {
            std::printf( "%d ", Geom_t::zone_nis[i] );
        }
        std::printf( "\n" );
    }
}

void Geom_t::ReadAndModifyGrid( const std::string &gridName )
{
    if ( Cmpi::pid == 0 )
    {
        std::vector<float> xcoor;
        ::cgns_read_and_modify_grid( xcoor, gridName );
        int ni = xcoor.size();
        Geom_t::ni_global = ni;
        Geom_t::ni_ghost = 2;
        Geom_t::ni_global_total = Geom_t::ni_global + Geom_t::ni_ghost;
        Geom_t::xcoor_global = new float[ Geom_t::ni_global_total ];
        ::CopyGrid( xcoor, xcoor_global );

        int nZones = Cmpi::nproc;
        Geom_t::zone_nis.resize( nZones );
        int grid_ni = ( Geom_t::ni_global + nZones - 1 ) / nZones;
        int ni_last = Geom_t::ni_global - ( nZones - 1 ) * ( grid_ni - 1 );

        for ( int i = 0; i < nZones - 1; ++ i )
        {
            Geom_t::zone_nis[i] = grid_ni;
        }
        Geom_t::zone_nis[nZones - 1] = ni_last;
        std::printf( "zone ni----------------------\n" );
        for ( int i = 0; i < nZones; ++ i )
        {
            std::printf( "%d ", Geom_t::zone_nis[i] );
        }
        std::printf( "\n" );
    }
}


void Geom_t::LoadGrid()
{
    if ( Geom_t::gridobj == 0 )
    {
        Geom_t::GenerateGrid();
    }
    else if ( Geom_t::gridobj == 1 )
    {
        Geom_t::ReadGrid( Geom_t::gridName );
    }
    else if ( Geom_t::gridobj == 2 )
    {
        Geom_t::ReadAndModifyGrid( Geom_t::gridName );
    }
}

Geom::Geom()
{
    this->xcoor = 0;
    this->ds = 0;
}

Geom::~Geom()
{
    delete [] this->xcoor;
    delete [] this->ds;
}

void Geom::Init()
{
    this->nZones = Cmpi::nproc;
    this->zoneId = Cmpi::pid;
    this->ni = Geom_t::zone_nis[ this->zoneId ];
    this->ni_total = this->ni + Geom_t::ni_ghost;
    this->xcoor = new float[ this->ni_total ];
    this->ds = new float[ this->ni_total ];

    this->bcSolver = new BoundarySolver{};
    this->bcSolver->zoneId = zoneId;
    this->bcSolver->Init( zoneId, nZones, this->ni );
}

void Geom::GenerateGrid()
{
    if ( Cmpi::pid == 0 )
    {
        std::cout << "Running on " << Cmpi::nproc << " nodes" << std::endl;

        for ( int i = 0; i < this->ni_total; ++ i )
        {
            this->xcoor[ i ] = Geom_t::xcoor_global[ i ];
        }
        int istart = 0;
        #ifdef PRJ_ENABLE_MPI
        for ( int ip = 1; ip < Cmpi::nproc; ++ ip )
        {
            int ni_tmp = Geom_t::zone_nis[ ip ];
            int ni_total_tmp = ni_tmp + Geom_t::ni_ghost;
            istart += ( ni_tmp - 1 );
            float * source_start = Geom_t::xcoor_global + istart;
            MPI_Send( source_start, ni_total_tmp, MPI_FLOAT, ip, 0, MPI_COMM_WORLD );
        }
        #endif
    }
    else
    {
        #ifdef PRJ_ENABLE_MPI
        MPI_Recv(this->xcoor, this->ni_total, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        #endif
    }

    std::printf("print xcoor: process id = %d Cmpi::nproc = %d\n", Cmpi::pid, Cmpi::nproc );
}

void Geom::ComputeGeom()
{
    for ( int i = 1; i < this->ni_total - 1; ++ i )
    {
        this->ds[ i ] = this->xcoor[ i ] - this->xcoor[ i - 1 ];
    }

    this->ds[ 0 ] = this->ds[ 1 ];
    this->ds[ this->ni + 1 ] = this->ds[ this->ni ];
}
