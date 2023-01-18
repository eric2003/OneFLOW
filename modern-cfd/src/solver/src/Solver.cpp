/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

#include "Solver.h"
#include "SolverDetail.h"
#include "Project.h"
#include "tools.h"
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "Cmpi.h"
#include "Grid.h"
#include "Geom.h"
#include "CfdPara.h"
#include "Visual.h"
#ifdef PRJ_ENABLE_CGNS
#include <cgnslib.h>
#endif

float SquareFun( float xm )
{
    if ( xm >= 0.5 && xm <= 1.0 )
    {
        return 2.0;
    }
    return 1.0;
}

void Theory( float time, float c, std::vector<float>& theory, std::vector<float>& xcoor )
{
    int ni = xcoor.size();
    float xs = c * time;
    for ( int i = 0; i < ni; ++ i )
    {
        float xm = xcoor[i];
        float xm_new = xm - xs;
        float fm = SquareFun( xm_new );
        theory[i] = fm;
    }
}

Solver::Solver()
{
    SolverInit();
}

Solver::~Solver()
{
    ;
}

void Solver::Init()
{
}

void Solver::Run( CfdPara * cfd_para, Geom * geom )
{
    this->CfdSolve( cfd_para, geom );
}

void Solver::AllocateField( Geom * geom )
{
    this->q = new float[ geom->ni_total ];
    this->qn = new float[ geom->ni_total ];
    this->timestep = new float[ geom->ni_total ];
}

void Solver::DeallocateField( Geom * geom )
{
    delete [] this->q;
    delete [] this->qn;
    delete [] this->timestep;
}

void Solver::InitField( CfdPara * cfd_para, Geom * geom )
{
    if ( cfd_para->irestart == 0 )
    {
        this->SetInflowField( cfd_para, geom );
    }
    else
    {
        this->ReadField( cfd_para, geom );
    }
}

void Solver::SetInflowField( CfdPara * cfd_para, Geom * geom )
{
    this->iterstep = 0;
    this->time_now = 0.0;
    for ( int i = 0; i < geom->ni_total; ++ i )
    {
        float fm = SquareFun( geom->xcoor[ i ] );
        this->q[ i ] = fm;
    }
}

void Solver::CfdSolve( CfdPara * cfd_para, Geom * geom )
{
    this->AllocateField( geom );
    this->InitField( cfd_para, geom );
    this->SolveField( cfd_para, geom );
    this->SaveField( cfd_para, geom );
    this->Visualize( cfd_para, geom );
    this->DeallocateField( geom );
}

void Solver::Timestep( CfdPara * cfd_para, Geom * geom )
{
    this->dt = std::numeric_limits<float>::max();
    for ( int i = 0; i < geom->ni_total; ++ i )
    {
        this->timestep[ i ] = geom->ds[ i ] * cfd_para->cfl / cfd_para->cspeed;
        this->dt = std::min( this->timestep[ i ], this->dt );
    }
    //std::cout << " this->dt = " << this->dt << std::endl;
}

void Solver::SolveField( CfdPara * cfd_para, Geom * geom )
{
    std::printf( " this->time_now = %f, cfd_para->simu_time = %f \n", this->time_now, cfd_para->simu_time );
    this->Timestep( cfd_para, geom );
    int ntimestep = cfd_para->simu_time / this->dt;
    while ( ( this->time_now < cfd_para->simu_time ) &&
            ( std::abs(this->time_now-cfd_para->simu_time) > std::numeric_limits<float>::epsilon() ) )
    {
        if ( Cmpi::IsServer() )
        {
            std::printf( " iStep = %d, nStep = %d \n", this->iterstep + 1, ntimestep );
            //std::printf( " this->time_now = %f, cfd_para->simu_time = %f \n", this->time_now, cfd_para->simu_time );
        }

        this->Boundary( q, geom );
        this->Timestep( cfd_para, geom );

        int nCpuThreads = 1;
        omp_set_num_threads( nCpuThreads );
        #pragma omp parallel
        {
            int cpu_thread_id = omp_get_thread_num();
            int num_cpu_threads = omp_get_num_threads();
            //int gpu_id = -1;
            SetDevice( cpu_thread_id );
            CfdCopyVector( qn, q, geom->ni_total );
        }
        CfdScalarUpdate(this->q, this->qn, cfd_para->cspeed, this->timestep, geom->ds, geom->ni );
        ++ this->iterstep;
        this->time_now += this->dt;
        if ( this->iterstep % 10 == 0 )
        {
            this->SaveField( cfd_para, geom );
            //std::cin.get();
        }
    }
}

void Solver::Boundary( float * q, Geom * geom )
{
    BoundarySolver * bcSolver = geom->bcSolver;
    //physical boundary
    int nBFace = bcSolver->GetNBFace();
    //std::printf(" Boundary zoneID = %d nBFace = %d\n", bcSolver->zoneId, nBFace);
    for ( int iface = 0; iface < nBFace; ++ iface )
    {
        int bctype = bcSolver->bctypes[ iface ];
        int ghostcell_id = bcSolver->bc_ghostcells[ iface ];
        int bc_faceid = bcSolver->bc_faceids[ iface ];
        if ( bctype == BCInterface ) continue;
        if ( bctype == BCInflow )
        {
            float xm = geom->xcoor[ ghostcell_id ];
            q[ ghostcell_id ] = SquareFun( xm );
        }
        else if ( bctype == BCOutflow )
        {
            q[ ghostcell_id ] = q[ bc_faceid ];
        }
    }

    this->BoundaryInterface( q, geom );
}

void Solver::BoundaryInterface( float * q, Geom * geom )
{
    BoundarySolver * bcSolver = geom->bcSolver;
    int nIFace = bcSolver->GetNIFace();
    //std::printf( " BoundaryInterface nIFace = %d\n", nIFace );
    InterfaceSolver * interfaceSolver = bcSolver->interfaceSolver;
    interfaceSolver->SwapData( q );
    for ( int iface = 0; iface < nIFace; ++ iface )
    {
        int ghostcell_id = interfaceSolver->interface_ghost_cells[ iface ];
        //interfaceSolver->ShowInfo( iface );
        float bcvalue = interfaceSolver->GetBcValue( iface );
        q[ ghostcell_id ] = bcvalue;
    }
}

void Solver::ReadField( CfdPara * cfd_para, Geom * geom )
{
    //for ( int i = 0; i < geom->ni_total; ++ i )
    //{
    //    float fm = SquareFun( geom->xcoor[ i ] );
    //    this->q[ i ] = fm;
    //}

    std::string filename = add_string( Project::prj_results_dir, "./flow", geom->zoneId, ".dat" );
    std::fstream file;
    file.open( filename.c_str(), std::fstream::in | std::fstream::binary);

    int ni_total = geom->ni_total;
    file.read( reinterpret_cast<char *>(&this->iterstep), sizeof(int) );
    file.read( reinterpret_cast<char *>(&this->time_now), sizeof(float) );
    std::printf( "ReadField: this->iterstep = %d, this->time_now = %f, cfd_para->simu_time = %f \n", this->iterstep, this->time_now, cfd_para->simu_time );
    file.read( reinterpret_cast<char *>(&ni_total), sizeof(int) );
    file.read( reinterpret_cast<char *>(this->q), ni_total * sizeof(float) );
    file.close();
}

void Solver::SaveField( CfdPara * cfd_para, Geom * geom )
{
    std::string filename = add_string( Project::prj_results_dir, "./flow", geom->zoneId, ".dat" );
    std::fstream file;
    file.open( filename.c_str(), std::fstream::out | std::fstream::binary);

    int ni_total = geom->ni_total;
    std::printf( "SaveField: this->iterstep = %d, this->time_now = %f, cfd_para->simu_time = %f \n", this->iterstep, this->time_now, cfd_para->simu_time );
    file.write( reinterpret_cast<char *>(&this->iterstep), sizeof(int) );
    file.write( reinterpret_cast<char *>(&this->time_now), sizeof(float) );
    file.write( reinterpret_cast<char *>(&ni_total), sizeof(int) );
    file.write( reinterpret_cast<char *>(this->q), ni_total * sizeof(float) );
    file.close();
}

void Solver::Visualize( CfdPara * cfd_para, Geom * geom )
{
    std::string filename = add_string( Project::prj_results_dir, "./cfd", geom->zoneId, ".png" );
    Visual( this->q, geom->xcoor, geom->ni, filename );

    std::vector<float> q_global;
    std::vector<float> x_global;
    int root = 0;
    int tag = 0;
    if ( geom->zoneId != 0 )
    {
#ifdef PRJ_ENABLE_MPI
        MPI_Send( this->q, geom->ni_total, MPI_FLOAT, root, tag, MPI_COMM_WORLD );
#endif
    }
    else
    {
        std::vector<std::vector<float>> qvec( Cmpi::nproc );
        for ( int ip = 1; ip < Cmpi::nproc; ++ ip )
        {
            int ni_tmp = Geom_t::zone_nis[ ip ];
            int ni_total_tmp = ni_tmp + Geom_t::ni_ghost;

            qvec[ ip ].resize( ni_total_tmp );
        }
        qvec[ 0 ].insert( qvec[ 0 ].end(), this->q, this->q + geom->ni_total );

        for ( int ip = 1; ip < Cmpi::nproc; ++ ip )
        {
            int ni_tmp = Geom_t::zone_nis[ ip ];
            int ni_total_tmp = ni_tmp + Geom_t::ni_ghost;
#ifdef PRJ_ENABLE_MPI
            MPI_Recv( qvec[ ip ].data(), ni_total_tmp, MPI_FLOAT, ip, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
#endif
        }

        for ( int ip = 0; ip < Cmpi::nproc; ++ ip )
        {
            if ( ip == 0 )
            {
                q_global.insert( q_global.end(), qvec[ ip ].begin() + 1, qvec[ ip ].end() - 1 );
            }
            else
            {
                q_global.insert( q_global.end(), qvec[ ip ].begin() + 2, qvec[ ip ].end() - 1 );
            }
        }
        x_global.insert( x_global.end(), Geom_t::xcoor_global + 1, Geom_t::xcoor_global + Geom_t::ni_global + 1 );
        std::vector<float> theory;
        theory.resize( x_global.size() );
        Theory( cfd_para->simu_time, cfd_para->cspeed, theory, x_global );
        //Visual( q_global, theory, x_global, "./cfd.png" );
        std::string visual_filename = add_string( Project::prj_results_dir, "./cfd.pdf" );
        Visual( q_global, theory, x_global, visual_filename );
    }
}

