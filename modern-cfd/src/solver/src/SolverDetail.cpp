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

#include "SolverDetail.h"
#include "SolverDetailCpu.h"
#ifdef PRJ_ENABLE_CUDA
#include "SolverDetailCuda.h"
#endif

void SolverInit()
{
#ifdef PRJ_ENABLE_CUDA
    SolverInitCuda();
#else
    SolverInitCpu();
#endif
}

void SetDevice( int cpu_thread_id )
{
#ifdef PRJ_ENABLE_CUDA
    SetDeviceCuda( cpu_thread_id );
#else

#endif
}

void CfdCopyVector( float * a, float * b, int ni )
{
#ifdef PRJ_ENABLE_CUDA
    CfdCopyVectorCuda( a,  b, ni );
#else
    CfdCopyVectorCpu( a,  b, ni );
#endif
}

void CfdScalarUpdate( float * q, float * qn, float c, float * timestep, float * ds, int ni )
{
#ifdef PRJ_ENABLE_CUDA
    CfdScalarUpdateCuda( q, qn, c, timestep, ds, ni );
#else
    CfdScalarUpdateCpu( q, qn, c, timestep, ds, ni );
#endif

}