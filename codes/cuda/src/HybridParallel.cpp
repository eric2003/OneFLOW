/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include <mpi.h>

#ifdef ENABLE_CUDA
#include "PrintDevice.h"
#include "cuda_sub.h"
#endif

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_OPENACC
#include "myopenacc.h"
#endif

#include "HybridParallel.h"
#include "jacobi.h"
#include <stdio.h>
#include <iostream>

#include <iostream>
#include <algorithm>

BeginNameSpace( ONEFLOW )

HybridParallel::HybridParallel()
{
}

HybridParallel::~HybridParallel()
{
}

void HybridParallel::Run()
{
	int argc = 0;
	char ** argv = 0;
	this->HybridRun( argc, argv );
}

void HybridParallel::HybridRun( int argc, char ** argv )
{
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int numprocs,namelen,rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);
	printf("Hello from %d on %s out of %d\n",(rank+1),processor_name,numprocs);

#ifdef ENABLE_CUDA
	int num_gpus = -1;
	GetCudaDeviceCount( num_gpus );
	printf("number of CUDA devices:\t%d\n", num_gpus);
#endif

#ifdef ENABLE_OPENACC
	test_open_acc();
#endif
	Jacobi_Test();

#ifdef ENABLE_OPENMP
	printf("number of host CPUs:\t%d\n", omp_get_num_procs());	
#pragma omp parallel num_threads(4)
	{
		int nt = omp_get_thread_num();
		printf("Hello OneFLOW CFD: Hybrid MPI+Cuda+OpenACC+OpenMP! Thread = %d on CPU Rank%d on Machine %s\n", nt, rank + 1, processor_name);
	}
#endif

	MPI_Finalize();
}

EndNameSpace
