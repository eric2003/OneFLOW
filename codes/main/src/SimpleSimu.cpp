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
#include "SimpleSimu.h"
#include "Scalar.h"
#include "FieldSolver.h"
#include "FieldSolverOpenMP.h"
#include "FieldSolverCuda.h"
#include "AccelRuntime.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

void ToyModelSimu()
{
    std::cout << "ToyModelSimu\n";
    if ( AccelRuntime::Instance().IsAccelerator() )
    {
        FieldSolver fieldSolver;
        fieldSolver.Run();
        return;
    }

    //FieldSolverOpenMP * fieldSolverOpenMP = new FieldSolverOpenMP();
    //fieldSolverOpenMP->Run();
    //delete fieldSolverOpenMP;

    FieldSolver fieldSolver;
    fieldSolver.Run();

}


//void HybridParallel::Run()
//{
//    std::cout << "Run HybridParallel test\n";
//    // your original test logic
//}

//void CgnsTest::Run()
//{
//    std::cout << "Run CgnsTest\n";
//}

//void JsonTest::Run()
//{
//    std::cout << "Run JsonTest\n";
//}

EndNameSpace

// register test cases, outside namespace
//REGISTER_TEST_CASE(HybridParallel, "hybrid_parallel");
//REGISTER_TEST_CASE(MpiTest, "mpi_test");
//REGISTER_TEST_CASE(CgnsTest, "cgns_test");
//REGISTER_TEST_CASE(JsonTest, "json_test");
