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

#include "FieldPara.h"
#include "DataBase.h"
#include <iostream>
#include <fstream>
#include <vector>


BeginNameSpace( ONEFLOW )

FieldPara::FieldPara()
{
    ;
}

FieldPara::~FieldPara()
{
    ;
}

//void FieldPara::Init()
//{
//    //this->nx = 41;
//    //this->len = 2.0;
//    //this->dx = len / ( nx - 1.0 );
//    ////this->nt = 25;
//    //this->nt = 25;
//    //this->dt = 0.025;
//    //this->c  = 1;
//
//    this->nx  = ONEFLOW::GetDataValue< int >("scalar_nx");
//    this->len = ONEFLOW::GetDataValue< int >("scalar_len");
//    this->nt  = ONEFLOW::GetDataValue< int >("scalar_nt");
//    this->dt  = ONEFLOW::GetDataValue< Real >("scalar_dt");
//    this->c   = ONEFLOW::GetDataValue< Real >("scalar_c");
//    this->dx  = this->len / ( this->nx - 1.0 );
//}

void FieldPara::Init()
{
    int nbase = 40;
    int base_nt = 25;

    this->nx  = ONEFLOW::GetDataValue< int >("scalar_nx");
    this->len = ONEFLOW::GetDataValue< int >("scalar_len");
    this->c   = ONEFLOW::GetDataValue< Real >("scalar_c");

    double base_dx =  this->len / nbase;
    double base_dt = 0.025;
    double cfl = base_dt / base_dx;

    double total_t = 25 * 0.025;

    this->dx = len / ( this->nx - 1.0 );
    int nratio = ( this->nx - 1 ) / nbase;
    this->dt = base_dt / nratio;
    this->nt  = base_nt * nratio;
    std::cout << " nratio = " << nratio << "\n";
    std::cout << " this->dt = " << this->dt << "\n";
    std::cout << " this->nt = " << this->nt << "\n";
    std::cout << " this->dx = " << this->dx << "\n";
    std::cout << "  dx = " << base_dx / nratio << "\n";
}


EndNameSpace
