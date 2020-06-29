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

#include "UNsSpectrum.h"
#include "NsCom.h"
#include "UNsCom.h"
#include "UCom.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "HXMath.h"
#include "DataBase.h"
#include "FieldBase.h"
#include "UsdData.h"
#include "Ctrl.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UNsSpectrum::UNsSpectrum()
{
}

UNsSpectrum::~UNsSpectrum()
{
}

void UNsSpectrum::CalcImplicitSpectrum()
{
    this->CalcUnsteadySpectrum();

    this->AddInvSpectrum();

    this->AddVisSpectrum();

}

void UNsSpectrum::CalcUnsteadySpectrum()
{
    if ( ctrl.idualtime == 0 )//单时间步，注意:是usd.sp2!
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            ( * unsf.impsr )[ 0 ][ cId ] = ( usd.sp2 /  ( * unsf.timestep )[ 0 ][ cId ] ) * ( * ug.cvol )[ cId ];
        }
    }
    else
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            ( * unsf.impsr )[ 0 ][ cId ] = ( usd.sp1 / ( * unsf.timestep )[ 0 ][ cId ] + usd.sp2 / ctrl.pdt1 ) * ( * ug.cvol )[ cId ];
        }
    }
}

void UNsSpectrum::AddInvSpectrum()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * unsf.impsr )[ 0 ][ cId ] += ( * unsf.invsr )[ 0 ][ cId ];
    }
}

void UNsSpectrum::AddVisSpectrum()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * unsf.impsr )[ 0 ][ cId ] += ( * unsf.vissr )[ 0 ][ cId ];
    }
}

EndNameSpace