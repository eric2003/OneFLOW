/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "UINsSpectrum.h"
#include "INsCom.h"
#include "UINsCom.h"
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


UINsSpectrum::UINsSpectrum()
{
}

UINsSpectrum::~UINsSpectrum()
{
}

void UINsSpectrum::CmpImplicitSpectrum()
{
    this->CmpUnsteadySpectrum();

    this->AddInvSpectrum();

    this->AddVisSpectrum();

}

void UINsSpectrum::CmpUnsteadySpectrum()
{
    if ( ctrl.idualtime == 0 )//单时间步，注意:是usd.sp2!
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            ( * uinsf.impsr )[ 0 ][ cId ] = ( usd.sp2 /  ( * uinsf.timestep )[ 0 ][ cId ] ) * ( * ug.cvol )[ cId ];
        }
    }
    else
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            ( * uinsf.impsr )[ 0 ][ cId ] = ( usd.sp1 / ( * uinsf.timestep )[ 0 ][ cId ] + usd.sp2 / ctrl.pdt1 ) * ( * ug.cvol )[ cId ];
        }
    }
}

void UINsSpectrum::AddInvSpectrum()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * uinsf.impsr )[ 0 ][ cId ] += ( * uinsf.invsr )[ 0 ][ cId ];
    }
}

void UINsSpectrum::AddVisSpectrum()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * uinsf.impsr )[ 0 ][ cId ] += ( * uinsf.vissr )[ 0 ][ cId ];
    }
}

EndNameSpace