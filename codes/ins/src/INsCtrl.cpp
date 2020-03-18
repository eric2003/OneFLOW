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
#include "INsCtrl.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )

INsCtrl ins_ctrl;

INsCtrl::INsCtrl()
{
    ;
}

INsCtrl::~INsCtrl()
{
    ;
}

void INsCtrl::Init()
{
    this->time_integral = GetDataValue< int >( "time_integral" );
    this->linearTwoStepMethods = GetDataValue< int >( "linearTwoStepMethods" );

    this->ieigenfix = GetDataValue< int >( "ieigenfix" );
    this->centropy1 = GetDataValue< Real >( "centropy1" );
    this->centropy2 = GetDataValue< Real >( "centropy2" );

    this->idump = GetDataValue< int >( "idump" );
    this->inflowType = GetDataValue< int >( "inflowType" );
    this->ireadwdst = GetDataValue< int >( "ireadwdst" );
    this->showfield = GetDataValue< int >( "showfield" );

    this->isowallbc = GetDataValue< int >( "isowallbc" );
    
    if ( inflowType == 3 )
    {
        int numpp = 4;
        initplane.resize( numpp );
        int nEqu = 5;
        initflow1.resize( nEqu );
        initflow2.resize( nEqu );

        CopyArray( initplane, "initplane" );
        CopyArray( initflow1, "initflow1" );
        CopyArray( initflow2, "initflow2" );

        int kkk = 1;
    }

    rk_stage = GetDataValue< int >( "rk_stage" );
    for ( int i = 0; i < rk_stage; ++ i )
    {
        Real coef = 1.0 / ( i + 1 );
        rk_coef.push_back( coef );
    }

    pdt = GetDataValue< Real >( "global_dt" );
    maxTime = GetDataValue< Real >( "maxTime" );
    iexitflag = GetDataValue< int >( "iexitflag" );
    currTime = 0.0;
    ilim = GetDataValue< int >( "ilim" );
    vencat_coef = GetDataValue< Real >( "vencat_coef" );

    nrokplus = 0;

    ivischeme = GetDataValue< int >( "ivischeme" );
}

EndNameSpace