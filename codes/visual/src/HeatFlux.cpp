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

#include "HeatFlux.h"
#include "Zone.h"
#include "ZoneState.h"

BeginNameSpace( ONEFLOW )

HeatFlux heat_flux;

SurfaceValue::SurfaceValue()
{
    var = new RealField();
}

SurfaceValue::~SurfaceValue()
{
    delete var;
}

HeatFlux::HeatFlux()
{
    init_flag = false;
}

HeatFlux::~HeatFlux()
{
    DeAllocate();
}

void HeatFlux::Init()
{
    InitGlobal();
    Allocate();

    SurfaceValue * heat_sur = heat_flux.heatflux[ ZoneState::zid ];
    heat_sur->var->resize( 0 );

    SurfaceValue * fric_sur = heat_flux.fricflux[ ZoneState::zid ];
    fric_sur->var->resize( 0 );
}

void HeatFlux::InitGlobal()
{
    if ( init_flag ) return;
    init_flag = true;
    this->heatflux.resize( ZoneState::nZones, 0 );
    this->fricflux.resize( ZoneState::nZones, 0 );
    this->flag.resize( ZoneState::nZones, 0 );
}

void HeatFlux::Allocate()
{
    int zId = ZoneState::zid;
    if ( ! this->flag[ zId ] )
    {
        this->flag[ zId ] = 1;
        SurfaceValue * heat = new SurfaceValue();
        SurfaceValue * fric = new SurfaceValue();
        this->heatflux[ zId ] = heat;
        this->fricflux[ zId ] = fric;
    }
}

void HeatFlux::DeAllocate()
{
    int nSize = this->heatflux.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete this->heatflux[ i ];
        delete this->fricflux[ i ];
    }
}


EndNameSpace