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

#include "HeatFlux.h"
#include "Zone.h"
#include "ZoneState.h"

BeginNameSpace( ONEFLOW )

HeatFlux heat_flux;

SurfaceValue::SurfaceValue()
{
    // FIX: Use std::make_unique for exception-safe allocation.
    var = std::make_unique<RealField>();
}

SurfaceValue::~SurfaceValue()
{
    // std::unique_ptr automatically cleans up the RealField.
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

    // FIX: Use .get() to access the raw pointer from unique_ptr for short-term observation.
    SurfaceValue * heat_sur = heat_flux.heatflux[ ZoneState::zid ].get();
    heat_sur->var->resize( 0 );

    SurfaceValue * fric_sur = heat_flux.fricflux[ ZoneState::zid ].get();
    fric_sur->var->resize( 0 );
}

void HeatFlux::InitGlobal()
{
    if ( init_flag ) return;
    init_flag = true;
    this->heatflux.resize( ZoneState::nZones );
    this->fricflux.resize( ZoneState::nZones );
    this->flag.resize( ZoneState::nZones, 0 );
}

void HeatFlux::Allocate()
{
    int zId = ZoneState::zid;
    if ( ! this->flag[ zId ] )
    {
        this->flag[ zId ] = 1;
        // FIX: Use std::make_unique instead of new.
        this->heatflux[ zId ] = std::make_unique<SurfaceValue>();
        this->fricflux[ zId ] = std::make_unique<SurfaceValue>();
    }
}

void HeatFlux::DeAllocate()
{
    // FIX: clear() automatically invokes the destructor of std::unique_ptr,
    // safely releasing all SurfaceValue and their internal RealField objects.
    // No manual delete loop is needed, preventing memory leaks on exceptions.
    this->heatflux.clear();
    this->fricflux.clear();
    this->flag.clear();
    this->init_flag = false;
}

EndNameSpace