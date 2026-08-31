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

#include "BgField.h"
#include "Zone.h"
#include "ZoneState.h"
#include "SolverState.h"
#include "GridState.h"
#include "Multigrid.h"
#include "SolverMap.h"
#include "FieldWrap.h"
#include <iostream>

BeginNameSpace( ONEFLOW )
BasicBgField::BasicBgField()
{
}

BasicBgField::~BasicBgField()
{
    this->Free();
}

void BasicBgField::Init()
{
    int numberOfFields  = 2; //FIELD_FLOW = 0, FIELD_RHS = 1
    this->data.resize( SolverState::nSolver );

    for ( int solverIndex = 0; solverIndex < SolverState::nSolver; ++ solverIndex )
    {
        this->data[ solverIndex ].resize( numberOfFields );

        SolverState::solverIndex = solverIndex;
        SolverState::SetSolverTypeBySolverIndex( solverIndex );

        for ( int fid = 0; fid < numberOfFields; ++ fid )
        {
            this->data[ solverIndex ][ fid ].resize( MG::nMulti );
        
            for ( int gl = 0; gl < MG::nMulti; ++ gl )
            {
                GridState::gridLevel = gl;
                
                FieldWrap * fieldWrap = FieldHome::CreateField();
                this->data[ solverIndex ][ fid ][ gl ] = fieldWrap;
            }
        }
    }

}

void BasicBgField::Free()
{
    int numberOfSolvers = this->data.size();

    for ( int solverIndex = 0; solverIndex < numberOfSolvers; ++ solverIndex )
    {
        int nFields = this->data[ solverIndex ].size();

        for ( int fid = 0; fid < nFields; ++ fid )
        {
            int nGrids = this->data[ solverIndex ][ fid ].size();

            for ( int gl = 0; gl < nGrids; ++ gl )
            {
                delete this->data[ solverIndex ][ fid ][ gl ];
            }
        }
    }
}

HXVector< BasicBgField * > BgField::data;
bool BgField::flag = false;

BgField::BgField()
{
}

BgField::~BgField()
{
}

void BgField::Init()
{
    if ( BgField::flag ) return;
    BgField::flag = true;

    BgField::data.resize( ZoneState::nZones );

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        ZoneState::zid = iZone;
        BasicBgField * bbgField = new BasicBgField();
        bbgField->Init();
        BgField::data[ iZone ] = bbgField;
    }
}

void BgField::Free()
{
    for ( int iZone = 0; iZone < BgField::data.size(); ++ iZone )
    {
        delete BgField::data[ iZone ];
    }
}

FieldWrap * BgField::GetFieldWrap( int zid, int solverIndex, int fid, int gl )
{
    return BgField::data[ zid ]->data[ solverIndex ][ fid ][ gl ];
}

EndNameSpace
