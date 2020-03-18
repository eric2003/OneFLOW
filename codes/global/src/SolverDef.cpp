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
#include "SolverDef.h"

BeginNameSpace( ONEFLOW )

map< string, int > * solverTypeMap = 0;
map< string, int > * interfaceMap = 0;
map< string, int > * sendRecvMap = 0;
map< string, int > * fieldIdMap = 0;

void CreateSysMap()
{
    solverTypeMap = new map< string, int >;
    interfaceMap = new map< string, int >;
    sendRecvMap = new map< string, int >;
    fieldIdMap = new map< string, int >;

    ( * solverTypeMap )[ "CFD_SOLVER"        ] = ONEFLOW::CFD_SOLVER;
    ( * solverTypeMap )[ "NS_SOLVER"         ] = ONEFLOW::NS_SOLVER;
    ( * solverTypeMap )[ "TURB_SOLVER"       ] = ONEFLOW::TURB_SOLVER;
    ( * solverTypeMap )[ "INC_NS_SOLVER"     ] = ONEFLOW::INC_NS_SOLVER;
    ( * solverTypeMap )[ "HIGH_ORDER_SOLVER" ] = ONEFLOW::HIGH_ORDER_SOLVER;
    ( * solverTypeMap )[ "SCALAR_SOLVER"     ] = ONEFLOW::SCALAR_SOLVER;
    ( * solverTypeMap )[ "GRID_SOLVER"       ] = ONEFLOW::GRID_SOLVER;
    ( * solverTypeMap )[ "SIMU_SOLVER"       ] = ONEFLOW::SIMU_SOLVER;

    ( * interfaceMap )[ "INTERFACE_DATA"          ] = ONEFLOW::INTERFACE_DATA;
    ( * interfaceMap )[ "INTERFACE_DQ_DATA"       ] = ONEFLOW::INTERFACE_DQ_DATA;
    ( * interfaceMap )[ "INTERFACE_GRADIENT_DATA" ] = ONEFLOW::INTERFACE_GRADIENT_DATA;
    ( * interfaceMap )[ "INTERFACE_OVERSET_DATA"  ] = ONEFLOW::INTERFACE_OVERSET_DATA;

    ( * sendRecvMap )[ "GREAT_SEND" ] = ONEFLOW::GREAT_SEND;
    ( * sendRecvMap )[ "GREAT_RECV" ] = ONEFLOW::GREAT_RECV;

    ( * fieldIdMap )[ "FIELD_FLOW" ] = ONEFLOW::FIELD_FLOW;
    ( * fieldIdMap )[ "FIELD_RHS"  ] = ONEFLOW::FIELD_RHS;
}

int GetOppositeSendRecv( int iSr )
{
    if ( iSr == GREAT_SEND )
    {
        return GREAT_RECV;
    }
    else
    {
        return GREAT_SEND;
    }
}

EndNameSpace