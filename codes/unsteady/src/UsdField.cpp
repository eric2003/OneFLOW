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

#include "UsdField.h"
#include "UsdPara.h"
#include "UsdFieldConfig.h"
#include "FieldManager.h"
#include "FieldWrap.h"
#include "DataBase.h"
#include "Zone.h"
#include "Fatal.h"
#include "UnsGrid.h"

BeginNameSpace( ONEFLOW )


UsdField::UsdField()
{
}

UsdField::~UsdField()
{
}

void UsdField::Init()
{
    ;
}

void UsdField::InitBasic( int solverType )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    UsdFieldConfig * config =
        UsdFieldConfigRegistry::GetConfig(
            solverType );

    StringField flow;
    StringField residual;

    if ( config != nullptr )
    {
        const UsdFieldNames & fieldNames =
            config->GetFieldNames();

        flow =
            fieldNames.flow;

        residual =
            fieldNames.residual;
    }
    else
    {
        FieldManager * fieldManager =
            FieldManagerRegistry::GetFieldManager(
                solverType );

        if ( fieldManager == nullptr )
        {
            Fatal(
                "FieldManager is not registered for solverType" );
        }

        UsdPara * usdPara =
            &fieldManager->GetUsdPara();

        flow =
            usdPara->flow;

        residual =
            usdPara->residual;
    }

    this->flow.resize(
        flow.size() );

    for ( std::size_t i = 0;
        i < flow.size();
        ++ i )
    {
        this->flow[ i ] =
            GetFieldPointer< MRField >(
                grid,
                flow[ i ] );
    }

    this->residual.resize(
        residual.size() );

    if ( this->flow.size() < 3 )
    {
        Fatal(
            "Unsteady flow fields require at least 3 time levels." );
    }

    if ( this->residual.size() < 3 )
    {
        Fatal(
            "Unsteady residual fields require at least 3 time levels." );
    }

    for ( std::size_t i = 0;
        i < residual.size();
        ++ i )
    {
        this->residual[ i ] =
            GetFieldPointer< MRField >(
                grid,
                residual[ i ] );
    }
}

MRField * UsdField::GetFlow( HistoryLevel level )
{
    return this->flow[ static_cast< std::size_t >( level ) ];
}

MRField * UsdField::GetResidual( HistoryLevel level )
{
    return this->residual[ static_cast< std::size_t >( level ) ];
}


EndNameSpace
