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
#include "UsdFieldConfig.h"
#include "Fatal.h"

BeginNameSpace( ONEFLOW )

UsdFieldConfig::UsdFieldConfig()
{
}

void UsdFieldConfig::Init(
    const UsdFieldNames & fieldNames )
{
    this->fieldNames = fieldNames;
}

const UsdFieldNames & UsdFieldConfig::GetFieldNames() const
{
    return this->fieldNames;
}

std::map<
    int,
    std::unique_ptr< UsdFieldConfig > >
    UsdFieldConfigRegistry::data;

void UsdFieldConfigRegistry::AddConfig(
    int solverType )
{
    if ( UsdFieldConfigRegistry::data.find(
        solverType ) !=
        UsdFieldConfigRegistry::data.end() )
    {
        return;
    }

    UsdFieldConfigRegistry::data[ solverType ] =
        std::make_unique< UsdFieldConfig >();
}

UsdFieldConfig * UsdFieldConfigRegistry::GetConfig(
    int solverType )
{
    auto iter =
        UsdFieldConfigRegistry::data.find(
            solverType );

    if ( iter ==
        UsdFieldConfigRegistry::data.end() )
    {
        return nullptr;
    }

    return iter->second.get();
}

void UsdFieldConfigRegistry::SetConfig(
    int solverType,
    const UsdFieldNames & fieldNames )
{
    AddConfig( solverType );

    UsdFieldConfig * config =
        GetConfig( solverType );

    config->Init( fieldNames );
}

void UsdFieldConfigRegistry::FreeConfig()
{
    UsdFieldConfigRegistry::data.clear();
}

EndNameSpace