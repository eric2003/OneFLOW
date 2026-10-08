/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#include "GridWorkflow.h"
#include "GridConversion.h"
#include "DomainInp.h"
#include "Partition.h"

BeginNameSpace( ONEFLOW )

void ExecuteGridConversion(
    const GridConfig & config,
    const std::string & caseDir )
{
    ConvertGrid( config, caseDir );
}

void ExecuteDomainInpWorkflow(
    const GridConfig & /*config*/,
    const std::string & /*caseDir*/ )
{
    DomainInp domainInp;
    domainInp.Run();
}

void ExecutePartitionWorkflow(
    const GridConfig & /*config*/,
    const std::string & /*caseDir*/ )
{
    Partition partition;
    partition.Run();
}

EndNameSpace
