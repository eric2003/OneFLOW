/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ParaFile.h"
#include <string>
#include <vector>

namespace ONEFLOW {

/**
 * @brief Reads the script files listed by the active project's control manifest.
 *
 * Manifest order is preserved. Each script is parsed and committed before
 * the next script is processed; this operation is not transactional.
 */
void ReadPrjScript()
{
    std::vector< std::string > scriptFileNameList;
    ONEFLOW::ReadScriptFileNameList( scriptFileNameList );
    ONEFLOW::ReadMultiScriptFiles( scriptFileNameList );
}

/**
 * @brief Reads a case's manifest and commits its scripts in manifest order.
 *
 * This function owns the composition of manifest discovery and ordered
 * multi-file ingestion. Process coordination, broadcast, and database
 * diagnostics are deliberately handled by ReadControlInfo().
 */
void ReadPrjScript( const std::string & caseDir )
{
    std::vector< std::string > scriptFileNameList;
    ONEFLOW::ReadScriptFileNameList( caseDir, scriptFileNameList );
    ONEFLOW::ReadMultiScriptFiles( scriptFileNameList );
}

} // namespace ONEFLOW
