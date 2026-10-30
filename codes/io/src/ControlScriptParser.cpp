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

#include "ParaFile.h"
#include "ConfigLoader.h"
#include "ConfigDatabaseAdapter.h"
#include "TextFileParser.h"
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

void ReadOneFLOWScriptFile( TextFileParser & textFileParser )
{
    ConfigLoader loader( GetParameterArraySize );
    loader.ParseFromParser( textFileParser );
    ConfigDatabaseAdapter::Commit( loader.Document() );
}

void ReadOneFLOWScriptFile( const std::string & fileName )
{
    ConfigLoader loader( GetParameterArraySize );
    loader.ParseFile( fileName );
    ConfigDatabaseAdapter::Commit( loader.Document() );
}

void ReadMultiScriptFiles( const std::vector< std::string > & scriptFileNameList )
{
    for ( const std::string & scriptFileName : scriptFileNameList )
    {
        ONEFLOW::ReadOneFLOWScriptFile( scriptFileName );
    }
}

EndNameSpace
