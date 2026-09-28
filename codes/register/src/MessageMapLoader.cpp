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

#include "MessageMapLoader.h"
#include "Message.h"
#include "CmxTaskNames.h"
#include "TextFileParser.h"
#include "Prj.h"
#include "Fatal.h"
#include <iostream>
#include <sstream>

BeginNameSpace( ONEFLOW )

namespace {

const char * const kCmxTaskNameTable[] = {
    kInitFlowFieldTaskName,
    kPostProcessTaskName,
    kStoreRhsTaskName,
    kRestrictAllQTaskName,
    kLoadQTaskName,
    kLoadResidualsTaskName,
    kUpdateResidualsTaskName,
    kRestrictDefectTaskName,
    kModifyCoarseGridTaskName,
    kModifyFineGridTaskName,
    kRecoverCoarseGridTaskName,
    kRecoverResidualsTaskName,
    kZeroResidualsTaskName,
    kCalcTimeStepTaskName,
    kCalcLhsTaskName,
    kUpdateFlowFieldTaskName,
    kCalcBoundaryTaskName,
    kZeroDqFieldTaskName,
    kInitLusgsTaskName,
    kLusgsLowerSweepTaskName,
    kExchangeInterfaceDqTaskName,
    kLusgsUpperSweepTaskName,
    kUpdateFlowFieldLusgsTaskName,
    kSolTurbTaskName,
    kSolHeatTaskName,
    kCalcUnsteadyCriterionTaskName,
    kCalcMetricsTaskName,
    kFillWallStructTaskName,
    kCalcWallDistTaskName,
    kWriteWallDistTaskName,
    kReadWallDistTaskName,
    kAllocateWallDistTaskName,
    kInitFirstTaskName,
    kInitRestartTaskName,
    kReadRestartTaskName,
    kInitInsRestartTaskName,
    kReadInsRestartTaskName,
    kInitFinalTaskName,
    kUploadInterfaceDataTaskName,
    kUpdateInterfaceDataTaskName,
    kDownloadInterfaceDataTaskName,
    kDumpResidualTaskName,
    kDumpAerodynamicTaskName,
    kDumpPressureCoeffTaskName,
    kDumpHeatfluxCoeffTaskName,
    kDumpRestartTaskName,
    kDumpLaminarPlateTaskName,
    kDumpTurbPlateTaskName,
    kVisualizationTaskName,
    kUpdateUnsteadyFlowTaskName,
};

} // namespace

StringField CollectMissingCmxTaskNames()
{
    StringField missing;
    const int n = static_cast< int >( sizeof( kCmxTaskNameTable ) / sizeof( kCmxTaskNameTable[ 0 ] ) );
    for ( int i = 0; i < n; ++ i )
    {
        const char * name = kCmxTaskNameTable[ i ];
        if ( ! MessageMap::Contains( name ) )
        {
            missing.push_back( name );
        }
    }
    return missing;
}

void RequireCmxTaskNamesRegistered()
{
    const StringField missing = CollectMissingCmxTaskNames();
    if ( missing.empty() )
    {
        return;
    }

    std::ostringstream oss;
    oss << "MessageMap is missing "
        << missing.size()
        << " name(s) required by CmxTaskNames.h (source string table):\n";
    for ( std::size_t i = 0; i < missing.size(); ++ i )
    {
        oss << "  - " << missing[ i ] << "\n";
    }
    oss << "Add them to system action message files, or fix CmxTaskNames.h.";
    Fatal( oss.str() );
}


void CreateMsgMap()
{
    StringField fileNameList;
    GetMsgFileNameList( fileNameList );

    MessageMap::Init();

    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        MessageMap::ReadFile( fileNameList[ iFile ] );
    }
    
    // Source (CmxTaskNames) must be covered by runtime MessageMap after load.
    RequireCmxTaskNamesRegistered();
}

void GetMsgFileNameList( StringField & fileNameList )
{
    //\t is the tab key
    std::string separator = " =\r\n\t#$,;\"()";
    std::string msgFileName = Prj::GetSystemFileName( "action/actionFileList.txt" );

    TextFileParser textFileParser;
    textFileParser.OpenFile( msgFileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    // FIX: same class of bug already fixed in MessageMapImp::ReadFile -
    // ReadNextNonEmptyLine() only skips blank lines, not comment lines,
    // so a "# comment" line in actionFileList.txt would have its first
    // token treated as a real file name and appended to fileNameList,
    // causing a spurious failed file open downstream. Driving the loop
    // by ReadNextMeaningfulLine()'s return value skips both blank and
    // comment lines, and correctly signals end-of-file without an extra
    // trailing iteration.
    while ( textFileParser.ReadNextMeaningfulLine() )
    {
        std::string fileName = textFileParser.ReadNextWord();
        if ( fileName.empty() )
        {
            continue; // defensive: skip malformed/empty lines
        }
        std::string fullPathFileName = Prj::GetSystemFileName( "action/" + fileName );
        fileNameList.push_back( fullPathFileName );
    }

    textFileParser.CloseFile();
}

EndNameSpace
