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

#include "MsgMapImp.h"
#include "Message.h"
#include "TextFileParser.h"
#include "Prj.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

void CreateMsgMap()
{
    StringField fileNameList;
    GetMsgFileNameList( fileNameList );

    MessageMap::Init();

    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        MessageMap::ReadFile( fileNameList[ iFile ] );
    }
}

void GetMsgFileNameList( StringField & fileNameList )
{
    //\t is the tab key
    std::string separator  = " =\r\n\t#$,;\"()";
    std::string msgFileName = Prj::system_root + "action/" + "actionFileList.txt";

    TextFileParser textFileParser;
    textFileParser.OpenFile( msgFileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    while ( ! textFileParser.ReachTheEndOfFile()  )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        std::string fileName = textFileParser.ReadNextWord();
        std::string fullPathFileName = Prj::system_root + "action/" + fileName;
        fileNameList.push_back( fullPathFileName );
    }

    textFileParser.CloseFile();
}


EndNameSpace
