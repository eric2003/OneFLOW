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

#include "MsgMapImp.h"
#include "Message.h"
#include "FileIO.h"
#include "FileUtil.h"
#include "SimuCtrl.h"
#include <iostream>
using namespace std;

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
    //\tÎªtab¼ü
    string separator  = " =\r\n\t#$,;\"()";
    string msgFileName = SimuCtrl::system_root + "action/" + "actionFileList.txt";

    FileIO ioFile;
    ioFile.OpenFile( msgFileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile()  )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        string fileName = ioFile.ReadNextWord();
        string fullPathFileName = SimuCtrl::system_root + "action/" + fileName;
        fileNameList.push_back( fullPathFileName );
    }

    ioFile.CloseFile();
}


EndNameSpace