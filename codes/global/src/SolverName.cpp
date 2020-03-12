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
#include "SolverName.h"
#include "SimuCtrl.h"
#include "OStream.h"
#include "FileIO.h"
#include "StrUtil.h"

BeginNameSpace( ONEFLOW )


void GetSolverFileNames( const string & solverName, StringField & fileNameList )
{
    //\tÎªtab¼ü
    string separator = " =\r\n\t#$,;\"()";

    OStream ostr;
    ostr.ClearAll();
    ostr << SimuCtrl::system_root << solverName << "/function/";
    string baseDir = ostr.str();
    ostr << "fileList.txt";
    string keyFileName = ostr.str();

    FileIO ioFile;
    ioFile.OpenFile( keyFileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile()  )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        string fileName = ioFile.ReadNextWord();

        fileName = AddString( baseDir, fileName );

        fileNameList.push_back( fileName );
    }

    ioFile.CloseFile();
}


EndNameSpace