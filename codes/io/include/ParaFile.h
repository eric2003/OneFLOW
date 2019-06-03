/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#pragma once
#include "AsciiFileIO.h"
#include "DataBase.h"
#include "DataBook.h"

BeginNameSpace( ONEFLOW )

bool IsArrayParameter( const string & lineOfName );
void ReadBasicData( AsciiFileRead & asciiFileRead );
void AnalysisArrayParameter( AsciiFileRead & asciiFileRead, int keyWordIndex );
int AnalysisScalarParameter( AsciiFileRead & asciiFileRead, int keyWordIndex );
int GetParameterArraySize( const string & word );
void ReadHXFile( const std::string & fileName );

void ReadControlInformation();
void ReadPrjBaseDir();
void ReadHXScript();
void ReadMultiFile();
void BroadcastControlParameterToAllProcessors();

int GetNumberOfParameterFiles();
std::string GetParameterFileName( int iFile = 0 );

void CompressData( DataBase * dataBase, DataBook *& dataBook );
void DecompressData( DataBase * dataBase, DataBook * dataBook );

void CompressData( DataBook *& dataBook );
void DecompressData( DataBook * dataBook );

EndNameSpace