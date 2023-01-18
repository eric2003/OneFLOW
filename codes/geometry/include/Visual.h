/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "HXDefine.h"
#include <fstream>


BeginNameSpace( ONEFLOW )

class NodeMesh;
class FaceMesh;
class CellMesh;
class Mesh;

class Visual
{
public:
    Visual();
    ~Visual();
public:
    static int numberOfWords;
public:
    static void SetNumberOfWords( int numberOfWords );
    static void Show( Mesh * mesh );
public:
    static void DumpTitle( std::fstream & file, Mesh * mesh );
    static void DumpCoordinate( std::fstream & file, Mesh * mesh );
    static void DumpCoordinate( std::fstream & file, RealField & coordinate );
    static void DumpFaceNodesLink( std::fstream & file, Mesh * mesh );
    static void DumpFaceElementLink( std::fstream & file, Mesh * mesh );
    static void DumpFaceElementLink( std::fstream & file, HXSize_t nFaces, HXSize_t numberOfElements, IntField & faceElementIndex );
};


EndNameSpace
