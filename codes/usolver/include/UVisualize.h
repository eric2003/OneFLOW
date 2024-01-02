/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "Visualize.h"
#include <sstream>


BeginNameSpace( ONEFLOW )

class DataBook;

class VisualTool
{
public:
    VisualTool();
    ~VisualTool();
public:
    StringField title;
    HXVector< MRField * > qNodeField;
public:
    void Init();
    void AddTitle( const std::string & varName );
    MRField * AddField( RealField & qCellField, const std::string & varName );
    MRField * AddField( const std::string & varName );
    MRField * CreateField( const std::string & varName, int nEqu = 1 );
};

class BcVisual
{
public:
    BcVisual();
    ~BcVisual();
public:
    void Calc( int bcType );
    void Calcf2n( int bcType );
    void ResolveElementEdge();
    void Dump( std::ostringstream & oss, VisualTool * visualTool, std::string & bcTitle );
    void DumpDebug( std::ostringstream & oss, VisualTool * visualTool, std::string & bcTitle );
    void DumpSeveralElement();
public:
    LinkField f2n;
    IntField l2g;

    LinkField e2n;
    IntField lcell;
    IntField rcell;
};

class UVisualize : public Visualize
{
public:
    UVisualize();
    ~UVisualize() override;
public:
    void Visual() override;
    bool NeedVisualField();
    void CalcNodeField( VisualTool * visualTool );
    void ShowField( std::ostringstream & oss, VisualTool * visualTool );
    void ShowBc( std::ostringstream & oss, VisualTool * visualTool );
    void ShowBcDebugTest( std::ostringstream & oss, VisualTool * visualTool );
    void ExtractLinkNum( LinkField & f2n, IntField & fnNumber );
    int  GetTotalNumFaceNodes( LinkField & f2n );
};

void CalcMach( MRField * r, MRField * u, MRField * v, MRField * w, MRField * p, MRField * gama, MRField * mach );

EndNameSpace
