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

#include "Su2Grid.h"
#include "GridTypes.h"
#include "HXStd.h"
#include "ElementHome.h"
#include "CgnsFactory.h"
#include "UnstructuredCgnsZoneWriter.h"
#include "CgnsZbase.h"
#include "CgnsZsection.h"
#include "CgnsZbc.h"
#include "CgnsZbcBoco.h"
#include "CgnsBcBoco.h"
#include "CgnsSection.h"
#include "CgnsZone.h"
#include "CgnsCoor.h"
#include "GridMediator.h"
#include "Boundary.h"
#include "Grid.h"
#include "StrGrid.h"
#include "UnsGrid.h"
#include "CalcGrid.h"
#include "NodeMesh.h"
#include "GridState.h"
#include "TextFileParser.h"
#include "Dimension.h"
#include "HXMath.h"
#include "BcRecord.h"
#include <iostream>



BeginNameSpace( ONEFLOW )

Su2Bc::Su2Bc()
{
    Init();
}

Su2Bc::~Su2Bc()
{
    ;
}

void Su2Bc::Init()
{
    bcList.insert("HEATFLUX");
    bcList.insert("FAR");
    using String2IntPair = std::pair< std::string, int >;
    bcNameToValueMap.insert(String2IntPair("HEATFLUX", BCWall));
    bcNameToValueMap.insert(String2IntPair("FAR", BCFarfield));
}

void Su2Bc::AddBc( const std::string & geoName, const std::string & bcName )
{
    using StringPair = std::pair< std::string, std::string >;
    bcMap.insert(StringPair(geoName, bcName));
}

void Su2Bc::Process( const StringField & markerBCNameList, const StringField & markerNameList )
{
    for ( int i = 0; i < markerBCNameList.size(); ++ i )
    {
        const std::string & bcName = markerBCNameList[ i ];
        if (bcList.find(bcName) != bcList.end())
        {
            const std::string & geoName = markerNameList[ i ];
            this->AddBc(geoName, bcName);
        }
    }
}

std::string Su2Bc::GetBcName( const std::string & geoName ) const
{
    std::map<std::string, std::string>::const_iterator iter;
    iter = bcMap.find(geoName);
    if (iter!= bcMap.end())
    {
        return iter->second;
    }
    return "";
}

int Su2Bc::GetCgnsBcType( const std::string & geoName ) const
{
    const std::string bcName = this->GetBcName(geoName);
    return bcNameToValueMap.find(bcName)->second;
}

Su2Grid::Su2Grid()
    : nZone(1)
{
    vtkmap.assign(MAX_VTK_TYPE, 0);
    vtkmap[VTK_TYPE::VERTEX] = 1;
    vtkmap[VTK_TYPE::LINE] = 2;
    vtkmap[VTK_TYPE::TRIANGLE] = 3;
    vtkmap[VTK_TYPE::QUADRILATERAL] = 4;
    vtkmap[VTK_TYPE::TETRAHEDRON] = 4;
    vtkmap[VTK_TYPE::HEXAHEDRON] = 8;
    vtkmap[VTK_TYPE::PRISM] = 6;
    vtkmap[VTK_TYPE::PYRAMID] = 5;
}

void Su2Grid::ReadSu2Grid( GridMediator * gridMediator )
{
}

void Su2Grid::ReadSu2GridAscii( const std::string & fileName, const std::string & caseDir )
{
    TextFileParser textFileParser;
    std::string separator  = " =\r\n\t#$,;";
    if ( ! caseDir.empty() )
    {
        textFileParser.OpenCaseFile(
            caseDir, fileName, std::ios_base::in );
    }
    else
    {
        textFileParser.OpenPrjFile( fileName, std::ios_base::in );
    }
    textFileParser.SetDefaultSeparator( separator );

    this->nZone = 1;

    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        while ( ! textFileParser.ReachTheEndOfFile() )
        {
            textFileParser.ReadNextNonEmptyLine();
            std::string word = textFileParser.ReadNextWord();

            if ( word == "NDIME" )
            {
                this->ndim = textFileParser.ReadNextDigit< int >();
                continue;
            }
            else if ( word == "NELEM" )
            {
                this->nElem = textFileParser.ReadNextDigit< int >();

                for ( int iElem = 0; iElem < this->nElem; ++ iElem )
                {
                    textFileParser.ReadNextNonEmptyLine();
                    int vtk_type = textFileParser.ReadNextDigit< int >();
                    int nVertex = vtkmap[ vtk_type ];
                    elemVTKType.push_back( vtk_type );
                    IntField elem;
                    for ( int iV = 0; iV < nVertex; ++ iV )
                    {
                        int ip = textFileParser.ReadNextDigit< int >();
                        elem.push_back( ip );
                    }
                    this->elems.push_back( elem );
                    //int id = textFileParser.ReadNextDigit< int >();
                    int id = iElem + 1;
                    elemId.push_back( id );
                }
                continue;
            }
            else if ( word == "NPOIN" )
            {
                this->nPoin = textFileParser.ReadNextDigit< int >();
                for ( int ip = 0; ip < this->nPoin; ++ ip )
                {
                    textFileParser.ReadNextNonEmptyLine();
                    Real xx = textFileParser.ReadNextDigit< Real >();
                    Real yy = textFileParser.ReadNextDigit< Real >();
                    Real zz = 0;
                    int id = textFileParser.ReadNextDigit< int >();
                    this->xN.push_back( xx );
                    this->yN.push_back( yy );
                    this->zN.push_back( zz );
                }
                continue;
            }
            else if ( word == "NMARK" )
            {
                mmark.nMarker = textFileParser.ReadNextDigit< int >();
                mmark.CreateMarkerList( mmark.nMarker );

                for ( int im = 0; im < this->mmark.nMarker; ++ im )
                {
                    textFileParser.ReadNextNonEmptyLine();
                    std::string tag = textFileParser.ReadNextWord();
                    std::string name = textFileParser.ReadNextWord();
                    Marker * marker = & mmark.markerList[ im ];
                    marker->name = name;
                    marker->bcName = su2Bc.GetBcName( name );
                    marker->cgns_bcType = su2Bc.GetCgnsBcType(name);
                    textFileParser.ReadNextNonEmptyLine();
                    std::string marker_elems = textFileParser.ReadNextWord();
                    marker->nElem = textFileParser.ReadNextDigit< int >();

                    for ( int ielem = 0; ielem < marker->nElem; ++ ielem )
                    {
                        textFileParser.ReadNextNonEmptyLine();

                        int eVtkType = textFileParser.ReadNextDigit< int >();

                        marker->eTypes.push_back( eVtkType );
                        
                        int eCgnsType = vtk_CgnsMap.vtk2Cgns[ eVtkType ];
                        int nFENode = GetElementNodeNumbers( eCgnsType );

                        IntField elem;
                        for ( int i = 0; i < nFENode; ++ i )
                        {
                            int ip = textFileParser.ReadNextDigit< int >();
                            elem.push_back( ip );
                        }
                        marker->elems.push_back( elem );
                    }
                }
                continue;
            }


        }
    }

    textFileParser.CloseFile();
}

void Su2Grid::MarkBoundary( const std::string & su2cfgFile, const std::string & caseDir )
{
    TextFileParser textFileParser;
    std::string separator = " =\r\n\t#$,;()";
    if ( ! caseDir.empty() )
    {
        textFileParser.OpenCaseFile(
            caseDir, su2cfgFile, std::ios_base::in );
    }
    else
    {
        textFileParser.OpenPrjFile( su2cfgFile, std::ios_base::in );
    }
    textFileParser.SetDefaultSeparator(separator);

    StringField su2Comment;
    su2Comment.push_back("%");
    textFileParser.ResetCommentString( su2Comment );

    StringField markerBCNameList;
    StringField markerNameList;

    while (!textFileParser.ReachTheEndOfFile())
    {
        textFileParser.ReadNextNonEmptyLine();
        std::string word = textFileParser.ReadNextWord();

        if (word.substr(0, 7) != "MARKER_") continue;
        std::string bcName = word.substr(7);
        markerBCNameList.push_back(bcName);
        word = textFileParser.ReadNextWord();
        markerNameList.push_back(word);

    }

    su2Bc.Process( markerBCNameList, markerNameList );

    textFileParser.CloseFile();
}

void Su2Grid::Su2ToOneFlowGrid( const GridConfig & config, const std::string & caseDir )
{
    this->MarkBoundary( config.bcFile, caseDir );
    this->ReadSu2GridAscii( config.sourceFile, caseDir );

    ONEFLOW::Su2ToOneFlowGrid( *this );
}

void Su2Grid::FillSU2CgnsZone(
    CgnsZone &cgnsZone,
    bool convertToInnerDataStandard)
{
    UnstructuredCgnsZoneWriter::Write(
        *this, cgnsZone, convertToInnerDataStandard);
}

void Su2ToOneFlowGrid( Su2Grid & su2Grid )
{
    CgnsFactory cgnsFactory;
    cgnsFactory.Su2ToOneFlowGrid( su2Grid );
}

EndNameSpace
