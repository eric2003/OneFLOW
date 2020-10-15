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

#include "Su2Grid.h"
#include "HXStd.h"
#include "ElementHome.h"
#include "CgnsFactory.h"
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
#include "GridPara.h"
#include "Grid.h"
#include "StrGrid.h"
#include "UnsGrid.h"
#include "BgGrid.h"
#include "CalcGrid.h"
#include "HXPointer.h"
#include "NodeMesh.h"
#include "BgGrid.h"
#include "GridState.h"
#include "FileIO.h"
#include "Dimension.h"
#include "HXMath.h"
#include "Zone.h"
#include "BcRecord.h"
#include "DataBase.h"
#include <iostream>

using namespace std;

BeginNameSpace( ONEFLOW )

VTK_CgnsMap vtk_CgnsMap;

int VTK_TYPE::VERTEX        = 1;
int VTK_TYPE::LINE          = 3;
int VTK_TYPE::TRIANGLE      = 5;
int VTK_TYPE::QUADRILATERAL = 9;
int VTK_TYPE::TETRAHEDRON   = 10;
int VTK_TYPE::HEXAHEDRON    = 12;
int VTK_TYPE::PRISM         = 13;
int VTK_TYPE::PYRAMID       = 14;

VTK_TYPE::VTK_TYPE()
{
    ;
}

VTK_TYPE::~VTK_TYPE()
{
    ;
}

VTK_CgnsMap::VTK_CgnsMap()
{
    this->Init();
}

VTK_CgnsMap::~VTK_CgnsMap()
{
    ;
}

void VTK_CgnsMap::Init()
{
    typedef pair< int, int > IntPair;

    vtk2Cgns.insert( IntPair( VTK_TYPE::VERTEX       , NODE    ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::LINE         , BAR_2   ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::TRIANGLE     , TRI_3   ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::QUADRILATERAL, QUAD_4  ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::TETRAHEDRON  , TETRA_4 ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::HEXAHEDRON   , HEXA_8  ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::PRISM        , PENTA_6 ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::PYRAMID      , PYRA_5  ) );
}

SecMarker::SecMarker()
{
    ;
}

SecMarker::~SecMarker()
{
    ;
}

SecMarkerManager::SecMarkerManager()
{
    ;
}

SecMarkerManager::~SecMarkerManager()
{
    int nType = data.size();
    for ( int i = 0; i < nType; ++ i )
    {
        delete data[ i ];
    }
}

void SecMarkerManager::Alloc( int nType )
{
    this->nType = nType;
    data.resize( nType );
    for ( int i = 0; i < nType; ++ i )
    {
        data[ i ] = new SecMarker();
    }
}

int SecMarkerManager::CalcTotalElem()
{
    int nType = data.size();
    int nElem = 0;
    for ( int i = 0; i < nType; ++ i )
    {
        SecMarker * sec = this->data[ i ];
        nElem += sec->nElem;
    }
    return nElem;
}

MarkerManager::MarkerManager()
{
    ;
}

MarkerManager::~MarkerManager()
{
    int nMarker = markerList.size();
    for ( int i = 0; i < nMarker; ++ i )
    {
        delete markerList[ i ];
    }
}

void MarkerManager::CreateMarkerList( int nMarker )
{
    markerList.resize( nMarker );
    for ( int i = 0; i < nMarker; ++ i )
    {
        markerList[ i ] = new Marker();
    }
}

void MarkerManager::CalcSecMarker( SecMarkerManager * secMarkerManager )
{
    int nMarker = this->markerList.size();
    IntSet typeSet;
    for ( int iMarker = 0; iMarker < nMarker; ++ iMarker )
    {
        Marker * marker = this->markerList[ iMarker ];
        for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
        {
            int type = marker->eTypes[ iElem ];
            typeSet.insert( type );
        }
    }

    ONEFLOW::Set2Array( typeSet, types );

    l2g.resize( nMarker );
    for ( int iMarker = 0; iMarker < nMarker; ++ iMarker )
    {
        Marker * marker = this->markerList[ iMarker ];
        l2g[ iMarker ].resize( marker->nElem );
    }

    int nType = types.size();

    secMarkerManager->Alloc( nType );

    int gid = 0;
    for ( int iType = 0; iType < nType; ++ iType )
    {
        int eType = types[ iType ];
        SecMarker * secMarker = secMarkerManager->data[ iType ];
        secMarker->vtk_type = eType;
        secMarker->cgns_type = vtk_CgnsMap.vtk2Cgns[ eType ];
        secMarker->name = ElementTypeName[ secMarker->cgns_type ];
        for ( int iMarker = 0; iMarker < nMarker; ++ iMarker )
        {
            Marker * marker = this->markerList[ iMarker ];
            for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
            {
                int type = marker->eTypes[ iElem ];
                if ( type == eType )
                {
                    secMarker->elems.push_back( marker->elems[ iElem ] );
                    this->l2g[ iMarker ][ iElem ] = gid ++;
                }
            }
        }
        secMarker->nElem = secMarker->elems.size();
    }
    int kkk = 1;

}

VolumeSecManager::VolumeSecManager()
{
    ;
}

VolumeSecManager::~VolumeSecManager()
{
    ;
}

void VolumeSecManager::CalcVolSec( Su2Grid* su2Grid, SecMarkerManager * secMarkerManager )
{
    IntSet typeSet;

    for ( int iElem = 0; iElem < su2Grid->nElem; ++ iElem )
    {
        int eVtkType = su2Grid->elemVTKType[ iElem ];
        typeSet.insert( eVtkType );
    }
    
    ONEFLOW::Set2Array( typeSet, types );

    int nType = types.size();
    l2g.resize( nType );

    secMarkerManager->Alloc( nType );

    for ( int iType = 0; iType < nType; ++ iType )
    {
        int eType = types[ iType ];
        SecMarker * secMarker = secMarkerManager->data[ iType ];
        secMarker->vtk_type = eType;
        secMarker->cgns_type = vtk_CgnsMap.vtk2Cgns[ eType ];
        secMarker->name = ElementTypeName[ secMarker->cgns_type ];
        int gid = 0;
        for ( int iElem = 0; iElem < su2Grid->nElem; ++ iElem )
        {
            int e_VtkType = su2Grid->elemVTKType[ iElem ];
            if ( eType == e_VtkType )
            {
                secMarker->elems.push_back( su2Grid->elems[ iElem ] );
                this->l2g[ iType ].push_back( iElem );
            }
        }
        secMarker->nElem = secMarker->elems.size();
    }
    int kkk = 1;
}

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
    typedef pair< string, int > String2IntPair;
    bcNameToValueMap.insert(String2IntPair("HEATFLUX", BCWall));
    bcNameToValueMap.insert(String2IntPair("FAR", BCFarfield));
}

void Su2Bc::AddBc(string& geoName, string& bcName)
{
    typedef pair< string, string > StringPair;
    bcMap.insert(StringPair(geoName, bcName));
}

void Su2Bc::Process( StringField &markerBCNameList, StringField& markerNameList)
{
    for ( int i = 0; i < markerBCNameList.size(); ++ i )
    {
        string &bcName = markerBCNameList[i];
        if (bcList.find(bcName) != bcList.end())
        {
            string& geoName = markerNameList[i];
            this->AddBc(geoName, bcName);
        }
    }
}

string Su2Bc::GetBcName(string& geoName)
{
    map<string, string>::iterator iter;
    iter = bcMap.find(geoName);
    if (iter!= bcMap.end())
    {
        return iter->second;
    }
    return "";
}

int Su2Bc::GetCgnsBcType(string& geoName)
{
    string bcName = this->GetBcName(geoName);
    return bcNameToValueMap.find(bcName)->second;
}

Su2Grid::Su2Grid()
{
    ndim = 0;
    nPoin = 0;
    nElem = 0;
    nZone = 1;
    vtkmap.resize( MAX_VTK_TYPE, 0 );
    vtkmap[ VTK_TYPE::VERTEX ] = 1;
    vtkmap[ VTK_TYPE::LINE ] = 2;
    vtkmap[ VTK_TYPE::TRIANGLE ] = 3;
    vtkmap[ VTK_TYPE::QUADRILATERAL ] = 4;
    vtkmap[ VTK_TYPE::TETRAHEDRON ] = 4;
    vtkmap[ VTK_TYPE::HEXAHEDRON ] = 8;
    vtkmap[ VTK_TYPE::PRISM ] = 6;
    vtkmap[ VTK_TYPE::PYRAMID ] = 5;
}

Su2Grid::~Su2Grid()
{
}

void Su2Grid::ReadSu2Grid( GridMediator * gridMediator )
{
}

void Su2Grid::ReadSu2GridAscii( string & fileName )
{
    FileIO ioFile;
    string separator  = " =\r\n\t#$,;";
    ioFile.OpenPrjFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    this->nZone = 1;

    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        while ( ! ioFile.ReachTheEndOfFile() )
        {
            ioFile.ReadNextNonEmptyLine();
            string word = ioFile.ReadNextWord();

            if ( word == "NDIME" )
            {
                this->ndim = ioFile.ReadNextDigit< int >();
                continue;
            }
            else if ( word == "NELEM" )
            {
                this->nElem = ioFile.ReadNextDigit< int >();
                for ( int iElem = 0; iElem < this->nElem; ++ iElem )
                {
                    ioFile.ReadNextNonEmptyLine();
                    int vtk_type = ioFile.ReadNextDigit< int >();
                    int nVertex = vtkmap[ vtk_type ];
                    elemVTKType.push_back( vtk_type );
                    IntField elem;
                    for ( int iV = 0; iV < nVertex; ++ iV )
                    {
                        int ip = ioFile.ReadNextDigit< int >();
                        elem.push_back( ip );
                    }
                    this->elems.push_back( elem );
                    int id = ioFile.ReadNextDigit< int >();
                    elemId.push_back( id );
                }
                continue;
            }
            else if ( word == "NPOIN" )
            {
                this->nPoin = ioFile.ReadNextDigit< int >();
                for ( int ip = 0; ip < this->nPoin; ++ ip )
                {
                    ioFile.ReadNextNonEmptyLine();
                    Real xx = ioFile.ReadNextDigit< Real >();
                    Real yy = ioFile.ReadNextDigit< Real >();
                    Real zz = 0;
                    int id = ioFile.ReadNextDigit< int >();
                    this->xN.push_back( xx );
                    this->yN.push_back( yy );
                    this->zN.push_back( zz );
                }
                continue;
            }
            else if ( word == "NMARK" )
            {
                mmark.nMarker = ioFile.ReadNextDigit< int >();
                mmark.CreateMarkerList( mmark.nMarker );

                for ( int im = 0; im < this->mmark.nMarker; ++ im )
                {
                    ioFile.ReadNextNonEmptyLine();
                    string tag = ioFile.ReadNextWord();
                    string name = ioFile.ReadNextWord();
                    Marker * marker = mmark.markerList[ im ];
                    marker->name = name;
                    marker->bcName = su2Bc.GetBcName( name );
                    marker->cgns_bcType = su2Bc.GetCgnsBcType(name);
                    ioFile.ReadNextNonEmptyLine();
                    string marker_elems = ioFile.ReadNextWord();
                    marker->nElem = ioFile.ReadNextDigit< int >();

                    for ( int ielem = 0; ielem < marker->nElem; ++ ielem )
                    {
                        ioFile.ReadNextNonEmptyLine();

                        int eVtkType = ioFile.ReadNextDigit< int >();

                        marker->eTypes.push_back( eVtkType );
                        
                        int eCgnsType = vtk_CgnsMap.vtk2Cgns[ eVtkType ];
                        int nFENode = GetElementNodeNumbers( eCgnsType );

                        IntField elem;
                        for ( int i = 0; i < nFENode; ++ i )
                        {
                            int ip = ioFile.ReadNextDigit< int >();
                            elem.push_back( ip );
                        }
                        marker->elems.push_back( elem );
                    }
                }
                continue;
            }

            int kkk = 1;
        }
    }

    ioFile.CloseFile();
}

void Su2Grid::MarkBoundary( string & su2cfgFile)
{
    FileIO ioFile;
    string separator = " =\r\n\t#$,;()";
    ioFile.OpenPrjFile(su2cfgFile, ios_base::in);
    ioFile.SetDefaultSeparator(separator);

    StringField su2Comment;
    su2Comment.push_back("%");
    ioFile.ResetCommentString( su2Comment );

    StringField markerBCNameList;
    StringField markerNameList;

    while (!ioFile.ReachTheEndOfFile())
    {
        ioFile.ReadNextNonEmptyLine();
        string word = ioFile.ReadNextWord();

        if (word.substr(0, 7) != "MARKER_") continue;
        string bcName = word.substr(7);
        markerBCNameList.push_back(bcName);
        word = ioFile.ReadNextWord();
        markerNameList.push_back(word);
        int kkk = 1;
    }

    su2Bc.Process(markerBCNameList, markerNameList);

    ioFile.CloseFile();
}

void Su2Grid::Su2ToOneFlowGrid()
{
    string gridFile = ONEFLOW::GetDataValue< string >( "sourceGridFileName" );
    string su2cfgFile = ONEFLOW::GetDataValue< string >("sourceGridBcName");
    this->MarkBoundary(su2cfgFile);
    this->ReadSu2GridAscii( gridFile );

    Grids grids( nZone );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsFactory * cgnsFactory = new CgnsFactory();
        int cgnsZoneId = iZone + 1;
        CgnsZone * cgnsZone = cgnsFactory->CreateOneUnsCgnsZone( cgnsZoneId );

        FillSU2CgnsZone( this, cgnsZone );

        cgnsFactory->CgnsToOneFlowGrid( grids[ iZone ], iZone );

        delete cgnsFactory;
    }

    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void FillSU2CgnsZone( Su2Grid* su2Grid, CgnsZone * cgnsZone )
{
    int nNode = su2Grid->xN.size();
    int nCell = su2Grid->nElem;

    cgnsZone->cgnsCoor->SetNNode( nNode );
    cgnsZone->cgnsCoor->SetNCell( nCell );

    NodeMesh * nodeMesh = cgnsZone->cgnsCoor->GetNodeMesh();

    nodeMesh->CreateNodes( nNode );
    nodeMesh->xN = su2Grid->xN;
    nodeMesh->yN = su2Grid->yN;
    nodeMesh->zN = su2Grid->zN;
    
    SecMarkerManager volSec;

    su2Grid->volSec.CalcVolSec( su2Grid, & volSec );
    SecMarkerManager bcSec;
    su2Grid->mmark.CalcSecMarker( &bcSec );

    int nVolSec = volSec.nType;
    int nBcSec = bcSec.nType;

    int nSection = nVolSec + nBcSec;

    CgnsZsection * cgnsZsection = cgnsZone->cgnsZsection;

    cgnsZsection->nSection = nSection;
    cgnsZsection->CreateCgnsSection();

    int nVolCell = volSec.CalcTotalElem();
 
    int sumElem = 0;
    for ( int iSection = 0; iSection < nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = cgnsZsection->GetCgnsSection( iSection );
        SecMarker * sec = 0;
        if ( iSection < nVolSec )
        {
            sec = volSec.data[ iSection ];
        }
        else
        {
            int jSection = iSection - nVolSec;
            sec = bcSec.data[ jSection ];
        }
            

        int nElem = sec->nElem;
        cgnsSection->SetSectionInfo( sec->name, sec->cgns_type, sumElem + 1, sumElem + nElem );
        cgnsSection->CreateConnList();
        sumElem += nElem;

        int pos = 0;
        for ( int iElem = 0; iElem < nElem; ++ iElem )
        {
            IntField & elem = sec->elems[ iElem ];
            int nNode = elem.size();
            for ( int i = 0; i < nNode; ++ i )
            {
                cgnsSection->connList[ pos ++ ]= elem[ i ] + 1;
            }
        }
        int kkk = 1;
    }

    for ( int iSection = 0; iSection < nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = cgnsZsection->GetCgnsSection( iSection );
        cgnsSection->SetElemPosition();
    }

    CgnsZbc * cgnsZbc = cgnsZone->cgnsZbc;
    cgnsZbc->cgnsZbcBoco->ReadZnboco( su2Grid->mmark.nMarker );
    cgnsZbc->cgnsZbcBoco->CreateCgnsZbc();

    for ( int iMarker = 0; iMarker < su2Grid->mmark.nMarker; ++ iMarker )
    {
        Marker * marker = su2Grid->mmark.markerList[ iMarker ];
        string & name = marker->name;
        string& bcName = marker->bcName;

        CgnsBcBoco * cgnsBcBoco = cgnsZbc->cgnsZbcBoco->GetCgnsBc( iMarker );
        cgnsBcBoco->name = name;
        cgnsBcBoco->gridLocation = CellCenter;
        cgnsBcBoco->nElements    = marker->nElem;
        cgnsBcBoco->bcType = static_cast< BCType_t >( marker->cgns_bcType );
        cgnsBcBoco->pointSetType = PointList;
        cgnsBcBoco->CreateCgnsBcConn();

        for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
        {
            int elemId = su2Grid->mmark.l2g[ iMarker ][ iElem ];
            cgnsBcBoco->connList[ iElem ] = elemId + 1 + nVolCell;
        }
        
        //string bcName = GetCgnsBcName( cgnsBcBoco->bcType );
    }

    cgnsZone->ConvertToInnerDataStandard();
    int kkk = 1;
}

EndNameSpace