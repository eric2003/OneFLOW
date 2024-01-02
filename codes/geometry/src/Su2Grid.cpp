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
    typedef std::pair< int, int > IntPair;

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
    typedef std::pair< std::string, int > String2IntPair;
    bcNameToValueMap.insert(String2IntPair("HEATFLUX", BCWall));
    bcNameToValueMap.insert(String2IntPair("FAR", BCFarfield));
}

void Su2Bc::AddBc(std::string& geoName, std::string& bcName)
{
    typedef std::pair< std::string, std::string > StringPair;
    bcMap.insert(StringPair(geoName, bcName));
}

void Su2Bc::Process( StringField &markerBCNameList, StringField& markerNameList)
{
    for ( int i = 0; i < markerBCNameList.size(); ++ i )
    {
        std::string &bcName = markerBCNameList[i];
        if (bcList.find(bcName) != bcList.end())
        {
            std::string& geoName = markerNameList[i];
            this->AddBc(geoName, bcName);
        }
    }
}

std::string Su2Bc::GetBcName(std::string& geoName)
{
    std::map<std::string, std::string>::iterator iter;
    iter = bcMap.find(geoName);
    if (iter!= bcMap.end())
    {
        return iter->second;
    }
    return "";
}

int Su2Bc::GetCgnsBcType(std::string& geoName)
{
    std::string bcName = this->GetBcName(geoName);
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

void Su2Grid::ReadSu2GridAscii( std::string & fileName )
{
    FileIO ioFile;
    std::string separator  = " =\r\n\t#$,;";
    ioFile.OpenPrjFile( fileName, std::ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    this->nZone = 1;

    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        while ( ! ioFile.ReachTheEndOfFile() )
        {
            ioFile.ReadNextNonEmptyLine();
            std::string word = ioFile.ReadNextWord();

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
                    //int id = ioFile.ReadNextDigit< int >();
                    int id = iElem + 1;
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
                    std::string tag = ioFile.ReadNextWord();
                    std::string name = ioFile.ReadNextWord();
                    Marker * marker = mmark.markerList[ im ];
                    marker->name = name;
                    marker->bcName = su2Bc.GetBcName( name );
                    marker->cgns_bcType = su2Bc.GetCgnsBcType(name);
                    ioFile.ReadNextNonEmptyLine();
                    std::string marker_elems = ioFile.ReadNextWord();
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

void Su2Grid::MarkBoundary( std::string & su2cfgFile)
{
    FileIO ioFile;
    std::string separator = " =\r\n\t#$,;()";
    ioFile.OpenPrjFile(su2cfgFile, std::ios_base::in);
    ioFile.SetDefaultSeparator(separator);

    StringField su2Comment;
    su2Comment.push_back("%");
    ioFile.ResetCommentString( su2Comment );

    StringField markerBCNameList;
    StringField markerNameList;

    while (!ioFile.ReachTheEndOfFile())
    {
        ioFile.ReadNextNonEmptyLine();
        std::string word = ioFile.ReadNextWord();

        if (word.substr(0, 7) != "MARKER_") continue;
        std::string bcName = word.substr(7);
        markerBCNameList.push_back(bcName);
        word = ioFile.ReadNextWord();
        markerNameList.push_back(word);
        int kkk = 1;
    }

    su2Bc.Process( markerBCNameList, markerNameList );

    ioFile.CloseFile();
}

void Su2Grid::Su2ToOneFlowGrid()
{
    std::string gridFile = ONEFLOW::GetDataValue< std::string >( "sourceGridFileName" );
    std::string su2cfgFile = ONEFLOW::GetDataValue< std::string >("sourceGridBcName");
    this->MarkBoundary(su2cfgFile);
    this->ReadSu2GridAscii( gridFile );

    ONEFLOW::Su2ToOneFlowGrid( this );
}

void Su2Grid::FillSU2CgnsZone( CgnsZone * cgnsZone )
{
    int nNodes = this->xN.size();
    int nCells = this->nElem;

    cgnsZone->cgnsCoor->SetNNode( nNodes );
    cgnsZone->cgnsCoor->SetNCell( nCells );

    NodeMesh * nodeMesh = cgnsZone->cgnsCoor->GetNodeMesh();

    nodeMesh->CreateNodes( nNodes );
    nodeMesh->xN = this->xN;
    nodeMesh->yN = this->yN;
    nodeMesh->zN = this->zN;
    
    SecMarkerManager volSec;

    this->volSec.CalcVolSec( this, & volSec );
    SecMarkerManager bcSec;
    this->mmark.CalcSecMarker( &bcSec );

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
            int nNodes = elem.size();
            for ( int i = 0; i < nNodes; ++ i )
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
    cgnsZbc->cgnsZbcBoco->ReadZnboco( this->mmark.nMarker );
    cgnsZbc->cgnsZbcBoco->CreateCgnsZbc();

    for ( int iMarker = 0; iMarker < this->mmark.nMarker; ++ iMarker )
    {
        Marker * marker = this->mmark.markerList[ iMarker ];
        std::string & name = marker->name;
        std::string& bcName = marker->bcName;

        CgnsBcBoco * cgnsBcBoco = cgnsZbc->cgnsZbcBoco->GetCgnsBc( iMarker );
        cgnsBcBoco->name = name;
        cgnsBcBoco->gridLocation = CellCenter;
        cgnsBcBoco->nElements    = marker->nElem;
        cgnsBcBoco->bcType = static_cast< BCType_t >( marker->cgns_bcType );
        cgnsBcBoco->pointSetType = PointList;
        cgnsBcBoco->CreateCgnsBcBoco();

        for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
        {
            int elemId = this->mmark.l2g[ iMarker ][ iElem ];
            cgnsBcBoco->connList[ iElem ] = elemId + 1 + nVolCell;
        }
        
        //string bcName = GetCgnsBcName( cgnsBcBoco->bcType );
    }

    cgnsZone->cgnsZoneType = Unstructured;

    cgnsZone->ConvertToInnerDataStandard();
    int kkk = 1;
}

void Su2ToOneFlowGrid( Su2Grid* su2Grid )
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    cgnsFactory->Su2ToOneFlowGrid( su2Grid );

    delete cgnsFactory;
}

EndNameSpace
