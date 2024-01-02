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

#include "UCom.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include "PIO.h"
#include "Stop.h"
#include "DataBase.h"
#include "OStream.h"
#include "HXMath.h"
#include "FileUtil.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

UGeom ug;

UGeom::UGeom()
{
    ;
}

UGeom::~UGeom()
{
    ;
}

void UGeom::Init()
{
    //UnsGrid * grid = Zone::GetUnsGrid();
    grid = Zone::GetUnsGrid();

    ug.nBFaces = grid->nBFaces;
    ug.nCells = grid->nCells;
    ug.nTCell = grid->nCells + grid->nBFaces;
    ug.nFaces = grid->nFaces;

    this->SetStEd( F_TOTAL );
    this->CreateBcTypeRegion();

    FaceTopo * faceTopo = grid->faceTopo;
    ug.lcf = & faceTopo->lCells;
    ug.rcf = & faceTopo->rCells;

    FaceMesh * faceMesh = grid->faceMesh;
    CellMesh * cellMesh = grid->cellMesh;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;

    ug.xfn = & faceMesh->xfn;
    ug.yfn = & faceMesh->yfn;
    ug.zfn = & faceMesh->zfn;
    ug.vfn = & faceMesh->vfn;
    ug.farea = & faceMesh->area;

    ug.vfx = & faceMesh->vfx;
    ug.vfy = & faceMesh->vfy;
    ug.vfz = & faceMesh->vfz;

    ug.xfc = & faceMesh->xfc;
    ug.yfc = & faceMesh->yfc;
    ug.zfc = & faceMesh->zfc;

    ug.xcc = & cellMesh->xcc;
    ug.ycc = & cellMesh->ycc;
    ug.zcc = & cellMesh->zcc;

    ug.cvol  = & cellMesh->vol;
    ug.cvol1 = & cellMesh->vol;
    ug.cvol2 = & cellMesh->vol;

    ug.blankf = & cellTopo->blank;

    cellTopo->CalcC2f( faceTopo );

    ug.c2f = & cellTopo->c2f;

    //ug.ireconface = 0;
    ug.ireconface = 1;
}

void UGeom::CreateBcTypeRegion()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
    bcRecord->CreateBcTypeRegion();

    ug.bcRecord = bcRecord;
}

void UGeom::SetStEd( int flag )
{
    if ( flag == F_INNER )
    {
        this->ist = 0;
        this->ied = ug.nCells;
    }
    else if ( flag == F_GHOST )
    {
        this->ist = ug.nCells;
        this->ied = ug.nTCell;
    }
    else if ( flag == F_TOTAL )
    {
        this->ist = 0;
        this->ied = ug.nTCell;
    }
}

void UGeom::DumpCellFace( int cId )
{
    int nFaces = ( * this->c2f )[ cId ].size();
    for ( int fId = 0; fId < nFaces; ++ fId )
    {
        std::cout << ( * this->c2f )[ cId ][ fId ] << " ";
    }
    std::cout << std::endl;
}

void AddF2CField( MRField * cellField, MRField * faceField )
{
    int nEqu = cellField->GetNEqu();
    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
        }
    }

    for ( int fId = ug.nBFaces; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
            ( * cellField )[ iEqu ][ ug.rc ] += ( * faceField )[ iEqu ][ ug.fId ];
        }
    }
}

void AddF2CFieldDebug( MRField * cellField, MRField * faceField )
{
    int nEqu = cellField->GetNEqu();
    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( ug.lc == 9 || ug.rc == 9 )
        {
            std::cout << " fId = " << fId << "\n";
            int iEqu = 0;
            std::cout << ( * faceField )[ iEqu ][ ug.fId ] << "\n";
        }

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
        }
        if ( ug.lc == 9 || ug.rc == 9 )
        {
            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                int cc = 9;
                std::cout << ( * cellField )[ iEqu ][ cc ] << "\n";
            }
        }
    }

    for ( int fId = ug.nBFaces; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];
        if ( ug.lc == 9 || ug.rc == 9 )
        {
            std::cout << " fId = " << fId << "\n";
            int iEqu = 0;
            std::cout << ( * faceField )[ iEqu ][ ug.fId ] << "\n";
        }

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
            ( * cellField )[ iEqu ][ ug.rc ] += ( * faceField )[ iEqu ][ ug.fId ];
        }
        if ( ug.lc == 9 || ug.rc == 9 )
        {
            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                int cc = 9;
                std::cout << ( * cellField )[ iEqu ][ cc ] << "\n";
            }
        }
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        int cc = 9;
        std::cout << ( * cellField )[ iEqu ][ cc ] << "\n";
    }
    int kkk = 1;
}

std::string HXDebug::fileName1;
std::string HXDebug::fileName2;

HXDebug::HXDebug()
{
    ;
}

HXDebug::~HXDebug()
{
    ;
}

std::string HXDebug::GetFullFileName( const std::string & fileName, int startStrategy )
{
    std::string newFileName = fileName;
    if ( startStrategy == 1 )
    {
        newFileName = AddSymbolToFileName( fileName, "Continue" );
    }

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << "debug/" << newFileName;

    std::string fullFileName = ONEFLOW::StrIO.str();
    return fullFileName;
}

void HXDebug::DumpResField( const std::string & fileName )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "res" );
    HXDebug::DumpField( fileName, res );
}

void HXDebug::DumpField( const std::string & fileName, MRField * field )
{
    std::fstream file;

    HXDebug::fileName1 = HXDebug::GetFullFileName( fileName, 0 );
    HXDebug::fileName2 = HXDebug::GetFullFileName( fileName, 1 );

    std::string fullFileName = HXDebug::fileName1;
    int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
    if ( startStrategy == 1 )
    {
        fullFileName = HXDebug::fileName2;
        std::cout << " HXDebug::fileName1 = " << HXDebug::fileName1 << " HXDebug::fileName2 = " << HXDebug::fileName2 << "\n";
    }

    PIO::OpenPrjFile( file, fullFileName, std::ios_base::out );

    int nEqu = field->GetNEqu();

    HXWrite( & file, nEqu );
    int nElems = ( * field )[ 0 ].size();
    HXWrite( & file, nElems );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXWrite( &file, ( * field )[ iEqu ] );
    }

    PIO::CloseFile( file );
}

MRField * HXDebug::ReadField( const std::string & fileName )
{
    std::fstream file;

    PIO::OpenPrjFile( file, fileName, std::ios_base::in );
    int nEqu = 0;
    HXRead( & file, nEqu );
    int nCells = 0;
    HXRead( & file, nCells );

    std::cout << " nEqu = " << nEqu << " nCells = " << nCells << "\n";

    MRField * field = new MRField( nEqu, nCells );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXRead( &file, ( * field )[ iEqu ] );
    }

    PIO::CloseFile( file );

    return field;
}

void HXDebug::CompareFile( Real mindiff, int idump )
{
    int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
    if ( startStrategy != 1 ) return;

    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * field1 = HXDebug::ReadField( HXDebug::fileName1 );
    MRField * field2 = HXDebug::ReadField( HXDebug::fileName2 );
    int nEqu = field1->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXVector< Real > & f1 = ( * field1 )[ iEqu ];
        HXVector< Real > & f2 = ( * field2 )[ iEqu ];
        int nCells = f1.size();
        for ( int iCell = 0; iCell < nCells; ++ iCell )
        {
            Real diff = f1[ iCell ] - f2[ iCell ];
            if ( ABS( diff ) > mindiff )
            {
                std::cout << " iEqu = " << iEqu << " iCell = " << iCell << " diff = " << diff << "\n";
                std::cout << " v1 = " << std::setiosflags( std::ios::fixed ) << std::setprecision( 10 ) << f1[ iCell ];
                std::cout << " v2 = " << std::setiosflags( std::ios::fixed ) << std::setprecision( 10 ) << f2[ iCell ] << "\n";
                if ( idump == 0 )
                {
                    HXDebug::DumpCellInfo( iCell );
                }
                break;
            }
        }
    }

    delete field1;
    delete field2;
}

void HXDebug::DumpCellInfo( int iCell )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nFaces = ( * ug.c2f )[ iCell ].size();
    int nBFaces = ( * ug.bcRecord ).bcType.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int fid = ( * ug.c2f )[ iCell ][ iFace ];
        int bctype = -1;
        if ( fid < nBFaces )
        {
            bctype = ( * ug.bcRecord ).bcType[ fid ];
        }
        std::cout << " fid = " << fid << " bctype = " << bctype << "\n";
    }
}

void HXDebug::CheckNANField( MRField * field )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    int nEqu = field->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        int nElems = ( * field )[ iEqu ].size();
        for ( int iElem = 0; iElem < nElems; ++ iElem )
        {
            Real value = ( * field )[ iEqu ][ iElem ];
            if ( NotANumber( value ) )
            {
                std::cout << " iEqu = " << iEqu << " iElem = " << iElem << " nElems = " << nElems << " value = " << value << "\n";
                Stop( "NotANumber" );
            }
        }
    }
}

EndNameSpace
