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

#include "Plot3D.h"
#include "CgnsFactory.h"
#include "GridMediator.h"
#include "Stop.h"
#include "Prj.h"
#include "StrGrid.h"
#include "NodeMesh.h"
#include "BgGrid.h"
#include "GridState.h"
#include "FileIO.h"
#include "FileO.h"

#include "Dimension.h"
#include "HXMath.h"
#include "Zone.h"
#include "ZoneState.h"
#include "BcRecord.h"
#include "DataBase.h"
#include "GridPara.h"
#include <iostream>



BeginNameSpace( ONEFLOW )

Plot3D::Plot3D()
{
}

Plot3D::~Plot3D()
{
}

void Plot3D::ReadPlot3D( GridMediator * gridMediator )
{
    Plot3D::ReadCoor( gridMediator );
    Plot3D::ReadBc( gridMediator );
}

void Plot3D::ReadCoor( GridMediator * gridMediator )
{
    int fileFormat = GetDataValue< int >( "fileFormat" );

    if ( fileFormat == BINARY )
    {
        Plot3D::ReadCoorBinary( gridMediator );
    }
    else
    {
        Plot3D::ReadCoorAscii( gridMediator );
    }
}

void Plot3D::DumpCoor( GridMediator * gridMediator )
{
    int fileFormat = GetDataValue< int >( "fileFormat" );

    if ( fileFormat == BINARY )
    {
        Plot3D::DumpCoorBinary( gridMediator );
    }
    else
    {
        Plot3D::DumpCoorAscii( gridMediator );
    }
}

void Plot3D::ReadCoorBinary( GridMediator * gridMediator )
{
    std::string & fileName = gridMediator->gridFile;
    std::string separator  = " =\r\n\t#$,;";

    std::fstream file;
    Prj::OpenPrjFile( file, fileName, std::ios_base::in|std::ios_base::binary );

    HXRead( & file, gridMediator->numberOfZones );
    gridMediator->gridVector.resize( gridMediator->numberOfZones );

    bool nkflag = GetPlot3D_NKFlag();

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        int ni, nj, nk = 1;
        HXRead( & file, ni );
        HXRead( & file, nj );

        if ( nkflag )
        {
            HXRead( & file, nk );
        }

        Grid * gridstr = ONEFLOW::CreateGrid( ONEFLOW::SMESH );
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        gridMediator->gridVector[ iZone ] = grid;
        grid->id = iZone;
        grid->ni = ni;
        grid->nj = nj;
        grid->nk = nk;
        grid->SetBasicDimension();
        grid->nodeMesh->CreateNodes( grid->nNodes );

        grid->SetLayout();

        int wordWidth = 8;

        std::cout << std::setiosflags( std::ios::right );
        std::cout << " iZone = " << std::setw( wordWidth ) << iZone;
        std::cout << " nZone = " << std::setw( wordWidth ) << gridMediator->numberOfZones;
        std::cout << " ni, nj, nk = ";
        std::cout << std::setw( wordWidth ) << ni;
        std::cout << std::setw( wordWidth ) << nj;
        std::cout << std::setw( wordWidth ) << nk;
        std::cout << std::endl;
    }

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );
        int numberOfNodes = grid->nNodes;

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        HXRead( & file, grid->nodeMesh->xN );
        HXRead( & file, grid->nodeMesh->yN );

        if ( nkflag )
        {
            HXRead( & file, grid->nodeMesh->zN );
        }
        else
        {
            grid->nodeMesh->zN = 0;
        }
    }

    Prj::CloseFile( file );
}

void Plot3D::DumpCoorBinary( GridMediator * gridMediator )
{
    std::string & fileName = gridMediator->gridFile;

    std::fstream file;
    Prj::OpenPrjFile( file, fileName, std::ios_base::out|std::ios_base::binary );

    int numberOfZones = gridMediator->numberOfZones;
    HXWrite( & file, numberOfZones );

    bool nkflag = GetPlot3D_NKFlag();

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        HXWrite( & file, ni );
        HXWrite( & file, nj );

        if ( nkflag )
        {
            HXWrite( & file, nk );
        }

        int wordWidth = 8;

        std::cout << std::setiosflags( std::ios::right );
        std::cout << " iZone = " << std::setw( wordWidth ) << iZone;
        std::cout << " nZone = " << std::setw( wordWidth ) << gridMediator->numberOfZones;
        std::cout << " ni, nj, nk = ";
        std::cout << std::setw( wordWidth ) << ni;
        std::cout << std::setw( wordWidth ) << nj;
        std::cout << std::setw( wordWidth ) << nk;
        std::cout << "\n";
    }

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        HXWrite( & file, grid->nodeMesh->xN );
        HXWrite( & file, grid->nodeMesh->yN );

        if ( nkflag )
        {
            HXWrite( & file, grid->nodeMesh->zN );
        }
    }

    Prj::CloseFile( file );
}

void Plot3D::ReadCoorAscii( GridMediator * gridMediator )
{
    std::string & fileName = gridMediator->gridFile;

    FileIO ioFile;
    std::string separator  = " =\r\n\t#$,;";
    ioFile.OpenPrjFile( fileName, std::ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    gridMediator->numberOfZones = ioFile.ReadNextDigit< int >();
    gridMediator->gridVector.resize( gridMediator->numberOfZones );

    bool nkflag = GetPlot3D_NKFlag();

    int zCount = 0;
    while ( zCount < gridMediator->numberOfZones )
    {
        int ni, nj, nk = 1;
        ni = ioFile.ReadNextDigit< int >();
        nj = ioFile.ReadNextDigit< int >();

        if ( nkflag )
        {
            nk = ioFile.ReadNextDigit< int >();
        }

        Grid * gridstr = ONEFLOW::CreateGrid( ONEFLOW::SMESH );
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        gridMediator->gridVector[ zCount ] = grid;
        grid->id = zCount;
        grid->ni = ni;
        grid->nj = nj;
        grid->nk = nk;
        grid->SetBasicDimension();
        grid->nodeMesh->CreateNodes( grid->nNodes );

        grid->SetLayout();

        int wordWidth = 8;

        std::cout << std::setiosflags( std::ios::right );
        std::cout << " iZone = " << std::setw( wordWidth ) << zCount;
        std::cout << " nZone = " << std::setw( wordWidth ) << gridMediator->numberOfZones;
        std::cout << " ni, nj, nk = ";
        std::cout << std::setw( wordWidth ) << ni;
        std::cout << std::setw( wordWidth ) << nj;
        std::cout << std::setw( wordWidth ) << nk;
        std::cout << std::endl;
        ++ zCount;
    }

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );
        int numberOfNodes = grid->nNodes;

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        int total_size = 2 * numberOfNodes;

        if ( nkflag )
        {
            total_size = 3 * numberOfNodes;
        }

        RealField coor;

        Plot3D::ReadCoor( & ioFile, coor, total_size );
        int pos = 0;
        for ( int i = 0; i < numberOfNodes; ++ i )
        {
            grid->nodeMesh->xN[ i ] = coor[ i + pos ];
        }
        pos += numberOfNodes;

        for ( int i = 0; i < numberOfNodes; ++ i )
        {
            grid->nodeMesh->yN[ i ] = coor[ i + pos ];
        }
        pos += numberOfNodes;

        if ( nkflag )
        {
            for ( int i = 0; i < numberOfNodes; ++ i )
            {
                grid->nodeMesh->zN[ i ] = coor[ i + pos ];
            }
            pos += numberOfNodes;
        }
        else
        {
            grid->nodeMesh->zN = 0;
        }
    }

    ioFile.CloseFile();
}

void Plot3D::DumpCoorAscii( GridMediator * gridMediator )
{
    std::string & fileName = gridMediator->gridFile;

    FileO fileO;
    fileO.OpenPrjFile( fileName, std::ios_base::out );

    int numberOfZones = gridMediator->numberOfZones;
    fileO.WriteLine( numberOfZones );

    bool nkflag = GetPlot3D_NKFlag();

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        fileO.Write( ni );
        fileO.Write( nj );
        if ( nkflag )
        {
            fileO.Write( nk );
        }
        fileO.WriteEndLine();

        int wordWidth = 8;
        std::cout << std::setiosflags( std::ios::right );
        std::cout << " iZone = " << std::setw( wordWidth ) << iZone;
        std::cout << " nZone = " << std::setw( wordWidth ) << gridMediator->numberOfZones;
        std::cout << " ni, nj, nk = ";
        std::cout << std::setw( wordWidth ) << ni;
        std::cout << std::setw( wordWidth ) << nj;
        std::cout << std::setw( wordWidth ) << nk;
        std::cout << "\n";
    }

    for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        fileO.DumpCoorAscii( grid->nodeMesh->xN );
        fileO.DumpCoorAscii( grid->nodeMesh->yN );

        if ( nkflag )
        {
            fileO.DumpCoorAscii( grid->nodeMesh->yN );
        }
    }

    fileO.CloseFile();
}

void Plot3D::DumpCoorAscii( std::fstream & file, RealField & coor )
{
    int wordWidth = 15;
    int nWord = 5;
    int nPoint = coor.size();
    int nCount = 0;
    for ( int i = 0; i < nPoint; ++ i )
    {
        file << std::setw( wordWidth ) << coor[ i ];
        nCount ++;
        if ( nCount % nWord == 0 )
        {
            file << "\n";
        }
    }
}

void Plot3D::ReadBc( GridMediator * gridMediator )
{
    std::string & bcName = gridMediator->bcFile;
    //\t is the tab key
    std::string separator = " =\r\n#$,;";

    FileIO ioFile;
    ioFile.OpenPrjFile( bcName, std::ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    ioFile.ReadNextNonEmptyLine();
    int flowSolverIndex = ioFile.ReadNextDigit< int >();
    std::cout << " flowSolverIndex = " << flowSolverIndex << std::endl;

    ioFile.ReadNextNonEmptyLine();
    int numberOfZones = ioFile.ReadNextDigit< int >();
    std::cout << " numberOfZones = " << numberOfZones << std::endl;

    if ( numberOfZones != gridMediator->numberOfZones )
    {
        Stop( "nzone in boundary is not consistent with nzone in grid!\n" );
    }

    bool readPid = false;

    ZoneState::pid.resize( numberOfZones );

    IntField zoneidlist;

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );
        ioFile.ReadNextNonEmptyLine();

        int ni = ioFile.ReadNextDigit< int >();
        int nj = ioFile.ReadNextDigit< int >();
        int nk = 1;
        if ( ONEFLOW::IsThreeD() )
        {
            nk = ioFile.ReadNextDigit< int >();
        }

        std::cout << " ni, nj, nk = " << ni << " " << nj << " " << nk << std::endl;
        int ndif = ONEFLOW::ABS( ni - grid->ni )
                 + ONEFLOW::ABS( nj - grid->nj )
                 + ONEFLOW::ABS( nk - grid->nk );
        if ( ndif != 0 )
        {
            std::cout << "The dimension in " << iZone + 1 << "block boundary file is not consistent with the dimension in grid!\n";
        }

        std::string word = ioFile.ReadNextWord();

        if ( word == "proc" )
        {
            readPid = true;
            int pid = ioFile.ReadNextDigit< int >();
            ZoneState::pid[ iZone ] = pid;
        }

        ioFile.ReadNextNonEmptyLine();
        std::string blockName = ioFile.ReadNextWord();
        ioFile.ReadNextNonEmptyLine();

        int nBcRegions = ioFile.ReadNextDigit< int >();

        grid->bcRegionGroup->Create( nBcRegions );
        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            int imin, imax, jmin, jmax, kmin, kmax;

            ioFile.ReadNextNonEmptyLine();
            imin = ioFile.ReadNextDigit< int >();
            imax = ioFile.ReadNextDigit< int >();

            jmin = ioFile.ReadNextDigit< int >();
            jmax = ioFile.ReadNextDigit< int >();

            if ( ONEFLOW::IsThreeD() )
            {
                kmin = ioFile.ReadNextDigit< int >();
                kmax = ioFile.ReadNextDigit< int >();
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            int bcType = ioFile.ReadNextDigit< int >();
            BcRegion * bcRegion = new BcRegion( iZone, ir );
            bcRegion->s->SetRegion( imin, imax, jmin, jmax, kmin, kmax );
            bcRegion->s->zid = iZone;
            bcRegion->bcType = bcType;
            bcRegionGroup->SetBcRegion( ir, bcRegion );

            if ( bcType == 3 )
            {
                zoneidlist.push_back( iZone );
            }

            if ( bcType < 0 )
            {
                ioFile.ReadNextNonEmptyLine();

                imin = ioFile.ReadNextDigit< int >();
                imax = ioFile.ReadNextDigit< int >();

                jmin = ioFile.ReadNextDigit< int >();
                jmax = ioFile.ReadNextDigit< int >();

                if ( ONEFLOW::IsThreeD() )
                {
                    kmin = ioFile.ReadNextDigit< int >();
                    kmax = ioFile.ReadNextDigit< int >();
                }
                else
                {
                    kmin = 1;
                    kmax = 1;
                }

                bcRegion->t->SetRegion( imin, imax, jmin, jmax, kmin, kmax );
                bcRegion->t->zid = ioFile.ReadNextDigit< int >();

            }
        }
    }

    int kkk = 1;

    ioFile.CloseFile();
}

void Plot3D::DumpBc( GridMediator * gridMediator )
{
    std::string & bcName = gridMediator->bcFile;

    std::fstream file;
    Prj::OpenPrjFile( file, bcName, std::ios_base::out );

    int flowSolverIndex = 1;
    file << flowSolverIndex << "\n";
    std::cout << " flowSolverIndex = " << flowSolverIndex << "\n";

    int numberOfZones = gridMediator->numberOfZones;

    file << numberOfZones << "\n";
    std::cout << " numberOfZones = " << numberOfZones << "\n";

    bool readPid = false;

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        file << ni << " " << nj;
        if ( ONEFLOW::IsThreeD() )
        {
            file << nk;
        }
        file << "\n";

        std::cout << " ni, nj, nk = " << ni << " " << nj << " " << nk << "\n";

        std::string blockName = grid->name;
        file << blockName << "\n";

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        int nBcRegions = bcRegionGroup->regions->size();

        file << nBcRegions << "\n";

        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = bcRegionGroup->GetBcRegion( ir );

            int imin, imax, jmin, jmax, kmin, kmax;
            BasicRegion * s = bcRegion->s;
            imin = s->start[ 0 ];
            imax = s->end[ 0 ];
            jmin = s->start[ 1 ];
            jmax = s->end[ 1 ];
            kmin = s->start[ 2 ];
            kmax = s->end[ 2 ];
            int width = 5;
            file << std::setiosflags( std::ios::right );
            file << std::setw( width ) << imin;
            file << std::setw( width ) << imax;
            file << std::setw( width ) << jmin;
            file << std::setw( width ) << jmax;

            if ( ONEFLOW::IsThreeD() )
            {
                file << std::setw( width ) << kmin;
                file << std::setw( width ) << kmax;
            }

            int bcType = bcRegion->bcType;
            file << std::setw( width ) << bcType;
            file << "\n";

            if ( bcType < 0 )
            {
                BasicRegion * t = bcRegion->t;
                imin = t->start[ 0 ];
                imax = t->end[ 0 ];
                jmin = t->start[ 1 ];
                jmax = t->end[ 1 ];
                kmin = t->start[ 2 ];
                kmax = t->end[ 2 ];
                file << std::setiosflags( std::ios::right );
                file << std::setw( width ) << imin;
                file << std::setw( width ) << imax;
                file << std::setw( width ) << jmin;
                file << std::setw( width ) << jmax;

                if ( ONEFLOW::IsThreeD() )
                {
                    file << std::setw( width ) << kmin;
                    file << std::setw( width ) << kmax;
                }
                file << std::setw( width ) << t->zid;
                file << "\n";
            }
        }
    }

    Prj::CloseFile( file );
}

void Plot3D::ReadCoor( FileIO * ioFile, RealField & coordinate )
{
    HXSize_t numberOfNodes = coordinate.size();
    HXSize_t i = 0;
    while ( i < numberOfNodes )
    {
        int num = 1;
        Real tmp = ioFile->ReadNextDigit< Real >( num );
        for ( HXSize_t j = 0; j < num; ++ j )
        {
            coordinate[ i ] = tmp;
            ++ i;
            if ( i >= numberOfNodes )
            {
                int left = num - j;
                std::cout << " unread elem = " << left << std::endl;
                break;
            }
        }
    }
}

void Plot3D::ReadCoor( FileIO * ioFile, RealField & coor, int total_size )
{
    while ( coor.size() < total_size )
    {
        int num = 1;
        Real tmp = ioFile->ReadNextDigit< Real >( num );

        for ( int i = 0; i < num; ++ i )
        {
            coor.push_back( tmp );
        }
    }
}

void Plot3D::Plot3DToCgns( ZgridMediator * zgridMediator )
{
    std::cout << "plot3d to cgns\n";
    GridMediator * gridMediator = new GridMediator();
    gridMediator->gridFile = grid_para.gridFile;
    gridMediator->bcFile = grid_para.bcFile;
    gridMediator->targetFile = grid_para.targetFile;

    gridMediator->gridType = grid_para.filetype;
    gridMediator->ReadGrid();
    gridMediator->AddDefaultName();

    zgridMediator->SetDeleteFlag( true );
    zgridMediator->AddGridMediator( gridMediator );
}


bool GetPlot3D_NKFlag()
{
    int iplot3d = GetDataValue< int >( "iplot3d" );

    bool is3d           = ( Dim::dimension == ONEFLOW::THREE_D );
    bool gridgen_plot3d = ( iplot3d == 0 );
    bool readnkflag     = is3d || gridgen_plot3d;

    if ( Dim::dimension == ONEFLOW::TWO_D )
    {
        if ( gridgen_plot3d )
        {
            std::cout << "The grid type is gridgen plot3d format\n";
        }
        else
        {
            std::cout << "The grid type is standard plot3d format\n";
        }
    }

    return readnkflag;
}
EndNameSpace
