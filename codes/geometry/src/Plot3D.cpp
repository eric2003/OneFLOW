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

#include "Plot3D.h"
#include "GridMediator.h"
#include "BasicIO.h"
#include "Prj.h"
#include "StrGrid.h"
#include "NodeMesh.h"
#include "BgGrid.h"
#include "GridState.h"
#include "AsciiFileIO.h"
#include "Dimension.h"
#include "HXMath.h"
#include "Zone.h"
#include "BcRecord.h"
#include "DataBase.h"
#include <iostream>

using namespace std;

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

void Plot3D::ReadCoorBinary( GridMediator * gridMediator )
{
    string & fileName = gridMediator->gridFile;
    string separator  = " =\r\n\t#$,;";

	fstream file;
    ONEFLOW::OpenPrjFile( file, fileName, ios_base::in|ios_base::binary );

    HXRead( & file, gridMediator->numberOfZones );
    gridMediator->gridVector.resize( gridMediator->numberOfZones );

    bool readnkflag = GetPlot3D_NKFlag();

	for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
	{
		int ni, nj, nk = 1;
        HXRead( & file, ni );
        HXRead( & file, nj );

        if ( readnkflag )
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
        grid->nodeMesh->CreateNodes( grid->nNode );

        grid->SetLayout();

		int wordWidth = 8;

		cout << setiosflags( ios::right );
		cout << " iZone = " << setw( wordWidth ) << iZone;
		cout << " nZone = " << setw( wordWidth ) << gridMediator->numberOfZones;
		cout << " ni, nj, nk = ";
		cout << setw( wordWidth ) << ni;
		cout << setw( wordWidth ) << nj;
		cout << setw( wordWidth ) << nk;
		cout << endl;
	}

	for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
	{
		StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );
		int numberOfNodes = grid->nNode;

		int ni = grid->ni;
		int nj = grid->nj;
		int nk = grid->nk;

        HXRead( & file, grid->nodeMesh->xN );
        HXRead( & file, grid->nodeMesh->yN );

		if ( readnkflag )
		{
			HXRead( & file, grid->nodeMesh->zN );
		}
		else
		{
			grid->nodeMesh->zN = 0;
		}
	}

	ONEFLOW::CloseFile( file );
}

void Plot3D::ReadCoorAscii( GridMediator * gridMediator )
{
    string & fileName = gridMediator->gridFile;

	AsciiFileRead ioFile;
    string separator  = " =\r\n\t#$,;";
	ioFile.OpenPrjFile( fileName, ios_base::in );
	ioFile.SetDefaultSeparator( separator );

    gridMediator->numberOfZones = ioFile.ReadNextDigit< int >();
    gridMediator->gridVector.resize( gridMediator->numberOfZones );

    bool readnkflag = GetPlot3D_NKFlag();

    int zCount = 0;
	while ( zCount < gridMediator->numberOfZones )
	{
		int ni, nj, nk = 1;
        ni = ioFile.ReadNextDigit< int >();
        nj = ioFile.ReadNextDigit< int >();

        if ( readnkflag )
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
        grid->nodeMesh->CreateNodes( grid->nNode );

        grid->SetLayout();

		int wordWidth = 8;

		cout << setiosflags( ios::right );
		cout << " iZone = " << setw( wordWidth ) << zCount;
		cout << " nZone = " << setw( wordWidth ) << gridMediator->numberOfZones;
		cout << " ni, nj, nk = ";
		cout << setw( wordWidth ) << ni;
		cout << setw( wordWidth ) << nj;
		cout << setw( wordWidth ) << nk;
		cout << endl;
        ++ zCount;
	}

	for ( int iZone = 0; iZone < gridMediator->numberOfZones; ++ iZone )
	{
		StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ iZone ] );
		int numberOfNodes = grid->nNode;

		int ni = grid->ni;
		int nj = grid->nj;
		int nk = grid->nk;

        int total_size = 2 * numberOfNodes;
        //if ( Dim::dimension == ONEFLOW::THREE_D )
		if ( readnkflag )
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

        //if ( Dim::dimension == ONEFLOW::THREE_D )
		if ( readnkflag )
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

void Plot3D::ReadBc( GridMediator * gridMediator )
{
    string & bcName = gridMediator->bcFile;
	//\tÎªtab¼ü
	string separator = " =\r\n#$,;";

	AsciiFileRead ioFile;
	ioFile.OpenPrjFile( bcName, ios_base::in );
	ioFile.SetDefaultSeparator( separator );

	ioFile.ReadNextNonEmptyLine();
	int flowSolverIndex = ioFile.ReadNextDigit< int >();
	cout << " flowSolverIndex = " << flowSolverIndex << endl;

	ioFile.ReadNextNonEmptyLine();
	int numberOfZones = ioFile.ReadNextDigit< int >();
	cout << " numberOfZones = " << numberOfZones << endl;

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

		cout << " ni, nj, nk = " << ni << " " << nj << " " << nk << endl;
		int ndif = ONEFLOW::ABS( ni - grid->ni )
                 + ONEFLOW::ABS( nj - grid->nj )
                 + ONEFLOW::ABS( nk - grid->nk );
		if ( ndif != 0 )
		{
			cout << "The dimension in " << iZone + 1 << "block boundary file is not consistent with the dimension in grid!\n";
		}

		string word = ioFile.ReadNextWord();

		if ( word == "proc" )
		{
			readPid = true;
			int pid = ioFile.ReadNextDigit< int >();
		    ZoneState::pid[ iZone ] = pid;
		}

		ioFile.ReadNextNonEmptyLine();
		string blockName = ioFile.ReadNextWord();
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

void Plot3D::ReadCoor( AsciiFileRead * ioFile, RealField & coordinate )
{
    Int numberOfNodes = coordinate.size();
    Int i = 0;
    while ( i < numberOfNodes )
    {
        int num = 1;
        Real tmp = ioFile->ReadNextDigit< Real >( num );
        for ( Int j = 0; j < num; ++ j )
        {
            coordinate[ i ] = tmp;
            ++ i;
            if ( i >= numberOfNodes )
            {
                int left = num - j;
                cout << " unread elem = " << left << endl;
                break;
            }
        }
    }
}

void Plot3D::ReadCoor( AsciiFileRead * ioFile, RealField & coor, int total_size )
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

bool GetPlot3D_NKFlag()
{
    int iplot3d = GetDataValue< int >( "iplot3d" );

    bool is3d           = ( Dim::dimension == ONEFLOW::THREE_D );
    bool gridgen_plot3d = ( iplot3d == 0 );
    bool readnkflag     = is3d || gridgen_plot3d;

    return readnkflag;
}
EndNameSpace