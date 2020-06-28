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

#include "Cylinder.h"
#include "Transfinite.h"
#include "CurveLine.h"
#include "CurveMesh.h"
#include "PointMachine.h"
#include "CurveMachine.h"
#include "GridMachine.h"
#include "LineMachine.h"
#include "Prj.h"
#include "DataBaseIO.h"
#include "Boundary.h"
#include "HXMath.h"
#include "DataBase.h"
#include "FileIO.h"
#include "FileUtil.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

DomainData::DomainData()
{
    ;
}

DomainData::~DomainData()
{
    ;
}

void DomainData::Alloc()
{
    ONEFLOW::AllocateVector(x, ni, nj);
    ONEFLOW::AllocateVector(y, ni, nj);
    ONEFLOW::AllocateVector(z, ni, nj);
}

void DomainData::Symmetry( DomainData * datain )
{
    this->ni = datain->ni;
    this->nj = datain->nj;
    this->Alloc();

    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            Real xx = datain->x[ i ][ j ];
            Real yy = datain->y[ i ][ j ];
            Real zz = datain->z[ i ][ j ];

            this->x[ i ][ j ] =   xx;
            this->y[ i ][ j ] = - yy;
            this->z[ i ][ j ] =   zz;
        }
    }
}

void DomainData::Join( DomainData * d1, DomainData * d2 )
{
    this->ni = 2 * d1->ni - 1;
    this->nj = d1->nj;
    this->Alloc();

    for ( int j = 0; j < d2->nj; ++ j )
    {
        for ( int i = 0; i < d2->ni; ++ i )
        {
            int ii = d2->ni - i - 1;
            Real xx = d2->x[ ii ][ j ];
            Real yy = d2->y[ ii ][ j ];
            Real zz = d2->z[ ii ][ j ];

            this->x[ i ][ j ] = xx;
            this->y[ i ][ j ] = yy;
            this->z[ i ][ j ] = zz;
        }
    }

    for ( int j = 0; j < d1->nj; ++ j )
    {
        for ( int i = 0; i < d1->ni; ++ i )
        {
            int ii = i;
            Real xx = d1->x[ ii ][ j ];
            Real yy = d1->y[ ii ][ j ];
            Real zz = d1->z[ ii ][ j ];

            int i1 = i + d2->ni - 1;

            this->x[ i1 ][ j ] = xx;
            this->y[ i1 ][ j ] = yy;
            this->z[ i1 ][ j ] = zz;
        }
    }
}

Cylinder::Cylinder()
{
    this->strCurveLoop = new StrCurveLoop();
}

Cylinder::~Cylinder()
{
    delete this->strCurveLoop;
}

void Cylinder::Run( int igene )
{
    if ( igene == 3 )
    {
        this->HalfCylinder();
        //this->QuarterCylinder();
    }
    else if ( igene == 4 )
    {
        this->GenePlate();
    }
}

void Cylinder::GenePlate()
{
    grid_Machine.ReadScript();
    grid_Machine.GeneGrid();


    int kkk = 1;
}

void Cylinder::HalfCylinder()
{
    strCurveLoop->ni = 61;
    strCurveLoop->nj = 81;
    domain_data.ni = strCurveLoop->ni;
    domain_data.nj = strCurveLoop->nj;
    domain_data.Alloc();

    this->SetBoundaryGrid();
    this->GeneDomain();
    symm_domain.Symmetry( & domain_data );
    final_domain.Join( & domain_data, & symm_domain );

    this->nZone = 1;

    IntField bcList;
    bcList.push_back( BC::SOLID_SURFACE );
    bcList.push_back( BC::INFLOW );
    bcList.push_back( BC::OUTFLOW );
    bcList.push_back( BC::OUTFLOW );

    DumpGrid( "/grid/halfcylinder.grd", & final_domain );
    DumpBcFile( "/grid/halfcylinder.inp", & final_domain, bcList );
    this->ToTecplot( "/grid/halfcylinder-tecplot.dat", & final_domain );
}

void Cylinder::QuarterCylinder()
{
    int ni = 61;
    int nj = 81;
    domain_data.ni = ni;
    domain_data.nj = nj;
    domain_data.Alloc();

    CurveLine s1, s2, s3, s4;

    s1.Alloc( ni );
    s2.Alloc( ni );

    s3.Alloc( nj );
    s4.Alloc( nj );

    this->SetBoundaryGrid();
    this->GeneDomain();

    this->nZone = 1;
    IntField bcList;
    bcList.push_back( BC::SOLID_SURFACE );
    bcList.push_back( BC::INFLOW );
    bcList.push_back( BC::SYMMETRY );
    bcList.push_back( BC::OUTFLOW );

    DumpGrid( "/grid/cylinder.grd", & domain_data );
    DumpBcFile( "/grid/cylinder.inp", & domain_data, bcList );
    this->ToTecplot( "/grid/cylinder-tecplot.dat", & domain_data );
}

void Cylinder::DumpGrid( const string & fileName, DomainData * domain )
{
    RealField xN;
    RealField yN;
    RealField zN;

    for ( int j = 0; j < domain->nj; ++ j )
    {
        for ( int i = 0; i < domain->ni; ++ i )
        {
            xN.push_back( domain->x[ i ][ j ] );
            yN.push_back( domain->y[ i ][ j ] );
            zN.push_back( domain->z[ i ][ j ] );
        }
    }

    fstream file;
    OpenPrjFile( file, fileName, ios_base::out | ios_base::binary );
    int nZone = 1;
    int nk = 1;
    HXWrite( & file, nZone );
    HXWrite( & file, domain->ni );
    HXWrite( & file, domain->nj );
    HXWrite( & file, nk );

    HXWrite( & file, xN );
    HXWrite( & file, yN );
    HXWrite( & file, zN );

    CloseFile( file );
}

void Cylinder::DumpBcFile( const string & fileName, DomainData * domain, IntField & bcList )
{
    fstream file;
    OpenPrjFile( file, fileName, ios_base::out );
    int solver = 1;
    string zName = "A";
    int nBc = 4;
    file << solver << endl;
    file << nZone << endl;
    file << domain->ni << " " << domain->nj << endl;
    file << zName << endl;
    file << nBc << endl;

    DumpBc( file, 1         , domain->ni, 1         , 1         , bcList[ 0 ] );
    DumpBc( file, 1         , domain->ni, domain->nj, domain->nj, bcList[ 1 ] );
    DumpBc( file, 1         , 1         , 1         , domain->nj, bcList[ 2 ] );
    DumpBc( file, domain->ni, domain->ni, 1         , domain->nj, bcList[ 3 ] );

    CloseFile( file );
}

void Cylinder::ToTecplot( const string & fileName, DomainData * domain )
{
    int ni = domain->ni;
    int nj = domain->nj;
    int nk = 1;

    fstream file;
    OpenPrjFile( file, fileName, ios_base::out );
    file << " VARIABLES = \"X\" \"Y\" \"Z\"" << "\n";
    file << "ZONE DATAPACKING = BLOCK, I = " << ni << ", J = " << nj << ", K = " << nk << "\n";

    ONEFLOW::ToTecplot( file, domain->x, ni, nj, nk );
    ONEFLOW::ToTecplot( file, domain->y, ni, nj, nk );
    ONEFLOW::ToTecplot( file, domain->z, ni, nj, nk );

    CloseFile( file );
}

void Cylinder::CalcCircleCenter( PointType & p1, PointType & p2, PointType & p0, PointType & pcenter )
{
    Real L = CalcDist( p1, p0 );
    Real H = CalcDist( p2, p0 );
    Real R = ( L * L + H * H ) / ( 2 * L );
    Real S = R - L;
    Real coef = R / L;
    pcenter.x = p1.x - coef * ( p1.x - p0.x );
    pcenter.y = p1.y - coef * ( p1.y - p0.y );
    pcenter.z = p1.z - coef * ( p1.z - p0.z );
}

void Cylinder::SetBoundaryGrid()
{
    string fileName = GetDataValue< string >( "gridLayoutFileName" );
    string separator = " =\r\n\t#$,;\"";

    FileIO ioFile;

    ioFile.OpenPrjFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    ioFile.ReadNextNonEmptyLine();

    Real x0 = ioFile.ReadNextDigit< Real >();
    Real y0 = ioFile.ReadNextDigit< Real >();
    Real z0 = ioFile.ReadNextDigit< Real >();

    ioFile.ReadNextNonEmptyLine();
    Real x1 = ioFile.ReadNextDigit< Real >();
    Real y1 = ioFile.ReadNextDigit< Real >();
    Real z1 = ioFile.ReadNextDigit< Real >();

    ioFile.ReadNextNonEmptyLine();
    Real x2 = ioFile.ReadNextDigit< Real >();
    Real y2 = ioFile.ReadNextDigit< Real >();
    Real z2 = ioFile.ReadNextDigit< Real >();

    ioFile.ReadNextNonEmptyLine();
    Real x3 = ioFile.ReadNextDigit< Real >();
    Real y3 = ioFile.ReadNextDigit< Real >();
    Real z3 = ioFile.ReadNextDigit< Real >();

    ioFile.ReadNextNonEmptyLine();
    Real x4 = ioFile.ReadNextDigit< Real >();
    Real y4 = ioFile.ReadNextDigit< Real >();
    Real z4 = ioFile.ReadNextDigit< Real >();

    ioFile.ReadNextNonEmptyLine();
    Real ds_start = ioFile.ReadNextDigit< Real >();
    Real ds_end   = ioFile.ReadNextDigit< Real >();
    this->beta = ioFile.ReadNextDigit< Real >();

    ioFile.CloseFile();

    point_Machine.AddPoint( x1, y1, z1 );
    point_Machine.AddPoint( x2, y2, z2 );
    point_Machine.AddPoint( x3, y3, z3 );
    point_Machine.AddPoint( x4, y4, z4 );
    point_Machine.AddPoint( x0, y0, z0 );

    PointType * p1 = point_Machine.GetPoint( 0 );
    PointType * p2 = point_Machine.GetPoint( 1 );
    PointType * p3 = point_Machine.GetPoint( 2 );
    PointType * p4 = point_Machine.GetPoint( 3 );
    PointType * p0 = point_Machine.GetPoint( 4 );

    curve_Machine.AddLine( p1->id, p2->id );
    curve_Machine.AddLine( p3->id, p4->id );
    curve_Machine.AddCircle( p1->id, p3->id, p0->id );
    curve_Machine.AddParabolic( p2->id, p4->id );

    this->strCurveLoop->AddCurve( 0 );
    this->strCurveLoop->AddCurve( 1 );
    this->strCurveLoop->AddCurve( 2 );
    this->strCurveLoop->AddCurve( 3 );

    this->strCurveLoop->SetDimension();

    CurveLine * s1 = this->strCurveLoop->GetCurve( 0 );
    CurveLine * s2 = this->strCurveLoop->GetCurve( 1 );
    CurveLine * s3 = this->strCurveLoop->GetCurve( 2 );
    CurveLine * s4 = this->strCurveLoop->GetCurve( 3 );

    s3->GenerateCurveLine();
    s4->GenerateCurveLine();
    s1->MakeCircle();
    s2->GenerateParabolicLine();
 }

void Cylinder::GeneDomain()
{
    CurveLine * s1 = this->strCurveLoop->GetCurve( 0 );
    CurveLine * s2 = this->strCurveLoop->GetCurve( 1 );
    CurveLine * s3 = this->strCurveLoop->GetCurve( 2 );
    CurveLine * s4 = this->strCurveLoop->GetCurve( 3 );

    int ni = s1->nNode;
    int nj = s3->nNode;

    int il = 0;
    int ir = ni - 1;

    int jl = 0;
    int jr = nj - 1;

    for ( int i = 0; i < ni; ++ i )
    {
        domain_data.x[ i ][ jl ] = s1->x[ i ];
        domain_data.y[ i ][ jl ] = s1->y[ i ];
        domain_data.z[ i ][ jl ] = s1->z[ i ];

        domain_data.x[ i ][ jr ] = s2->x[ i ];
        domain_data.y[ i ][ jr ] = s2->y[ i ];
        domain_data.z[ i ][ jr ] = s2->z[ i ];
    }

    for ( int j = 0; j < nj; ++ j )
    {
        domain_data.x[ il ][ j ] = s3->x[ j ];
        domain_data.y[ il ][ j ] = s3->y[ j ];
        domain_data.z[ il ][ j ] = s3->z[ j ];

        domain_data.x[ ir ][ j ] = s4->x[ j ];
        domain_data.y[ ir ][ j ] = s4->y[ j ];
        domain_data.z[ ir ][ j ] = s4->z[ j ];
    }

    TransfiniteInterpolation( domain_data.x, ni, nj );
    TransfiniteInterpolation( domain_data.y, ni, nj );
    TransfiniteInterpolation( domain_data.z, ni, nj );

    this->beta = 1.002;
    RealField nbx, nby, nbz;
    s1->CalcNormal( nbx, nby, nbz );
    nbx[ 0 ] = - 1.0;
    nby[ 0 ] = 0.0;
    nbz[ 0 ] = 0.0;

    nbx[ ni - 1 ] = 0.0;
    nby[ ni - 1 ] = 1.0;
    nbz[ ni - 1 ] = 0.0;
    AlgebraInterpolation( domain_data.x, domain_data.y, domain_data.z, ni, nj, nbx, nby, nbz, this->beta );
}

void ToTecplot( fstream & file, RealField2D & coor, int ni, int nj, int nk )
{
    int numberOfWords = 5;

    int icount = 0;
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            for ( int i = 0; i < ni; ++ i )
            {
                file << coor[ i ][ j ] << " ";
                if ( ( icount + 1 ) % numberOfWords == 0 ) file << "\n";
                icount ++;
            }
        }
    }
    if ( ( icount + 1 ) % numberOfWords != 0 ) file << "\n";
}

void DumpBc( fstream &file, int imin, int imax, int jmin, int jmax, int bcType )
{
    int width = 5;
    file << setw( width ) << imin << setw( width ) << imax;
    file << setw( width ) << jmin << setw( width ) << jmax;
    file << setw( width ) << bcType << endl;
}


EndNameSpace