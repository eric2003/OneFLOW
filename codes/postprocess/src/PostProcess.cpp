/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "PostProcess.h"
#include "PIO.h"
#include "FileIO.h"
#include "HXMath.h"
#include "StrUtil.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

Post::Post()
{
}

Post::~Post()
{
}

void Post::Run()
{
    cout << " Post::Run()\n";

    fstream file;

    string fileName = "grid/wallaero.dat";

    //\t is the tab key
    string separator = " =\r\n\t#$,;\"";

    FileIO ioFile;
    ioFile.OpenPrjFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );
    int count = 0;
    //vector< Real > xList, yList, zList;
    MakeCurveClass makeCurve;
    while ( ! ioFile.ReachTheEndOfFile() )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        count ++;
        //ioFile.DumpLineContentToScreen();

        if ( count <= 8 )
        {
            cout << "header = ";
            ioFile.DumpLineContentToScreen();
            continue;
        }

        Real x = ioFile.ReadNextDigit< Real >();
        Real y = ioFile.ReadNextDigit< Real >();
        Real z = ioFile.ReadNextDigit< Real >();
        Real cp = ioFile.ReadNextDigit< Real >();
        Real cf = ioFile.ReadNextDigit< Real >();
        //cout << x << " " << y << " " << z << "\n";
        //makeCurve.AddPoint( x, y, z );
        makeCurve.AddPointValue( x, y, z, cp, cf );
    }
    makeCurve.Run();
    ioFile.CloseFile();
}

void PostSimu()
{
    Post * post = new Post();
    post->Run();
    delete post;
}

VectDir::VectDir()
{
    ;
}

VectDir::~VectDir()
{
    ;
}

void VectDir::Normalize()
{
    Real ds = DIST( dx, dy );
    dx /= ( ds + SMALL );
    dy /= ( ds + SMALL );
}

Real VectDir::Angle( VectDir * vd )
{
    Real cp = dx * vd->dx + dy * vd->dy;
    Real val = 180.0 / PI;
    Real ret = acos( cp ) * val;
    return ret;
}

void GetCurveTitle( StringField & title )
{
    title.resize( 0 );
    title.push_back("title=\"THE FLOW FIELD OF ONEFLOW\"");
    title.push_back("variables=");
    title.push_back("\"x\"");
    title.push_back("\"y\"");
    title.push_back("\"z\"");
}

void GetCurveValueTitle( StringField & title )
{
    title.resize( 0 );
    title.push_back("title=\"THE FLOW FIELD OF ONEFLOW\"");
    title.push_back("variables=");
    title.push_back("\"x\"");
    title.push_back("\"y\"");
    title.push_back("\"z\"");
    title.push_back("\"-cp\"");
    title.push_back("\"cf\"");
}

PointEdgeClass::PointEdgeClass()
{
}

PointEdgeClass::~PointEdgeClass()
{
}

void PointEdgeClass::EdgeToPoint( vector< IntField > & edgeList )
{
    for ( int i = 0; i < edgeList.size(); ++ i )
    {
        IntField &pq = edgeList[ i ];
        this->idset.insert( pq[ 0 ] );
        this->idset.insert( pq[ 1 ] );
    }

    for ( set<int>::iterator iter = idset.begin(); iter != idset.end(); ++ iter )
    {
        this->idlist.push_back( * iter );
    }
}

CurveData::CurveData()
{
    file_prestr = "";
}

CurveData::~CurveData()
{
    ;
}

void CurveData::SetFilePreStr( const string & file_prestr )
{
    this->file_prestr = file_prestr;
}

void CurveData::AddPoint( Real x, Real y, Real z )
{
    this->xList.push_back( x );
    this->yList.push_back( y );
    this->zList.push_back( z );
}

void CurveData::AddValue( Real cp, Real cf )
{
    this->cpList.push_back( cp );
    this->cfList.push_back( cf );
}

void CurveData::AddPointValue( Real x, Real y, Real z, Real cp, Real cf )
{
    this->AddPoint( x, y, z );
    this->AddValue( cp, cf );
}

void CurveData::FindExtremePoint()
{
    int N = xList.size();
    for ( int i = 0; i < N; ++ i )
    {
        extremeList.push_back( 1 );
    }

    for ( int p = 0; p < N; ++ p )
    {
        cout << " p = " << p << " N = " << N << "\n";
        for ( int q = p + 1; q < N; ++ q )
        {
            for ( int r = q + 1; r < N; ++ r )
            {
                for ( int s = 0; s < N; ++ s )
                {
                    if ( s == p || s == q || s == r || extremeList[ s ] == 0 )
                    {
                        continue;
                    }

                    if ( InTriangle( p, q, r, s ) )
                    {
                        extremeList[ s ] = 0;
                    }
                }
            }
        }
    }
    int count = 0;
    for ( int i = 0; i < N; ++ i )
    {
        if ( extremeList[ i ] == 1 )
        {
            count ++;
        }
    }
    cout << "number of extreme points = " << count << " total points = " << N << "\n";
}

bool CurveData::InTriangle( int p, int q, int r, int s )
{
    bool pqLeft = ToLeft( p, q, s );
    bool qrLeft = ToLeft( q, r, s );
    bool rpLeft = ToLeft( r, p, s );

    return ( pqLeft == qrLeft ) && ( qrLeft == rpLeft );
}

bool CurveData::ToLeft( int p, int q, int s )
{
    Real xp = xList[ p ];
    Real yp = yList[ p ];
    Real xq = xList[ q ];
    Real yq = yList[ q ];
    Real xs = xList[ s ];
    Real ys = yList[ s ];

    Real s2 = xp * yq - yp * xq +
        xq * ys - yq * xs +
        xs * yp - ys * xp;

    if ( s2 > 0.0 ) return true;
    return false;
}

bool CurveData::CheckEdge( int p, int q )
{
    bool lEmpty = true;
    bool rEmpty = true;
    int N = xList.size();
    for ( int k = 0; k < N && ( lEmpty || rEmpty ); ++ k )
    {
        if ( k != p && k != q )
        {
            this->ToLeft( p, q, k ) ?
                lEmpty = false :
                rEmpty = false;
        }
    }
    if ( lEmpty || rEmpty )
    {
        this->extremeList[ p ] = 1;
        this->extremeList[ q ] = 1;
        return true;
    }
    return false;
}

void CurveData::AddExtremeEdge( int p, int q )
{
    IntField pq;
    pq.push_back( p );
    pq.push_back( q );

    std::sort( pq.begin(), pq.end() );

    p2p[ p ].insert( q );
    p2p[ q ].insert( p );

    int id = searchEdgeList.size();

    HXSort< IntField > edge( pq, id );

    set < HXSort< IntField > >::iterator iter = this->searchEdgeList.find( edge );

    if ( iter == this->searchEdgeList.end() )
    {
        this->extremeEdgeList.push_back( pq );
        this->searchEdgeList.insert( edge );
    }
}

bool CurveData::FindEdge( int p, int q )
{
    IntField pq;
    pq.push_back( p );
    pq.push_back( q );

    std::sort( pq.begin(), pq.end() );

    int id = searchEdgeList.size();

    HXSort< IntField > edge( pq, id );

    set < HXSort< IntField > >::iterator iter = this->searchEdgeList.find( edge );

    return iter != this->searchEdgeList.end();
}

void CurveData::InitP2p()
{
    int N = xList.size();
    p2p.resize( N );
}

void CurveData::FindExtremeEdge()
{
    this->InitP2p();
    int N = xList.size();
    for ( int i = 0; i < N; ++ i )
    {
        extremeList.push_back( 0 );
    }

    for ( int p = 0; p < N; ++ p ) //test
    {
        cout << " p = " << p << " N = " << N << "\n";
        for ( int q = p + 1; q < N; ++ q ) //each
        {
            bool flag = CheckEdge( p, q ); //directed edge pq

            if ( flag )
            {
                AddExtremeEdge( p, q );
            }
        }
    }
}

int CurveData::FindXminIndex()
{
    int index = 0;
    Real xmin = xList[ index ];

    int N = xList.size();
    
    for ( int i = 1; i < N; ++ i )
    {
        if ( xList[ i ] < xmin )
        {
            xmin = xList[ i ];
            index = i;
        }
    }
    return index;
}

bool CurveData::FindNearestPoint( int p0, int &pNearest )
{
    if ( idset.empty() ) return false;
    Real minds = LARGE;
    pNearest = -1;
    for ( set< int >::iterator iter = idset.begin(); iter != idset.end(); ++ iter )
    {
        int p1 = ( * iter );
        Real dx = xList[ p1 ] - xList[ p0 ];
        Real dy = yList[ p1 ] - yList[ p0 ];
        Real ds = DIST( dx, dy );
        if ( ds < minds )
        {
            minds = ds;
            pNearest = p1;
        }
    }
    return true;
}

void CurveData::AddVector( int p )
{
    newidlist.push_back( p );
}

void CurveData::RemovePoint( int p )
{
    idset.erase( p );
}

int CurveData::FindNextExtremeEdgePoint( int p1, int p2 )
{
    if ( p2p[ p2 ].size() != 1 )
    {
        cout << " p2 = " << p2 << "p2p[ p2 ].size()=" << p2p[ p2 ].size() << "\n";
    }
    for ( set<int>::iterator iter = p2p[ p2 ].begin(); iter != p2p[ p2 ].end(); ++ iter )
    {
        int p3 = ( * iter );
        if ( p3 != p1 ) return p3;
    }
    return -1;
}

void CurveData::SetExtremePoint()
{
    pec.EdgeToPoint( extremeEdgeList );
}

IntField CurveData::GetFistEdge()
{
    IntField &pq = extremeEdgeList[ 0 ];
    return pq;
}

bool CurveData::FindNextPoint( int p1, int & pNext, VectDir * vd )
{
    if ( idset.empty() ) return false;
    Real minds = LARGE;
    pNext = -1;
    for ( set< int >::iterator iter = idset.begin(); iter != idset.end(); ++ iter )
    {
        int p2 = ( * iter );
        Real dx = xList[ p2 ] - xList[ p1 ];
        Real dy = yList[ p2 ] - yList[ p1 ];
        Real ds = DIST( dx, dy );
        if ( ds < minds )
        {
            minds = ds;
            pNext = p2;
        }
    }

    VectDir vecd;
    vecd.dx = xList[ pNext ] - xList[ p1 ];
    vecd.dy = yList[ pNext ] - yList[ p1 ];
    vecd.Normalize();

    Real angle = vecd.Angle( vd );
    if ( angle < 30 ) return true;

    int pCandidate = pNext;
    Real dCandidate = minds;

    cout << " angle = " << angle << "\n";

    minds = LARGE;
    pNext = -1;
    for ( set< int >::iterator iter = idset.begin(); iter != idset.end(); ++ iter )
    {
        int p2 = ( * iter );
        vecd.dx = xList[ p2 ] - xList[ p1 ];
        vecd.dy = yList[ p2 ] - yList[ p1 ];
        Real ds = DIST( vecd.dx, vecd.dy );
        vecd.Normalize();
        angle = vecd.Angle( vd );
        cout << " angle = " << angle << " ds = " << ds << "\n";
        if ( angle >= 30 ) continue;

        if ( ds < minds )
        {
            minds = ds;
            pNext = p2;
        }
    }

    cout << " pNext = " << pNext << "\n";

    if ( pNext == - 1 )
    {
        pNext = pCandidate;
    }

    return true;
}

void CurveData::ConstructFinalList()
{
    IntField &pq = this->GetFistEdge();

    int p1 = pq[ 0 ];
    int p2 = pq[ 1 ];

    //Init idset
    idset.clear();
    int N = this->xList.size();
    for ( int i = 0; i < N; ++ i )
    {
        idset.insert( i );
    }

    this->RemovePoint( p1 );
    this->AddVector( p1 );
    this->RemovePoint( p2 );
    this->AddVector( p2 );

    //find next point
    //1.same direction, 0 < angle < 90
    //2.opposite direction angle > 90
    //The minimum distance is the first principle, and then there are special cases. 
    //In very special cases, some points with the minimum distance are not needed.
    //It is necessary to consider whether they are going along the curve
    int count = 0;
    while ( true )
    {
        VectDir vd;
        vd.dx = xList[ p2 ] - xList[ p1 ];
        vd.dy = yList[ p2 ] - yList[ p1 ];
        vd.Normalize();

        p1 = p2;
        bool flag = FindNextPoint( p1, p2, & vd );
        if ( ! flag ) break;
        this->RemovePoint( p2 );
        this->AddVector( p2 );
        count ++;
        if ( count >= 30000 ) break;
    }
}

void CurveData::CreateXSortList()
{
    this->FindExtremeEdge();
    this->SetExtremePoint();
    this->VisualCreateList( & this->pec );

    this->ConstructFinalList();
    this->VisualCurveValue();
}

void CurveData::VisualCreateList()
{
    string fileName = AddString( "grid/", this->file_prestr, "createdwallaero.dat");
    VisualCreateList( this->newidlist, fileName );
}


void CurveData::VisualCurveValue()
{
    string fileName = AddString( "grid/", this->file_prestr, "createdwallaerovalue.dat");
    VisualCurveValue( this->newidlist, fileName );
}

void CurveData::VisualCreateList( PointEdgeClass * pec )
{
    string fileName = AddString( "grid/", this->file_prestr, "created_pec_wallaero.dat");
    VisualCreateList( pec->idlist, fileName );
}

void CurveData::VisualCreateList( vector<int> & newidlist, const string& fileName )
{
    vector< Real > x, y, z;
    for ( int i = 0; i < newidlist.size(); ++ i )
    {
        int id = newidlist[ i ];
        Real xm = this->xList[ id ];
        Real ym = this->yList[ id ];
        Real zm = this->zList[ id ];
        x.push_back( xm );
        y.push_back( ym );
        z.push_back( zm );
    }

    this->Visual( fileName, x, y, z );
}

void CurveData::VisualCurveValue( vector<int> & newidlist, const string& fileName )
{
    vector< Real > x, y, z, cp, cf;
    for ( int i = 0; i < newidlist.size(); ++ i )
    {
        int id = newidlist[ i ];
        Real xm = this->xList[ id ];
        Real ym = this->yList[ id ];
        Real zm = this->zList[ id ];
        Real cpm = this->cpList[ id ];
        Real cfm = this->cfList[ id ];
        x.push_back( xm );
        y.push_back( ym );
        z.push_back( zm );
        cp.push_back( cpm );
        cf.push_back( cfm );
    }

    this->VisualCurveValue( fileName, x, y, z, cp, cf );
}

void CurveData::Visual( const string & fileName )
{
    this->Visual( fileName, this->xList, this->yList, this->zList );
}

void CurveData::Visual( const string & fileName, vector< Real > & x, vector< Real > & y, vector< Real > & z )
{
    fstream file;

    PIO::ParallelOpenPrj( file, fileName, ios_base::out );

    StringField title;
    GetCurveTitle( title );

    for ( size_t i = 0; i < title.size(); ++ i )
    {
        file << title[i] << "\n";
    }

    file << "Zone i = " << x.size() << " F = POINT " << "\n";

    for ( size_t i = 0; i < x.size(); ++ i )
    {
        file << x[i] << " ";
        file << y[i] << " ";
        file << z[i] << "\n";
    }

    PIO::Close( file );
}

void CurveData::VisualCurveValue( const string & fileName, vector< Real > & x, vector< Real > & y, vector< Real > & z, vector< Real > & cp, vector< Real > & cf )
{
    fstream file;

    PIO::ParallelOpenPrj( file, fileName, ios_base::out );

    StringField title;
    GetCurveValueTitle( title );

    for ( size_t i = 0; i < title.size(); ++ i )
    {
        file << title[i] << "\n";
    }

    int N = x.size();

    file << "Zone i = " << N << " F = POINT " << "\n";

    for ( size_t i = 0; i < N; ++ i )
    {
        file << x[i] << " ";
        file << y[i] << " ";
        file << z[i] << " ";
        file << cp[i] << " ";
        file << cf[i] << "\n";
    }

    PIO::Close( file );
}

MakeCurveClass::MakeCurveClass()
{
    ;
}

MakeCurveClass::~MakeCurveClass()
{
    ;
}

void MakeCurveClass::AddPoint( Real x, Real y, Real z )
{
    totalPart.AddPoint( x, y, z );
}

void MakeCurveClass::AddPointValue( Real x, Real y, Real z, Real cp, Real cf )
{
    totalPart.AddPointValue( x, y, z, cp, cf );
}

void MakeCurveClass::Run()
{
    this->SplitCurve();
    this->NLR7301MainSplit();
    this->Visual();
}

void MakeCurveClass::NLR7301MainSplit()
{
    this->mainPart.SetFilePreStr( "main" );
    this->mainPart.CreateXSortList();
    this->slapPart.SetFilePreStr( "slap" );
    this->slapPart.CreateXSortList();
}

void MakeCurveClass::SplitCurve()
{
    int N = this->totalPart.xList.size();
    Real x1 = 0.89;
    Real x2 = 0.95;
    Real y1 = 0.01;
    int flag = 0;
    for ( int i = 0; i < N; ++ i )
    {
        Real xm = this->totalPart.xList[ i ];
        Real ym = this->totalPart.yList[ i ];
        Real zm = this->totalPart.zList[ i ];
        Real cp = this->totalPart.cpList[ i ];
        Real cf = this->totalPart.cfList[ i ];
        if ( xm <= x1 )
        {
            mainPart.AddPoint( xm, ym, zm );
            mainPart.AddValue( cp, cf );
        }
        else if ( xm >= x2 )
        {
            slapPart.AddPoint( xm, ym, zm );
            slapPart.AddValue( cp, cf );
        }
        else
        {
            if ( ym >= y1 )
            {
                mainPart.AddPoint( xm, ym, zm );
                mainPart.AddValue( cp, cf );
            }
            else
            {
                slapPart.AddPoint( xm, ym, zm );
                slapPart.AddValue( cp, cf );
            }
        }
    }
}

void MakeCurveClass::Visual()
{
    string fileName = "grid/modifiedwallaero.dat";
    this->totalPart.Visual( fileName );

    string fileName1 = "grid/mainwallaero.dat";
    this->mainPart.Visual( fileName1 );

    string fileName2 = "grid/slapwallaero.dat";
    this->mainPart.Visual( fileName2 );
}

EndNameSpace