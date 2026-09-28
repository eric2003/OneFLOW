#include "Grid.h"
#include <iostream>
#include <vector>

Region::Region()
{
}

Region::~Region()
{
}

Region & Region::operator = ( const Region & rhs )
{
    if ( this == & rhs ) return * this;

    this->start = rhs.start;
    this->end = rhs.end;

    return * this;
}

void Region::SetRegion( std::vector<int> & pnts )
{
    int index_dim = pnts.size() / 2;
    this->start.resize( index_dim );
    this->end.resize( index_dim );
    for ( int m = 0; m < index_dim; ++ m )
    {
        this->start[ m ] = pnts[ m ];
        this->end[ m ] = pnts[ index_dim + m ];
    }
}

void Region::Print()
{
    int nSize = this->start.size();
    std::cout << "start:(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->start[ m ];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
    std::cout << "end  :(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->end[ m ];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
}

Coor::Coor()
{
    ;
}

Coor::~Coor()
{
}

void Coor::DumpCoor()
{
    double * xd = reinterpret_cast<double *>( const_cast<char *>( coord.data() ) );
    for ( int i = 0; i < this->nNodes; ++ i )
    {
        //std::cout << coord[i] << " ";
        std::cout << xd[ i ] << " ";
        if ( ( i + 1 ) % 5 == 0 ) std::cout << "\n";
    }
    std::cout << "\n";
}

void Coor::DumpCoorX( std::vector<double> &x )
{
    double * xd = reinterpret_cast<double *>( const_cast<char *>( coord.data() ) );
    for ( int i = 0; i < this->nNodes; ++ i )
    {
        x[ i ] = xd[ i ];
    }
}

ZoneBc::ZoneBc()
{
    ;
}

ZoneBc::~ZoneBc()
{
}

ZoneBc1To1::ZoneBc1To1()
{
    ;
}

ZoneBc1To1::~ZoneBc1To1()
{
}

Zone::Zone()
{
    ;
}

Zone::~Zone()
{
    for ( int i = 0; i < bccos.size(); ++ i )
    {
        delete bccos[ i ];
    }

    for ( int i = 0; i < coors.size(); ++ i )
    {
        delete coors[ i ];
    }
}

int Trans::M[ 3 ][ 3 ];
std::vector<int> Trans::transform;

int Trans::sgn( int x )
{
    if ( x >= 0 )
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

int Trans::del( int x, int y )
{
    if ( std::abs( x ) == std::abs( y ) )
    {
        return 1;
    }
    return 0;
}

void Trans::ZeroMatrix()
{
    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            Trans::M[ i ][ j ] = 0;
        }
    }
}

void Trans::CalcTransformMatrix()
{
    int dim = Trans::transform.size();
    if ( dim == 1 )
    {
        int a = Trans::transform[ 0 ];
        int sgna = Trans::sgn( a );
        int a1 = Trans::del( a, 1 );
        Trans::M[ 0 ][ 0 ] = sgna * a1;
    }
    else if ( dim == 2 )
    {
        int a = Trans::transform[ 0 ];
        int b = Trans::transform[ 1 ];
        int sgna = Trans::sgn( a );
        int sgnb = Trans::sgn( b );
        int a1 = Trans::del( a, 1 );
        int a2 = Trans::del( a, 2 );
        int b1 = Trans::del( b, 1 );
        int b2 = Trans::del( b, 2 );
        Trans::M[ 0 ][ 0 ] = sgna * a1;
        Trans::M[ 1 ][ 0 ] = sgna * a2;
        Trans::M[ 0 ][ 1 ] = sgnb * b1;
        Trans::M[ 1 ][ 1 ] = sgnb * b2;
    }
    else if ( dim == 3 )
    {
        int a = Trans::transform[ 0 ];
        int b = Trans::transform[ 1 ];
        int c = Trans::transform[ 2 ];
        int sgna = Trans::sgn( a );
        int sgnb = Trans::sgn( b );
        int sgnc = Trans::sgn( c );
        int a1 = Trans::del( a, 1 );
        int a2 = Trans::del( a, 2 );
        int a3 = Trans::del( a, 3 );
        int b1 = Trans::del( b, 1 );
        int b2 = Trans::del( b, 2 );
        int b3 = Trans::del( b, 3 );
        int c1 = Trans::del( c, 1 );
        int c2 = Trans::del( c, 2 );
        int c3 = Trans::del( c, 3 );
        Trans::M[ 0 ][ 0 ] = sgna * a1;
        Trans::M[ 1 ][ 0 ] = sgna * a2;
        Trans::M[ 2 ][ 0 ] = sgna * a3;
        Trans::M[ 0 ][ 1 ] = sgnb * b1;
        Trans::M[ 1 ][ 1 ] = sgnb * b2;
        Trans::M[ 2 ][ 1 ] = sgnb * b3;
        Trans::M[ 0 ][ 2 ] = sgnc * c1;
        Trans::M[ 1 ][ 2 ] = sgnc * c2;
        Trans::M[ 2 ][ 2 ] = sgnc * c3;
    }
}

Transform::Transform()
{
    //int dim = Dim::dim;
    int dim = 1;
    this->diff.resize( dim );
    this->mul.resize( dim );
}

Transform::~Transform()
{
    ;
}

void Transform::Init()
{
    Trans::ZeroMatrix();
    Trans::transform = this->transform;
    Trans::CalcTransformMatrix();

    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            this->Mt[ i ][ j ] = Trans::M[ i ][ j ];
        }
    }
}

void Transform::MapIndex( std::vector<int> & index1, std::vector<int> & index2 )
{
    int dim = index1.size();
    for ( int m = 0; m < dim; ++ m )
    {
        this->diff[ m ] = index1[ m ] - this->begin1[ m ];
    }

    this->Multiply( diff, this->mul );

    for ( int m = 0; m < dim; ++ m )
    {
        index2[ m ] = this->mul[ m ] + this->begin2[ m ];
    }

}

void Transform::Multiply( std::vector<int> & a, std::vector<int> & b )
{
    int dim = a.size();
    for ( int i = 0; i < dim; ++ i )
    {
        b[ i ] = 0;
        for ( int j = 0; j < dim; ++ j )
        {
            b[ i ] += this->Mt[ i ][ j ] * a[ j ];
        }
    }
}
