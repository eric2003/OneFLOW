#pragma once
#include <vector>
#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <fstream>

template<typename... Args>  
auto SQR(Args... args) {  
    return (... + (args * args));
}

class vec1d
{
public:
    std::vector<double> data;
    int ist = 0;
public:
    void Allocate( int ist, int ied, double value = 0 )
    {
        int nelement = ied - ist + 1;
        this->data.resize( nelement, value );
        this->ist = ist;
    }
    std::size_t size()
    {
        return this->data.size();
    }

    double operator [] ( int i ) const
    {
        return data[ i - ist ];
    }

    double & operator [] ( int i )
    {
        return data[ i - ist ];
    }

    vec1d & operator = ( const vec1d & rhs )
    {
        if ( this == & rhs ) return * this;
        this->data = rhs.data;

        return * this;
    }
};


void DumpCsvFile( const std::string & filename, vec1d & x, std::vector<vec1d> & u );
double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int N, vec1d & u, vec1d & f );
void wenoR( int N, vec1d & u, vec1d & f );
void rhs( int N, double dx, vec1d & u, vec1d & r );
void boundary( int N, vec1d & u );
void numerical();

void runge_kutta( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt );
void runge_kutta_stage1( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt );
void runge_kutta_stage2( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt );
void runge_kutta_stage3( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt );

class Solver
{
public:
    Solver();
    ~Solver();
public:
    void Solve();
};
