#pragma once
#include "Vec1d.h"
#include "Field.h"
#include <string>

template<typename... Args>  
auto SQR(Args... args) {  
    return ( ... + ( args * args ) );
}

class WenoField : public Field
{
public:
    int ni;
    int nt;
    double dx, dt, total_time;
    int ns;
public:
    void Init( Grid * grid );
public:
    void PhysicalBoundary( Zone * zone );
    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );
public:
    void RungeKutta( Zone * zone, int istage );
    void RungeKutta3Stage0( Zone * zone );
    void RungeKutta3Stage1( Zone * zone );
    void RungeKutta3Stage2( Zone * zone );
    void Rhs( Vec1d & u, Vec1d & res );
public:
    void UpdateOldField();
    void PostProcess( Grid * grid );
    void DumpField( const std::string & filename, Vec1d & x, Vec1d & u );
};

void DumpCsvFile( const std::string & filename, Vec1d & x, std::vector<Vec1d> & u );
double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int N, Vec1d & u, Vec1d & f );
void wenoR( int N, Vec1d & u, Vec1d & f );
void rhs( int N, double dx, Vec1d & u, Vec1d & r );
void boundary( int N, Vec1d & u );

void runge_kutta( int N, Vec1d & un, Vec1d & u, Vec1d & res, double dx, double dt );
void runge_kutta_stage1( int N, Vec1d & un, Vec1d & u, Vec1d & res, double dx, double dt );
void runge_kutta_stage2( int N, Vec1d & un, Vec1d & u, Vec1d & res, double dx, double dt );
void runge_kutta_stage3( int N, Vec1d & un, Vec1d & u, Vec1d & res, double dx, double dt );

void numerical();
