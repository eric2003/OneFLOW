#pragma once
#include "Vec1d.h"
#include <string>

template<typename... Args>  
auto SQR(Args... args) {  
    return (... + (args * args));
}

void DumpCsvFile( const std::string & filename, Vec1d & x, std::vector<Vec1d> & u );
double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int N, Vec1d & u, Vec1d & f );
void wenoR( int N, Vec1d & u, Vec1d & f );
void rhs( int N, double dx, Vec1d & u, Vec1d & r );
void boundary( int N, Vec1d & u );

void runge_kutta( int N, Vec1d & un, Vec1d & ut, Vec1d & res, double dx, double dt );
void runge_kutta_stage1( int N, Vec1d & un, Vec1d & ut, Vec1d & res, double dx, double dt );
void runge_kutta_stage2( int N, Vec1d & un, Vec1d & ut, Vec1d & res, double dx, double dt );
void runge_kutta_stage3( int N, Vec1d & un, Vec1d & ut, Vec1d & res, double dx, double dt );

void numerical();
