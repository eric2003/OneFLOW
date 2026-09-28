#pragma once
#include <vector>
#include "hxmath.h"
#include "vec1d.h"

void crabL( double w1, double w2, double w3,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 );
void crabR( double w1, double w2, double w3,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 );
void crwcL( double v1, double v2, double v3, double v4, double v5,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 );
void crwcR( double v1, double v2, double v3, double v4, double v5,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 );
void crwenoL( int ni, Vec1d & u, std::vector<double> & f );
void crwenoR( int ni, Vec1d & u, std::vector<double> & f );
void rhs_crweno( int nx, double dx, Vec1d & u, Vec1d & r );
void numerical( RhsPtr rhs );

