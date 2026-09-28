#pragma once
#include <vector>
#include <string>

class Vec1d;

double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int ni, Vec1d & u, std::vector<double> & f );
void wenoR( int ni, Vec1d & u, std::vector<double> & f );
void rhs( int ni, double dx, Vec1d & u, Vec1d & r );
