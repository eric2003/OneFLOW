#pragma once
#include <vector>
#include <string>

double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int nx, std::vector<double> & u, std::vector<double> & f );
void wenoR( int nx, std::vector<double> & u, std::vector<double> & f );
void rhs( int nx, double dx, std::vector<double> & u, std::vector<double> & r );
