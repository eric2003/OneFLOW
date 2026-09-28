#pragma once
#include "Vec1d.h"
#include "Field.h"
#include <string>

class WenoField : public Field
{
public:
    int nt;
    double dx, total_time;
    int ns;
public:
    void Init( Grid * grid );
public:
    void PhysicalBoundary( Zone * zone );
    void InterfaceBoundary( Zone * zone );
    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );
    void InterfaceBc( Region & region );
public:
    void Rhs( Vec1d & u, Vec1d & res );
public:
    void UpdateOldField();
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, Vec1d & u );
};

double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void wenoL( int N, Vec1d & u, Vec1d & f );
void wenoR( int N, Vec1d & u, Vec1d & f );

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

void crwenoL( int ni, Vec1d & u, Vec1d & f );
void crwenoR( int ni, Vec1d & u, Vec1d & f );

