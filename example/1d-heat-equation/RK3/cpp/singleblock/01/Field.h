#pragma once
#include <vector>
#include <iostream>
#include <iomanip> 
#include <numbers>
#include "global.h"
#include "cgnslib.h"
#include "CgnsGrid.h"

class Zone;
class Grid;

class Field
{
public:
    std::vector<double> u_e;
    std::vector<double> u, un, u1;
    std::vector<double> r;
    std::vector<double> error;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void Init( Grid * grid );
    void Solve( Zone * zone );
public:
    void FTCS( Zone * zone );
    void RungeKutta3( Zone * zone );
    void Rhs( std::vector<double> & u, std::vector<double> & r );

    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );

    void PhysicalBoundary( Zone * zone );

    void PostProcess( Grid * grid );
    void UpdateOldField();
    void UpdateRungeKuttaOldField( int istage = 0 );
public:
    void RungeKutta( Zone * zone, int istage );
    void RungeKutta3Stage0( Zone * zone );
    void RungeKutta3Stage1( Zone * zone );
    void RungeKutta3Stage2( Zone * zone );
};
