#pragma once
#include <vector>
#include <iostream>
#include <iomanip> 
#include <numbers>
#include "global.h"
#include "cgnslib.h"

class Field
{
public:
    std::vector<double> u_e;
    std::vector<double> u, un;
    std::vector<double> error;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void Init( Grid * grid );
    void Solve( Grid * grid );

    void Boundary( Grid * grid );
    void PhysicalBoundary( Grid * grid );
    void InterfaceBoundary( Grid * grid );

    void PostProcess( Grid * grid );
    void AddData( PointFactory & ptfactory, Grid * grid, std::vector<double> & global_x, std::vector<double> & global_ue, std::vector<double> & global_un );
    void update( std::vector<double> & un, std::vector<double> & u );
};