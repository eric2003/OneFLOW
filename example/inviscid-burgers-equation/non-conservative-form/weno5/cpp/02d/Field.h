#pragma once
#include <vector>
#include "Vec1d.h"

class Zone;
class Grid;
class Region;

class Field
{
public:
    std::vector<double> u_e;
    Vec1d u, un;
    Vec1d r;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void Init( Grid * grid );
public:
    void FTCS( Zone * zone );
    void CN( Zone * zone );
    void ICP( Zone * zone );
    void RungeKutta( Zone * zone, int istage );
    void RungeKutta3Stage0( Zone * zone );
    void RungeKutta3Stage1( Zone * zone );
    void RungeKutta3Stage2( Zone * zone );
public:
    void Rhs( Vec1d & u, Vec1d & r );

    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );

    void PhysicalBoundary( Zone * zone );
    void UpdateOldField();
};
