#pragma once
#include <vector>
#include "Vec1d.h"

class Zone;
class Grid;
class Region;

class Field
{
public:
    Field() {}
    virtual ~Field() {};
public:
    virtual void Init( Grid * grid ) {}
    virtual void UpdateOldField() {}
    virtual void FTCS( Zone * zone ) {}
    virtual void CrankNicolsonSeries( Zone * zone );
    virtual void CN( Zone * zone ) {};
    virtual void ICP( Zone * zone ) {}
    virtual void DumpField( Grid * grid ) {}
    virtual void PostProcess( Grid * grid ) {}
    virtual void Rhs( Vec1d & u, Vec1d & r ) {};
public:
    void RungeKutta( Zone * zone, int nStage, int istage );
    void RungeKutta1( Zone * zone, int istage );
    void RungeKutta3( Zone * zone, int istage );
    void RungeKutta3Stage0( Zone * zone );
    void RungeKutta3Stage1( Zone * zone );
    void RungeKutta3Stage2( Zone * zone );
public:
    void PhysicalBoundary( Zone * zone );
    void InterfaceBoundary( Zone * zone );
    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );
    void InterfaceBc( Region & region );
public:
    Vec1d u_e;
    Vec1d u, un;
    Vec1d res;
    int ni ,nic; //nnode, ncell;
    int nx;
    double dt;
};

