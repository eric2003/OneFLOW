#pragma once
#include <vector>
#include <fstream>
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
    virtual void Init( std::fstream & file, Grid * grid ) {}
    virtual void InitFieldCommon( Grid * grid ) {}
    virtual void InitFieldAsRestart( Grid * grid ) {}
    virtual void ReadFlowField( std::fstream & file, Grid * grid ) {}
    virtual void UpdateOldField() {}
    virtual void FTCS( Zone * zone ) {}
    virtual void CrankNicolsonSeries( Zone * zone );
    virtual void CN( Zone * zone ) {};
    virtual void ICP( Zone * zone ) {}
    virtual void DumpField( Grid * grid ) {}
    virtual void PostProcess( Grid * grid ) {}
    virtual void Rhs( Vec1d & u, Vec1d & r ) {};
    virtual void Rhs( VecWrap & u, VecWrap & r ) {};
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
    virtual void Boundary( Region & region, int bcType );
    virtual void InflowBc( Region & region );
    virtual void OutflowBc( Region & region );
    virtual void DirichletBc( Region & region );
    virtual void ExtrapolateBc( Region & region );
    virtual void InterfaceBc( Region & region );
public:
    Grid * grid;
    VecWrap u, un;
    VecWrap res;
    int nequ = 1;
    int ni ,nic; //nnode, ncell;
    int nx;
    double dt;
};

