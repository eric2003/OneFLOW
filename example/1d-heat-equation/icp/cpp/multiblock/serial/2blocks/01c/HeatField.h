#pragma once
#include <vector>
#include "Vec1d.h"
#include "Field.h"

class HeatField :  public Field
{
public:
    int nt;
    double dx, total_time;
    double alpha, beta;
public:
    void Init( Grid * grid );
    void InitHeatEquation( Grid * grid );
public:
    void FTCS( Zone * zone );
    void CN_Old( Zone * zone );
    void CN( Zone * zone );
    void ICP( Zone * zone );
public:
    void Rhs( Vec1d & u, Vec1d & r );

    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );

    void PhysicalBoundary( Zone * zone );
    void UpdateOldField();
public:
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, Vec1d & u );
};

