#pragma once
#include "Vec1d.h"
#include "Field.h"
#include <string>

class BurgersField : public Field
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
    void InviscidResidual( Vec1d & u, Vec1d & res );
    void ViscousResidual( Vec1d & u, Vec1d & res );
public:
    void UpdateOldField();
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, Vec1d & u );
};


