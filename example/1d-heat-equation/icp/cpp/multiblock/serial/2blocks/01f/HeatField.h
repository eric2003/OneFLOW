#pragma once
#include <vector>
#include "Vec1d.h"
#include "Field.h"

class HeatField :  public Field
{
public:
    int nt;
    double dx;
    double alpha, beta;
public:
    void Init( Grid * grid );
public:
    void FTCS( Zone * zone );
    void CN( Zone * zone );
    void ICP( Zone * zone );
public:
    void Rhs( Vec1d & u, Vec1d & res );
    void InviscidResidual( Vec1d & u, Vec1d & res );
    void ViscousResidual( Vec1d & u, Vec1d & res );
    void UpdateOldField();
public:
    void DumpField( Grid * grid );
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, Vec1d & u );
};

