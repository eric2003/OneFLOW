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
    void Init( std::fstream & file, Grid * grid );
public:
    void FTCS( Zone * zone );
    void CN( Zone * zone );
    void ICP( Zone * zone );
public:
    void Rhs( VecWrap & u, VecWrap & res );
    void InviscidResidual( VecWrap & u, VecWrap & res );
    void ViscousResidual( VecWrap & u, VecWrap & res );
    void UpdateOldField();
public:
    void DumpField( Grid * grid );
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, VecWrap & u );
};

