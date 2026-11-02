#pragma once
#include <vector>
#include "Vec1d.h"
#include "Field.h"

class ConvectionField :  public Field
{
public:
    int nt;
    double dx;
    double c;
    double alpha, beta;
public:
    void InitFieldCommon( Grid * grid ) override;
    void InitFieldAsRestart( Grid * grid ) override;
    void ReadFlowField( std::fstream & file, Grid * grid ) override;
    void ReadFlowField( std::fstream & file, Vec1d & x, VecWrap & u );
public:
    void FTCS( Zone * zone );
    void CN( Zone * zone );
public:
    void Rhs( VecWrap & u, VecWrap & res );
    void InviscidResidual( VecWrap & u, VecWrap & res );
    void InviscidNonConservative( VecWrap & u, VecWrap & res );
    void InviscidConservative( VecWrap & u, VecWrap & res );
    void ViscousResidual( VecWrap & u, VecWrap & res );
    void UpdateOldField();
public:
    void DumpField( Grid * grid );
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, VecWrap & u );
};

