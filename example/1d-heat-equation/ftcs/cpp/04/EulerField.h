#pragma once
#include "Vec1d.h"
#include "Field.h"
#include <string>

class EulerField : public Field
{
public:
    int nt;
    double dx;
public:
    void Init( Grid * grid );
public:
    void Rhs( VecWrap & u, VecWrap & res );
    void InviscidResidual( VecWrap & u, VecWrap & res );
    void ViscousResidual( VecWrap & u, VecWrap & res );
    void InviscidConservative( VecWrap & u, VecWrap & res );
    void WaveSpeed( VecWrap & um, VecWrap & psm );
public:
    void LaxFriedrichs( VecWrap & u, VecWrap & res );
    void Rusanov( VecWrap & u, VecWrap & res );
    void burgers_fluxes( int ist, int ied, VecWrap & u, VecWrap & f );
public:
    void UpdateOldField();
    void DumpField( Grid * grid );
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, VecWrap & u );
};


