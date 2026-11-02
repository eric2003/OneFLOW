#pragma once
#include "Vec1d.h"
#include "Field.h"
#include <string>

class EulerField : public Field
{
public:
    int nt;
    double dx;
    double gamma;
public:
    void Init( Grid * grid );
    void InitSodShockTube( Grid * grid );
public:
    void Rhs( VecWrap & u, VecWrap & res );
    void InviscidResidual( VecWrap & u, VecWrap & res );
    void ViscousResidual( VecWrap & u, VecWrap & res );
    void InviscidConservative( VecWrap & u, VecWrap & res );
    void WaveSpeed( VecWrap & qL, VecWrap & qR, Vec1d & ps );
public:
    void LaxFriedrichs( VecWrap & u, VecWrap & res );
    void Rusanov( VecWrap & u, VecWrap & res );
    void burgers_fluxes( int ist, int ied, VecWrap & u, VecWrap & f );
    void euler_fluxes( int ist, int ied, VecWrap & u, VecWrap & f );
    void rusanov_flux( VecWrap & qL, VecWrap & qR, VecWrap & fL, VecWrap & fR, VecWrap & f );
public:
    void UpdateOldField();
    void DumpField( Grid * grid );
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, VecWrap & u );
public:
    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );
};


