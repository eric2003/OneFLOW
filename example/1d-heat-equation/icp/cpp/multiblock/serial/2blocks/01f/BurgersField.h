#pragma once
#include "Vec1d.h"
#include "Field.h"
#include <string>

class BurgersField : public Field
{
public:
    int nt;
    double dx;
    int ns;
public:
    void Init( Grid * grid );
public:
    void Rhs( Vec1d & u, Vec1d & res );
    void InviscidResidual( Vec1d & u, Vec1d & res );
    void ViscousResidual( Vec1d & u, Vec1d & res );
public:
    void UpdateOldField();
    void DumpField( Grid * grid );
    void PostProcess( Grid * grid );
    void DumpField( Vec1d & x, Vec1d & u );
};


