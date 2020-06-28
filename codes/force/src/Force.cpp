/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "Force.h"
#include "Zone.h"
#include "UnsGrid.h"
#include "BcRecord.h"
#include "FaceTopo.h"
#include "Boundary.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "DataBase.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Grad.h"
#include "VisGrad.h"
#include "NsCom.h"
#include "Parallel.h"
#include "Iteration.h"
#include "FileIO.h"
#include "INsIdx.h"
#include "INsCom.h"
#include <sstream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

AeroForceInfo aeroForceInfo;
AeroCom aeroCom;

Force::Force()
{
    ;
}

Force::~Force()
{
    ;
}

Force::Force( const Force & rhs )
{
    * this = rhs;
}

Force::Force( const Real & value )
{
    * this = value;
}

Force & Force::operator += ( const Force & rhs )
{
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;

    return *this;
}

Force & Force::operator /= ( const Real & value )
{
    this->x /= value;
    this->y /= value;
    this->z /= value;

    return *this;
}

Force & Force::operator = ( const Real & value )
{
    this->x = value;
    this->y = value;
    this->z = value;

    return *this;
}

Force & Force::operator = ( const Force & rhs )
{
    this->x = rhs.x;
    this->y = rhs.y;
    this->z = rhs.z;

    return *this;
}

Force operator + ( const Force & f1, const Force & f2 )
{
    Force f = f1;
    f += f2;
    return f;
}

Force operator / ( const Force & f, const Real & value )
{
    Force ff = f;
    ff /= value;
    return ff;
}

AeroCom::AeroCom()
{
    ;
}

AeroCom::~AeroCom()
{
    ;
}

void AeroCom::Init()
{
    this->aref = GetDataValue< Real >( "aref" );
    this->lref = GetDataValue< Real >( "lref" );
    this->xref = GetDataValue< Real >( "xref" );
    this->yref = GetDataValue< Real >( "yref" );
    this->zref = GetDataValue< Real >( "zref" );
    this->aoa  = GetDataValue< Real >( "aoa_degree" ) * PI / 180.0;
    this->aos  = GetDataValue< Real >( "sideslip_degree" ) * PI / 180.0;
    this->sina = sin( aoa );
    this->cosa = cos( aoa );
    this->sinb = sin( aos );
    this->cosb = cos( aos );
    this->cForce  = aref;
    this->cMoment = aref * lref;
}

Real AeroCom::CalcCL( Force * f )
{
    Real cl = - sina * f->x  + cosa * f->y;
    return cl;
}

Real AeroCom::CalcCD( Force * f )
{
    Real cd = cosa * cosb * f->x + sina * cosb * f->y + sinb * f->z;
    return cd;
}

Real AeroCom::CalcCF( Force * f, Real area )
{
    Real dir = cosa * f->x + sina * f->y;
    Real cf = DIST( f->x, f->y, f->z ) / area * SIGN( 1.0, dir );
    return cf;
}

AeroForce::AeroForce()
{
    ;
}

AeroForce::~AeroForce()
{
    ;
}

void AeroForce::Init()
{
    total = 0;
    pres = 0;
    vis = 0;
    power = 0;
}

void AeroForce::SumForce()
{
    total = pres + vis;
}

void AeroForce::CalcPower()
{
    power = - 1.0 * ( vfx * total.x + vfy * total.y + vfz * total.z );
}

void AeroForce::CalcMoment( Real xc, Real yc, Real zc )
{
    mom.x = ( yc - aeroCom.yref ) * total.z - ( zc - aeroCom.zref ) * total.y;
    mom.y = ( zc - aeroCom.zref ) * total.x - ( xc - aeroCom.xref ) * total.z;
    mom.z = ( xc - aeroCom.xref ) * total.y - ( yc - aeroCom.yref ) * total.x;
}

void AeroForce::AddForce( AeroForce * rhs )
{
    this->total += rhs->total;
    this->pres  += rhs->pres;
    this->vis   += rhs->vis;
    this->power += rhs->power;
}

AeroForceInfo::AeroForceInfo()
{
    ;
}

AeroForceInfo::~AeroForceInfo()
{
    ;
}

void AeroForceInfo::Init()
{
    aeroCom.Init();
    totalForce.Init();
}

void AeroForceInfo::CollectForce()
{
    this->Sum( & totalForce.total );
    this->Sum( & totalForce.mom );
    this->Sum( & totalForce.pres );
    this->Sum( & totalForce.power );
}

void AeroForceInfo::Sum( Force * force )
{
    Force tf = 0.0;

    HXReduceReal( & force->x, & tf.x, 1, PL_SUM );
    HXReduceReal( & force->y, & tf.y, 1, PL_SUM );
    HXReduceReal( & force->z, & tf.z, 1, PL_SUM );

    * force = tf;
}

void AeroForceInfo::Sum( Real * var )
{
    Real sumVar = 0.0;

    HXReduceReal( var, & sumVar, 1, PL_SUM );

    * var = sumVar;
}

void AeroForceInfo::CalcCoef()
{
    cf = totalForce.total / aeroCom.cForce;
    cpres = totalForce.pres / aeroCom.cForce;
    cmom = totalForce.mom / aeroCom.cMoment;
    cpower = totalForce.power / aeroCom.cForce;

    const Real AR = 9.5;
    cl      = aeroCom.CalcCL( & cf );
    cd      = aeroCom.CalcCD( & cf );
    cd_pres = aeroCom.CalcCD( & cpres );
    cd_vis  = cd - cd_pres;
    cdl     = cd - SQR( cl ) / ( PI * AR );

    pres_center = SIGN( 1.0, cf.y ) * cmom.z / ( ABS( cf.y ) + SMALL );
}


EndNameSpace