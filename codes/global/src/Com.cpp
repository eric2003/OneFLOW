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
#include "Com.h"
#include "UCom.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

GCom gcom;

GCom::GCom()
{
    ;
}

GCom::~GCom()
{
    ;
}

void GCom::SetGeometry()
{
    this->xfn   = ( * ug.xfn   )[ ug.fId ];
    this->yfn   = ( * ug.yfn   )[ ug.fId ];
    this->zfn   = ( * ug.zfn   )[ ug.fId ];
    this->vfn   = ( * ug.vfn   )[ ug.fId ];
    this->farea = ( * ug.farea )[ ug.fId ];

    this->swapflag = false;

    if ( ug.rc == ug.cId )
    {
        this->swapflag = true;

        SWAP( ug.lc, ug.rc );
        this->Reverse();
    }

    this->xcc2   = ( * ug.xcc )[ ug.rc ];
    this->ycc2   = ( * ug.ycc )[ ug.rc ];
    this->zcc2   = ( * ug.zcc )[ ug.rc ];

    this->xcc1   = ( * ug.xcc )[ ug.lc ];
    this->ycc1   = ( * ug.ycc )[ ug.lc ];
    this->zcc1   = ( * ug.zcc )[ ug.lc ];
}

void GCom::Reverse()
{
    this->xfn = - this->xfn;
    this->yfn = - this->yfn;
    this->zfn = - this->zfn;
    this->vfn = - this->vfn;
    //faceArea unchanged
}

void GCom::CalcTangent()
{
    // Get first tangential
    this->idegenerate = false;
    if ( this->xfn != 0.0 )
    {
        this->t1x =   this->yfn;
        this->t1y = - this->xfn;
        this->t1z =   0.0;
    }
    else if ( this->yfn != 0.0 )
    {
        this->t1x = - this->yfn;
        this->t1y =   this->xfn;
        this->t1z =   0.0;
    }
    else if ( this->zfn != 0.0 )
    {
        this->t1x =   0.0;
        this->t1y = - this->zfn;
        this->t1z =   this->yfn;
    }
    else
    {
        this->t1x = 0.0;
        this->t1y = 0.0;
        this->t1z = 0.0;

        this->t2x = 0.0;
        this->t2y = 0.0;
        this->t2z = 0.0;
        this->idegenerate = true;
        return;
    }

    //normalize the tangential vector
    Real dtmp = 1.0 / DIST( this->t1x, this->t1y, this->t1z );
    this->t1x *= dtmp;
    this->t1y *= dtmp;
    this->t1z *= dtmp;

    //Get second tangential vector by cross dot t1 to normal
    gcom.t2x = gcom.yfn * gcom.t1z - gcom.zfn * gcom.t1y;
    gcom.t2y = gcom.zfn * gcom.t1x - gcom.xfn * gcom.t1z;
    gcom.t2z = gcom.xfn * gcom.t1y - gcom.yfn * gcom.t1x;
}

EndNameSpace