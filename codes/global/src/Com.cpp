/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
	Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
	this->fnx   = ( * ug.fnx   )[ ug.fId ];
	this->fny   = ( * ug.fny   )[ ug.fId ];
	this->fnz   = ( * ug.fnz   )[ ug.fId ];
	this->fvn   = ( * ug.fvn   )[ ug.fId ];
	this->farea = ( * ug.farea )[ ug.fId ];

    this->swapflag = false;

	if ( ug.rc == ug.cId )
	{
        this->swapflag = true;

		SWAP( ug.lc, ug.rc );
		this->Reverse();
	}

	this->ccx2   = ( * ug.ccx )[ ug.rc ];
    this->ccy2   = ( * ug.ccy )[ ug.rc ];
    this->ccz2   = ( * ug.ccz )[ ug.rc ];

	this->ccx1   = ( * ug.ccx )[ ug.lc ];
    this->ccy1   = ( * ug.ccy )[ ug.lc ];
    this->ccz1   = ( * ug.ccz )[ ug.lc ];
}

void GCom::Reverse()
{
    this->fnx = - this->fnx;
    this->fny = - this->fny;
    this->fnz = - this->fnz;
    this->fvn = - this->fvn;
	//faceArea unchanged
}

void GCom::CmpTangent()
{
    // Get first tangential
    this->idegenerate = false;
    if ( this->fnx != 0.0 )
    {
        this->t1x =   this->fny;
        this->t1y = - this->fnx;
        this->t1z =   0.0;
    }
    else if ( this->fny != 0.0 )
    {
        this->t1x = - this->fny;
        this->t1y =   this->fnx;
        this->t1z =   0.0;
    }
    else if ( this->fnz != 0.0 )
    {
        this->t1x =   0.0;
        this->t1y = - this->fnz;
        this->t1z =   this->fny;
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
    gcom.t2x = gcom.fny * gcom.t1z - gcom.fnz * gcom.t1y;
    gcom.t2y = gcom.fnz * gcom.t1x - gcom.fnx * gcom.t1z;
    gcom.t2z = gcom.fnx * gcom.t1y - gcom.fny * gcom.t1x;
}

EndNameSpace