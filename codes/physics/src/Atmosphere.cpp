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

#include "Atmosphere.h"
#include "SimuCtrl.h"
#include "HXMath.h"
#include "FileIO.h"
#include <cmath>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Atmosphere atmosphere;

Atmosphere::Atmosphere()
{
    flag = false;
}

Atmosphere::~Atmosphere()
{
    ;
}

void Atmosphere::Init()
{
    if ( flag ) return;
    flag = true;
    FileIO ioFile;
    string fileName = SimuCtrl::system_root +"physics/atmosphere.txt";
    ioFile.OpenFile( fileName, ios_base::in );

    //\tÎªtab¼ü
    string keyWordSeparator = " ()\r\n\t#$,;\"";
    ioFile.SetDefaultSeparator( keyWordSeparator );

    while ( ! ioFile.ReachTheEndOfFile() )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        Real hm = ioFile.ReadNextDigit< Real >();
        Real pm = ioFile.ReadNextDigit< Real >();
        Real rm = ioFile.ReadNextDigit< Real >();
        Real tm = ioFile.ReadNextDigit< Real >();
        hList.push_back( hm );
        pList.push_back( pm );
        rList.push_back( rm );
        tList.push_back( tm );
    }

    ioFile.CloseFile();
}

void Atmosphere::GetAirPara( const Real & hKm )
{
    GetAirPara( hKm, tm, pm, rm, cm );
}

void Atmosphere::GetAirPara( const Real & hKm, Real & temperature, Real & pressure, Real & density, Real & soundSpeed )
{
    Real hMeter = hKm * 1000;

    Real rGas = 287.053;
    Real g0 = 9.80665;
    Real rp = 6.37111e6;
    Real g = SQR( rp / ( rp + hMeter ) ) * g0;

    Real gama = 1.4;

    int idx = FindIndex( hMeter );

    Real hst = hList[ idx ];
    Real tst = tList[ idx ];
    Real hed = hList[ idx + 1 ];
    Real ted = tList[ idx + 1 ];

    Real rst = rList[ idx ];
    Real pst = pList[ idx ];

    if ( idx == 1 ||
        idx == 4 ||
        idx == 7 )
    {
        temperature = tst;
        pressure = pst * exp( - g * ( hMeter - hst ) / ( rGas * tst ) );
        density = rst * exp( - g * ( hMeter - hst ) / ( rGas * tst ) );
    }
    else
    {
        Real tmp = ( ted - tst ) / ( hed - hst );
        temperature = tst + tmp * hMeter;
        pressure = pst * pow( temperature / tst, -g / ( rGas * tmp ) );
        density = rst * pow( temperature / tst, -1.0 - g / ( rGas * tmp ) );
    }

    soundSpeed = sqrt( gama * rGas * temperature );
    pressure = pressure * 100.0;

    cout << pressure << " " << temperature << " " << density << " " << soundSpeed << "\n";
}

int Atmosphere::FindIndex( Real hMeter )
{
    int n = hList.size();
    for ( int i = 0; i < n - 1; ++ i )
    {
        Real hst = hList[ i ];
        Real hed = hList[ i + 1 ];
        if ( hst <= hMeter && hMeter < hed )
        {
            return i;
        }
    }
    return n - 1;
}


EndNameSpace
