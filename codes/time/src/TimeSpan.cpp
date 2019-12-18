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

#include "TimeSpan.h"
#include "HXDefine.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

TimeSpan::TimeSpan()
{
    ResetTime();
}

TimeSpan::~TimeSpan()
{
    ;
}

void TimeSpan::ShowTimeSpan( const string & title )
{
    t = clock();
    clock_t timeSpan = t - t_old;
    t_old = t;
    cout << title << " Time elapsed : " << static_cast< Real >( timeSpan ) / CLOCKS_PER_SEC << " seconds" << "\n";
}

void TimeSpan::ResetTime()
{
    t_old = clock();
    t = t_old;
}

EndNameSpace