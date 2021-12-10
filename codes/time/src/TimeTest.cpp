/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "TimeTest.h"
#include <iostream>
#include <thread>

BeginNameSpace( ONEFLOW )

TimeTest::TimeTest()
{
    this->Start();
}

TimeTest::~TimeTest()
{
    ;
}

void TimeTest::Start()
{
    this->startTime = std::chrono::system_clock::now();
    this->bRunning = true;
}

void TimeTest::Stop()
{
    this->endTime = std::chrono::system_clock::now();
    this->bRunning = false;
}

double TimeTest::ElapsedMilliseconds()
{
    std::chrono::time_point<std::chrono::system_clock> endTime;

    if ( bRunning )
    {
        this->endTime = std::chrono::system_clock::now();
    }

    return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
}

double TimeTest::ElapsedSeconds()
{
    return ElapsedMilliseconds() / 1000.0;
}

void TimeTest::RunTest()
{
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time in nanoseconds: "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << " ns" << std::endl;

    std::cout << "Elapsed time in microseconds: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " mus" << std::endl;

    std::cout << "Elapsed time in milliseconds: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;

    std::cout << "Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " sec";
}


EndNameSpace
