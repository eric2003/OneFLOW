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

#include "Sod.h"
#include "CurveLine.h"
#include "CurveMesh.h"
#include "FileUtil.h"
#include "Prj.h"
#include "DataBaseIO.h"
#include "Boundary.h"
#include "HXMath.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

Sod::Sod()
{
    ;
}

Sod::~Sod()
{
    ;
}

void Sod::Run()
{
    //this->SodGrid();
    this->Theory();
}

void Sod::SodGrid()
{
    int ni = 101;
    int nj = 51;
    int nNodes = ni * nj;

    RealField2D x;
    RealField2D y;

    AllocateVector( x, ni, nj );
    AllocateVector( y, ni, nj );

    Real xl = 0.0;
    Real xr = 1.0;
    Real yl = 0.0;
    Real yr = 0.1;

    Real dx = ( xr - xl ) / ( ni - 1 );
    Real dy = ( yr - yl ) / ( nj - 1 );

    RealField xN;
    RealField yN;
    RealField zN;
    //Generate an evenly spaced, planar grid.
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            x[ i ][ j ] = xl + i * dx;
            y[ i ][ j ] = yl + j * dy;
            xN.push_back( x[ i ][ j ] );
            yN.push_back( y[ i ][ j ] );
            zN.push_back( 0 );
        }
    }

    std::fstream file;
    Prj::OpenPrjFile( file, "/grid/sod2d.grd", std::ios_base::out|std::ios_base::binary );
    int nZone = 1;
    int nk = 1;
    HXWrite( & file, nZone );
    HXWrite( & file, ni );
    HXWrite( & file, nj );
    HXWrite( & file, nk );

    HXWrite( & file, xN );
    HXWrite( & file, yN );
    HXWrite( & file, zN );

    CloseFile( file );

    Prj::OpenPrjFile( file, "/grid/sod2d.inp", std::ios_base::out );
    int solver = 1;
    std::string zName = "A";
    int nBc = 4;
    file << solver << std::endl;
    file << nZone << std::endl;
    file << ni << " " << nj << std::endl;
    file << zName << std::endl;
    file << nBc << std::endl;
    file << 1 << " " << ni  << " " << 1  << " " << 1  << " " << BC::SYMMETRY << std::endl;
    file << 1  << " " << ni  << " " << nj  << " " << nj  << " " << BC::SYMMETRY << std::endl;
    file << 1  << " " << 1  << " " << 1  << " " << nj  << " " << BC::OUTFLOW << std::endl;
    file << ni  << " " << ni  << " " << 1  << " " << nj  << " " << BC::OUTFLOW << std::endl;
    CloseFile( file );
}

void Sod::Theory()
{
//  Calcs the exact solution for the Sod's shock-tube problem.
//  See J.D. Anderson, Modern Compressible Flow (1984) for details.
//
//     --------------------------------------------
//     |                     |                    |
//     |    p4, r4, u4       |     p1, r1, u1     |   tm = 0.0
//     |                     |                    |
//     --------------------------------------------
//     xl                    xd                   xr
//                p4>p1
//
//
//     --------------------------------------------
//     |      |      |           | up       | W   |
//     |   4  |<-----|    3    --|->  2   --|-->  |
//     |      |      |           |          |     |
//     --------------------------------------------
//     xl    expansion          slip     shock    xr
//
//-----------------------------------------------------------------

    //Set constants.
    Real gam = 1.4;
    Real gm1 = gam - 1.0;
    Real gp1 = gam + 1.0;

    //Set initial states (non-dimensional).
    Real p4 = 1.0;
    Real r4 = 1.0;
    Real u4 = 0.0;

    Real p1 = 0.1;
    Real r1 = 0.125;
    Real u1 = 0.0;

    Real tol = 1.0E-05;
    Real tm = 0.20;

    //Set dimensions of shocktube.
    Real xl = 0.0;
    Real xr = 1.0;
    Real xd = 0.5;

    //Calc acoustic velocities.
    Real a1 = sqrt( gam * p1 / r1 );
    Real a4 = sqrt( gam * p4 / r4 );

    //Use a Newton-secant iteration to compute p2p1.
    Real p2p1;
    int iterr;
    this->sp2p1( gam, p1, a1, p4, a4, p2p1, iterr, tol );

    std::cout << "p2p1 = " << p2p1 << std::endl;
    Real t2t1 = p2p1 * ( gp1 / gm1 + p2p1 ) / ( 1.0 + gp1 * p2p1 / gm1 );

    Real r2r1 = ( 1.0 + gp1 * p2p1 / gm1 ) / ( gp1 / gm1 + p2p1 );

    //shock-wave speed.
    Real wsp = a1 * sqrt( gp1 * ( p2p1 - 1.0 ) / ( 2.0 * gam ) + 1.0 );

    //Shock location.
    xs = xd + wsp * tm;

    //State 2.
    Real p2 = p2p1 * p1;
    Real r2 = r2r1 * r1;

    //State 3.
    Real p3 = p2;

    //Isentropic between 3 and 4.
    Real r3 = r4 * pow( p3 / p4, 1.0 / gam );
    Real a3 = sqrt( gam * p3 / r3 );

    //Speed of contact discontinuity.
    Real up = 2.0 * a4 * ( 1.0 - pow( p2 / p4, 0.5 * gm1 / gam ) )/ gm1;
    Real u2 = up;
    Real u3 = up;

    //Mach numbers.
    Real rmach1 = u1 / sqrt( gam * p1 / r1 );
    Real rmach2 = u2 / sqrt( gam * p2 / r2 );
    Real rmach3 = u3 / sqrt( gam * p3 / r3 );
    Real rmach4 = u4 / sqrt( gam * p4 / r4 );

    //Location of contact discontinuity.
    xc = xd + up * tm;

    //Location of expansion region.
    xhead = xd + ( u4 - a4 ) * tm;
    xtail = xd + ( u3 - a3 ) * tm;

    //Write out some data.
    std::cout << std::endl;
    std::cout << "gamma             = " << gam << std::endl;
    std::cout << "diaphram location = " << xd << std::endl;
    std::cout << "time              = " << tm << std::endl;
    std::cout << std::endl;
    std::cout << "(1) p1 = " << p1 << std::endl;
    std::cout << "    r1 = " << r1 << std::endl;
    std::cout << "    u1 = " << u1 << std::endl;
    std::cout << "    m1 = " << rmach1 << std::endl;
    std::cout << std::endl;
    std::cout << "Shock speed    = " << wsp << std::endl;
    std::cout << "Shock location = " << xs << std::endl;
    std::cout << std::endl;
    std::cout << "(2) p2 = " << p2 << std::endl;
    std::cout << "    r2 = " << r2 << std::endl;
    std::cout << "    u2 = " << u2 << std::endl;
    std::cout << "    m2 = " << rmach2 << std::endl;
    std::cout << std::endl;
    std::cout << "Contact discontinuity speed    = " << up << std::endl;
    std::cout << "Contact discontinuity location = " << xc << std::endl;
    std::cout << std::endl;
    std::cout << "(3) p3 = " << p3;
    std::cout << "    r3 = " << r3;
    std::cout << "    u3 = " << u3;
    std::cout << "    m3 = " << rmach3;
    std::cout << std::endl;
    std::cout << "Expansion region head = " << xhead << std::endl;
    std::cout << "Expansion region tail = " << xtail << std::endl;
    std::cout << std::endl;
    std::cout << "(4) p4 = " << p4 << std::endl;
    std::cout << "    r4 = " << r4 << std::endl;
    std::cout << "    u4 = " << u4 << std::endl;
    std::cout << "    m4 = " << rmach4 << std::endl;

    //Write out to files.
    std::fstream file;
    Prj::OpenPrjFile( file, "/grid/sod_theory.dat", std::ios_base::out );
    StringField title;
    title.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    title.push_back( "variables=" );
    title.push_back( "\"x\"" );
    title.push_back( "\"r\"" );
    title.push_back( "\"u\"" );
    title.push_back( "\"m\"" );
    title.push_back( "\"p\"" );

    for ( UInt i = 0; i < title.size(); ++ i )
    {
        file << title[ i ] << std::endl;
    }

    int nxp = 21;
    int nNodes = nxp + 8;

    file << " zone  i = " << nNodes << std::endl;

    file << xl << " " << r4 << " " << u4 << " " << rmach4 << " " << p4 << std::endl;
    file << xhead << " " << r4 << " " << u4 << " " << rmach4 << " " << p4 << std::endl;

    for ( int n = 1; n <= nxp; ++ n )
    {
        Real xx = xhead + ( xtail - xhead )  * n / ( nxp + 1.0 );
        Real ux = u4 + u3 * ( xx - xhead ) / ( xtail - xhead );
        Real px = p4 * pow( 1.0 - 0.5 * gm1 * ( ux / a4 ), 2.0 * gam / gm1 );
        Real rx = r4 * pow( 1.0 - 0.5 * gm1 * ( ux / a4 ), 2.0 / gm1 );
        Real mx = ux / sqrt( gam * px / rx );
        file << xx << " " << rx << " " << ux << " " << mx << " " << px << std::endl;
    }

    file << xtail << " " << r3 << " " << u3 << " " << rmach3 << " " << p3 << std::endl;
    file << xc << " " << r3 << " " << u3 << " " << rmach3 << " " << p3 << std::endl;
    file << xc << " " << r2 << " " << u2 << " " << rmach2 << " " << p2 << std::endl;
    file << xs << " " << r2 << " " << u2 << " " << rmach2 << " " << p2 << std::endl;
    file << xs << " " << r1 << " " << u1 << " " << rmach1 << " " << p1 << std::endl;
    file << xr << " " << r1 << " " << u1 << " " << rmach1 << " " << p1 << std::endl;
    CloseFile( file );
}

void Sod::sp2p1( Real gam, Real p1, Real a1, Real p4, Real a4, Real & p2p1, int & iterr, Real tol )
{
    //Uses Newton-secant method to iterate on eqn 7.94 (Anderson, 
    //1984) to fine p2p1 across moving shock wave.

    Real gm1 = gam - 1.0;
    Real gp1 = gam + 1.0;

    //Initialize p2p1 for starting guess
    Real p2p1m = 0.9 * p4 / p1;
    Real t1 = - 2.0 * gam / gm1;

    Real t2 = gm1 * ( a1 / a4 ) * ( p2p1m - 1.0 );
    Real t3 = 2.0 * gam * ( 2.0 * gam + gp1 * ( p2p1m - 1.0 ) );
    Real fm = p4 / p1 - p2p1m * pow( 1.0 - t2 / sqrt( t3 ), t1 );

    //Perturb p2p1
    p2p1 = 0.95 * p2p1m;

    //Begin iteration
    int iter = 0;
    int itmax = 20;

    while( true )
    {
        iter = iter + 1;

        t2 = gm1 * ( a1 / a4 ) * ( p2p1 - 1.0 );
        t3 = 2.0 * gam * ( 2.0 * gam + gp1 * ( p2p1 - 1.0 ) );

        Real f = p4 / p1 - p2p1 * pow( 1.0 - t2 / sqrt( t3 ), t1 );

        std::cout << "iter, p2p1, f: " << iter << " " << p2p1 << " " << f << std::endl;

        if ( ABS( f ) <= tol || iter >= itmax ) break;

        Real p2p1n = p2p1 - f * ( p2p1 - p2p1m ) / ( f - fm );
        p2p1m = p2p1;
        fm = f;
        p2p1 = p2p1n;
    }

    //c...Check to see if maximum iterations reached
    iterr = 0;
    if ( iter < itmax )  iterr = 1;
}

EndNameSpace
