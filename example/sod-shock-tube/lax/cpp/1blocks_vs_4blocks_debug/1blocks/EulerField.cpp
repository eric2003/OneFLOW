#include "EulerField.h"
#include "Weno.h"
#include "hxmath.h"
#include "Global.h"
#include "Grid.h"
#include "cgnslib.h"
#include <format>
#include <numbers>
#include <print>
#include <fstream>
#include <iostream>

void EulerField::InitFieldCommon( Grid * grid )
{
    this->grid = grid;
    this->nequ = Global::nequ;
    this->ni = grid->ni;
    this->nic = grid->nic;
    grid->CalcMetrics();
    if ( Global::ifinite_volume == 1 )
    {
        this->nx = this->nic;
    }
    else
    {
        // finite difference
        this->nx = this->ni;
    }

    Vec1d & x = grid->x;
    this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    this->dt = Global::dt;
    this->nt = std::round( Global::total_time / dt );

    std::print( "ni={}\n", ni );
    std::print( "dt={}\n", dt );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );
    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "nt * dt = " << nt * dt << "\n";

    Global::nt = nt;

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    this->u.Allocate( this->nequ, ist, ied );
    this->un.Allocate( this->nequ, ist, ied );
    this->res.Allocate( this->nequ, 0, this->nx ); //N+1

    gamma = 1.4; // specific gas ratio
}

void EulerField::InitFieldAsRestart( Grid * grid )
{
    this->InitSodShockTube( grid );
}

void EulerField::InitSodShockTube( Grid * grid )
{
    //Sod's Riemann problem
    // Left side
    double rhoL = 1.0;
    double uL = 0.0;
    double pL = 1.0;
    // Right side
    double rhoR = 0.125;
    double uR = 0.0;
    double pR = 0.1;

    double xc = 0.5; //seperator location

    if ( Global::ifinite_volume == 0 )
    {
    }
    else
    {
        //cell center
        Vec1d & xcc = grid->xcc;
        Vec1d & q0 = this->u.vec( 0 );
        Vec1d & q1 = this->u.vec( 1 );
        Vec1d & q2 = this->u.vec( 2 );

        double rho, u, p, e;
        //i=0,1,...,nx-1
        for ( int i = 0; i < nic; ++ i )
        {
            if ( xcc[ i ] > xc )
            {
                rho = rhoR;
                u = uR;
                p = pR;
            }
            else
            {
                rho = rhoL;
                u = uL;
                p = pL;
            }
            e = p / ( rho * ( gamma - 1.0 ) ) + 0.5 * u * u;

            //conservative variables
            q0[ i ] = rho;
            q1[ i ] = rho * u;
            q2[ i ] = rho * e;
        }
    }

    int kkk = 1;
}

void EulerField::ReadFlowField( std::fstream & file, Grid * grid )
{
    if ( Global::ifinite_volume == 1 )
    {
        this->ReadFlowField( file, grid->xcc, u );
    }
    else
    {
        this->ReadFlowField( file, grid->x, u );
    }
}

void EulerField::ReadFlowField( std::fstream & file, Vec1d & x, VecWrap & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::string line; 
        std::getline( file, line );
        std::stringstream ss( line );
        std::string item;
        std::vector<std::string> row;
        while ( std::getline(ss, item, ' ') )
        {  
            row.push_back( item );
        }
        double rho  = std::atof( row[ 1 ].data() );
        double rhou = std::atof( row[ 2 ].data() );
        double rhoe = std::atof( row[ 3 ].data() );
        this->u[ 0 ][ i ] = rho;
        this->u[ 1 ][ i ] = rhou;
        this->u[ 2 ][ i ] = rhoe;
        if ( std::isnan( rho ) )
        {
            int kkk = 1;
        }
    }
    int kkk = 1;
}

void EulerField::InviscidResidual( VecWrap & u, VecWrap & res )
{
    this->InviscidConservative( u, res );
}

void EulerField::WaveSpeed( VecWrap & qL, VecWrap & qR, Vec1d & ps )
{
    //spectral radius of Jacobian
    double gm1 = gamma - 1.0;
    for ( int i = 0; i <= nx; ++ i )
    {
        // left state
        double rhoL = qL[ 0 ][ i ];
        double uL = qL[ 1 ][ i ] / rhoL;
        double eL = qL[ 2 ][ i ] / rhoL;
        double pL = gm1 * ( rhoL * eL - 0.5 * rhoL * ( uL * uL ) );
        double hL = eL + pL / rhoL;

        // Right state;
        double rhoR = qR[ 0 ][ i ];
        double uR = qR[ 1 ][ i ] / rhoR;
        double eR = qR[ 2 ][ i ] / rhoR;
        double pR = gm1 * ( rhoR * eR - 0.5 * rhoR * ( uR * uR ) );
        double hR = eR + pR / rhoR;

        double alpha = 1.0 / ( std::sqrt( std::abs( rhoL ) ) + std::sqrt( std::abs( rhoR ) ) );

        double ubar = ( std::sqrt( std::abs( rhoL ) ) * uL + std::sqrt( std::abs( rhoR ) ) * uR ) * alpha;
        double hbar = ( std::sqrt( std::abs( rhoL ) ) * hL + std::sqrt( std::abs( rhoR ) ) * hR ) * alpha;
        double cbar = std::sqrt( std::abs( gm1 * ( hbar - 0.5 * ubar * ubar ) ) );

        ps[ i ] = std::abs( cbar + ubar );
    }

    for ( int i = ps.ist; i < 0; ++ i )
    {
        ps[ i ] = ps[ 0 ];
    }

    for ( int i = nx + 1; i <= ps.ied; ++ i )
    {
        ps[ i ] = ps[ nx ];
    }
}

//void EulerField::LaxWaveSpeed( VecWrap & q, Vec1d & ps )
//{
//    //spectral radius of Jacobian
//    double gm1 = gamma - 1.0;
//    for ( int i = 0; i <= nx; ++ i )
//    {
//        // left state
//        double rhom = q[ 0 ][ i ];
//        double um = q[ 1 ][ i ] / rhom;
//        double em = q[ 2 ][ i ] / rhom;
//        double pm = gm1 * ( rhom * em - 0.5 * rhom * ( um * um ) );
//        double hm = em + pm / rhom;
//
//        double ubar = std::sqrt( std::abs( rhom ) ) * um;
//        double hbar = std::sqrt( std::abs( rhom ) ) * hm;
//        double cbar = std::sqrt( std::abs( gm1 * ( hbar - 0.5 * ubar * ubar ) ) );
//
//        ps[ i ] = std::abs( cbar + ubar );
//    }
//
//    for ( int i = ps.ist; i < 0; ++ i )
//    {
//        ps[ i ] = ps[ 0 ];
//    }
//
//    for ( int i = nx + 1; i <= ps.ied; ++ i )
//    {
//        ps[ i ] = ps[ nx ];
//    }
//}

void EulerField::LaxWaveSpeed( VecWrap & q, Vec1d & ps )
{
    //spectral radius of Jacobian
    double gm1 = gamma - 1.0;
    int ist = ps.ist;
    int ied = ps.ied;
    for ( int i = ist; i <= ied; ++ i )
    {
        // left state
        double rhom = q[ 0 ][ i ];
        double um = q[ 1 ][ i ] / rhom;
        double em = q[ 2 ][ i ] / rhom;
        double pm = gm1 * ( rhom * em - 0.5 * rhom * ( um * um ) );
        double hm = em + pm / rhom;

        double ubar = std::sqrt( std::abs( rhom ) ) * um;
        double hbar = std::sqrt( std::abs( rhom ) ) * hm;
        double cbar = std::sqrt( std::abs( gm1 * ( hbar - 0.5 * ubar * ubar ) ) );

        ps[ i ] = std::abs( cbar + ubar );
    }
}

void EulerField::LaxFriedrichs( VecWrap & u, VecWrap & res )
{
    auto lambda_value = [&](double xm, bool &flag, int &im ) { 
        im = -1;
        flag = false;
        for ( int i = 0; i < nic; ++ i )
        {
            if ( std::abs( grid->xcc[ i ] - xm ) < 1.0e-5 )
            {
                im = i;
                flag = true;
                int kkk = 1;
            }
        }
        };  
    bool flag1 = false;
    int im1 = -1;
    bool flag2 = false;
    int im2 = -1;

    lambda_value( 0.748046875, flag1, im1 );
    lambda_value( 0.751953125, flag2, im2 );

    if ( flag2  && ( Global::iter + 1) == 1441 )
    {
        int kkk = 1;
    }

    if ( flag1 && ( Global::iter + 1) == 1441 )
    {
        int kkk = 1;
    }
    VecWrap fL, fR;
    fL.Allocate( this->nequ, 0, nx, 0 );
    fR.Allocate( this->nequ, 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    VecWrap f, fP, fN;

    f.Allocate( this->nequ, ist, ied, 0 );
    fP.Allocate( this->nequ, ist, ied, 0 );
    fN.Allocate( this->nequ, ist, ied, 0 );

    euler_fluxes( ist, ied, u, f );

    Vec1d ps;
    ps.Allocate( ist, ied, 0 );
    LaxWaveSpeed( u, ps );

    // left and right side fluxes at the interface
    for ( int m = 0; m < nequ; ++ m )
    {
        for ( int i = ist; i <= ied; ++ i )
        {
            fP[ m ][ i ] = 0.5 * ( f[ m ][ i ] + ps[ i ] * u[ m ][ i ] );
            fN[ m ][ i ] = 0.5 * ( f[ m ][ i ] - ps[ i ] * u[ m ][ i ] );
        }
    }

    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, fP, fL );
        crwenoR( nx, fN, fR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, fP, fL );
        wenoR( nx, fN, fR );
    }

    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & FL = fL.vec( m );
        Vec1d & FR = fR.vec( m );
        Vec1d & res = this->res.vec( m );
        for ( int i = 0; i < nx; ++ i )
        {
            res[ i ] -= ( FL[ i + 1 ] - FL[ i ] ) / dx + ( FR[ i + 1 ] - FR[ i ] ) / dx;
        }
    }

    if ( flag1 && ( Global::iter + 1) == 1441 )
    {
        int kkk = 1;
    }
    if ( flag2 && ( Global::iter + 1) == 1441 )
    {
        int kkk = 1;
    }
}

void EulerField::burgers_fluxes( int ist, int ied, VecWrap & u, VecWrap & f )
{
    Vec1d & F = f.vec();
    Vec1d & U = u.vec();
    for ( int i = ist; i <= ied; ++ i )
    {
        F[ i ] = 0.5 * U[ i ] * U[ i ];
    }
}

//Calculate fluxes
void EulerField::euler_fluxes( int ist, int ied, VecWrap & q, VecWrap & f )
{
    //i=0,1,...,nx
    for ( int i = ist; i <= ied; ++ i )
    {
        double rho = q[ 0 ][ i ];
        double rhou = q[ 1 ][ i ];
        double rhoe = q[ 2 ][ i ];
        double p = ( gamma - 1.0 ) * ( rhoe - 0.5 * SQR( rhou ) / rho );
        f[ 0 ][ i ] = rhou;
        f[ 1 ][ i ] = rhou * rhou / rho + p;
        f[ 2 ][ i ] = rhou * rhoe / rho + p * rhou / rho;
        if ( std::isnan( rho ) )
        {
            int kkk = 1;
        }

    }
}

void EulerField::rusanov_flux( VecWrap & qL, VecWrap & qR, VecWrap & fL, VecWrap & fR, VecWrap & f )
{
    Vec1d ps;
    ps.Allocate( 0, nx, 0 );

    WaveSpeed( qL, qR, ps );

    for ( int m = 0; m < nequ; ++ m )
    {
        for ( int i = 0; i <= nx; ++ i )
        {
            //Interface fluxes (Rusanov)
            f[ m ][ i ] = 0.5 * ( fR[ m ][ i ] + fL[ m ][ i ] ) - 0.5 * ps[ i ] * ( qR[ m ][ i ] - qL[ m ][ i ] );
        }
    }
}

void EulerField::hllc_flux( VecWrap & qL, VecWrap & qR, VecWrap & fL, VecWrap & fR, VecWrap & f )
{
    std::vector<double> Ds( 3, 0 );
    double gm = gamma - 1.0;

    // i=0,1,...,nx
    for ( int i = 0; i <= nx; ++ i )
    {
        //left state
        double rhoL = qL[ 0 ][ i ];
        double uL = qL[ 1 ][ i ] / rhoL;
        double eL = qL[ 2 ][ i ] / rhoL;
        double pL = gm * ( rhoL * eL - 0.5 * rhoL * ( uL * uL ) );
        double hL = eL + pL / rhoL;
        double cL = std::sqrt( std::abs( gamma * pL / rhoL ) );

        //right state
        double rhoR = qR[ 0 ][ i ];
        double uR = qR[ 1 ][ i ] / rhoR;
        double eR = qR[ 2 ][ i ] / rhoR;
        double pR = gm * ( rhoR * eR - 0.5 * rhoR * ( uR * uR ) );
        double hR = eR + pR / rhoR;
        double cR = std::sqrt( std::abs( gamma * pR / rhoR ) );

        //compute SL and Sr
        double SL = std::min( uL, uR ) - std::max( cL, cR );
        double SR = std::max( uL, uR ) + std::max( cL, cR );

        //compute compound speed
        double term1 = pR - pL + rhoL * uL * ( SL - uL ) - rhoR * uR * ( SR - uR );
        double term2 = rhoL * ( SL - uL ) - rhoR * ( SR - uR );
        double SP = term1 / term2; //never get zero;

        //compute compound pressure
        double PLR = 0.5 * ( pL + pR + rhoL * ( SL - uL ) * ( SP - uL ) + rhoR * ( SR - uR ) * ( SP - uR ) );

        //compute D
        Ds[ 0 ] = 0.0;
        Ds[ 1 ] = 1.0;
        Ds[ 2 ] = SP;

        if ( std::isnan( fL[ 0 ][ i ] ) || std::isnan( fR[ 0 ][ i ] ) )
        {
            int kkk = 1;
        }

        if ( SL >= 0.0 )
        {
            for ( int m = 0; m < nequ; ++ m )
            {
                f[ m ][ i ] = fL[ m ][ i ];
            }
        }
        else if ( SR <= 0.0 )
        {
            for ( int m = 0; m < nequ; ++ m )
            {
                f[ m ][ i ] = fR[ m ][ i ];
            }
        }
        else if ( ( SP >= 0.0 ) && ( SL <= 0.0 ) )
        {
            for ( int m = 0; m < nequ; ++ m )
            {
                f[ m ][ i ] = ( SP * ( SL * qL[ m ][ i ] - fL[ m ][ i ] ) + SL * PLR * Ds[ m ] ) / ( SL - SP );
            }
        }
        else if ( ( SP <= 0.0 ) && ( SR >= 0.0 ) )
        {
            for ( int m = 0; m < nequ; ++ m )
            {
                f[ m ][ i ] = ( SP * ( SR * qR[ m ][ i ] - fR[ m ][ i ] ) + SR * PLR * Ds[ m ] ) / ( SR - SP );
            }
        }
    }
}

void EulerField::roe_flux( VecWrap & qL, VecWrap & qR, VecWrap & fL, VecWrap & fR, VecWrap & f )
{
    std::vector<double> dd( 3, 0 );
    std::vector<double> dF( 3, 0 );
    std::vector<double> dQ( 3, 0 );
    double gm = gamma - 1.0;

    // i=0,1,...,nx
    for ( int i = 0; i <= nx; ++ i )
    {
        //Left and right states:
        double rhoL = qL[ 0 ][ i ];
        double uL = qL[ 1 ][ i ] / rhoL;
        double eL = qL[ 2 ][ i ] / rhoL;
        double pL = gm * ( rhoL * eL - 0.5 * rhoL * ( uL * uL ) );
        double hL = eL + pL / rhoL;

        double rhoR = qR[ 0 ][ i ];
        double uR = qR[ 1 ][ i ] / rhoR;
        double eR = qR[ 2 ][ i ] / rhoR;
        double pR = gm * ( rhoR * eR - 0.5 * rhoR * ( uR * uR ) );
        double hR = eR + pR / rhoR;

        double alpha = 1.0 / ( std::sqrt( std::abs( rhoL ) ) + std::sqrt( std::abs( rhoR ) ) );

        double uu = ( std::sqrt( std::abs( rhoL ) ) * uL + std::sqrt( std::abs( rhoR ) ) * uR ) * alpha;
        double hh = ( std::sqrt( std::abs( rhoL ) ) * hL + std::sqrt( std::abs( rhoR ) ) * hR ) * alpha;
        double aa = std::sqrt( std::abs( gm * ( hh - 0.5 * uu * uu ) ) );

        double D11 = std::abs( uu );
        double D22 = std::abs( uu + aa );
        double D33 = std::abs( uu - aa );

        double beta = 0.5 / ( aa * aa );
        double phi2 = 0.5 * gm * uu * uu;

        //Right eigenvector matrix
        double R11 = 1.0;
        double R21 = uu;
        double R31 = phi2 / gm;
        double R12 = beta;
        double R22 = beta * ( uu + aa );
        double R32 = beta * ( hh + uu * aa );
        double R13 = beta;
        double R23 = beta * ( uu - aa );
        double R33 = beta * ( hh - uu * aa );

        //Left eigenvector matrix
        double L11 = 1.0 - phi2 / ( aa * aa );
        double L12 = gm * uu / ( aa * aa );
        double L13 = -gm / ( aa * aa );

        double L21 = phi2 - uu * aa;
        double L22 = - gm * uu + aa;
        double L23 = gm;

        double L31 = phi2 + uu * aa;
        double L32 = - gm * uu - aa;
        double L33 = gm;

        for ( int m = 0; m < nequ; ++ m )
        {
            dQ[ m ] = qR[ m ][ i ] - qL[ m ][ i ];
        }

        dd[ 0 ] = D11 * ( L11 * dQ[ 0 ] + L12 * dQ[ 1 ] + L13 * dQ[ 2 ] );
        dd[ 1 ] = D22 * ( L21 * dQ[ 0 ] + L22 * dQ[ 1 ] + L23 * dQ[ 2 ] );
        dd[ 2 ] = D33 * ( L31 * dQ[ 0 ] + L32 * dQ[ 1 ] + L33 * dQ[ 2 ] );

        dF[ 0 ] = R11 * dd[ 0 ] + R12 * dd[ 1 ] + R13 * dd[ 2 ];
        dF[ 1 ] = R21 * dd[ 0 ] + R22 * dd[ 1 ] + R23 * dd[ 2 ];
        dF[ 2 ] = R31 * dd[ 0 ] + R32 * dd[ 1 ] + R33 * dd[ 2 ];

        for ( int m = 0; m < nequ; ++ m )
        {
            //Interface fluxes (Roe)
            f[ m ][ i ] = 0.5 * (  fR[ m ][ i ] + fL[ m ][ i ] ) - 0.5 * dF[ m ];
        }
    }
}

void EulerField::Hllc( VecWrap & u, VecWrap & res )
{
    VecWrap uL, uR;
    uL.Allocate( this->nequ, 0, nx, 0 );
    uR.Allocate( this->nequ, 0, nx, 0 );

    //WENO Reconstruction
    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, u, uL );
        crwenoR( nx, u, uR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, u, uL );
        wenoR( nx, u, uR );
    }

    //left and right side fluxes at the interface
    VecWrap fL, fR;
    fL.Allocate( this->nequ, 0, nx, 0 );
    fR.Allocate( this->nequ, 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    //Computing fluxes
    euler_fluxes( 0, nx, uL, fL );
    euler_fluxes( 0, nx, uR, fR );

    //fluxes at the interface   
    VecWrap f;
    f.Allocate( this->nequ, 0, nx, 0 );

    //compute Riemann solver using Roe scheme(flux at interface)
    hllc_flux( uL, uR, fL, fR, f );

    //Interface fluxes (Rusanov)
    for ( int m = 0; m < nequ; ++ m )
    {
        for ( int i = 0; i < nx; ++ i )
        {
            res[ m ][ i ] -= ( f[ m ][ i + 1 ] - f[ m ][ i ] ) / dx;
            if ( std::isnan( res[ m ][ i ] ) )
            {
                int kkk = 1;
            }
        }
    }
}

void EulerField::Roe( VecWrap & u, VecWrap & res )
{
    VecWrap uL, uR;
    uL.Allocate( this->nequ, 0, nx, 0 );
    uR.Allocate( this->nequ, 0, nx, 0 );

    //WENO Reconstruction
    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, u, uL );
        crwenoR( nx, u, uR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, u, uL );
        wenoR( nx, u, uR );
    }

    //left and right side fluxes at the interface
    VecWrap fL, fR;
    fL.Allocate( this->nequ, 0, nx, 0 );
    fR.Allocate( this->nequ, 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    //Computing fluxes
    euler_fluxes( 0, nx, uL, fL );
    euler_fluxes( 0, nx, uR, fR );

    //fluxes at the interface   
    VecWrap f;
    f.Allocate( this->nequ, 0, nx, 0 );

    //compute Riemann solver using Roe scheme(flux at interface)
    roe_flux( uL, uR, fL, fR, f );

    //Interface fluxes (Rusanov)
    for ( int m = 0; m < nequ; ++ m )
    {
        for ( int i = 0; i < nx; ++ i )
        {
            res[ m ][ i ] -= ( f[ m ][ i + 1 ] - f[ m ][ i ] ) / dx;
        }
    }

}

void EulerField::Rusanov( VecWrap & u, VecWrap & res )
{
    VecWrap uL, uR;
    uL.Allocate( this->nequ, 0, nx, 0 );
    uR.Allocate( this->nequ, 0, nx, 0 );

    //WENO Reconstruction
    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, u, uL );
        crwenoR( nx, u, uR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, u, uL );
        wenoR( nx, u, uR );
    }

    //left and right side fluxes at the interface
    VecWrap fL, fR;
    fL.Allocate( this->nequ, 0, nx, 0 );
    fR.Allocate( this->nequ, 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    //Computing fluxes
    euler_fluxes( 0, nx, uL, fL );
    euler_fluxes( 0, nx, uR, fR );

    //fluxes at the interface   
    VecWrap f;
    f.Allocate( this->nequ, 0, nx, 0 );

    //compute Riemann solver using HLLC scheme
    rusanov_flux( uL, uR, fL, fR, f );

    //Interface fluxes (Rusanov)
    for ( int m = 0; m < nequ; ++ m )
    {
        for ( int i = 0; i < nx; ++ i )
        {
            res[ m ][ i ] -= ( f[ m ][ i + 1 ] - f[ m ][ i ] ) / dx;
        }
    }

}

void EulerField::InviscidConservative( VecWrap & u, VecWrap & res )
{
    if ( Global::scheme.inviscid == to_int( BasicScheme::HLLC ) )
    {
        this->Hllc( u, res );
    }
    else if ( Global::scheme.inviscid == to_int( BasicScheme::LAX ) )
    {
        this->LaxFriedrichs( u, res );
    }
    else if ( Global::scheme.inviscid == to_int( BasicScheme::Roe ) )
    {
        this->Roe( u, res );
    }
    else if ( Global::scheme.inviscid == to_int( BasicScheme::Rusanov ) )
    {
        this->Rusanov( u, res );
    }
}

void EulerField::ViscousResidual( VecWrap & u, VecWrap & res )
{
    ;
}

void EulerField::Rhs( VecWrap & u, VecWrap & res )
{
    res = 0;
    InviscidResidual( u, res );
    if ( Global::iviscous > 0 )
    {
        ViscousResidual( u, res );
    }
}

void EulerField::UpdateOldField()
{
    this->un = this->u;
}

void EulerField::DumpField( Grid * grid )
{
    if ( Global::ifinite_volume == 1 )
    {
        this->DumpField( grid->xcc, u );
    }
    else
    {
        this->DumpField( grid->x, u );
    }
}

void EulerField::PostProcess( Grid * grid )
{
    if ( Global::ifinite_volume == 1 )
    {
        this->DumpField( grid->xcc, u );
    }
    else
    {
        this->DumpField( grid->x, u );
    }
}

//void EulerField::DumpField( Vec1d & x, VecWrap & u )
//{
//    for ( int i = 0; i < x.size(); ++ i )
//    {
//        Global::file_string += std::format( "{:.25f}", x[ i ] );
//        for ( int m = 0; m < nequ; ++ m )
//        {
//            Vec1d & u = this->u.vec( m );
//            Global::file_string += std::format( " {:.25f}", u[ i ] );
//        }
//        Global::file_string += std::format( "\n" );
//    }
//}

void EulerField::DumpField( Vec1d & x, VecWrap & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.25f}", x[ i ] );
        double rho  = u[ 0 ][ i ];
        double rhou = u[ 1 ][ i ];
        double rhoe = u[ 2 ][ i ];
        double um = rhou / rho;
        Global::file_string += std::format( " {:.25f}", rho );
        Global::file_string += std::format( " {:.25f}", rhou );
        Global::file_string += std::format( " {:.25f}", rhoe );
        Global::file_string += std::format( " {:.25f}", um );

        Global::file_string += std::format( "\n" );
    }
}

void EulerField::Boundary( Region &region, int bcType )
{
    if ( bcType == BCInflow )
    {
        this->InflowBc( region );
    }
    else if ( bcType == BCExtrapolate || bcType == BCOutflow )
    {
        this->OutflowBc( region );
    }
}

void EulerField::InflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];

    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        int ib = i - 1; //index from 0
        if ( i == 1 )
        {
            idir = -1;
        }
        else
        {
            ib -= Global::ifinite_volume;
        }
        int in = ib;
        int ig1 = ib + idir;

        for ( int m = 0; m < nequ; ++ m )
        {
            u[ m ][ ig1 ] = u[ m ][ in ];
        }

        if ( Global::nghost >= 2 )
        {
            int ig2 = ig1 + idir;
            for ( int m = 0; m < nequ; ++ m )
            {
                u[ m ][ ig2 ] = u[ m ][ in ];
            }

            if ( Global::nghost >= 3 )
            {
                int ig3 = ig2 + idir;
                for ( int m = 0; m < nequ; ++ m )
                {
                    u[ m ][ ig3 ] = u[ m ][ in ];
                }
            }
        }
    }
}

void EulerField::OutflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];
    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        int ib = i - 1; //index from 0
        if ( i == 1 )
        {
            idir = -1;
        }
        else
        {
            ib -= Global::ifinite_volume;
        }
        int in = ib;
        int ig1 = ib + idir;

        for ( int m = 0; m < nequ; ++ m )
        {
            u[ m ][ ig1 ] = u[ m ][ in ];
        }

        if ( Global::nghost >= 2 )
        {
            int ig2 = ig1 + idir;
            for ( int m = 0; m < nequ; ++ m )
            {
                u[ m ][ ig2 ] = u[ m ][ in ];
            }

            if ( Global::nghost >= 3 )
            {
                int ig3 = ig2 + idir;
                for ( int m = 0; m < nequ; ++ m )
                {
                    u[ m ][ ig3 ] = u[ m ][ in ];
                }
            }
        }
    }
}