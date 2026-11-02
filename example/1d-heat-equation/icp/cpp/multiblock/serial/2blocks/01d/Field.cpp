#include "Field.h"
#include "hxmath.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>

void Field::RungeKutta( Zone * zone, int nStage, int istage )
{
    if ( nStage == 1 )
    {
        this->RungeKutta1( zone, istage );
    }
    else if ( nStage == 3 )
    {
        this->RungeKutta3( zone, istage );
    }
}

void Field::RungeKutta1( Zone * zone, int istage )
{
    this->RungeKutta3Stage0( zone );
}

void Field::RungeKutta3( Zone * zone, int istage )
{
    if ( istage == 0 )
    {
        this->RungeKutta3Stage0( zone );
        return;
    }

    if ( istage == 1 )
    {
        this->RungeKutta3Stage1( zone );
        return;
    }

    if ( istage == 2 )
    {
        this->RungeKutta3Stage2( zone );
        return;
    }
}

void Field::RungeKutta3Stage0( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = u[ i ] + dt * res[ i ];
    }
}

void Field::RungeKutta3Stage1( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = 0.75 * un[ i ] + 0.25 * u[ i ] + 0.25 * dt * res[ i ];
    }
}

void Field::RungeKutta3Stage2( Zone * zone )
{
    this->Rhs( this->u, this->res );

    double c1 = 1.0 / 3.0;
    double c2 = 2.0 / 3.0;
    double c3 = 2.0 / 3.0;

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = c1 * un[ i ] + c2 * u[ i ] + c3 * dt * res[ i ];
    }
}
