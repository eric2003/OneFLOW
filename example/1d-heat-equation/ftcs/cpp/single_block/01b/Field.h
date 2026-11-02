#pragma once
#include <vector>
#include <iostream>
#include <iomanip> 
#include <numbers>
#include "global.h"
#include "cgnslib.h"
#include "CgnsGrid.h"

class Zone;
class Grid;

class Para
{
public:
    static int nt;
    static double dx, dt, t;
    static double alpha, beta;
public:
    static void Init();
};

class Field
{
public:
    std::vector<double> u_e;
    std::vector<double> u, un;
    std::vector<double> error;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void ModifyGrid( Grid * grid );
    void Init( Grid * grid );

    void Solve( Zone * zone );
    void Boundary( Region & region, int bcType );
    void InflowBc( Region & region );
    void OutflowBc( Region & region );

    void PhysicalBoundary( Zone * zone );

    void PostProcess( Grid * grid );
    void Update( std::vector<double> & un, std::vector<double> & u );
};