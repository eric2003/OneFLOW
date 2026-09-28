#pragma once
#include <vector>
#include "Global.h"
#include "Grid.h"
#include "Field.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );

enum class Scheme
{
    FTCS = 0,
    CN, //Crank¨CNicolson
    ICP, //Implicit Compact Pade (ICP) Scheme
    WENO,
    RungeKutta
};

class Solver
{
public:
    Solver();
    ~Solver();
public:
    Scheme scheme;
    int nghost;
public:
    void Run();
    void ReadGrid();
    void InitFields();
    void InitTopo();
    void SolveFields();
    void DumpInitialFields();
public:
    void Boundary();
    void UpdateOldField();
    void UploadInterfaceField();
    void UpdateInterfaceField();
    void DownloadInterfaceField();
    void ExchangeInterfaceField();
    void PrintField( std::vector<double> & f );
public:
    void PostProcess();
public:
    void FTCS();
    void CN();
    void ICP();
    void RungeKutta();
    void RungeKutta( int istage );
};
