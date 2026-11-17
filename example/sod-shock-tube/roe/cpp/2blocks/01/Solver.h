#pragma once
#include <vector>
#include "Global.h"
#include "Grid.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );

class Solver
{
public:
    Solver();
    ~Solver();
public:
    std::string gridfile;
public:
    void Init();
    void Run();
    void ReadGrid();
    void InitFields();
    void InitFieldCommon();
    void InitFieldAsRestart();
    void ReadFlowField();
    void CreateField();
    void InitTopo();
    void SolveFields();
    void TimeIntegral();
    void DumpInitialFields();
    void Read_iter();
    void Dump_iter();
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
    void DumpField();
    void DumpField( const std::string & filename );
public:
    void FTCS();
    void CrankNicolsonSeries();
    void ICP();
public:
    void RungeKutta( int nStage );
    void RungeKutta( int nStage, int istage );
};
