#include "Solver.h"
#include <fstream>

double compute_l2norm( int ni, std::vector<double> & r )
{
    double rms = 0.0;
    for ( int i = 1; i < ni - 1; ++ i )
    {
        rms += r[ i ] * r[ i ];
    }
    rms = std::sqrt( rms / ( ni - 2 ) );
    return rms;
}

double compute_max_error( int ni, std::vector<double> & u_error )
{
    double val_max = -1;
    int ipos = -1;
    for ( int i = 1; i < ni - 1; ++ i )
    {
        //val_max = std::max( val_max, std::abs( u_error[ i ] ) );
        if ( val_max < std::abs( u_error[ i ] ) )
        {
            ipos = i;
            val_max = std::abs( u_error[ i ] );
        }
    }
    std::cout << " ipos = " << ipos << "\n";
    return val_max;
}

void DumpErrorDetails( std::vector<double> &u_error )
{
    int ni = u_error.size();
    double rms_error = compute_l2norm( ni, u_error );
    double max_error = compute_max_error( ni, u_error );
    std::cout << "max_error = " << std::setprecision(15) << max_error << "\n";
    //create output file for L2-norm
    std::fstream file;
    file.open("output.txt", std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "Error details: \n");
    std::format_to(std::ostream_iterator<char>(file), "L-2 Norm = {0}\n", rms_error);
    std::format_to(std::ostream_iterator<char>(file), "Maximum Norm = {0}\n", max_error);
    file.close();
}

void DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<double> &ue, std::vector<double> &un, std::vector<double> &uerror )
{
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "x ue un uerror\n");
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::format_to(std::ostream_iterator<char>(file), "{:.16f} {:.16f} {:.16f} {:.16f}\n", x[i], ue[i], un[i], uerror[i] );
    }
    file.close();
}

void Solver::Run()
{
    this->ReadGrid();
    this->InitFields();
    this->SolveMultiZones();
    this->PostProcess();
}

void Solver::ReadGrid()
{
    std::string fileName = "../heat1d2blocksv3.cgns";
    ReadCgnsGridBaseZone( fileName );
    int kkk = 1;
    ReadCgnsGrid( fileName );
    int nZones = Global::grids.size();
    PointFactory ptfactory;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        grid->zoneIndex = iZone;
        int ni = grid->x.size();
        for ( int i = 0; i < ni; ++ i )
        {
            int pid = ptfactory.AddPoint( Point( grid->x[ i ] ) );
        }
        int kkk = 1;
    }
}

void Solver::InitFields()
{
    int nZones = Global::grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Field * field = new Field();
        Global::fields.push_back( field );
        field->Init( grid );
    }
}

void Solver::SolveMultiZones()
{
    for ( int it = 0; it < Global::nt; ++ it )
    {
        ExchangeInterfaceValue();
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            Field * field = Global::fields[ iZone ];
            field->Solve( grid );
        }
    }
}

void Solver::ExchangeInterfaceValue()
{
    int nZones = Global::grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        InterFaceZone * interfacezone = Global::interfacezones[ iZone ];
        int nInterFace = interfacezone->left_zones.size();
        interfacezone->left_u.resize( nInterFace );
        interfacezone->right_u.resize( nInterFace );
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];

        InterFaceZone * interfacezone = Global::interfacezones[ iZone ];
        int nInterFace = interfacezone->left_zones.size();

        for ( int i = 0; i < nInterFace; ++ i )
        {
            int left_cell = interfacezone->left_cells[ i ];
            int donor_zoneid = interfacezone->right_zones[ i ];
            int donor_cell = interfacezone->right_cells[ i ];
            Field * donor_field = Global::fields[ donor_zoneid ];

            interfacezone->left_u[ i ] = field->u[ left_cell ];
            interfacezone->right_u[ i ] = donor_field->u[ donor_cell ];
        }
    }
}

void Solver::PostProcess()
{
    int nZones = Global::grids.size();
    std::vector<double> u_e;
    std::vector<double> un;
    std::vector<double> x;

    PointFactory ptfactory;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Field * field = Global::fields[ iZone ];
        field->AddData( ptfactory, grid, x, u_e, un );
    }

    std::cout << " x.size() = " << x.size() << "\n";

    this->PostProcess( x, u_e, un );
}

void Solver::PostProcess( std::vector<double> &x, std::vector<double> &u_e, std::vector<double> &un )
{
    int ni = x.size();
    //compute L2 norm of the error
    std::vector<double> u_error( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        u_error[ i ] = un[ i ] - u_e[ i ];
    }

    ::DumpErrorDetails( u_error );

    std::string csvname = "field_final.csv";
    ::DumpCsvFile( csvname, x, u_e, un, u_error );
}

void Solver::PrintField( std::vector<double> &f )
{
    int icount = 0;
    for ( int i = 0; i < f.size(); ++ i )
    {
        std::cout << std::setprecision(15) << f[ i ] << " ";
        icount ++;
        if ( icount % 5 == 0 )
        {
            std::cout << "\n";
        }
    }
    std::cout << "\n";
    std::cout << "\n";
}