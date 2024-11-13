#include "cgnslib.h"
#include "global.h"
#include "CgnsGrid.h"
#include <iostream>
#include <vector>
#include <numbers>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include <string>
#include <map>
#include <algorithm>

class Grid;
double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );
void DumpErrorDetails( std::vector<double> & u_error );
void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );

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

class Grid;
class Field;


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
    void Init( Grid * grid )
    {
        this->ni = grid->x.size();
        std::cout << "ni = " << ni << "\n";

        std::vector<double> & x = grid->x;
        this->dx = x[ 1 ] - x[ 0 ];
        this->dt = dx / 10.0;
        this->t = 1.0;
        this->nt = static_cast<int>( t / dt ) + 1;
        std::cout << "nt = " << nt << "\n";

        Global::nt = nt;

        this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
        this->beta = this->alpha * dt / ( dx * dx );
        std::cout << "alpha = " << std::setprecision( 15 ) << this->alpha << "\n";
        std::cout << "beta = " << std::setprecision( 15 ) << this->beta << "\n";

        int nghost = 2;
        int ni_total = ni + nghost;

        u_e.resize( ni_total );
        u.resize( ni_total );
        un.resize( ni_total );

        int ist = 0;
        int ied = ni - 1;

        // 0(bc) 1(ist) 2 3 ... ni

        //int ist = 0 + 1;
        //int ied = ni;

        for ( int i = ist; i <= ied; ++ i )
        {
            u_e[ i ] = - std::exp( -t ) * std::sin( std::numbers::pi * x[i] ); //theory solution
            un[ i ] = - std::sin( std::numbers::pi * x[ i ] ); //initial condition @ t=0
        }
        //un[ ist - 1 ] = 0.0;
        //un[ ied + 1 ] = 0.0;
        un[ 0 ] = 0.0;
        un[ ni - 1 ] = 0.0;
    }

    void Solve( Grid * grid )
    {
        Boundary( grid );
        for ( int i = 1; i < ni - 1; ++ i )
        {
            u[ i ] = un[ i ] + beta * ( un[ i + 1 ] - 2.0 * un[ i ] + un[ i - 1 ] );
        }
        this->update( un, u );
    }

    void Boundary( Grid * grid )
    {
        this->InterfaceBoundary( grid );
        this->PhysicalBoundary( grid );
    }

    void PhysicalBoundary( Grid * grid )
    {
        BC * bc = Global::bcs[ grid->zoneIndex ];

        int nBFace = bc->bctypes.size();
        for ( int i = 0; i < nBFace; ++ i )
        {
            int bctype = bc->bctypes[ i ];
            int faceId = bc->faceids[ i ];
            if ( bctype == BCInflow )
            {
                u[ faceId ] = 0.0;
            }
            else if ( bctype == BCOutflow )
            {
                u[ faceId ] = 0.0;
            }
        }
    }

    void InterfaceBoundary( Grid * grid )
    {
        InterFaceZone * interfacezone = Global::interfacezones[ grid->zoneIndex ];
        int nInterfaces = interfacezone->left_u.size();
        for ( int i = 0; i < nInterfaces; ++ i )
        {
            double ul = interfacezone->left_u[ i ];
            double ur = interfacezone->right_u[ i ];
            double um = 0.5 * ( ul + ur );

            int faceId = interfacezone->face_ids[ i ];

            u[ faceId ] = um;
            int kkk = 1;
        }
    }

    void PostProcess( Grid * grid )
    {
        std::vector<double> & x = grid->x;
        //compute L2 norm of the error
        std::vector<double> u_error( ni );
        for ( int i = 0; i < ni; ++ i )
        {
            u_error[ i ] = un[ i ] - u_e[ i ];
        }
    }

    void AddData( PointFactory &ptfactory, Grid * grid, std::vector<double> &global_x, std::vector<double> &global_ue, std::vector<double> &global_un )
    {
        int ni = grid->x.size();
        for ( int i = 0; i < ni; ++ i )
        {
            double xm = grid->x[ i ];
            bool flag = ptfactory.FindPoint( Point( xm ) );
            if ( i == 0 || ( i == ni - 1 ) )
            {
                if ( ! flag )
                {
                    ptfactory.AddPoint( Point( xm ) );
                    global_x.push_back( grid->x[ i ] );
                    global_ue.push_back( this->u_e[ i ] );
                    global_un.push_back( this->un[ i ] );
                }
            }
            else
            {
                global_x.push_back( grid->x[ i ] );
                global_ue.push_back( this->u_e[ i ] );
                global_un.push_back( this->un[ i ] );
            }
        }
    }

    void update( std::vector<double> &un, std::vector<double> &u )
    {
        for ( int i = 0; i < u.size(); ++ i )
        {
            un[ i ] = u[ i ];
        }
    }
};

class Solver
{
public:
    void Run()
    {
        this->ReadGrid();
        this->InitFields();
        this->SolveMultiZones();
        this->PostProcess();
    }

    void ReadGrid()
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
                int pid = ptfactory.AddPoint( Point(grid->x[i]) );
            }
            int kkk = 1;
        }
    }

    void InitFields()
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

    void SolveMultiZones()
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

    void ExchangeInterfaceValue()
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

    void PostProcess()
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

    void PostProcess( std::vector<double> &x, std::vector<double> &u_e, std::vector<double> &un )
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

    void PrintField( std::vector<double> &f )
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
};

int main( int argc, char ** argv )
{
    Solver solver;
    solver.Run();
    return 0;
}
