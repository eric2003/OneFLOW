#pragma once
#include <vector>
#include "global.h"
#include "CgnsGrid.h"
#include "Field.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );
void DumpErrorDetails( std::vector<double> & u_error );
void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );

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