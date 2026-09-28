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
        if ( val_max < std::abs( u_error[ i ] ) )
        {
            ipos = i;
            val_max = std::abs( u_error[ i ] );
        }
    }
    std::cout << " ipos = " << ipos << "\n";
    return val_max;
}

void Solver::Run()
{
    this->ReadGrid();
    this->InitTopo();
    this->InitFields();
    this->SolveMultiZones();
    this->PostProcess();
}

void Solver::ReadGrid()
{
    std::string fileName = "../heat1d2blocksv3.cgns";
    ReadCgnsGridBaseZone( fileName );
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

void Solver::InitTopo()
{
    int nZones = Global::grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Interface * interface = new Interface();
        Global::interfaces.push_back( interface );

        IJK * ijk = new IJK();
        Global::ijks.push_back( ijk );
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Zone * zone = Global::zones[ iZone ];
        Interface * interface = Global::interfaces[ iZone ];

        int nbc1to1s = zone->bc1to1s.size();
        for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
        {
            ZoneBc1To1 * bc1to1 = zone->bc1to1s[ ibc1to1 ];
            int zoneid = bc1to1->zoneid;
            int donor_zoneid = bc1to1->donor_zoneid;
            Region region;
            region.SetRegion( bc1to1->pnts );

            Region donor_region;
            donor_region.SetRegion( bc1to1->donor_pnts );

            Transform transform;
            transform.begin1 = region.start;
            transform.begin2 = donor_region.start;
            transform.transform = bc1to1->transform;
            transform.Init();

            interface->CalcInterface( &transform, region.start, region.end, donor_zoneid );
        }
    }
}

void Solver::InitFields()
{
    Para::Init();
    int nZones = Global::grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Field * field = new Field();
        Global::fields.push_back( field );
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Field * field = Global::fields[ iZone ];
        field->Init( grid );
        Zone * zone = Global::zones[ iZone ];
        field->PhysicalBoundary( zone );
    }

    ExchangeInterfaceField();

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        field->Update( field->un, field->u );
    }
}

void Solver::SolveMultiZones()
{
    for ( int it = 0; it < Global::nt; ++ it )
    {
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Field * field = Global::fields[ iZone ];
            Zone * zone = Global::zones[ iZone ];
            zone->zoneIndex = iZone;
            field->Solve( zone );
        }
        Boundary();
    }
}

void Solver::Boundary()
{
    int nZones = Global::grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->PhysicalBoundary( zone );
    }
    ExchangeInterfaceField();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->UpdateOldField();
    }
}

void Solver::ExchangeInterfaceField()
{
    int nZones = Global::interfaces.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];

        Interface * interface = Global::interfaces[ iZone ];
        int nInterFaces = interface->zoneList.size();
        int index_dim = 1;
        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            int ijkpos = index_dim * iFace;
            int i_ghost_cell = interface->ijk_ghosts[ ijkpos + 0 ];
            int i_donor_cell = interface->ijk_donors[ ijkpos + 0 ];
            
            Field * donor_field = Global::fields[ donor_zoneid ];

            double donor_value = donor_field->u[ i_donor_cell ];
            field->u[ i_ghost_cell ] = donor_value;
        }
    }
}

void Solver::PostProcess()
{
    Post post;
    post.Process();
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


void Post::Process()
{
    this->ReorderZones();
    this->GatherField();
    this->DumpField();
}

void Post::ReorderZones()
{
    int nZones = Global::grids.size();
    std::vector<double> xmin_list;
    std::vector<double> xmax_list;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        int ni = grid->x.size();
        double xmin = grid->x[ 0 ];
        double xmax = grid->x[ 0 ];
        int imin = 0;
        int imax = 0;
        for ( int i = 0; i < ni; ++ i )
        {
            if ( xmin > grid->x[ i ] )
            {
                xmin = grid->x[ i ];
                imin = i;
            }
            if ( xmax < grid->x[ i ] )
            {
                xmax = grid->x[ i ];
                imax = i;
            }
        }
        xmin_list.push_back( xmin );
        xmax_list.push_back( xmax );
    }
    std::vector<std::pair<double, int>> pairs;
    int nSize = xmin_list.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        pairs.push_back( std::make_pair( xmin_list[i], i ) );
    }

    std::sort( pairs.begin(), pairs.end() );

    for ( int i = 0; i < nSize; ++ i )
    {
        this->zoneids.push_back( pairs[i].second );
    }

    int kkk = 1;
}

void Post::GatherField()
{
    int nZones = Global::grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        int zoneid = this->zoneids[ iZone ];
        Grid * grid = Global::grids[ zoneid ];
        Field * field = Global::fields[ zoneid ];

        int ni = grid->x.size();

        int ist = 1;
        int ied = ni;

        int dir = 1;
        if ( grid->x[ 0 ] > grid->x[ 1 ] ) dir = -1;

        if ( dir == 1 )
        {
            if ( iZone == nZones - 1 )
            {
                for ( int i = ist; i <= ied; ++ i )
                {
                    double xm = grid->x[ i - ist ];
                    this->x.push_back( xm );
                    this->u_e.push_back( field->u_e[ i ] );
                    this->un.push_back( field->un[ i ] );
                    this->u.push_back( field->u[ i ] );
                }

            }
            else
            {
                for ( int i = ist; i <= ied - 1; ++ i )
                {
                    double xm = grid->x[ i - ist ];
                    this->x.push_back( xm );
                    this->u_e.push_back( field->u_e[ i ] );
                    this->un.push_back( field->un[ i ] );
                    this->u.push_back( field->u[ i ] );
                }
            }
        }
        else
        {
            if ( iZone == nZones - 1 )
            {
                for ( int i = ied; i >= ist; -- i )
                {
                    double xm = grid->x[ i - ist ];
                    this->x.push_back( xm );
                    this->u_e.push_back( field->u_e[ i ] );
                    this->un.push_back( field->un[ i ] );
                    this->u.push_back( field->u[ i ] );
                }
            }
            else
            {
                for ( int i = ied; i >= ist + 1; -- i )
                {
                    double xm = grid->x[ i - ist ];
                    this->x.push_back( xm );
                    this->u_e.push_back( field->u_e[ i ] );
                    this->un.push_back( field->un[ i ] );
                    this->u.push_back( field->u[ i ] );
                }
            }
        }
    }
}

void Post::DumpField()
{
    int ni_total = this->x.size();
    std::cout << " DumpField x.size() = " << x.size() << "\n";
    //compute L2 norm of the error
    std::vector<double> u_error( ni_total );
    for ( int i = 0; i < ni_total; ++ i )
    {
        u_error[ i ] = un[ i ] - u_e[ i ];
    }

    this->DumpErrorDetails( u_error );

    //std::string csvname = "field_final.csv";
    //this->DumpCsvFile( csvname, x, u_e, un, u_error );
    //std::string csvname_test = "field_final_test.csv";
    //this->DumpCsvFile( csvname_test, x, u_e, un, u, u_error );

    std::string csvname = "field_final0.csv";
    this->DumpCsvFile( csvname, x, u_e, un, u_error );
    std::string csvname_test = "field_final_test0.csv";
    this->DumpCsvFile( csvname_test, x, u_e, un, u, u_error );
}

void Post::DumpErrorDetails( std::vector<double> &u_error )
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

void Post::DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<double> &ue, std::vector<double> &un, std::vector<double> &uerror )
{
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "x ue un uerror\n");
    for ( int i = 0; i < x.size(); ++ i )
    {
        //std::format_to(std::ostream_iterator<char>(file), "{:.16f} {:.16f} {:.16f} {:.16f}\n", x[i], ue[i], un[i], uerror[i] );
        std::format_to(std::ostream_iterator<char>(file), "{:.20f} {:.20f} {:.20f} {:.20f}\n", x[i], ue[i], un[i], uerror[i] );
    }
    file.close();
}

void Post::DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<double> &ue, std::vector<double> &un, std::vector<double> &u, std::vector<double> &uerror )
{
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "x ue un u uerror\n");
    for ( int i = 0; i < x.size(); ++ i )
    {
        //std::format_to(std::ostream_iterator<char>(file), "{:.16f} {:.16f} {:.16f} {:.16f}\n", x[i], ue[i], un[i], uerror[i] );
        std::format_to(std::ostream_iterator<char>(file), "{:.20f} {:.20f} {:.20f} {:.20f} {:.20f}\n", x[i], ue[i], un[i], u[i], uerror[i] );
    }
    file.close();
}