#include "Post.h"
#include "Parallel.h"
#include "ZoneState.h"
#include "Global.h"
#include "Field.h"
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <set>
#include <unordered_map>

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

void Post::Process()
{
    this->ReorderZones();
    this->GatherField();
    this->DumpField();
}

void Post::ReorderZones()
{
    std::vector<double> xmin_list;
    std::vector<double> xmax_list;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( !ZoneState::IsValid( iZone ) && !Parallel::IsServer() ) continue;

        double xmin = 0.0;
        double xmax = 0.0;

        if ( ZoneState::IsValid( iZone ) )
        {
            int local_zoneid = ZoneState::g2lzoneids[ iZone ];
            Grid * grid = Global::grids[ local_zoneid ];
            int ni = grid->x.size();
            xmin = grid->x[ 0 ];
            xmax = grid->x[ 0 ];
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
        }

        HXSendRecvData( &xmin, 1, ZoneState::pids[ iZone ], Parallel::serverid );
        HXSendRecvData( &xmax, 1, ZoneState::pids[ iZone ], Parallel::serverid );

        if ( Parallel::IsServer() )
        {
            xmin_list.push_back( xmin );
            xmax_list.push_back( xmax );
        }
    }

    if ( Parallel::IsServer() )
    {
        std::vector<std::pair<double, int>> pairs;
        int nSize = xmin_list.size();
        for ( int i = 0; i < nSize; ++ i )
        {
            pairs.push_back( std::make_pair( xmin_list[ i ], i ) );
        }

        std::sort( pairs.begin(), pairs.end() );

        for ( int i = 0; i < nSize; ++ i )
        {
            this->zoneids.push_back( pairs[ i ].second );
        }
    }

    int kkk = 1;
}

bool IsLastZone( int iZone );
bool IsLastZone( int iZone )
{
    return iZone == ZoneState::nZones - 1;
}

void Post::AddVector( std::vector<double> & a, std::vector<double> & b )
{
    for ( int i = 0; i < b.size(); ++ i )
    {
        a.push_back( b[ i ] );
    }
}

void Post::GatherField()
{
    this->zoneids.resize( ZoneState::nZones );
    HXBcastData( this->zoneids.data(), this->zoneids.size(), Parallel::serverid );
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "Post::GatherField() this->zoneids = \n";
    for ( int i = 0; i < this->zoneids.size(); ++ i )
    {
        std::cout << this->zoneids[ i ] << " ";
    }
    std::cout << "\n";

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int zoneid = this->zoneids[ iZone ];

        if ( !ZoneState::IsValid( zoneid ) && !Parallel::IsServer() ) continue;

        int nSize = -1;
        std::vector<double> local_x;
        std::vector<double> local_u_e;
        std::vector<double> local_un;
        std::vector<double> local_u;

        if ( Parallel::pid == ZoneState::pids[ zoneid ] )
        {
            int local_zoneid = ZoneState::g2lzoneids[ zoneid ];
            Grid * grid = Global::grids[ local_zoneid ];
            Field * field = Global::fields[ local_zoneid ];
            int ni = grid->x.size();

            int dir = 1;
            if ( grid->x[ 0 ] > grid->x[ 1 ] ) dir = -1;

            if ( dir == 1 )
            {
                if ( IsLastZone( iZone ) )
                {
                    for ( int i = 0; i < ni; ++ i )
                    {
                        double xm = grid->x[ i ];
                        local_x.push_back( xm );
                        local_u_e.push_back( field->u_e[ i ] );
                        local_un.push_back( field->un[ i ] );
                        local_u.push_back( field->u[ i ] );
                    }
                }
                else
                {
                    for ( int i = 0; i < ni - 1; ++ i )
                    {
                        double xm = grid->x[ i ];
                        local_x.push_back( xm );
                        local_u_e.push_back( field->u_e[ i ] );
                        local_un.push_back( field->un[ i ] );
                        local_u.push_back( field->u[ i ] );
                    }
                }
            }
            else
            {
                if ( IsLastZone( iZone ) )
                {
                    for ( int i = ni - 1; i >= 0; -- i )
                    {
                        double xm = grid->x[ i ];
                        local_x.push_back( xm );
                        local_u_e.push_back( field->u_e[ i ] );
                        local_un.push_back( field->un[ i ] );
                        local_u.push_back( field->u[ i ] );
                    }
                }
                else
                {
                    for ( int i = ni - 1; i >= 1; -- i )
                    {
                        double xm = grid->x[ i ];
                        local_x.push_back( xm );
                        local_u_e.push_back( field->u_e[ i ] );
                        local_un.push_back( field->un[ i ] );
                        local_u.push_back( field->u[ i ] );
                    }
                }
            }
            nSize = local_x.size();
        }

        int send_pid = ZoneState::pids[ zoneid ];

        HXSendRecvData( &nSize, 1, send_pid, Parallel::serverid );

        if ( Parallel::IsServer() )
        {
            local_x.resize( nSize );
            local_u_e.resize( nSize );
            local_un.resize( nSize );
            local_u.resize( nSize );
        }

        HXSendRecvData( local_x.data(), local_x.size(), send_pid, Parallel::serverid );
        HXSendRecvData( local_u_e.data(), local_u_e.size(), send_pid, Parallel::serverid );
        HXSendRecvData( local_un.data(), local_un.size(), send_pid, Parallel::serverid );
        HXSendRecvData( local_u.data(), local_u.size(), send_pid, Parallel::serverid );

        if ( Parallel::IsServer() )
        {
            this->AddVector( this->x, local_x );
            this->AddVector( this->u_e, local_u_e );
            this->AddVector( this->un, local_un );
            this->AddVector( this->u, local_u );
        }
    }
}

void Post::DumpField()
{
    if ( !Parallel::IsServer() ) return;
    int ni = this->x.size();
    std::cout << " DumpField x.size() = " << x.size() << "\n";
    //compute L2 norm of the error
    std::vector<double> u_error( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        u_error[ i ] = un[ i ] - u_e[ i ];
    }

    this->DumpErrorDetails( u_error );

    std::string csvname = "field_final.csv";
    this->DumpCsvFile( csvname, x, u_e, un, u_error );
}

void Post::DumpErrorDetails( std::vector<double> &u_error )
{
    int ni = u_error.size();
    double rms_error = compute_l2norm( ni, u_error );
    double max_error = compute_max_error( ni, u_error );
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
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
        std::format_to(std::ostream_iterator<char>(file), "{:.20f} {:.20f} {:.20f} {:.20f}\n", x[i], ue[i], un[i], uerror[i] );
    }
    file.close();
}
