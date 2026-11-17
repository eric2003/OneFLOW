#include "Solver.h"
#include <fstream>
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
    std::string fileName = "../heat1d4blocksv1.cgns";
    ReadCgnsGridBaseZone( fileName );
    ReadCgnsGrid( fileName );
}

void Solver::InitTopo()
{
    int nZones = Global::grids.size();
    Global::donor_zone_sets.resize( nZones );
    Global::donor_zones.resize( nZones );
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Interface * interface = new Interface();
        interface->zoneid = iZone;
        Global::interfaces.push_back( interface );
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
        int nInterfaces = interface->zoneList.size();
        interface->data_recv.resize( nInterfaces );
        interface->data_send.resize( nInterfaces );
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];
        int nInterFaces = interface->zoneList.size();
        std::set<int> &donor_zoneSet = Global::donor_zone_sets[ iZone ];
        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            donor_zoneSet.insert( donor_zoneid );
        }

        std::vector<int> &donor_zones = Global::donor_zones[ iZone ];
        for ( std::set<int>::iterator iter = donor_zoneSet.begin(); iter != donor_zoneSet.end(); ++ iter )
        {
            donor_zones.push_back( *iter );
        }

        interface->neighbor_donor_zones = donor_zones;

        std::unordered_map<int, int> donor_zonelocal;

        for ( int idonor = 0; idonor < donor_zones.size(); ++ idonor )
        {
            int donor_zone = donor_zones[ idonor ];
            donor_zonelocal.insert( std::make_pair( donor_zone, idonor ) );
        }
        int index_dim = 1;
        int ndonors = donor_zones.size();
        std::vector<std::vector<int>> & neighbor_donorfaces = interface->neighbor_donorfaces;
        neighbor_donorfaces.resize( ndonors );

        std::vector<std::vector<int>> & sub_local_faceids = interface->sub_local_faceids;
        sub_local_faceids.resize( ndonors );

        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            int ineighbor = donor_zonelocal[ donor_zoneid ];
            std::vector<int> &donorfaces = neighbor_donorfaces[ ineighbor ];
            int global_faceid = interface->global_faceids[ iFace ];
            donorfaces.push_back( global_faceid );
            int local_faceid = interface->local_faceids[ iFace ];

            std::vector<int> & sub_local_faces = sub_local_faceids[ ineighbor ];
            sub_local_faces.push_back( local_faceid );
        }
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> &neighbor_donor_zones = interface->neighbor_donor_zones;
        int ndonor_zones = neighbor_donor_zones.size();

        std::vector<std::vector<int>> & neighbor_donorfaces = interface->neighbor_donorfaces;

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            //std::cout << "donorzone = " << donorzone << " donorfaces = ";
            std::vector<int> & neighbor_donorface = neighbor_donorfaces[ iNei ];

            Global::SendGeom( iZone, donorzone, neighbor_donorface );
        }
        int kkk = 1;
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];
        std::vector<std::vector<int>> & neighbor_donorfaces = interface->neighbor_donorfaces;

        std::vector<int> &neighbor_donor_zones = interface->neighbor_donor_zones;
        int ndonor_zones = neighbor_donor_zones.size();
        std::cout << "zone = " << iZone << " donorzone = ";
        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            std::cout << donorzone;
            if( iNei != ndonor_zones - 1 )
            {
                std::cout << ",";
            }
        }
        std::cout << "\n";

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            std::cout << "donorzone = " << donorzone << " donorfaces = ";
            std::vector<int> & neighbor_donorface = neighbor_donorfaces[ iNei ];
            int ndonorfaces = neighbor_donorface.size();
            for ( int idonorface = 0; idonorface < ndonorfaces; ++ idonorface )
            {
                std::cout <<  neighbor_donorface[ idonorface ];
                if( idonorface != ndonorfaces - 1 )
                {
                    std::cout << ",";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "zone = " << iZone << "\n";

        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> & send_to_zones = interface->send_to_zones;
        int nsend_zones = send_to_zones.size();

        std::vector<std::vector<int>> & donorfaces_for_send = interface->donorfaces_for_send;

        for ( int iSend = 0; iSend < nsend_zones; ++ iSend )
        {
            int zone_to_send = send_to_zones[ iSend ];
            std::cout << "zone_to_send = " << zone_to_send << " donorfaces_for_send = ";
            std::vector<int> & donorfaces = donorfaces_for_send[ iSend ];
            int nFace = donorfaces.size();
            for ( int iFace = 0; iFace < nFace; ++ iFace )
            {
                std::cout <<  donorfaces[ iFace ];
                if( iFace != nFace - 1 )
                {
                    std::cout << ",";
                }
            }
            std::cout << "\n";
        }
        int kkk = 1;
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "zone = " << iZone << "\n";

        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> & send_to_zones = interface->send_to_zones;
        int nsend_zones = send_to_zones.size();

        std::vector<std::vector<int>> & donorfaces_for_send = interface->donorfaces_for_send;

        for ( int iSend = 0; iSend < nsend_zones; ++ iSend )
        {
            int zone_to_send = send_to_zones[ iSend ];
            std::vector<int> & donorfaces = donorfaces_for_send[ iSend ];
            Global::ReSendGeomTest( iZone, zone_to_send, donorfaces );
        }
        int kkk = 1;
    }

    std::cout << "\n";
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "zone = " << iZone << "\n";

        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> & neighbor_donor_zones = interface->neighbor_donor_zones;
        int ndonor_zones = neighbor_donor_zones.size();

        std::vector<std::vector<int>> & donorfaces_for_send = interface->donorfaces_for_send;

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            int tmp_recv_donor_zone = interface->tmp_recv_donor_zones[ iNei ];

            std::cout << "donorzone = " << donorzone << " tmp_recv_donor_zone = " << tmp_recv_donor_zone << "\n";
            std::vector<int> & neighbor_donorfaces = interface->neighbor_donorfaces[ iNei ];
            std::cout << "neighbor_donorfaces\n";
            for ( int i = 0; i < neighbor_donorfaces.size(); ++ i )
            {
                std::cout << neighbor_donorfaces[ i ] << " ";
            }
            std::cout << "\n";
            std::vector<int> & tmp_recv_donorfaces = interface->tmp_recv_donorfaces[ iNei ];
            std::cout << "tmp_recv_donorfaces\n";
            for ( int i = 0; i < tmp_recv_donorfaces.size(); ++ i )
            {
                std::cout << tmp_recv_donorfaces[ i ] << " ";
            }
            std::cout << "\n";
        }
    }

    int kkk = 1;
}

void Solver::InitFields()
{
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
    }

    Boundary();
}

void Solver::SolveMultiZones()
{
    for ( int it = 0; it < Global::nt; ++ it )
    {
        int nZones = Global::zones.size();
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

void Solver::ExchangeInterfaceFieldOldVersion()
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

void Solver::ExchangeInterfaceFieldOldVersion1()
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
            int i_donor_cell = interface->ijk_donors[ ijkpos + 0 ];

            Field * donor_field = Global::fields[ donor_zoneid ];
            double donor_value = donor_field->u[ i_donor_cell ];
            interface->data_recv[ iFace ] = donor_value;
        }

        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            int ijkpos = index_dim * iFace;
            int i_ghost_cell = interface->ijk_ghosts[ ijkpos + 0 ];

            double donor_value = interface->data_recv[ iFace ];
            field->u[ i_ghost_cell ] = donor_value;
        }
    }
}

void Solver::ExchangeInterfaceField()
{
    int nZones = Global::interfaces.size();

    //zone0:
    //donor zone2 donorfaces(0)
    //zone0->send_geom(zone=2,facelist(0)) : 
    //  zone2->recv_geom(zone=0,facelist(0));
    //  zone2->get_value(facelist(0),valuelist(200));
    //  zone2->send_value(zone=0,valuelist(200));
    //zone0->recv_value(zone=2,valuelist(200));
    
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> & neighbor_donor_zones = interface->neighbor_donor_zones;
        int ndonor_zones = neighbor_donor_zones.size();

        std::cout << "zone = " << iZone << "\n";

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            std::cout << "donorzone = " << donorzone << "\n";
            std::vector<int> & neighbor_donorfaces = interface->neighbor_donorfaces[ iNei ];
            std::vector<int> & sub_local_faceids = interface->sub_local_faceids[ iNei ];
            for ( int i = 0; i < neighbor_donorfaces.size(); ++ i )
            {
                int global_faceid = neighbor_donorfaces[ i ];
                int local_faceid = interface->global_local_face_map[ global_faceid ];
                int mylocal = sub_local_faceids[ i ];
                std::cout << " global_faceid = " << global_faceid << " local_faceid = " << local_faceid 
                    << " mylocal = " << mylocal << "\n";
            }
        }
    }
    int kkk = 1;

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
            int i_donor_cell = interface->ijk_donors[ ijkpos + 0 ];

            Field * donor_field = Global::fields[ donor_zoneid ];
            double donor_value = donor_field->u[ i_donor_cell ];
            interface->data_recv[ iFace ] = donor_value;
        }

        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int donor_zoneid = interface->zoneList[ iFace ];
            int ijkpos = index_dim * iFace;
            int i_ghost_cell = interface->ijk_ghosts[ ijkpos + 0 ];

            double donor_value = interface->data_recv[ iFace ];
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

    std::string csvname = "field_final.csv";
    this->DumpCsvFile( csvname, x, u_e, un, u_error );
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
        std::format_to(std::ostream_iterator<char>(file), "{:.20f} {:.20f} {:.20f} {:.20f}\n", x[i], ue[i], un[i], uerror[i] );
    }
    file.close();
}
