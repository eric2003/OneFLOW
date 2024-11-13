#include "Solver.h"
#include "Parallel.h"
#include "ZoneState.h"
#include "global.h"
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

Solver::Solver()
{
    Parallel::Init();
}

Solver::~Solver()
{
    Parallel::Finalize();
}


void Solver::Run()
{
    this->ReadGrid();
    this->InitTopo();
    this->InitFields();
    ////if ( Parallel::pid == Parallel::serverid )
    ////{
        this->SolveMultiZones();
        this->PostProcess();
    ////}
}

void Solver::ReadGrid()
{
    std::string fileName = "../heat1d4blocksv1.cgns";
    ReadCgnsGridBaseZone( fileName );
    ReadCgnsGrid( fileName );
}

void Solver::InitTopo()
{
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "Solver::InitTopo() " << "\n";

    Global::donor_zone_sets.resize( LocalZone::nZones );
    Global::donor_zones.resize( LocalZone::nZones );
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        int global_zoneid = LocalZone::global_zoneids[ iZone ];
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "iZone = " << iZone << " global_zoneid = " << global_zoneid << "\n";

        Interface * interface = new Interface();
        interface->zoneid = iZone;
        Global::interfaces.push_back( interface );
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Zone * zone = Global::zones[ iZone ];

        Interface * interface = Global::interfaces[ iZone ];

        int nbc1to1s = zone->bc1to1s.size();
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout <<  "nbc1to1s = " << nbc1to1s << "\n";

        for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
        {
            ZoneBc1To1 * bc1to1 = zone->bc1to1s[ ibc1to1 ];
            int zoneid = bc1to1->zoneid;
            int donor_zoneid = bc1to1->donor_zoneid;
            std::cout << "Parallel::pid = " << Parallel::pid << " ";
            std::cout << "zoneid = " << zoneid << " donor_zoneid = " << donor_zoneid << "\n";

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

    for ( int iProc = 0; iProc < Parallel::nProc; ++ iProc )
    {
        int nSize = -1;
        if ( iProc == Parallel::pid )
        {
            nSize = Global::facePairList.size();
        }
        HXBcastData( &nSize, 1, iProc );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "nSize = " << nSize << "\n";
        std::vector<FacePair> tmp;
        if ( iProc == Parallel::pid )
        {
            tmp = Global::facePairList;
        }
        else
        {
            tmp.resize( nSize );
        }

        HXBcastData( tmp.data(), tmp.size(), iProc );
        Global::AddFacePairList( Global::mpi_facePairList, tmp );
    }

    for ( int i = 0; i < Global::mpi_facePairList.size(); ++ i )
    {
        FacePair &facePair = Global::mpi_facePairList[ i ];
        Global::InsertFacePairMap( facePair );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        facePair.Print();
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Zone * zone = Global::zones[ iZone ];

        Interface * interface = Global::interfaces[ iZone ];
        int nInterfaces = interface->local_faceids.size();
        for ( int iInterface = 0; iInterface < nInterfaces; ++ iInterface )
        {
            int local_faceid = interface->local_faceids[ iInterface ];
            int proc_global_faceid = interface->proc_global_faceids[ iInterface ];
            FacePair & facePair = Global::facePairList[ proc_global_faceid ];
            int global_faceid = Global::InsertFacePairMap( facePair );
            interface->global_faceids.push_back( global_faceid );
            interface->global_local_face_map.insert( std::make_pair( global_faceid, local_faceid ) );
        }
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
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

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> &neighbor_donor_zones = interface->neighbor_donor_zones;
        int ndonor_zones = neighbor_donor_zones.size();

        std::vector<std::vector<int>> & neighbor_donorfaces = interface->neighbor_donorfaces;

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            std::vector<int> & neighbor_donorface = neighbor_donorfaces[ iNei ];

            Interface * interface = Global::interfaces[ donorzone ];
            interface->SendGeom( iZone, neighbor_donorface );

            //Global::SendGeom( iZone, donorzone, neighbor_donorface );
        }
        int kkk = 1;
    }
    return;
}

void Solver::InitFields()
{
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "Solver::InitFields() ZoneState::nZones = " << ZoneState::nZones <<  "\n";
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValid( iZone ) ) continue;
        Field * field = new Field();
        Global::fields.push_back( field );
    }

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValid( iZone ) ) continue;
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
        for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
        {
            if ( ! ZoneState::IsValid( iZone ) ) continue;
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
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValid( iZone ) ) continue;
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->PhysicalBoundary( zone );
    }
    ExchangeInterfaceField();
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValid( iZone ) ) continue;
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->UpdateOldField();
    }
}

void Solver::ExchangeInterfaceField()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];

        std::vector<int> & send_to_zones = interface->send_to_zones;
        int nsend_zones = send_to_zones.size();

        std::vector<std::vector<int>> & donorfaces_for_send = interface->donorfaces_for_send;
        std::vector<std::vector<int>> & donorijk_for_send = interface->donorijk_for_send;
        std::vector<std::vector<double>> & donordata_for_send = interface->donordata_for_send;

        Field * field = Global::fields[ iZone ];

        for ( int iSend = 0; iSend < nsend_zones; ++ iSend )
        {
            int zone_to_send = send_to_zones[ iSend ];
            std::vector<int> & donorfaces = donorfaces_for_send[ iSend ];
            std::vector<int> & donorijks = donorijk_for_send[ iSend ];
            std::vector<double> & donordatas = donordata_for_send[ iSend ];

            int nface = donorfaces.size();
            int index_dim = 1;
            for ( int i = 0; i < nface; ++ i )
            {
                int ijkpos = index_dim * i;
                int i_donor_cell = donorijks[ ijkpos + 0 ];
                donordatas[ i ] = field->u[ i_donor_cell ];
            }
        }
    }
    int kkk = 1;

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];

        Interface * interface = Global::interfaces[ iZone ];
        int nInterFaces = interface->zoneList.size();

        std::vector<int> & neighbor_donor_zones = interface->neighbor_donor_zones;
        int ndonor_zones = neighbor_donor_zones.size();

        //std::vector<std::vector<double>> & recv_donor_fields = interface->recv_donor_fields;

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donorzone = neighbor_donor_zones[ iNei ];
            Interface * donor_interface = Global::interfaces[ donorzone ];
            int nSize = donor_interface->send_to_zones.size();
            int ipos = -1;
            for ( int i = 0; i < nSize; ++ i )
            {
                int szone = donor_interface->send_to_zones[ i ];
                if ( szone == iZone )
                {
                    ipos = i;
                    break;
                }
            }
            std::vector<double> & donordata = donor_interface->donordata_for_send[ ipos ];

            std::vector<int> & neighbor_donorfaces = interface->neighbor_donorfaces[ iNei ];
            std::vector<int> & sub_local_faceids = interface->sub_local_faceids[ iNei ];
            for ( int i = 0; i < neighbor_donorfaces.size(); ++ i )
            {
                int local_faceid = sub_local_faceids[ i ];
                double donor_value = donordata[ i ];
                interface->data_recv[ local_faceid ] = donor_value;
            }
        }

        int index_dim = 1;
        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
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
