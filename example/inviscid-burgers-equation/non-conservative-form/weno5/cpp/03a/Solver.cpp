#include "Solver.h"
#include "CgnsUtil.h"
#include "Parallel.h"
#include "Post.h"
#include "Weno.h"
#include "ZoneState.h"
#include "global.h"
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <print> 
#include <set>
#include <unordered_map>

Solver::Solver()
{
    Parallel::Init();
    //this->nghost = 1;
    //this->scheme = Scheme::FTCS;
    //this->scheme = Scheme::RungeKutta;
    //this->scheme = Scheme::CN;
    //this->scheme = Scheme::ICP;
    this->scheme = Scheme::WENO;
    this->nghost = 3;
    Global::nghost = this->nghost;
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
    this->SolveFields();
    //this->PostProcess();
}

void Solver::ReadGrid()
{
    //std::string fileName = "../burgers1d1blocksv1.cgns";
    std::string fileName = "../burgers1d2blocks.cgns";
    ReadCgnsGridBaseZone( fileName );
    ReadCgnsGrid( fileName );
}

void Solver::InitFields()
{
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "Solver::InitFields() ZoneState::nZones = " << ZoneState::nZones << "\n";
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = nullptr;
        if ( scheme == Scheme::WENO ) {
            field = new WenoField();
        }
        else {
            field = new FieldSub();
        }
        Global::fields.push_back( field );
    }

    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Field * field = Global::fields[ iZone ];
        field->Init( grid );
    }

    this->Boundary();
    this->UpdateOldField();
}

void Solver::DumpInitialFields()
{
    this->PostProcess();
}

void Solver::SolveFields()
{
    this->DumpInitialFields();
    for ( int it = 0; it < Global::nt; ++ it )
    {
        Global::iter = it;
        switch (scheme) {  
        case Scheme::FTCS:
            this->FTCS();
            break;
        case Scheme::CN:
            this->CN();
            break;
        case Scheme::ICP:
            this->ICP();
            break;
        case Scheme::WENO:
            this->RungeKutta();
            break;
        case Scheme::RungeKutta:  
            this->RungeKutta();
            break;  
        default:
            this->FTCS();
        }
        if ( ( Global::iter + 1 ) % 250 == 0 )
        {
            std::print( "it = {} nt = {}\n", Global::iter + 1, Global::nt );
            this->PostProcess();
        }

    }
}

void Solver::FTCS()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->FTCS( zone );
    }
    this->Boundary();
    this->UpdateOldField();
}

void Solver::CN()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->CN( zone );
    }
    this->Boundary();
    this->UpdateOldField();
}

void Solver::ICP()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->ICP( zone );
    }
    this->Boundary();
    this->UpdateOldField();
}

void Solver::RungeKutta()
{
    for ( int istage = 0; istage < 3; ++ istage )
    {
        this->RungeKutta( istage );
    }
    this->UpdateOldField();
}

void Solver::RungeKutta( int istage )
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        zone->zoneIndex = iZone;
        field->RungeKutta( zone, istage );
    }
    this->Boundary();
}

void Solver::Boundary()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->PhysicalBoundary( zone );
    }
    ExchangeInterfaceField();
}

void Solver::UpdateOldField()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Field * field = Global::fields[ iZone ];
        Zone * zone = Global::zones[ iZone ];
        field->UpdateOldField();
    }
}

void Solver::UploadInterfaceField()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];

        Field * field = Global::fields[ iZone ];

        int nsend_zones = interface->send_to_zones.size();
        for ( int iSend = 0; iSend < nsend_zones; ++ iSend )
        {
            int zone_to_send = interface->send_to_zones[ iSend ];
            std::vector<int> & donorfaces_for_send = interface->donorfaces_for_send[ iSend ];
            std::vector<int> & donorijk_for_send = interface->donorijk_for_send[ iSend ];
            std::vector<double> & donordata_for_send = interface->donordata_for_send[ iSend ];

            int nInterFaces = donorfaces_for_send.size();
            int index_dim = 1;
            for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
            {
                int ijkpos = index_dim * iFace * nghost;
                int data_pos = iFace * nghost;

                for ( int ig = 0; ig < nghost; ++ ig )
                {
                    int id_cell = donorijk_for_send[ ijkpos + ig ] - 1;
                    donordata_for_send[ data_pos + ig ] = field->u[ id_cell ];
                }
            }
        }
    }
}

void Solver::UpdateInterfaceField()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int ndonor_zones = Global::interfaceTopo.linkmap[ iZone ].size();
        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donor_zone = Global::interfaceTopo.linkmap[ iZone ][ iNei ];
            int send_pid = ZoneState::pids[ iZone ];
            int recv_pid = ZoneState::pids[ donor_zone ];
            int nsend = -1;
            std::vector<double> donordata;
            if ( Parallel::pid != send_pid && Parallel::pid != recv_pid ) continue;
            if ( Parallel::pid == send_pid )
            {
                int local_zoneid = ZoneState::g2lzoneids[ iZone ];
                Interface * interface = Global::interfaces[ local_zoneid ];
                donordata = interface->donordata_for_send[ iNei ];
                nsend = donordata.size();
            }
            HXSendRecvData( &nsend, 1, send_pid, recv_pid );

            if ( Parallel::pid == recv_pid )
            {
                donordata.resize( nsend );
            }
            HXSendRecvData( donordata.data(), donordata.size(), send_pid, recv_pid );

            if ( Parallel::pid == recv_pid )
            {
                int local_donor_zoneid = ZoneState::g2lzoneids[ donor_zone ];
                Interface * donor_interface = Global::interfaces[ local_donor_zoneid ];
                int nSize = donor_interface->neighbor_donor_zones.size();
                int ipos = -1;
                for ( int i = 0; i < nSize; ++ i )
                {
                    int nei_zone = donor_interface->neighbor_donor_zones[ i ];
                    if ( nei_zone == iZone )
                    {
                        ipos = i;
                        break;
                    }
                }

                std::vector<int> & neighbor_donorfaces = donor_interface->neighbor_donorfaces[ ipos ];
                std::vector<int> & sub_local_faceids = donor_interface->sub_local_faceids[ ipos ];
                for ( int i = 0; i < neighbor_donorfaces.size(); ++ i )
                {
                    int local_faceid = sub_local_faceids[ i ];

                    int index_dim = 1;

                    int ijkpos = index_dim * local_faceid * nghost;
                    int data_pos = local_faceid * nghost;
                    int donor_data_pos = i * nghost;

                    for ( int ig = 0; ig < nghost; ++ ig )
                    {
                        double donor_value = donordata[ donor_data_pos + ig ];
                        donor_interface->data_recv[ data_pos + ig ] = donor_value;
                    }
                }
            }
        }
    }
}

void Solver::DownloadInterfaceField()
{
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Interface * interface = Global::interfaces[ iZone ];
        Field * field = Global::fields[ iZone ];

        int nInterFaces = interface->zoneList.size();

        int index_dim = 1;
        for ( int iFace = 0; iFace < nInterFaces; ++ iFace )
        {
            int ijkpos = index_dim * iFace * nghost;
            int data_pos = iFace * nghost;
            for ( int ig = 0; ig < nghost; ++ ig )
            {
                int ig_cell = interface->ijk_ghosts[ ijkpos + ig ] - 1;

                double donor_value = interface->data_recv[ data_pos + ig ];
                field->u[ ig_cell ] = donor_value;
            }
        }
    }
}

void Solver::ExchangeInterfaceField()
{
    this->UploadInterfaceField();
    this->UpdateInterfaceField();
    this->DownloadInterfaceField();
}

void Solver::PostProcess()
{
    Global::file_string = {};
    std::string filename = std::format( "field_final{}.csv", Global::iter+1 );
    std::fstream file;
    file.open( filename.c_str(), std::fstream::out );
    for ( int iZone = 0; iZone < LocalZone::nZones; ++ iZone )
    {
        Grid * grid = Global::grids[ iZone ];
        Field * field = Global::fields[ iZone ];
        field->PostProcess( grid );
    }
    std::format_to(std::ostream_iterator<char>(file), "{}", Global::file_string );
    file.close();
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
        interface->zoneid = global_zoneid;
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

            interface->CalcInterface( &transform, region.start, region.end, donor_zoneid, this->nghost );
        }
        int nInterfaces = interface->zoneList.size();
        int nData = nInterfaces * Global::nghost;
        interface->data_recv.resize( nData );
        interface->data_send.resize( nData );
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

    Global::interfaceTopo.InitNeighborInfo();
    Global::interfaceTopo.SwapNeighborInfo();

    int kkk = 1;
}