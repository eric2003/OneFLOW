#include "Partition.h"
#include "CgnsGrid.h"
#include "SplitZone.h"
#include "cgnslib.h"
#include <string>

bool CmpSplitZone::operator() ( SplitZone * a, SplitZone * b ) const
{
    return a->zoneIndex < b->zoneIndex;
};


double CalSumOfProcSpeed( std::vector<double> & procSpeed );
double CalSumOfProcSpeed( std::vector<double> & procSpeed )
{
    int nProc = procSpeed.size();
    double sumOfProcSpeed = 0;
    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        sumOfProcSpeed += procSpeed[ iProc ];
    }
    return sumOfProcSpeed;
}

Partition::Partition()
{
    this->inputName = "../heat2d1blockv2.cgns";
    this->outName =  "../heat2dpart.cgns";
}

Partition::~Partition()
{
    ;
}

int Partition::GetMinValueIndex( double * data, int numberOfData )
{
    double minValue = data[ 0 ];
    int minimumValueIndex = 0;

    for ( int iData = 1; iData < numberOfData; ++ iData )
    {
        if ( minValue > data[ iData ] )
        {
            minValue      = data[ iData ];
            minimumValueIndex = iData;
        }
    }
    return minimumValueIndex;
}

int Partition::GetMinLoadProc()
{
    return this->GetMinValueIndex( this->timeProc.data(), this->nProc );
}

SplitZone * Partition::GetLargestZone()
{
    SplitZone * largestZone = 0;

    double nMaxCell = 0.0;

    for ( std::set< SplitZone * >::iterator iter = unsignedZoneGroup.begin(); iter != unsignedZoneGroup.end(); ++ iter )
    {
        double nSize = ( * iter )->GetNCell();
        if ( nMaxCell < nSize )
        {
            nMaxCell    = nSize;
            largestZone = * iter;
        }
    }

    return largestZone;
}

void Partition::Solve()
{
    std::cout << "Partition::Split()\n";
    //this->nProc = 2;
    //this->nProc = 3;
    this->nProc = 16;
    this->ReadGrid();
    this->Split();
    this->DumpPartitionedGridFile();
    this->FreeAllZones();
}

void Partition::ReadGrid()
{
    ReadCgnsGridBaseZone( this->inputName );
    ReadCgnsGrid( this->inputName );

    this->nOriZone = ::Global::zones.size();
    Dim::dim = Global::cell_dim;

    std::cout << "numberOfBlocks = " << nOriZone << "\n";

    this->refZone.resize( nOriZone );

    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        SplitZone * splitZone = new SplitZone();
        refZone[ iZone ] = splitZone;

        splitZone->SetZoneIndex( iZone );
        splitZone->SetProcId( 0 );
        splitZone->SetParentAndChild( 0 );
    }

    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        SplitZone * splitZone = refZone[ iZone ];

        Zone * zone = ::Global::zones[ iZone ];

        splitZone->dimInfo.dimList = zone->nijk;

        int nbccos = zone->bccos.size();
        std::cout << "nbccos = " << nbccos << "\n";

        for ( int ibcco = 0; ibcco < nbccos; ++ ibcco )
        {
            ZoneBc * zoneBc = zone->bccos[ ibcco ];

            PhysicalSplitBc * physicalSplitBc = new PhysicalSplitBc();
            splitZone->physicalSplitBcList.push_back( physicalSplitBc );
            physicalSplitBc->bcType = zoneBc->bcType;
            physicalSplitBc->SetRegion( zoneBc->pnts );
            physicalSplitBc->splitZone = splitZone;
        }

        int nbc1to1s = zone->bc1to1s.size();

        for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
        {
            ZoneBc1To1 * zoneBc1To1 = zone->bc1to1s[ ibc1to1 ];

            InterfaceSplitBc * interfaceSplitBc = new InterfaceSplitBc();
            splitZone->interfaceSplitBcList.push_back( interfaceSplitBc );
            interfaceSplitBc->zone = splitZone;
            interfaceSplitBc->region.SetRegion( zoneBc1To1->pnts );

            interfaceSplitBc->transform = zoneBc1To1->transform;
            interfaceSplitBc->CalcTransformMatrix();

            interfaceSplitBc->donor_zone = refZone[ zoneBc1To1->donor_zoneid ];
            interfaceSplitBc->donor_region.SetRegion( zoneBc1To1->donor_pnts );
        }
    }

    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        SplitZone * splitZone = refZone[ iZone ];

        int nbc1to1s = splitZone->interfaceSplitBcList.size();
        for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
        {
            InterfaceSplitBc * interfaceSplitBc = splitZone->interfaceSplitBcList[ ibc1to1 ];
            InterfaceInfo interfaceInfo;
            interfaceSplitBc->donor_zone->FindDonorInterface( interfaceSplitBc, interfaceInfo );
            interfaceSplitBc->donor_zone = interfaceInfo.donor_zone;
            interfaceSplitBc->donor_bc = interfaceInfo.donor_interface;
        }
    }
}

void Partition::Split()
{
    this->InitSplit();
    this->BinarySplit();
    this->ModifyZoneIndex();
    int kkk = 1;
}

void Partition::InitSplit()
{
    //nZone is the number of zones in the original grid
    //nProc is not the number of blocks in the end, but the number of processes that need to be divided.
    //This difference occurs because the structure grid can have many blocks per process.

    //this->nProc    = nProc;
    //this->nPhyBlocks = referenceBlock.size();
    //this->nZone  = this->nPhyBlocks;

    //tj Computation time on processor j
    //Nk Number of internal grid-cells for zone j
    //aveFloatOp:  Number of Real operations per gridcell
    //Gj The set of blocks allocated to process j
    //Pj Real operations per second for process j

    //The Greedy Load Balancing Algorithm
    //1.Compute an "ideal" uniprocessor time, Tuni = sum ( Nk * aveFloatOp )/ sum( Pj ).
    //2.Set the group of unassigned blocks, Gu, to contain all blocks. The number of blocks in Nu=M, where M is the total number of blocks in the mesh.
    //3.Start all processors with a very small fictitious zone to get a start tj=delt aveFloatOp/Pj, where delt<<1 is the size of the very small fictious zone.
    //This is used to find the fastest processor in a heterogeneous system, tj could otherwise be set to 0 in the homogeneous case.
    this->aveFloatOp = 1.0;
    this->delta      = 1.0e-5;
    this->procSpeed.resize( nProc );
    this->timeProc.resize( nProc );

    this->nTZone = this->nOriZone;

    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        procSpeed[ iProc ] = 1.0;
    }

    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        timeProc[ iProc ] = delta * aveFloatOp / procSpeed[ iProc ];
    }

    double sumProcSpeed = CalSumOfProcSpeed( procSpeed );

    double sumProcOp = 0.0;
    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        sumProcOp += refZone[ iZone ]->GetNCell() * aveFloatOp;
    }

    this->idealUniProcTime = sumProcOp / sumProcSpeed;

    this->eps = 0.05;

    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        unsignedZoneGroup.insert( refZone[ iZone ] );
    }

    //Init zoneStorage
    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        zoneStorage.insert( refZone[ iZone ] );
    }

    //Assign procSplitZone to each process
    procSplitZone.resize( nProc );

    this->minNumberCell = 1;
    for ( int m = 0; m < Dim::dim; ++ m )
    {
        this->minNumberCell *= 1;
    }

    int kkk = 1;

}

void Partition::FreeAllZones()
{
    int iCount = 0;
    for ( std::set < SplitZone * >::iterator iter = zoneStorage.begin(); iter != zoneStorage.end(); ++ iter )
    {
        //cout << "iCount = " << iCount << "\n";
        ++ iCount;
        delete * iter;
    }
}

int Partition::GetNZones( int iProc )
{
    return this->procSplitZone[ iProc ].size();
}

SplitZone * Partition::GetSplitZone( int iProc, int zoneId )
{
    return this->procSplitZone[ iProc ][ zoneId ];
}

void Partition::AddZoneToUnsignedGroup( SplitZone * zone )
{
    this->unsignedZoneGroup.insert( zone );
}

void Partition::RemoveZoneFromUnsignedGroup( SplitZone * zone )
{
    this->unsignedZoneGroup.erase( zone );
}

void Partition::AddZoneToProcess( int iProc, SplitZone * zone )
{
    this->procSplitZone[ iProc ].push_back( zone );
}

void Partition::AddNewZone( SplitZone * zone )
{
    zone->SetZoneIndex( nTZone ++ );

    this->AddZoneToUnsignedGroup( zone );

    this->zoneStorage.insert( zone );
}

SplitZone * Partition::Split( SplitZone * zone, double nCell )
{
    SplitZone * zoneL = new SplitZone();
    zoneL->SetParentAndChild( zone );

    SplitZone * zoneR = new SplitZone();
    zoneR->SetParentAndChild( zone );

    //Remove the original block from the UnsignedGroup
    this->RemoveZoneFromUnsignedGroup( zone );

    this->AddNewZone( zoneL );
    this->AddNewZone( zoneR );

    zone->Split( zoneL, zoneR, nCell );

    return zoneL;
}

void Partition::BinarySplit()
{
    while ( ! unsignedZoneGroup.empty() )
    {
        std::cout << "number of unsigned zone = " << unsignedZoneGroup.size() << "\n";
        int iProc = this->GetMinLoadProc();

        //Get the block with the largest number of grids and its size
        SplitZone * largestZone = this->GetLargestZone();

        double nMaxCell = largestZone->GetNCell();

        double aveFloatOp = this->aveFloatOp;
        double tau = timeProc[ iProc ] + nMaxCell * aveFloatOp / procSpeed[ iProc ];

        //It seems that EPS is still necessary
        //The idea here is as follows. For the candidate blocks which can be divided into blocks, 
        //if the difference between the target blocks and the candidate blocks is too wide, 
        //the binarySplit method is adopted. The binarySplit blocks do not directly join the process,
        //but enter the next round of screening.

        if ( tau > idealUniProcTime * ( 1.0 + eps ) )
        {
            double nRest = ( tau - idealUniProcTime ) * procSpeed[ iProc ] / aveFloatOp;

            double nTarget = nMaxCell - nRest;
            if ( nRest > this->minNumberCell && nTarget > this->minNumberCell )
            {
                //It is necessary to divide blocks only when a certain number of points are satisfied.
                double nHalf = nMaxCell / 2.0;
                nMaxCell = nTarget;
                double ratio = nMaxCell / nRest;
                double oRatio = 1.0 / ratio;

                if ( std::max( ratio, oRatio ) < 3.0 )
                {
                    int kkk = 1;
                    largestZone = this->Split( largestZone, nMaxCell );
                    timeProc[ iProc ] += nMaxCell * aveFloatOp / procSpeed[ iProc ];
                    this->RemoveZoneFromUnsignedGroup( largestZone );
                    this->AddZoneToProcess( iProc, largestZone );
                }
                else
                {
                    this->Split( largestZone, nHalf );
                }
            }
            else
            {
                //If the block is of the right size, join the process directly, and there is no need to split it.
                timeProc[ iProc ] += nMaxCell * aveFloatOp / procSpeed[ iProc ];
                this->RemoveZoneFromUnsignedGroup( largestZone );
                this->AddZoneToProcess( iProc, largestZone );
            }
        }
        else
        {
            //If the block is of the right size, join the process directly, and there is no need to split it.
            timeProc[ iProc ] += nMaxCell * aveFloatOp / procSpeed[ iProc ];
            this->RemoveZoneFromUnsignedGroup( largestZone );
            this->AddZoneToProcess( iProc, largestZone );
        }
    }
}

void Partition::ComputeCoor( SplitZone * splitZone, int icoord, std::vector<char> & coord )
{
    SplitZone * rootSplitZone = splitZone->GetRootZone();
    Zone * zone = Global::zones[ rootSplitZone->zoneIndex ];

    int nNodes = splitZone->dimInfo.ComputeNumNodes();
    std::cout << "nNodes = " << nNodes << "\n";

    coord.resize( nNodes * sizeof(double) );

    int index_dim = splitZone->dimInfo.dimList.size();
    Coor * coor = zone->coors[ icoord ];

    double * xlocal = reinterpret_cast<double *>( const_cast<char *>( coord.data() ) );
    double * xglobal = reinterpret_cast<double *>( const_cast<char *>( coor->coord.data() ) );
    if ( index_dim == 1 )
    {
        int icount = 0;
        int imin = splitZone->dimInfo.rmin[ 0 ] - 1;
        int imax = splitZone->dimInfo.rmax[ 0 ];
        for ( int i = imin; i < imax; ++ i )
        {
            xlocal[ icount++ ] = xglobal[ i ];
        }
    }
    else if ( index_dim == 2 )
    {
        int icount = 0;
        int imin = splitZone->dimInfo.rmin[ 0 ] - 1;
        int imax = splitZone->dimInfo.rmax[ 0 ];
        int jmin = splitZone->dimInfo.rmin[ 1 ] - 1;
        int jmax = splitZone->dimInfo.rmax[ 1 ];
        int ni = rootSplitZone->dimInfo.dimList[ 0 ];
        int nj = rootSplitZone->dimInfo.dimList[ 1 ];
        for ( int j = jmin; j < jmax; ++ j )
        {
            int jstride = j * ni;
            for ( int i = imin; i < imax; ++ i )
            {
                xlocal[ icount++ ] = xglobal[ jstride + i ];
            }
        }
    }
}

void Partition::DumpBc( SplitZone * splitZone )
{
}

void Partition::DumpPartitionedGridFile()
{
    std::cout << "\n";
    std::cout << "DumpPartitionedGridFile -------------------------------" << "\n";
    std::cout << "\n";
    int fileId = -1;
    std::string filename = this->outName;
    int ierr = cg_open( filename.c_str(), CG_MODE_WRITE, &fileId );
    std::cout << "fileId = " << fileId << "\n";

    std::vector<std::string> coordnames;
    coordnames.push_back( "CoordinateX" );
    coordnames.push_back( "CoordinateY" );
    coordnames.push_back( "CoordinateZ" );

    std::string basename = "Base";
    int baseId = -1;
    int icelldim = Global::cell_dim;
    int iphysdim = Global::phys_dim;
    cg_base_write( fileId, basename.c_str(), icelldim, iphysdim, &baseId );

    std::cout << "baseId = " << baseId << "\n";

    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        int nZones = this->GetNZones( iProc );

        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            ZoneType_t zoneType = Structured;

            SplitZone * splitZone = this->GetSplitZone( iProc, iZone );

            int index_dim = splitZone->dimInfo.dimList.size();
            std::cout << "index_dim = " << index_dim << "\n";

            splitZone->dimInfo.ComputeISize();

            int isize_dim = index_dim * 3;

            std::vector<cgsize_t> isize( index_dim * 3, 0 );
            for ( int m = 0; m < isize_dim; ++ m )
            {
                isize[ m ] = splitZone->dimInfo.isize[ m ];
            }

            int zoneId = -1;

            std::cout << "splitZone->newIndex = " << splitZone->newIndex << "\n";
            std::cout << "splitZone->zoneIndex = " << splitZone->zoneIndex << "\n";

            std::string zonename = "Zone" + std::to_string( splitZone->newIndex );

            cg_zone_write( fileId, baseId, zonename.c_str(), isize.data(), zoneType, &zoneId );

            std::cout << "zoneId = " << zoneId << "\n";
            std::cout << "zonename = " << zonename << "\n";
                
            splitZone->ComputeRminmax();

            SplitZone * rootSplitZone = splitZone->GetRootZone();

            Zone * zone = Global::zones[ rootSplitZone->zoneIndex ];

            int ncoords = zone->coors.size();
            std::cout << "ncoords = " << ncoords << "\n";

            for ( int icoord = 0; icoord < ncoords; ++ icoord )
            {
                int coorId = -1;
                DataType_t dataType = RealDouble;
                std::string coordname = coordnames[ icoord ];
                Coor * coor = zone->coors[ icoord ];

                int nNodes = splitZone->dimInfo.ComputeNumNodes();
                std::cout << "nNodes = " << nNodes << "\n";

                std::vector<char> coord;
                ComputeCoor( splitZone, icoord, coord );

                cg_coord_write( fileId, baseId, zoneId, dataType, coordname.c_str(), coord.data(), &coorId );

                std::cout << "fileId = " << fileId  << " baseId = " << baseId << " zoneId = " << zoneId  << " coorId = " << coorId << "\n";
            }

            std::vector< PhysicalSplitBc * > &bccoList = splitZone->physicalSplitBcList;

            int nbocos = bccoList.size();
            std::cout << "nbocos = " << nbocos << "\n";
            for ( int iboco = 0; iboco < nbocos; ++ iboco )
            {
                int bccoId = -1;
                GridLocation_t location = Vertex;

                std::string boconame = "BC" + std::to_string( iboco + 1 );

                PhysicalSplitBc * physicalSplitBc = bccoList[ iboco ];
                BCType_t bocotype = static_cast<BCType_t>( physicalSplitBc->bcType );
                PointSetType_t ptset_type = PointRange;
                cgsize_t npnts = 2;
                std::vector<cgsize_t> pnts;
                int nSize = physicalSplitBc->region.start.size();
                for ( int i = 0; i < nSize; ++ i )
                {
                    pnts.push_back( physicalSplitBc->region.start[ i ] );
                }
                for ( int i = 0; i < nSize; ++ i )
                {
                    pnts.push_back( physicalSplitBc->region.end[ i ] );
                }

                cg_boco_write( fileId, baseId, zoneId, boconame.c_str(), bocotype, ptset_type, npnts, pnts.data(), & bccoId );
            }

            //std::vector< InterfaceSplitBc * > & bc1to1List = splitZone->interfaceSplitBcList;
            std::vector< InterfaceSplitBc * > bc1to1List;
            splitZone->GetLeafInterfaceBc( bc1to1List );

            int n1to1 = bc1to1List.size();
            std::cout << "n1to1 = " << n1to1 << "\n";

            std::vector<int> itranfrm( index_dim );
            for ( int m = 0; m < index_dim; ++ m )
            {
                itranfrm[ m ] = m + 1;
            }
            for ( int i1to1 = 0; i1to1 < n1to1; ++ i1to1 )
            {
                int index_conn = -1;

                InterfaceSplitBc * interfaceSplitBc = bc1to1List[ i1to1 ];

                std::vector<cgsize_t> pnts;
                std::vector<cgsize_t> pntsdonor;
                int nSize = interfaceSplitBc->region.start.size();
                for ( int i = 0; i < nSize; ++ i )
                {
                    pnts.push_back( interfaceSplitBc->region.start[ i ] );
                }
                for ( int i = 0; i < nSize; ++ i )
                {
                    pnts.push_back( interfaceSplitBc->region.end[ i ] );
                }

                int nPatchSize = interfaceSplitBc->donor_region.start.size();
                for ( int i = 0; i < nPatchSize; ++ i )
                {
                    pntsdonor.push_back( interfaceSplitBc->donor_region.start[ i ] );
                }
                for ( int i = 0; i < nPatchSize; ++ i )
                {
                    pntsdonor.push_back( interfaceSplitBc->donor_region.end[ i ] );
                }
                int id = interfaceSplitBc->donor_zone->newIndex;
                int oldZoneIndex = interfaceSplitBc->donor_zone->zoneIndex;
                std::cout << "splitZone->newIndex = " << splitZone->newIndex << "\n";
                std::cout << "splitZone->zoneIndex = " << splitZone->zoneIndex << "\n";
                std::cout << "donor zone id = " << id << "\n";
                std::cout << "donor zone oldZoneIndex = " << oldZoneIndex << "\n";
                std::string donorname = "Zone" + std::to_string( id );
                std::string connectname = "Interface" + std::to_string( i1to1 );

                std::cout << "donorname = " << donorname << "\n";
                std::cout << "connectname = " << connectname << "\n";

                std::vector<int> & itranfrm = interfaceSplitBc->transform;

                cg_1to1_write( fileId, baseId, zoneId, connectname.c_str(), donorname.c_str(), pnts.data(), pntsdonor.data(), itranfrm.data(), &index_conn );
                std::cout << "index_conn = " << index_conn << "\n";
            }
        }
    }
    cg_close( fileId );
}

void Partition::ModifyZoneIndex()
{
    std::cout << "ModifyZoneIndex +++++++++++++++++++++++++ " << "\n";
    int globalZoneId = 0;
    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        std::cout << "iProc = " << iProc << "\n";
        int nZones = this->GetNZones( iProc );
        std::cout << "nZones in Proc " << iProc << "  = " << nZones << "\n";

        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            SplitZone * splitZone = this->GetSplitZone( iProc, iZone );
            //splitZone->SetZoneIndex( globalZoneId++ );
            splitZone->newIndex = ( globalZoneId++ );
            std::cout << "splitZone->newIndex = " << splitZone->newIndex << "\n";
            std::cout << "splitZone->zoneIndex = " << splitZone->zoneIndex << "\n";
        }
    }
}