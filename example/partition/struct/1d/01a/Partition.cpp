#include "Partition.h"
#include "CgnsGrid.h"
#include "cgnslib.h"
#include <string>

DimInfo::DimInfo()
{
    this->dim = 1;
    dimList.resize( this->dim, 1 );
    oriPoint.resize( this->dim, 0 );
}

DimInfo::~DimInfo()
{
}

DimInfo & DimInfo::operator = ( const DimInfo & rhs )
{
    if ( this == & rhs ) return * this;

    this->dim = rhs.dim;
    this->dimList = rhs.dimList;
    this->oriPoint = rhs.oriPoint;

    return * this;
}

int DimInfo::ComputeNumNodes()
{
    int nNodes = 1;
    for ( int m = 0; m < this->dimList.size(); ++ m )
    {
        nNodes *= this->dimList[ m ];
    }
    return nNodes;
}

void DimInfo::ComputeISize()
{
    int index_dim = this->dimList.size();
    this->isize.resize( index_dim * 3, 0 );
    for ( int m = 0; m < index_dim; ++ m )
    {
        isize[ m ] = this->dimList[ m ];
    }

    for ( int m = 0; m < index_dim; ++ m )
    {
        int start = index_dim * 1;
        isize[ start + m ] = this->dimList[ m ] - 1;
    }

    for ( int m = 0; m < index_dim; ++ m )
    {
        int start = index_dim * 2;
        isize[ start + m ] = 0;
    }

    for ( int i = 0; i < isize.size(); ++ i )
    {
        std::cout << "i = " << i << " isize = " << isize[i] << "\n";
    }
}

void DimInfo::ComputeRminmax( std::vector<int> & ref )
{
    int index_dim = this->dimList.size();
    this->rmin.resize( index_dim );
    this->rmax.resize( index_dim );

    for ( int m = 0; m < index_dim; ++ m )
    {
        rmin[ m ] = ref[ m ] + 1;
        rmax[ m ] = ref[ m ] + this->dimList[ m ];
    }

    std::cout << "rmin, rmax = \n";
    for ( int m = 0; m < index_dim; ++ m )
    {
        std::cout << rmin[ m ] << " " << rmax[ m ] << "\n";
    }
}

void DimInfo::Init( std::vector<int> & dimList )
{
    this->dimList = dimList;
    std::vector<int> zero( dimList.size(), 0 );
    this->oriPoint = zero;
}

void DimInfo::SetNewDimension( int axis, int newStart, int newWidth )
{
    this->oriPoint[ axis ] = newStart;
    this->dimList [ axis ] = newWidth;
}

SplitZone::SplitZone()
{
    this->oriZoneId = -1;
    this->parent  = 0;
    this->zoneIndex = 0;
    this->iProc = 0;
}

SplitZone::~SplitZone()
{
    ;
}

void SplitZone::ComputeRminmax()
{
    int index_dim = this->dimInfo.dimList.size();
    std::vector<int> ref( index_dim, 0 );

    SplitZone * rootSplitZone = 0;
    this->GetRootInfo( rootSplitZone, ref );

    this->oriZoneId = rootSplitZone->oriZoneId;
    std::cout << "this->oriZoneId = " << this->oriZoneId << "\n";

    this->dimInfo.ComputeRminmax( ref );
}

void SplitZone::GetRootInfo( SplitZone *& root, std::vector<int> & ref )
{
    for ( int i = 0; i < this->dimInfo.dimList.size(); ++ i )
    {
        ref[ i ] += this->dimInfo.oriPoint[ i ];
    }

    if ( this->parent )
    {
        this->parent->GetRootInfo( root, ref );
    }
    else
    {
        root = this;
    }
}


int SplitZone::GetNCell()
{
    int numberOfCells = 1;
    for ( int i = 0; i < this->dimInfo.dim; ++ i )
    {
        numberOfCells *= std::max( 1, ( this->dimInfo.dimList[ i ] - 1 ) );
    }

    return numberOfCells; 
}

void SplitZone::SetLeftDimension( std::vector<int> & dimList, int axisNeedSplit, int iSplit )
{
    int leftStart = 0;
    int leftWidth = iSplit - 0;

    this->dimInfo.Init( dimList );
    this->dimInfo.SetNewDimension( axisNeedSplit, leftStart , leftWidth  );
}

void SplitZone::SetRightDimension( std::vector<int> & dimList, int axisNeedSplit, int iSplit )
{
    int rightStart = iSplit - 1;
    int rightWidth = dimList[ axisNeedSplit ] - ( iSplit - 1 );

    this->dimInfo.Init( dimList );
    this->dimInfo.SetNewDimension( axisNeedSplit, rightStart, rightWidth );
}

void SplitZone::Split( SplitZone *& zoneL, SplitZone *& zoneR, double nCell )
{
    //Find the longest axis
    int axisNeedSplit = this->FindAxisNeedSplit();
    int iSplit        = this->FindSplitPosition( axisNeedSplit, nCell );

    zoneL = new SplitZone();
    zoneL->SetParent( this );

    zoneR = new SplitZone();
    zoneR->SetParent( this );

    zoneL->SetLeftDimension( this->dimInfo.dimList, axisNeedSplit, iSplit );
    zoneR->SetRightDimension( this->dimInfo.dimList, axisNeedSplit, iSplit );

    //zoneL->CreateBc();
    //zoneR->CreateBc();

    //CreateInterfaceBc( zoneL, zoneR );
}

int SplitZone::FindAxisNeedSplit()
{
    //Find the longest axis
    int axisNeedSplit = 0;
    int maxN = this->dimInfo.dimList[ 0 ];

    for ( int m = 1; m < this->dimInfo.dimList.size(); ++ m )
    {
        if ( this->IsPoleBc( m ) )
        {
            continue;
        }
        int Nm = this->dimInfo.dimList[ m ];
        if ( maxN < Nm )
        {
            maxN = Nm;
            axisNeedSplit = m;
        }
    }

    return axisNeedSplit;
}

int SplitZone::GetCrossSection( int axis )
{
    int nDim = this->dimInfo.dimList.size();
    int crossSection = 1;
    int N = nDim - 1;
    for ( int m = 0; m < N; ++ m )
    {
         int mm = ( axis + 1 ) % nDim;
         int nLineElements = std::max( 1, ( this->dimInfo.dimList[ mm ] - 1 ) );
         crossSection *= nLineElements;
    }
    return crossSection;
}

int SplitZone::FindSplitPosition( int axisNeedSplit, double nCell )
{
    int nDim = this->dimInfo.dimList.size();
    std::vector<int> crossSections( nDim );
    for ( int m = 0; m < nDim; ++ m )
    {
        crossSections[ m ] = this->GetCrossSection( m );
    }

    //guess the value
    double dSplit = nCell / crossSections[ axisNeedSplit ];
    int  iSplit = 1 + static_cast< int >( std::floor( dSplit + 0.5 ) );
    return iSplit;
}


bool SplitZone::IsPoleBc( int splitDirection )
{
    return false;
}


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
    ;
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
    this->nProc = 2;
    this->ReadGrid();
    this->Split();
    this->DumpPartitionedGridFile();
    this->FreeAllZones();
}

void Partition::ReadGrid()
{
    ReadCgnsGrid( "../heat1dni81.cgns" );

    this->nOriZone = ::Global::zones.size();

    std::cout << "numberOfBlocks = " << nOriZone << "\n";

    this->refZone.resize( nOriZone );

    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        refZone[ iZone ] = new SplitZone();
        refZone[ iZone ]->SetZoneIndex( iZone );
        refZone[ iZone ]->oriZoneId = iZone;
        refZone[ iZone ]->SetProcId( 0 );
    }

    for ( int iZone = 0; iZone < nOriZone; ++ iZone )
    {
        SplitZone * splitZone = refZone[ iZone ];

        Zone * zone = ::Global::zones[ iZone ];

        splitZone->dimInfo.dimList = zone->nijk;

        int nBcRegion = zone->bccos.size();
        std::cout << "nBcRegion = " << nBcRegion << "\n";

        //HXVector< ZoneBc * > & bcList = splitZone->GetZoneBcList();
        //bcList.resize( nBcRegion );

        for ( int iBcRegion = 0; iBcRegion < nBcRegion; ++ iBcRegion )
        {
        }

        //splitZone->NormalizeBc();
    }
}

void Partition::Split()
{
    this->InitSplit();
    this->BinarySplit();
    //this->SetZoneIndex();
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

    this->minNode = 3 * 3 * 2;
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
    SplitZone * zoneL = 0;
    SplitZone * zoneR = 0;

    zone->Split( zoneL, zoneR, nCell );

    this->AddNewZone( zoneL );
    this->AddNewZone( zoneR );
    this->RemoveZoneFromUnsignedGroup( zone );

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

        int kkk = 1;

        //It seems that EPS is still necessary
        //The idea here is as follows. For the candidate blocks which can be divided into blocks, 
        //if the difference between the target blocks and the candidate blocks is too wide, 
        //the binarySplit method is adopted. The binarySplit blocks do not directly join the process,
        //but enter the next round of screening.

        if ( tau > idealUniProcTime * ( 1.0 + eps ) )
        {
            double nRest = ( tau - idealUniProcTime ) * procSpeed[ iProc ] / aveFloatOp;

            double nTarget = nMaxCell - nRest;
            if ( nRest > this->minNode && nTarget > this->minNode )
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

void Partition::DumpPartitionedGridFile()
{
    int fileId = -1;
    std::string filename = "split1d.cgns";
    cg_open( filename.c_str(), CG_MODE_WRITE, &fileId);
    std::cout << "fileId = " << fileId << "\n";

    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        int icelldim = 1;
        int iphysdim = 1;
        std::string basename = "Base" + std::to_string( iProc );

        int baseId = -1;
        cg_base_write( fileId, basename.c_str(), icelldim, iphysdim, &baseId );

        std::cout << "baseId = " << baseId << "\n";

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

            std::string zonename = "B" + std::to_string( iProc ) + "Zone" + std::to_string( iZone );
            cg_zone_write( fileId, baseId, zonename.c_str(), isize.data(), zoneType, &zoneId );

            std::cout << "zoneId = " << zoneId << "\n";
                
            std::vector<cgsize_t> irmin( index_dim );
            std::vector<cgsize_t> irmax( index_dim );
            int nNodes = 1;
            for ( int m = 0; m < index_dim; ++ m )
            {
                /* lower range index */
                irmin[ m ] = 1;
                /* upper range index of vertices */
                irmax[ m ] = isize[ m ];
                nNodes *= irmax[ m ];
            }
            std::cout << "nNodes = " << nNodes << "\n";

            std::cout << "irmin, irmax = \n";
            for ( int m = 0; m < index_dim; ++ m )
            {
                std::cout << irmin[ m ] << " " << irmax[ m ] << "\n";
            }

            //std::vector<int> ref( index_dim, 0 );

            //SplitZone * rootSplitZone = 0;
            //splitZone->GetRootInfo( rootSplitZone, ref );

            //int oriZoneId = rootSplitZone->oriZoneId;
            //std::cout << "oriZoneId = " << oriZoneId << "\n";
            splitZone->ComputeRminmax();

            Zone * zone = Global::zones[ splitZone->oriZoneId ];

            //std::vector<cgsize_t> rmin( index_dim );
            //std::vector<cgsize_t> rmax( index_dim );

            //for ( int m = 0; m < index_dim; ++ m )
            //{
            //    rmin[ m ] = ref[ m ] + 1;
            //    rmax[ m ] = ref[ m ] + splitZone->dimInfo.dimList[ m ];
            //}

            //std::cout << "rmin, rmax = \n";
            //for ( int m = 0; m < index_dim; ++ m )
            //{
            //    std::cout << rmin[ m ] << " " << rmax[ m ] << "\n";
            //}

            int ncoords = 1;
            std::cout << "ncoords = " << ncoords << "\n";

            for ( int icoord = 0; icoord < ncoords; ++ icoord )
            {
                int coorId = -1;
                DataType_t dataType = RealDouble;
                std::string coordname = "X";
                Coor * coor = zone->coors[ icoord ];

                int nNodes = splitZone->dimInfo.ComputeNumNodes();
                std::cout << "nNodes = " << nNodes << "\n";

                std::vector<char> coord( nNodes * sizeof(double) );
                {
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
                        int ni = splitZone->dimInfo.dimList[ 0 ];
                        int nj = splitZone->dimInfo.dimList[ 1 ];
                        for ( int j = jmin; j < jmax; ++ j )
                        {
                            int jstride = j * ni;
                            for ( int i = imin; i < imax; ++ i )
                            {
                                xlocal[ icount++ ] = xglobal[ jstride + i ];
                            }
                        }
                    }
                    //for ( int i = 0; i < nNodes; ++ i )
                    //{
                    //    int j = i + ref[ 0 ];
                    //    xlocal[ i ] = xglobal[j];
                    //    std::cout << xlocal[ i ] << " ";
                    //    if ( ( i + 1 ) % 5 == 0 ) std::cout << "\n";
                    //}
                    //std::cout << "\n";
                }

                cg_coord_write( fileId, baseId, zoneId, dataType, coordname.c_str(), coord.data(), &coorId );

                std::cout << "cg_get_error() = " << cg_get_error() << "\n";
                std::cout << "fileId = " << fileId  << " baseId = " << baseId << " zoneId = " << zoneId  << " coorId = " << coorId << "\n";
            }

            //int nbocos = -1;

            //cg_nbocos( fileId, baseId, zoneId, &nbocos );
            //std::cout << "nbocos = " << nbocos << "\n";
            //for ( int iboco = 0; iboco < nbocos; ++ iboco )
            //{
            //    int bccoId = iboco + 1;
            //    GridLocation_t location;
            //    cg_boco_gridlocation_read( fileId, baseId, zoneId, bccoId, &location );
            //    std::cout << "iboco = " << iboco <<  " location = " << location << " GridLocationName = " << GridLocationName[location] << "\n";

            //    ZoneBc * zonebc = new ZoneBc();
            //    zone->bccos.push_back( zonebc );

            //    char boconame[ 33 ];
            //    BCType_t bocotype;
            //    PointSetType_t ptset_type;
            //    cgsize_t npnts = 0;
            //    std::vector<int> normalIndex(index_dim,-1);
            //    cgsize_t normalListSize = 0;
            //    DataType_t normalDataType;
            //    int ndataset = -1;

            //    cg_boco_info( fileId, baseId, zoneId, bccoId, boconame, &bocotype, &ptset_type, 
            //        &npnts, normalIndex.data(), &normalListSize, &normalDataType, &ndataset);
            //    std::cout << "boconame = " << boconame <<  " bocotype = " << bocotype << " BCTypeName = " << BCTypeName[bocotype] << "\n";
            //    std::cout << "ptset_type = " << ptset_type <<  " PointSetTypeName = " << PointSetTypeName[ptset_type] << "\n";
            //    std::cout << "npnts = " << npnts << "\n";
            //    std::cout << "normalIndex = ";
            //    for ( int i = 0; i < index_dim; ++ i )
            //    {
            //        std::cout << normalIndex[i] << " ";
            //    }
            //    std::cout << "\n";
            //    std::cout << "normalListSize = " << normalListSize << "\n";
            //    std::cout << "normalDataType = " << normalDataType << " DataTypeName = " << DataTypeName[normalDataType] << "\n";
            //    std::cout << "ndataset = " << ndataset << "\n";

            //    std::vector<char> normalList(nNodes*iphysdim*sizeof(double));

            //    std::vector<cgsize_t> pnts(npnts*index_dim);
            //    cg_boco_read( fileId, baseId, zoneId, bccoId, pnts.data(), normalList.data() );
            //    std::cout << "pnts = ";
            //    for ( int i = 0; i < pnts.size(); ++ i )
            //    {
            //        std::cout << pnts[ i ] << " ";
            //    }
            //    std::cout << "\n";

            //    double * normal_d = reinterpret_cast<double*>(const_cast<char *>(normalList.data()));

            //    //std::cout << "normalList = ";
            //    //for ( int i = 0; i < nNodes*iphysdim; ++ i )
            //    //{
            //    //    std::cout << normal_d[ i ] << " ";
            //    //}
            //    //std::cout << "\n";
            //}
        }
    }
    cg_close(fileId);
}