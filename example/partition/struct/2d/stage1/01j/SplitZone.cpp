#include "SplitZone.h"
#include "cgnslib.h"
#include <string>

int Dim::dim = -1;

void GetDomainRegion( RangeRegion & zoneBox, int iDomain, RangeRegion & domainBox  )
{
    domainBox.start = zoneBox.start;
    domainBox.end = zoneBox.end;

    int axis = iDomain / 2;
    int direction = iDomain % 2;

    if ( direction == 0 )
    {
        domainBox.end[ axis ] = domainBox.start[ axis ];
    }
    else if ( direction == 1 )
    {
        domainBox.start[ axis ] = domainBox.end[ axis ];
    }
}

DimInfo::DimInfo()
{
    this->Init( 1 );
}

DimInfo::~DimInfo()
{
}

void DimInfo::Init( int dim )
{
    this->dim = dim;
    dimList.resize( this->dim, 1 );
    oriPoint.resize( this->dim, 0 );
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
    this->dim = dimList.size();
    this->dimList = dimList;
    std::vector<int> zero( dimList.size(), 0 );
    this->oriPoint = zero;
}

void DimInfo::SetNewDimension( int axis, int newStart, int newWidth )
{
    this->oriPoint[ axis ] = newStart;
    this->dimList [ axis ] = newWidth;
}

void DimInfo::GetZoneGlobalBox( RangeRegion &zoneBox )
{
    int index_dim = this->dimList.size();
    zoneBox.start.resize( index_dim );
    zoneBox.end.resize( index_dim );
    for ( int m = 0; m < index_dim; ++ m )
    {
        zoneBox.start[ m ] = this->oriPoint[ m ] + 1;
        zoneBox.end[ m ] = this->oriPoint[ m ] + this->dimList[ m ];
    }
}

void DimInfo::GetZoneLocalBox( RangeRegion & globalBox, RangeRegion & localBox )
{
    int index_dim = this->dimList.size();
    localBox.start.resize( index_dim );
    localBox.end.resize( index_dim );

    for ( int m = 0; m < index_dim; ++ m )
    {
        localBox.start[ m ] = globalBox.start[ m ] - this->oriPoint[ m ];
        localBox.end[ m ] = globalBox.end[ m ] - this->oriPoint[ m ];
    }
}

RangeRegion::RangeRegion( int nSize )
{
    this->Init( nSize );
}

RangeRegion::~RangeRegion()
{
}

void RangeRegion::Init( int nSize )
{
    this->start.resize( nSize );
    this->end.resize( nSize );
}

RangeRegion & RangeRegion::operator = ( const RangeRegion & rhs )
{
    if ( this == & rhs ) return * this;

    this->start = rhs.start;
    this->end = rhs.end;

    return * this;
}

bool RangeRegion::operator == ( const RangeRegion & rhs ) const
{
    int nSize = this->start.size();
    for ( int m = 0; m < nSize; ++ m )
    {
        if ( this->start[ m ] != rhs.start[ m ] || 
            this->end  [ m ] != rhs.end  [ m ] )
        {
            return false;
        }
    }
    return true;
}

void RangeRegion::SetRegion( std::vector<int> & pnts )
{
    int index_dim = pnts.size() / 2;
    for ( int m = 0; m < index_dim; ++ m )
    {
        this->start[ m ] = pnts[ m ];
        this->end[ m ] = pnts[ index_dim + m ];
    }
}

bool RangeRegion::InRegion( const RangeRegion & region )
{
    int nSize = this->start.size();
    for ( int m = 0; m < nSize; ++ m )
    {
        if ( ( this->start[ m ] < region.start[ m ] ) || ( this->end[ m ] > region.end[ m ] ) )
        {
            return false;
        }
    }
    return true;
}

void RangeRegion::Normalize()
{
    int nSize = this->start.size();
    for ( int m = 0; m < nSize; ++ m )
    {
        int minv = std::min( this->start[ m ], this->end[ m ] );
        int maxv = std::max( this->start[ m ], this->end[ m ] );
        this->start[ m ] = minv;
        this->end[ m ] = maxv;
    }
}

void RangeRegion::Localize( std::vector<int> & oriPoint )
{
    int nSize = this->start.size();
    for ( int m = 0; m < nSize; ++ m )
    {
        this->start[ m ] -= oriPoint[ m ];
        this->end[ m ] -= oriPoint[ m ];
    }
}

void RangeRegion::AddShift( std::vector<int> & oriPoint )
{
    int nSize = this->start.size();
    for ( int m = 0; m < nSize; ++ m )
    {
        this->start[ m ] += oriPoint[ m ];
        this->end[ m ] += oriPoint[ m ];
    }
}

void RangeRegion::PrintGlobal( std::vector<int> & oriPoint )
{
    int nSize = this->start.size();
    std::cout << "global start:(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->start[ m ] + oriPoint[m];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
    std::cout << "global end  :(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->end[ m ] + oriPoint[m];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
}

void RangeRegion::Print()
{
    int nSize = this->start.size();
    std::cout << "start:(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->start[ m ];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
    std::cout << "end  :(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->end[ m ];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
}

bool CheckOverLapping( RangeRegion & overlapRegion );

bool ComputeOverlapRegion( RangeRegion & region1, RangeRegion & region2, RangeRegion & overlapRegion )
{
    int index_dim = region1.start.size();
    overlapRegion.Init( index_dim );
    for ( int m = 0; m < index_dim; ++ m )
    {
        int ijkmin1 = std::min( region1.start[ m ], region1.end[ m ] );
        int ijkmax1 = std::max( region1.start[ m ], region1.end[ m ] );
        int ijkmin2 = std::min( region2.start[ m ], region2.end[ m ] );
        int ijkmax2 = std::max( region2.start[ m ], region2.end[ m ] );
        if ( ( ijkmax1 < ijkmin2 ) || ( ijkmin1 > ijkmax2 ) )
        {
            return false;
        }

        int ijkmin = std::max( ijkmin1, ijkmin2 );
        int ijkmax = std::min( ijkmax1, ijkmax2 );

        overlapRegion.start[ m ] = ijkmin;
        overlapRegion.end  [ m ] = ijkmax;
    }

    return CheckOverLapping( overlapRegion );
}

bool CheckOverLapping( RangeRegion & overlapRegion )
{
    int index_dim = overlapRegion.start.size();
    int iCount = 0;
    for ( int m = 0; m < index_dim; ++ m )
    {
        if ( overlapRegion.start[ m ] == overlapRegion.end[ m ] )
        {
            ++ iCount;
        }
    }

    return ( iCount == 1 );
}

bool GetInterfaceRegion( RangeRegion & zoneBoxL, RangeRegion & zoneBoxR, RangeRegion & interfaceRegion )
{
    int index_dim = zoneBoxL.start.size();
    int nDomains = 2 * index_dim;

    for ( int iDomain = 0; iDomain < nDomains; ++ iDomain )
    {
        RangeRegion domainRegionL;
        GetDomainRegion( zoneBoxL, iDomain, domainRegionL );

        for ( int jDomain = 0; jDomain < nDomains; ++ jDomain )
        {
            RangeRegion domainRegionR;
            GetDomainRegion( zoneBoxR, jDomain, domainRegionR );

            if ( domainRegionL == domainRegionR )
            {
                interfaceRegion = domainRegionL;
                return true;
            }
        }
    }
    return false;
}

void CreateInterfaceBc( SplitZone * zoneL, SplitZone * zoneR )
{
    RangeRegion globalBoxL, globalBoxR, globalInterfaceRegion;

    zoneL->GetZoneGlobalBox( globalBoxL );
    zoneR->GetZoneGlobalBox( globalBoxR );

    if ( ! GetInterfaceRegion( globalBoxL, globalBoxR, globalInterfaceRegion ) )
    {
        return;
    }

    RangeRegion localInterfaceRegionL, localInterfaceRegionR;

    zoneL->GetZoneLocalBox( globalInterfaceRegion, localInterfaceRegionL );
    zoneR->GetZoneLocalBox( globalInterfaceRegion, localInterfaceRegionR );

    int kkk = 1;

    InterfaceSplitBc * interfaceSplitBcL = new InterfaceSplitBc();
    zoneL->interfaceSplitBcList.push_back( interfaceSplitBcL );
    interfaceSplitBcL->zone = zoneL;
    interfaceSplitBcL->region = localInterfaceRegionL;

    int nSize = localInterfaceRegionL.start.size();

    std::vector<int> transform;

    for ( int m = 0; m < nSize; ++ m )
    {
        transform.push_back( m + 1 );
    }

    interfaceSplitBcL->transform = transform;
    interfaceSplitBcL->CalcTransformMatrix();

    interfaceSplitBcL->donor_zone = zoneR;
    interfaceSplitBcL->donor_region = localInterfaceRegionR;

    InterfaceSplitBc * interfaceSplitBcR = new InterfaceSplitBc();
    zoneR->interfaceSplitBcList.push_back( interfaceSplitBcR );

    interfaceSplitBcR->zone = zoneR;
    interfaceSplitBcR->region = localInterfaceRegionR;

    interfaceSplitBcR->transform = transform;
    interfaceSplitBcR->CalcTransformMatrix();

    interfaceSplitBcR->donor_zone = zoneL;
    interfaceSplitBcR->donor_region = localInterfaceRegionL;
}


PhysicalSplitBc::PhysicalSplitBc()
{
    this->region.Init( Dim::dim );
}

PhysicalSplitBc::~PhysicalSplitBc()
{
}

void PhysicalSplitBc::SetRegion( std::vector<int> & pnts )
{
    this->region.SetRegion( pnts );
}

void PhysicalSplitBc::SetRegion( RangeRegion & region )
{
    this->region = region;
}

void PhysicalSplitBc::ChangeRegionToLocalCoordinate()
{
    int index_dim = this->region.start.size();
    std::vector<int> & oriPoint = splitZone->dimInfo.oriPoint;
    for ( int m = 0; m < index_dim; ++ m )
    {
        this->region.start[ m ] -= oriPoint[ m ];
        this->region.end[ m ] -= oriPoint[ m ];
    }
}

BasicSplitBc::BasicSplitBc()
{
    this->region.Init( Dim::dim );
}

BasicSplitBc::~BasicSplitBc()
{
    ;
}

void BasicSplitBc::SetRegion( std::vector<int> & pnts )
{
    this->region.SetRegion( pnts );
}

void BasicSplitBc::SetRegion( RangeRegion & region )
{
    this->region = region;
}

int trans::M[ 3 ][ 3 ];
std::vector<int> trans::transform;

int trans::sgn( int x )
{
    if ( x >= 0 )
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

int trans::del( int x, int y )
{
    if ( std::abs( x ) == std::abs( y ) )
    {
        return 1;
    }
    return 0;
}

void trans::ZeroMatrix()
{
    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            trans::M[ i ][ j ] = 0;
        }
    }
}

void trans::CalcTransformMatrix()
{
    int dim = trans::transform.size();
    if ( dim == 1 )
    {
        int a = trans::transform[ 0 ];
        int sgna = trans::sgn( a );
        int a1 = trans::del( a, 1 );
        trans::M[ 0 ][ 0 ] = sgna * a1;
    }
    else if ( dim == 2 )
    {
        int a = trans::transform[ 0 ];
        int b = trans::transform[ 1 ];
        int sgna = trans::sgn( a );
        int sgnb = trans::sgn( b );
        int a1 = trans::del( a, 1 );
        int a2 = trans::del( a, 2 );
        int b1 = trans::del( b, 1 );
        int b2 = trans::del( b, 2 );
        trans::M[ 0 ][ 0 ] = sgna * a1;
        trans::M[ 1 ][ 0 ] = sgna * a2;
        trans::M[ 0 ][ 1 ] = sgnb * b1;
        trans::M[ 1 ][ 1 ] = sgnb * b2;
    }
    else if ( dim == 3 )
    {
        int a = trans::transform[ 0 ];
        int b = trans::transform[ 1 ];
        int c = trans::transform[ 2 ];
        int sgna = trans::sgn( a );
        int sgnb = trans::sgn( b );
        int sgnc = trans::sgn( c );
        int a1 = trans::del( a, 1 );
        int a2 = trans::del( a, 2 );
        int a3 = trans::del( a, 3 );
        int b1 = trans::del( b, 1 );
        int b2 = trans::del( b, 2 );
        int b3 = trans::del( b, 3 );
        int c1 = trans::del( c, 1 );
        int c2 = trans::del( c, 2 );
        int c3 = trans::del( c, 3 );
        trans::M[ 0 ][ 0 ] = sgna * a1;
        trans::M[ 1 ][ 0 ] = sgna * a2;
        trans::M[ 2 ][ 0 ] = sgna * a3;
        trans::M[ 0 ][ 1 ] = sgnb * b1;
        trans::M[ 1 ][ 1 ] = sgnb * b2;
        trans::M[ 2 ][ 1 ] = sgnb * b3;
        trans::M[ 0 ][ 2 ] = sgnc * c1;
        trans::M[ 1 ][ 2 ] = sgnc * c2;
        trans::M[ 2 ][ 2 ] = sgnc * c3;
    }
}

InterfaceInfo::InterfaceInfo()
{
    this->donor_region.Init( Dim::dim );
}

InterfaceInfo::~InterfaceInfo()
{
}

InterfaceSplitBc::InterfaceSplitBc()
{
    this->region.Init( Dim::dim );
    this->donor_region.Init( Dim::dim );
}

InterfaceSplitBc::~InterfaceSplitBc()
{
}

void InterfaceSplitBc::CopyMatrix( InterfaceSplitBc * interfaceSplitBc )
{
    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            this->Mt[ i ][ j ] = interfaceSplitBc->Mt[ i ][ j ];
        }
    }
}

void InterfaceSplitBc::CalcTransformMatrix()
{
    trans::ZeroMatrix();
    trans::transform = this->transform;
    trans::CalcTransformMatrix();

    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            this->Mt[ i ][ j ] = trans::M[ i ][ j ];
        }
    }

}

void InterfaceSplitBc::mapindex( std::vector<int> & begin1, std::vector<int> & begin2,  std::vector<int> & index1, std::vector<int> & index2 )
{
    int dim = Dim::dim;
    std::vector<int> diff( dim );
    for ( int m = 0; m < dim; ++ m )
    {
        diff[ m ] = index1[ m ] - begin1[ m ];
    }

    std::vector<int> mul( dim );

    Multiply( Mt, diff, mul );

    for ( int m = 0; m < dim; ++ m )
    {
        index2[ m ] = mul[ m ] + begin2[ m ];
    }
}

void InterfaceSplitBc::Multiply( int Mt[][3], std::vector<int> & a, std::vector<int> & b )
{
    int dim = Dim::dim;
    for ( int i = 0; i < dim; ++ i )
    {
        b[ i ] = 0;
        for ( int j = 0; j < dim; ++ j )
        {
            b[ i ] += Mt[ i ][ j ] * a[ j ];
        }
    }
}

void InterfaceSplitBc::ChangeRegionToLocalCoordinate()
{
    int index_dim = this->region.start.size();
    std::vector<int> & oriPoint = this->zone->dimInfo.oriPoint;
    for ( int m = 0; m < index_dim; ++ m )
    {
        this->region.start[ m ] -= oriPoint[ m ];
        this->region.end[ m ] -= oriPoint[ m ];
    }
}

void InterfaceSplitBc::CalcSubDonorRegion( RangeRegion &subRegion, RangeRegion &subDonorRegion )
{
    std::vector<int> begin1 = this->region.start;
    std::vector<int> end1 = this->region.end;
    std::vector<int> begin2 = this->donor_region.start;
    std::vector<int> end2 = this->donor_region.end;

    std::vector<int> new_begin1 = subRegion.start;
    std::vector<int> new_end1 = subRegion.end;

    int dim = Dim::dim;
    std::vector<int> new_begin2( dim );
    std::vector<int> new_end2( dim );

    mapindex( begin1, begin2, new_begin1, new_begin2 );
    mapindex( begin1, begin2, new_end1, new_end2 );

    subDonorRegion.start = new_begin2;
    subDonorRegion.end = new_end2;
}

SplitZone::SplitZone()
{
    this->parent  = 0;
    this->zoneIndex = 0;
    this->iProc = 0;
    dimInfo.Init( Dim::dim );
    int kk = 1;
}

SplitZone::~SplitZone()
{
    ;
}

void SplitZone::GetZoneGlobalBox( RangeRegion & zoneBox )
{
    this->dimInfo.GetZoneGlobalBox( zoneBox );
}

void SplitZone::GetZoneLocalBox( RangeRegion & globalBox, RangeRegion & localBox )
{
    this->dimInfo.GetZoneLocalBox( globalBox, localBox );
}


void SplitZone::ComputeRminmax()
{
    int index_dim = this->dimInfo.dimList.size();
    std::vector<int> ref( index_dim, 0 );

    SplitZone * rootSplitZone = 0;
    this->GetRootInfo( rootSplitZone, ref );

    this->dimInfo.ComputeRminmax( ref );
}

SplitZone * SplitZone::GetRootZone()
{
    SplitZone * root = this;

    while ( true )
    {
        if ( root->parent )
        {
            root = root->parent;
        }
        else
        {
            return root;
        }
    };
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

void SplitZone::SetParentAndChild( SplitZone * parent )
{
    this->parent = parent;
    if ( parent )
    {
        parent->child.push_back( this );
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

    zoneL->SetLeftDimension( this->dimInfo.dimList, axisNeedSplit, iSplit );
    zoneR->SetRightDimension( this->dimInfo.dimList, axisNeedSplit, iSplit );

    zoneL->CreateBcFromParent();
    zoneR->CreateBcFromParent();
    
    CreateInterfaceBc( zoneL, zoneR );
}

void SplitZone::CreateBcFromParent()
{
    this->CreatePhysicalBcFromParent();
    this->CreateInterfaceBcFromParent();
}

void SplitZone::CreatePhysicalBcFromParent()
{
    int nParentPhysicalBcs = this->parent->physicalSplitBcList.size();

    for ( int iParentPhysicalBcs = 0; iParentPhysicalBcs < nParentPhysicalBcs; ++ iParentPhysicalBcs )
    {
        PhysicalSplitBc * parentPhysicalSplitBc = this->parent->physicalSplitBcList[ iParentPhysicalBcs ];

        this->CreatePhysicalBcFromParent( parentPhysicalSplitBc );
    }
}

void SplitZone::CreatePhysicalBcFromParent( PhysicalSplitBc * parentPhysicalSplitBc )
{
    int index_dim = this->dimInfo.dimList.size();
    int nDomains = 2 * index_dim;
    RangeRegion zoneBoxRegion( index_dim );
    this->dimInfo.GetZoneGlobalBox( zoneBoxRegion );

    for ( int iDomain = 0; iDomain < nDomains; ++ iDomain )
    {
        RangeRegion domainRegion( index_dim );
        GetDomainRegion( zoneBoxRegion, iDomain, domainRegion );

        RangeRegion overlapRegion;
        bool overlap_flag = ComputeOverlapRegion( parentPhysicalSplitBc->region, domainRegion, overlapRegion );
        if ( ! overlap_flag ) continue;

        this->CreatePhysicalBc( overlapRegion, parentPhysicalSplitBc->bcType );

        int kkk = 1;
    }
}

void SplitZone::CreateInterfaceBcFromParent()
{
    int nInterfaceSplitBcs = this->parent->interfaceSplitBcList.size();

    for ( int iInterfaceSplitBcs = 0; iInterfaceSplitBcs < nInterfaceSplitBcs; ++ iInterfaceSplitBcs )
    {
        InterfaceSplitBc * parentInterfaceSplitBc = this->parent->interfaceSplitBcList[ iInterfaceSplitBcs ];

        this->CreateInterfaceBcFromParent( parentInterfaceSplitBc );
    }
}

void SplitZone::CreateInterfaceBcFromParent( InterfaceSplitBc * parentInterfaceSplitBc )
{
    int index_dim = this->dimInfo.dimList.size();
    int nDomains = 2 * index_dim;
    RangeRegion zoneBoxRegion( index_dim );
    this->dimInfo.GetZoneGlobalBox( zoneBoxRegion );

    for ( int iDomain = 0; iDomain < nDomains; ++ iDomain )
    {
        RangeRegion domainRegion( index_dim );
        GetDomainRegion( zoneBoxRegion, iDomain, domainRegion );

        RangeRegion overlapRegion;
        bool overlap_flag = ComputeOverlapRegion( parentInterfaceSplitBc->region, domainRegion, overlapRegion );
        if ( ! overlap_flag ) continue;

        this->CreateInterfaceBcFromParent( parentInterfaceSplitBc, overlapRegion );

        int kkk = 1;
    }
}

void SplitZone::FindDonorInterface( InterfaceSplitBc * interfaceSplitBc, InterfaceInfo & interfaceInfo )
{
    RangeRegion donor_region = interfaceSplitBc->donor_region;
    donor_region.Normalize();
    interfaceInfo.donor_zone = this;
    for ( int iInterface = 0; iInterface < interfaceSplitBcList.size(); ++ iInterface )
    {
        InterfaceSplitBc * candidateDonorBc = interfaceSplitBcList[ iInterface ];
        bool flag = donor_region.InRegion( candidateDonorBc->region );
        if ( flag )
        {
            if ( donor_region == candidateDonorBc->region )
            {
                interfaceInfo.donor_interface = candidateDonorBc;
                return;
            }
            else
            {
                InterfaceSplitBc * donorBc = new InterfaceSplitBc();
                candidateDonorBc->child.push_back( donorBc );
                donorBc->region = donor_region;
                donorBc->zone = candidateDonorBc->zone;
                donorBc->donor_zone = interfaceSplitBc->zone;
                donorBc->transform = candidateDonorBc->transform;
                donorBc->CopyMatrix( candidateDonorBc );
                RangeRegion subDonorRegion( Dim::dim );
                candidateDonorBc->CalcSubDonorRegion( donorBc->region, subDonorRegion );
                std::vector<int> & oriPoint = interfaceSplitBc->zone->dimInfo.oriPoint;
                subDonorRegion.Localize( oriPoint );
                donorBc->donor_region = subDonorRegion;
                donorBc->donor_bc = interfaceSplitBc;
                interfaceInfo.donor_interface = donorBc;
                return;
            }
        }
    }
    return;
}

void SplitZone::GetChildDonorRegion( RangeRegion & subDonorRegion, InterfaceInfo & interfaceInfo )
{
    std::cout << "subDonorRegion = \n";
    subDonorRegion.Print();
    int nChild = this->child.size();
    for ( int ichild = 0; ichild < nChild; ++ ichild )
    {
        std::cout << " iChild = " << ichild << " nChild = " << nChild << "\n";
        SplitZone * sub_zone = this->child[ ichild ];
        std::cout << " zoneIndex = " << sub_zone->zoneIndex << "\n";
        int nInters = sub_zone->interfaceSplitBcList.size();
        for ( int iInter = 0; iInter < nInters; ++ iInter )
        {
            std::cout << " iInter = " << iInter << "\n";
            InterfaceSplitBc * interfaceSplitBc = sub_zone->interfaceSplitBcList[ iInter ];
            std::cout << "interfaceSplitBc local = \n";
            interfaceSplitBc->region.Print();
            std::cout << "interfaceSplitBc global = \n";
            interfaceSplitBc->region.PrintGlobal( sub_zone->dimInfo.oriPoint );
            RangeRegion g_region = interfaceSplitBc->region;
            g_region.AddShift( sub_zone->dimInfo.oriPoint );
            if ( g_region == subDonorRegion )
            {
                interfaceInfo.donor_zone = sub_zone;
                interfaceInfo.donor_interface = interfaceSplitBc;
                interfaceInfo.donor_region = interfaceSplitBc->region;
                return;
            }
        }
    }
}

void SplitZone::CreateInterfaceBcFromParent( InterfaceSplitBc * parentInterfaceSplitBc, RangeRegion & overlapRegion )
{
    InterfaceSplitBc * interfaceSplitBc = new InterfaceSplitBc();
    this->interfaceSplitBcList.push_back( interfaceSplitBc );

    interfaceSplitBc->region = overlapRegion;
    interfaceSplitBc->zone = this;
    interfaceSplitBc->transform = parentInterfaceSplitBc->transform;
    interfaceSplitBc->CopyMatrix( parentInterfaceSplitBc );
    RangeRegion subDonorRegion( Dim::dim );
    parentInterfaceSplitBc->CalcSubDonorRegion( overlapRegion, subDonorRegion );
    interfaceSplitBc->ChangeRegionToLocalCoordinate();
    interfaceSplitBc->donor_region = subDonorRegion;
    SplitZone * parent_donor_zone = parentInterfaceSplitBc->donor_zone;
    std::cout << " parentInterfaceSplitBc->zone->zoneIndex " << parentInterfaceSplitBc->zone->zoneIndex << "\n";
    std::cout << " parentInterfaceSplitBc->donor_zone->zoneIndex " << parent_donor_zone->zoneIndex << "\n";
    if ( parent_donor_zone->child.size() > 0 )
    {
        InterfaceInfo interfaceInfo;
        parent_donor_zone->GetChildDonorRegion( subDonorRegion, interfaceInfo );
        interfaceSplitBc->donor_zone = interfaceInfo.donor_zone;
        interfaceSplitBc->donor_bc = interfaceInfo.donor_interface;
        interfaceSplitBc->donor_region = interfaceInfo.donor_region;
        //compatibility
        InterfaceSplitBc * donor_bc = interfaceSplitBc->donor_bc;
        SplitZone * donor_zone = interfaceSplitBc->donor_zone;
        //if ( donor_bc->donor_zone->zoneIndex == interfaceSplitBc->zone->zoneIndex )
        if ( interfaceSplitBc->zone != donor_bc->donor_zone )
        {
            donor_bc->donor_zone = interfaceSplitBc->zone;
            donor_bc->donor_region = interfaceSplitBc->region;
           
            int kkk = 1;
        }
        int kkk = 1;
    }
    else
    {
        InterfaceInfo interfaceInfo;
        parent_donor_zone->FindDonorInterface( interfaceSplitBc, interfaceInfo );
        interfaceSplitBc->donor_zone = interfaceInfo.donor_zone;
        interfaceSplitBc->donor_bc = interfaceInfo.donor_interface;
    }
    
    int kkk = 1;
}

void SplitZone::CreatePhysicalBc( RangeRegion & overlapRegion, int bcType )
{
    PhysicalSplitBc * physicalSplitBc = new PhysicalSplitBc();
    this->physicalSplitBcList.push_back( physicalSplitBc );
    physicalSplitBc->SetRegion( overlapRegion );
    physicalSplitBc->bcType = bcType;
    physicalSplitBc->splitZone = this;
    physicalSplitBc->ChangeRegionToLocalCoordinate();
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

