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
        int kkk = 1;
        return;
    }

    RangeRegion localInterfaceRegionL, localInterfaceRegionR;

    zoneL->GetZoneLocalBox( globalInterfaceRegion, localInterfaceRegionL );
    zoneR->GetZoneLocalBox( globalInterfaceRegion, localInterfaceRegionR );

    int kkk = 1;

    SplitBc * splitBcL = new SplitBc();
    zoneL->AddSplitBc( splitBcL );

    SplitBc * splitBcR = new SplitBc();
    zoneR->AddSplitBc( splitBcR );

    splitBcL->splitZone = zoneL;
    splitBcL->SetRegion( localInterfaceRegionL );
    splitBcL->bcType = BCTypeUserDefined;

    splitBcR->splitZone = zoneR;
    splitBcR->SetRegion( localInterfaceRegionR );
    splitBcR->bcType = BCTypeUserDefined;

    //localInterfaceRegionL.MarkPatchBc();
    //localInterfaceRegionR.MarkPatchBc();

    SplitBcPatch * splitBcPatchL = new SplitBcPatch();
    splitBcL->patch = splitBcPatchL;
    splitBcPatchL->region = localInterfaceRegionR;
    splitBcPatchL->splitZone = zoneR;

    SplitBcPatch * splitBcPatchR = new SplitBcPatch();
    splitBcR->patch = splitBcPatchR;
    splitBcPatchR->region = localInterfaceRegionL;
    splitBcPatchR->splitZone = zoneL;

    //splitBcL->Normalize();
    //splitBcR->Normalize();
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

SplitBc::SplitBc()
{
    this->region.Init( Dim::dim );
}

SplitBc::~SplitBc()
{
    delete this->patch;
}

void SplitBc::SetRegion( std::vector<int> & pnts )
{
    this->region.SetRegion( pnts );
}

void SplitBc::SetRegion( RangeRegion & region )
{
    this->region = region;
}

void SplitBc::ChangeRegionToLocalCoordinate()
{
    int index_dim = this->region.start.size();
    std::vector<int> & oriPoint = splitZone->dimInfo.oriPoint;
    for ( int m = 0; m < index_dim; ++ m )
    {
        this->region.start[ m ] -= oriPoint[ m ];
        this->region.end[ m ] -= oriPoint[ m ];
    }
}

void SplitBc::GetLeafBcs( std::vector< SplitBc * > & leafBcs )
{
    if ( ! this->child )
    {
        leafBcs.push_back( this );
        return;
    }

    for ( int i = 0; i < this->child->size(); ++ i )
    {
        std::vector< SplitBc * > child_leafBcs;
        ( * child )[ i ]->GetLeafBcs( child_leafBcs );
        for ( int j = 0; j < child_leafBcs.size(); ++ j )
        {
            leafBcs.push_back( child_leafBcs[ j ] );
        }
    }

    return;
}

void SplitBc::CreateInterfaceBcFromParentBc( SplitBc * parentSplitBc )
{
    if ( parentSplitBc->bcType !=  BCTypeUserDefined ) return;

    this->patch = new SplitBcPatch();
    this->patch->splitZone = parentSplitBc->patch->splitZone;
    this->patch->SetRegion( parentSplitBc->patch->region );
    int kkk = 1;

    //Find the block corresponding to this docking boundary
    //From patch->splitZone, find the corresponding docking boundary according to the size of overlapt
    //The boundary conditions on the opposite side must change due to the split of this block.
    //Here, the relevant boundaries are added to the sub-boundaries.

    RangeRegion targetRegionBox;
    targetRegionBox = parentSplitBc->patch->region;

    SplitZone * targetZone = this->patch->splitZone;

    targetZone->CreateChildInterfaceBc( this );

    //this->patch->faceAxis = parentSplitBc->patch->faceAxis;

    //( * this->patch->vertexMap ) = ( * parentBc->patch->vertexMap );

    ////然后进行修正
    ////block1->map 1, 33 -> block1 369, 337
    ////block2->global -> map 1, 33 -> block1 369, 337
    ////block2->local  -> map 1, 33 -> block1 369, 337
    ////refOld + M * stold           = block1Startold
    ////refOld + M * ( stnew + ref ) = block1Startold
    ////refNew = refOld + M * ( ref )

    //for ( int m = 0; m < 3; ++ m )
    //{
    //    for ( int n = 0; n < 3; ++ n )
    //    {
    //        this->patch->vertexMap->pt[ m ] += this->patch->vertexMap->vm[ m ][ n ] * ( this->splitZone->sdim->oriPoint )[ n ];
    //    }
    //}

    //RangeBox tRangeBox;

    //MappingPatchPoint( this->patch, this->rangeBox, tRangeBox );
    //tRangeBox.SetOrder();

    ////找到和本对接边界相对应的块
    ////从patch->splitZone按照overlap的大小找到对应的对接边界，这样至少可以找到
    ////对面的对接边界条件由于本块的split必然要发生改变，这里将相关边界加入子边界
    //this->patch->splitZone->CreateChildInterfaceBc( tRangeBox, this->splitZone );
}

SplitBcPatch::SplitBcPatch()
{
    this->region.Init( Dim::dim );
}

SplitBcPatch::~SplitBcPatch()
{
}

void SplitBcPatch::SetRegion( std::vector<int> & pnts )
{
    this->region.SetRegion( pnts );
}

void SplitBcPatch::SetRegion( RangeRegion & region )
{
    this->region = region;
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

InterfaceSplitBc::InterfaceSplitBc()
{
    //this->left = new BasicSplitBc();
    //this->right = new BasicSplitBc();
    this->region.Init( Dim::dim );
    this->donor_region.Init( Dim::dim );
}

InterfaceSplitBc::~InterfaceSplitBc()
{
    //delete this->left;
    //delete this->right;
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

void InterfaceSplitBc::RecalcDonorRegion()
{
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


int SplitZone::GetNCell()
{
    int numberOfCells = 1;
    for ( int i = 0; i < this->dimInfo.dim; ++ i )
    {
        numberOfCells *= std::max( 1, ( this->dimInfo.dimList[ i ] - 1 ) );
    }

    return numberOfCells; 
}

void SplitZone::GetLeafBcList( std::vector< SplitBc * > & leafBcList )
{
    int nSplitBcs = this->splitBcList.size();

    for ( int iSplitBcs = 0; iSplitBcs < nSplitBcs; ++ iSplitBcs )
    {
        SplitBc * splitBc = splitBcList[ iSplitBcs ];

        std::vector< SplitBc * > leafBcs;
        splitBc->GetLeafBcs( leafBcs );
        for ( int j = 0; j < leafBcs.size(); ++ j )
        {
            leafBcList.push_back( leafBcs[ j ] );
        }
    }
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

void SplitZone::CreatePhysicalBcFromParent(  PhysicalSplitBc * parentPhysicalSplitBc )
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

void SplitZone::CreateInterfaceBcFromParent(  InterfaceSplitBc * parentInterfaceSplitBc )
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

void SplitZone::CreateInterfaceBcFromParent( InterfaceSplitBc * parentInterfaceSplitBc, RangeRegion & overlapRegion )
{
    InterfaceSplitBc * interfaceSplitBc = new InterfaceSplitBc();
    interfaceSplitBc->region = overlapRegion;
    interfaceSplitBc->zone = this;
    interfaceSplitBc->donor_zone = parentInterfaceSplitBc->donor_zone;
    interfaceSplitBc->CopyMatrix( parentInterfaceSplitBc );
    RangeRegion subDonorRegion( Dim::dim );
    parentInterfaceSplitBc->CalcSubDonorRegion( overlapRegion, subDonorRegion );
    interfaceSplitBc->donor_region = subDonorRegion;
    SplitZone * splitZone = interfaceSplitBc->donor_zone;
    //splitZone->
}

void SplitZone::CreateChildInterfaceBc( SplitBc * targetInterfaceBc )
{
    std::vector< SplitBc * > & splitBcList = this->GetSplitBcList();

    int nSplitBcs = splitBcList.size();

    for ( int iSplitBc = 0; iSplitBc < nSplitBcs; ++ iSplitBc )
    {
        SplitBc * splitBc = splitBcList[ iSplitBc ];

        //Only consider the interface boundaries
        if ( splitBc->bcType !=  BCTypeUserDefined ) continue;
        SplitBc * candidateInterfaceBc = splitBc;

        if ( this->CreateChildInterfaceBc( candidateInterfaceBc, targetInterfaceBc ) )
        {
            return;
        }
    }
}

bool SplitZone::CreateChildInterfaceBc( SplitBc * candidateInterfaceBc, SplitBc * targetInterfaceBc )
{
    RangeRegion overlapRegion;
    bool overlap_flag = ComputeOverlapRegion( candidateInterfaceBc->region, targetInterfaceBc->patch->region, overlapRegion );
    if ( ! overlap_flag ) return false;

    //1. 如果splitBc和interfaceBc有交集，但是interfaceBc没有child，则此splitBc作为interfaceBc的一个child
    //2. 如果splitBc和interfaceBc有交集，且interfaceBc有child，如果splitBc和所有child没有交集，
    //   则splitBc加入interfaceBc的一个child
    //3. 如果splitBc和interfaceBc有交集，且interfaceBc有child，如果splitBc和child有交集，则此child作为新的interfaceBc

    std::vector< SplitBc * > * child = candidateInterfaceBc->child;

    if ( candidateInterfaceBc->child )
    {
        for ( int iChild = 0; iChild < candidateInterfaceBc->child->size(); ++ iChild )
        {
            SplitBc * childSplitBc = ( * candidateInterfaceBc->child )[ iChild ];
            if ( this->CreateChildInterfaceBc( childSplitBc, targetInterfaceBc ) ) return true;
        }
    }
    else
    {
        candidateInterfaceBc->child = new std::vector< SplitBc * >();
    }

    SplitBc * childInterfaceBc = new SplitBc();

    childInterfaceBc->SetRegion( overlapRegion );
    childInterfaceBc->bcType = candidateInterfaceBc->bcType;
    childInterfaceBc->splitZone = candidateInterfaceBc->splitZone;

    candidateInterfaceBc->child->push_back( childInterfaceBc );

    //vector< int > & parentFaceDirectionContainer = interfaceBc->GetFaceDirectionContainer();
    //vector< int > & faceDirectionContainer  = childInterfaceBc->GetFaceDirectionContainer();

    //for ( int m = 0; m < 3; ++ m )
    //{
    //    faceDirectionContainer[ m ] = parentFaceDirectionContainer[ m ];
    //}

    //将此新生成的对接边界加入本块里面
    SplitBcPatch * patch = new SplitBcPatch();
    childInterfaceBc->patch = patch;
    patch->splitZone = targetInterfaceBc->splitZone;
    patch->region = targetInterfaceBc->region;

    //patch->CopyVertexMappingContent( interfaceBc->GetZoneBoundaryConditionPatch()->GetVertexMapping() );

    //vector< int > & pNextFaceDirectionContainer = interfaceBc->GetZoneBoundaryConditionPatch()->GetFaceDirectionContainer();
    //vector< int > & patchFaceDirectionContainer = patch->GetFaceDirectionContainer();

    //for ( int m = 0; m < 3; ++ m )
    //{
    //    patchFaceDirectionContainer[ m ] = pNextFaceDirectionContainer[ m ];
    //}

    //原对接边界满足
    //t0 = tref + M * s0
    //新的对接边界满足
    //t0newblock + refNewblock = tref + M * s0
    //t0newblock = ( tref - refNewblock ) + M * s0

    //然后进行修正
    //for ( int m = 0; m < 3; ++ m )
    //{
    //    ( patch->GetVertexMapping() )[ m ][ 3 ] -= ( sourceZone->GetOriginalPointIndexContainer() )[ m ];
    //}

    return true;
}

void SplitZone::AddSplitBc( SplitBc * splitBc )
{
    if ( ! splitBc ) return;
    this->splitBcList.push_back( splitBc );
}

void SplitZone::CreateSplitBc( SplitBc * parentSplitBc )
{
    if ( ! parentSplitBc->child )
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
            bool overlap_flag = ComputeOverlapRegion( parentSplitBc->region, domainRegion, overlapRegion );
            if ( ! overlap_flag ) continue;

            this->CreateSplitBcFromOverlapRegion( parentSplitBc, overlapRegion );

            int kkk = 1;
        }
    }
    else
    {
        int nChilds = parentSplitBc->child->size();
        for ( int iChild = 0; iChild < nChilds; ++ iChild )
        {
            this->CreateSplitBc( ( * parentSplitBc->child )[ iChild ] );
        }
    }
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

void SplitZone::CreateSplitBcFromOverlapRegion( SplitBc * parentSplitBc, RangeRegion & overlapRegion )
{
    SplitBc * splitBc = new SplitBc();
    this->AddSplitBc( splitBc );

    splitBc->SetRegion( overlapRegion );
    splitBc->bcType = parentSplitBc->bcType;
    splitBc->splitZone = this;
    splitBc->ChangeRegionToLocalCoordinate();
    splitBc->CreateInterfaceBcFromParentBc( parentSplitBc );
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

