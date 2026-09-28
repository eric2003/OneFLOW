#include "Plot.h"
#include "CgnsBase.h"

Edge::Edge() : Edge(0, 0)
{
}

Edge::Edge(std::size_t p1, std::size_t p2)
{
    this->p1 = p1;
    this->p2 = p2;
}

Edge::~Edge()
{
    ;
}

bool Edge::operator < ( const Edge & rhs ) const
{
    std::size_t sp1 = std::min( p1, p2 );
    std::size_t sp2 = std::max( p1, p2 );

    std::size_t tp1 = std::min( rhs.p1, rhs.p2 );
    std::size_t tp2 = std::max( rhs.p1, rhs.p2 );

    if ( sp1 != tp1 )
    {
        return sp1 < tp1;
    }

    if ( sp2 != tp2 )
    {
        return sp2 < tp2;
    }
    return false;
}


Edges::Edges()
{
    ;
}

Edges::~Edges()
{
    ;
}

Nodes::Nodes()
{
    ;
}

Nodes::~Nodes()
{
    ;
}

Uns2D::Uns2D()
{
    ;
}

Uns2D::~Uns2D()
{
    ;
}

void Uns2D::CalcUnsFaceInfo( Section * section )
{
    std::vector<int> faceNodes;
    int num_edge = 0;
    int num_node = 0;
    for ( int iElem = 0; iElem < section->nElements; ++ iElem )
    {
        int st = section->conn_offsets[ iElem ];
        int ed = section->conn_offsets[ iElem + 1 ];
        int nNode = ed - st;
        faceNodes.resize( 0 );
        for ( int i = st; i < ed; ++ i )
        {
            int node = section->conn[ i ]-1;
            auto it = this->nodes.datamap.find( node );
            if ( it == this->nodes.datamap.end() )
            {
                int node_id = num_node;
                this->nodes.datamap.insert(std::pair<int, int>{node,node_id});
                this->nodes.index.push_back(node);
                faceNodes.push_back( node_id );
                num_node ++;
            }
            else
            {
                faceNodes.push_back( it->second );
            }
        }

        for ( int i = 0; i < nNode; ++ i )
        {
            int p1 = faceNodes[ i ];
            int p2 = faceNodes[ ( i + 1 ) % nNode ];
            Edge edge( p1, p2 );
            auto it = this->edges.data.find( edge );
            if ( it == this->edges.data.end() )
            {
                int edge_id = num_edge;
                this->edges.data.insert(std::pair<Edge, int>{edge,edge_id});
                this->edges.edgeList.push_back( edge );
                num_edge ++;

                this->left_elements.resize( num_edge );
                this->right_elements.resize( num_edge );
                this->left_elements[ edge_id ] = iElem;
                this->right_elements[ edge_id ] = -1;
            }
            else
            {
                this->right_elements[ it->second ] = iElem;
            }
        }
    }

    this->nNodes = this->nodes.index.size();
    this->nFaces = this->edges.GetNumberOfEdges();
    this->nCells = section->nElements;
    int totalNumFaceNodes = nFaces * 2;
}

void Uns2D::PlotSectionMesh( std::ostringstream & oss, Section * section )
{
    int totalNumFaceNodes = nFaces * 2;

    // output for Tecplot
    std::string zoneTitle = "ZONE T=\"" + std::string{ section->sectionName } + "\"";
    oss << zoneTitle << "\n";
    oss << "ZONETYPE = FEPolygon\n";
    oss << "DATAPACKING = BLOCK\n";
    oss << "Nodes    = " << this->nNodes << std::endl;
    oss << "Faces    = " << this->nFaces << std::endl;
    oss << "Elements = " << this->nCells << std::endl;
    //oss << "TotalNumFaceNodes = " << totalNumFaceNodes << std::endl;
    oss << "NumConnectedBoundaryFaces = 0\n";
    oss << "TotalNumBoundaryConnections = 0\n";
    this->DumpSectionCoor( oss, section );
}

void Uns2D::DumpSectionCoor( std::ostringstream & oss, Section * section )
{
    int nWords = 5;
    Zone * zone = section->zone;
    for ( int iCoord = 0; iCoord < zone->nCoords; ++ iCoord )
    {
        Coor * coor = zone->coors[iCoord];
        if ( coor->dataType == RealSingle )
        {
            float * field = static_cast<float *>( coor->data );
            ::DumpCoor( oss, field, this->nodes.index, nWords );
        }
        else
        {
            double * field = static_cast<double *>( coor->data );
            ::DumpCoor( oss, field, this->nodes.index, nWords );
        }
    }
}

void Uns2D::DumpFaceElementLink( std::ostringstream & oss, Section * section )
{
    oss << "# left elements\n";
    this->DumpFaceElementLink( oss, section, this->left_elements );
    oss << "# right elements\n";
    this->DumpFaceElementLink( oss, section, this->right_elements );
}

void Uns2D::DumpFaceElementLink( std::ostringstream & oss, Section * section, std::vector<int> &elements )
{
    int icount = 0;
    int nWords = 10;
    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        int eId = elements[ iFace ] + 1;
        oss << eId << " ";
        icount ++;
        if ( icount % nWords == 0 ) oss << std::endl;
    }
    if ( icount % nWords != 0 ) oss << std::endl;
}

void Uns2D::DumpFaceNodeLink( std::ostringstream & oss, Section * section )
{
    int icount = 0;
    int nWords = 10;

    oss << "# face nodes\n";

    if ( section->elementType == NGON_n )
    {
        int icount = 0;
        int nWords = 5;

        for ( int i = 0; i < this->edges.edgeList.size(); ++ i )
        {
            Edge &edge = this->edges.edgeList[ i ];
            oss <<  edge.p1 + 1 << " " << edge.p2 + 1 << " ";
            icount ++;
            if ( icount % nWords == 0 ) oss << std::endl;
        }
        if ( icount % nWords != 0 ) oss << std::endl;
    }
}

Plot::Plot()
{
    ;
}

Plot::~Plot()
{
    ;
}

void Plot::PlotZoneMesh( std::ostringstream & oss, Zone * zone )
{
    DumpTecplotHeader( oss );

    int nNodes = zone->nNodes;
    int nFaces = zone->nFaces;
    int nCells = zone->nCells;
    int totalNumFaceNodes = zone->totalNumFaceNodes;

    // output for Tecplot
    std::string zoneTitle = "ZONE T=\"" + std::string{ zone->zoneName } + "\"";
    oss << zoneTitle << "\n";
    oss << "ZONETYPE = FEPolyhedron\n";
    oss << "DATAPACKING = BLOCK\n";
    oss << "Nodes    = " << nNodes << std::endl;
    oss << "Faces    = " << nFaces << std::endl;
    oss << "Elements = " << nCells << std::endl;
    oss << "TotalNumFaceNodes = " << totalNumFaceNodes << std::endl;
    oss << "NumConnectedBoundaryFaces = 0\n";
    oss << "TotalNumBoundaryConnections = 0\n";
}

void DumpTecplotHeader( std::ostringstream & oss )
{
    std::vector< std::string > titles;

    titles.push_back( "TITLE=\"The Mesh\"" );
    titles.push_back( "VARIABLES=" );
    titles.push_back( "\"X\"" );
    titles.push_back( "\"Y\"" );
    titles.push_back( "\"Z\"" );

    for ( std::size_t i = 0; i < titles.size(); ++ i )
    {
        oss << titles[ i ] << std::endl;
    }
}


