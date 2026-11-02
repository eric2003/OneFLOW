#include <iostream>
#include <string>
#include <vector>
#include <sstream>      // std::ostringstream
#include <iomanip>      //std::setw
#include <set>
#include <map>

const int ONE_D = 1;
const int TWO_D = 2;
const int THREE_D = 3;

class Zone;
class Section;

class Edge
{
public:
    Edge();
    Edge(std::size_t p1, std::size_t p2);
    ~Edge();
    bool operator < ( const Edge & rhs ) const;
public:
    std::size_t p1, p2;
};

class Edges
{
public:
    Edges();
    ~Edges();
public:
    std::size_t GetNumberOfEdges() { return this->data.size(); }
public:
    std::map<Edge,int> data;
    std::vector<Edge> edgeList;
};

class Nodes
{
public:
    Nodes();
    ~Nodes();
public:
    std::map<int,int> datamap;
    std::vector<int> index;
};


class Section;

class Uns2D
{
public:
    Uns2D();
    ~Uns2D();
public:
    Nodes nodes;
    Edges edges;
    std::size_t nNodes;
    std::size_t nFaces;
    std::size_t nCells;

    std::vector<int> left_elements;
    std::vector<int> right_elements;
public:
    void CalcUnsFaceInfo( Section * section );
    void PlotSectionMesh( std::ostringstream & oss, Section * section );
    void DumpSectionCoor( std::ostringstream & oss, Section * section );
    void DumpFaceElementLink( std::ostringstream & oss, Section * section );
    void DumpFaceElementLink( std::ostringstream & oss, Section * section, std::vector<int> & elements );
    void DumpFaceNodeLink( std::ostringstream & oss, Section * section );
};

class Plot
{
public:
    Plot();
    ~Plot();
public:
    void PlotZoneMesh( std::ostringstream & oss, Zone * zone );
};

template< typename T >
void DumpCoor( std::ostringstream & oss, T *field, int nNodes, int nWords )
{
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        oss <<  std::setw(16) << std::scientific << field[ iNode ] << " ";
        if ( ( iNode + 1 ) % nWords == 0 ) oss << std::endl;
    }
    if ( nNodes % nWords != 0 ) oss << std::endl;
}

template< typename T >
void DumpCoor( std::ostringstream & oss, T *field, std::vector<int> &nodes, int nWords )
{
    int icount = 0;
    for ( int i = 0; i < nodes.size(); ++ i )
    {
        int node_id = nodes[ i ];
        oss <<  std::setw(16) << std::scientific << field[ node_id ] << " ";
        icount ++;
        if ( icount % nWords == 0 ) oss << std::endl;
    }
    if ( icount % nWords != 0 ) oss << std::endl;
}

void DumpTecplotHeader( std::ostringstream & oss );



