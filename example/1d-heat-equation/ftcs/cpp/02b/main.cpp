#include "cgnslib.h"
#include <iostream>
#include <vector>
#include <numbers>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include <string>
#include <map>
#include <algorithm>

class Grid;
void ReadCgnsGrid( const std::string & filename );
void ReadCgnsGridBaseZone( const std::string & filename );
double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );
void DumpErrorDetails( std::vector<double> & u_error );
void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );

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
        //val_max = std::max( val_max, std::abs( u_error[ i ] ) );
        if ( val_max < std::abs( u_error[ i ] ) )
        {
            ipos = i;
            val_max = std::abs( u_error[ i ] );
        }
    }
    std::cout << " ipos = " << ipos << "\n";
    return val_max;
}

void DumpErrorDetails( std::vector<double> &u_error )
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

void DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<double> &ue, std::vector<double> &un, std::vector<double> &uerror )
{
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "x ue un uerror\n");
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::format_to(std::ostream_iterator<char>(file), "{:.16f} {:.16f} {:.16f} {:.16f}\n", x[i], ue[i], un[i], uerror[i] );
    }
    file.close();
}

class Face
{
public:
    std::vector<int> nodes;
    std::vector<int> sorted_nodes;
public:
    void AddFace( std::vector<int> & nodes )
    {
        this->nodes = nodes;
        this->sorted_nodes = nodes;
        std::sort( this->sorted_nodes.begin(), this->sorted_nodes.end() );
    }
    bool operator < ( const Face & rhs ) const
    {
        return this->sorted_nodes < rhs.sorted_nodes;
    } 
};

class FaceList
{
public:
    std::map<Face, int> face_map;
public:
    void AddFace( const Face & face )
    {
        std::map<Face, int>::iterator iter;
        iter = face_map.find( face );
        if ( iter == face_map.end() )
        {
            int id = face_map.size();
            face_map.insert( std::make_pair( face, id ) );
        }
    }
};

class BaseZone
{
public:
    //int baseId;
    std::string zone_name;
    bool operator < ( const BaseZone & rhs ) const
    {
        //if ( this->baseId != rhs.baseId )
        //{
        //    return this->baseId < rhs.baseId;
        //}

        return this->zone_name < rhs.zone_name;
    } 
};

class BaseZoneList
{
public:
    std::map<BaseZone, int> basezone_map;
public:
    void AddBaseZone( const BaseZone & baseZone )
    {
        std::map<BaseZone, int>::iterator iter;
        iter = basezone_map.find( baseZone );
        if ( iter == basezone_map.end() )
        {
            int id = basezone_map.size();
            basezone_map.insert( std::make_pair( baseZone, id ) );
        }
    }

    int FindBaseZone( const BaseZone & baseZone )
    {
        std::map<BaseZone, int>::iterator iter = basezone_map.find( baseZone );
        if ( iter != basezone_map.end() )
        {
            return iter->second;
        }
        return -1;
    }
};

class Point
{
public:
    typedef Point point_type;
    using value_type = double;
public:
    Point()
    {
        this->x = 0;
    }
    Point(const value_type &x)
    {
        this->x = x;
    }
public:
    value_type x;
public:
    bool operator < ( const Point & rhs ) const
    {
        value_type dx = x - rhs.x;
        value_type diff = 1.0e-10;

        if ( std::abs( dx ) > diff ) return x < rhs.x;

        return false;
    }
};

class PointFactory
{
public:
    using value_type = double;
    using point_type = Point;
    using iterator = std::map< point_type, int >::iterator;
protected:
    std::map< point_type, int > pointList;
    std::vector<point_type> pointArray;
public:
    int GetNPoints() { return pointArray.size(); }
    bool FindPoint( const point_type & point, iterator & iter )
    {
        iter = pointList.find( point );
        if ( iter == pointList.end() )
        {
            return false;
        }
        return true;
    }

    int AddPoint( const point_type & point )
    {
        PointFactory::iterator iter;
        if ( FindPoint( point, iter ) )
        {
            return iter->second;
        }
        else
        {
            int index = this->GetNPoints();

            pointList.insert( std::make_pair(point, index) );
            pointArray.push_back( point );

            return index;
        }
    }
};

class InterFaceZone
{
public:
    std::vector<int> face_ids;
    std::vector<int> left_zones;
    std::vector<int> right_zones;
    std::vector<int> left_cells;
    std::vector<int> right_cells;
    std::vector<double> left_u;
    std::vector<double> right_u;
};

class BC
{
public:
    std::vector<int> bctypes;
    std::vector<int> faceids;
};

class Grid;
class Field;

class Global
{
public:
    static std::vector<Grid *> grids;
    static std::vector<Field *> fields;
    static std::vector<BC *> bcs;
    static std::vector<InterFaceZone *> interfacezones;
    static PointFactory boundary_points;
    static FaceList boundary_faces;
    static FaceList inter_faces;
    static BaseZoneList zone_names;
    static int nt;
};

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;
std::vector<BC *> Global::bcs;
std::vector<InterFaceZone *> Global::interfacezones;
PointFactory Global::boundary_points;
FaceList Global::boundary_faces;
FaceList Global::inter_faces;
BaseZoneList Global::zone_names;
int Global::nt = -1;

class Grid
{
public:
    int zoneIndex;
    std::vector<double> x;
};

class Field
{
public:
    std::vector<double> u_e;
    std::vector<double> u, un;
    std::vector<double> error;
public:
    int ni;
    int nt;
    double dx, dt, t;
    double alpha, beta;
public:
    void Init( Grid * grid )
    {
        this->ni = grid->x.size();
        std::cout << "ni = " << ni << "\n";

        std::vector<double> & x = grid->x;
        this->dx = x[ 1 ] - x[ 0 ];
        this->dt = dx / 10.0;
        this->t = 1.0;
        this->nt = static_cast<int>( t / dt ) + 1;
        std::cout << "nt = " << nt << "\n";

        Global::nt = nt;

        this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
        this->beta = this->alpha * dt / ( dx * dx );
        std::cout << "alpha = " << std::setprecision( 15 ) << this->alpha << "\n";
        std::cout << "beta = " << std::setprecision( 15 ) << this->beta << "\n";

        int nghost = 2;
        int ni_total = ni + nghost;

        u_e.resize( ni_total );
        u.resize( ni_total );
        un.resize( ni_total );

        int ist = 0;
        int ied = ni - 1;

        // 0(bc) 1(ist) 2 3 ... ni

        //int ist = 0 + 1;
        //int ied = ni;

        for ( int i = ist; i <= ied; ++ i )
        {
            u_e[ i ] = - std::exp( -t ) * std::sin( std::numbers::pi * x[i] ); //theory solution
            un[ i ] = - std::sin( std::numbers::pi * x[ i ] ); //initial condition @ t=0
        }
        //un[ ist - 1 ] = 0.0;
        //un[ ied + 1 ] = 0.0;
        un[ 0 ] = 0.0;
        un[ ni - 1 ] = 0.0;
    }

    void Solve( Grid * grid )
    {
        Boundary( grid );
        for ( int i = 1; i < ni - 1; ++ i )
        {
            u[ i ] = un[ i ] + beta * ( un[ i + 1 ] - 2.0 * un[ i ] + un[ i - 1 ] );
        }
        this->update( un, u );
    }

    void Boundary( Grid * grid )
    {
        this->InterfaceBoundary( grid );
        this->PhysicalBoundary( grid );
    }

    void PhysicalBoundary( Grid * grid )
    {
        BC * bc = Global::bcs[ grid->zoneIndex ];

        int nBFace = bc->bctypes.size();
        for ( int i = 0; i < nBFace; ++ i )
        {
            int bctype = bc->bctypes[ i ];
            int faceId = bc->faceids[ i ];
            if ( bctype == BCInflow )
            {
                u[ faceId ] = 0.0;
            }
            else if ( bctype == BCOutflow )
            {
                u[ faceId ] = 0.0;
            }
        }
    }

    void InterfaceBoundary( Grid * grid )
    {
        InterFaceZone * interfacezone = Global::interfacezones[ grid->zoneIndex ];
        int nInterfaces = interfacezone->left_u.size();
        for ( int i = 0; i < nInterfaces; ++ i )
        {
            double ul = interfacezone->left_u[ i ];
            double ur = interfacezone->right_u[ i ];
            double um = 0.5 * ( ul + ur );

            int faceId = interfacezone->face_ids[ i ];

            u[ faceId ] = um;
            int kkk = 1;
        }
    }

    void PostProcess( Grid * grid )
    {
        std::vector<double> & x = grid->x;
        //compute L2 norm of the error
        std::vector<double> u_error( ni );
        for ( int i = 0; i < ni; ++ i )
        {
            u_error[ i ] = un[ i ] - u_e[ i ];
        }
    }

    void AddData( Grid * grid, std::vector<double> &global_x, std::vector<double> &global_ue, std::vector<double> &global_un )
    {
        int ni = grid->x.size();
        for ( int i = 0; i < ni; ++ i )
        {
            global_x.push_back( grid->x[ i ] );
            global_ue.push_back( this->u_e[ i ] );
            global_un.push_back( this->un[ i ] );
        }
    }

    void update( std::vector<double> &un, std::vector<double> &u )
    {
        for ( int i = 0; i < u.size(); ++ i )
        {
            un[ i ] = u[ i ];
        }
    }
};

class Solver
{
public:
    void Run()
    {
        this->ReadGrid();
        this->InitFields();
        this->SolveMultiZones();
        this->PostProcess();
    }

    void ReadGrid()
    {
        ReadCgnsGridBaseZone( "../heat1d2blocks.cgns" );
        int kkk = 1;
        ReadCgnsGrid( "../heat1d2blocks.cgns" );
        int nZones = Global::grids.size();
        PointFactory ptfactory;
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            grid->zoneIndex = iZone;
            int ni = grid->x.size();
            for ( int i = 0; i < ni; ++ i )
            {
                int pid = ptfactory.AddPoint( Point(grid->x[i]) );
                //face.nodes.push_back( pid );
            }
            int kkk = 1;
        }
    }

    void InitFields()
    {
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            Field * field = new Field();
            Global::fields.push_back( field );
            field->Init( grid );
        }
    }

    void SolveMultiZones()
    {
        for ( int it = 0; it < Global::nt; ++ it )
        {
            ExchangeInterfaceValue();
            int nZones = Global::grids.size();
            for ( int iZone = 0; iZone < nZones; ++ iZone )
            {
                Grid * grid = Global::grids[ iZone ];
                Field * field = Global::fields[ iZone ];
                field->Solve( grid );
            }
        }
    }

    void ExchangeInterfaceValue()
    {
        int nZones = Global::grids.size();
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            InterFaceZone * interfacezone = Global::interfacezones[ iZone ];
            int nInterFace = interfacezone->left_zones.size();
            interfacezone->left_u.resize( nInterFace );
            interfacezone->right_u.resize( nInterFace );
        }

        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Field * field = Global::fields[ iZone ];

            InterFaceZone * interfacezone = Global::interfacezones[ iZone ];
            int nInterFace = interfacezone->left_zones.size();

            for ( int i = 0; i < nInterFace; ++ i )
            {
                int left_cell = interfacezone->left_cells[ i ];
                int donor_zoneid = interfacezone->right_zones[ i ];
                int donor_cell = interfacezone->right_cells[ i ];
                Field * donor_field = Global::fields[ donor_zoneid ];

                interfacezone->left_u[ i ] = field->u[ left_cell ];
                interfacezone->right_u[ i ] = donor_field->u[ donor_cell ];
            }
        }
    }

    void PostProcess()
    {
        int nZones = Global::grids.size();
        std::vector<double> u_e;
        std::vector<double> un;
        std::vector<double> x;

        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            Grid * grid = Global::grids[ iZone ];
            Field * field = Global::fields[ iZone ];
            field->AddData( grid, x, u_e, un );
        }

        std::cout << " x.size() = " << x.size() << "\n";

        this->PostProcess( x, u_e, un );
    }

    void PostProcess( std::vector<double> &x, std::vector<double> &u_e, std::vector<double> &un )
    {
        int ni = x.size();
        //compute L2 norm of the error
        std::vector<double> u_error( ni );
        for ( int i = 0; i < ni; ++ i )
        {
            u_error[ i ] = un[ i ] - u_e[ i ];
        }

        ::DumpErrorDetails( u_error );

        std::string csvname = "field_final.csv";
        ::DumpCsvFile( csvname, x, u_e, un, u_error );
    }

    void PrintField( std::vector<double> &f )
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
};

void ReadCgnsGridBaseZone( const std::string & filename )
{
    int fileId = -1;
    cg_open( filename.c_str(), CG_MODE_READ, &fileId);
    std::cout << "fileId = " << fileId << "\n";
    int nbases = -1;
    cg_nbases( fileId, &nbases );
    std::cout << "nbases = " << nbases << "\n";

    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";
        int nzones = -1;
        cg_nzones( fileId, baseId, &nzones );
        std::cout << "nzones = " << nzones << "\n";
        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            int index_dim = -1;
            cg_index_dim( fileId, baseId, zoneId, &index_dim );
            std::cout << "index_dim = " << index_dim << "\n";

            std::vector<cgsize_t> isize( index_dim * 3 );

            char zonename[33];
            cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            for ( int i = 0; i < isize.size(); ++ i )
            {
                std::cout << "i = " << i << " isize = " << isize[i] << "\n";
            }

            BaseZone baseZone;
            //baseZone.baseId = baseId;
            baseZone.zone_name = zonename;

            Global::zone_names.AddBaseZone( baseZone );
        }
    }
    cg_close( fileId );
}

void ReadCgnsGrid( const std::string & filename )
{
    int fileId = -1;
    cg_open( filename.c_str(), CG_MODE_READ, &fileId);
    std::cout << "fileId = " << fileId << "\n";
    int nbases = -1;
    cg_nbases( fileId, &nbases );
    std::cout << "nbases = " << nbases << "\n";

    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";
        int nzones = -1;
        cg_nzones(fileId, baseId, &nzones);
        std::cout << "nzones = " << nzones << "\n";
        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            int index_dim = -1;
            cg_index_dim( fileId, baseId, zoneId, &index_dim );
            std::cout << "index_dim = " << index_dim << "\n";

            std::vector<cgsize_t> isize( index_dim * 3 );

            char zonename[33];
            cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            for ( int i = 0; i < isize.size(); ++ i )
            {
                std::cout << "i = " << i << " isize = " << isize[i] << "\n";
            }

            std::vector<cgsize_t> irmin(index_dim);
            std::vector<cgsize_t> irmax(index_dim);
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

            ZoneType_t zoneType;
            cg_zone_type( fileId, baseId, zoneId, &zoneType );
            std::cout << "zoneType = " << zoneType << " ZoneTypeName = " << ZoneTypeName[ zoneType ] << "\n";
            int ncoords = -1;
            cg_ncoords( fileId, baseId, zoneId, &ncoords );
            std::cout << "ncoords = " << ncoords << "\n";

            Grid * grid = new Grid();
            Global::grids.push_back( grid );

            BC * bc = new BC();
            Global::bcs.push_back( bc );

            InterFaceZone * interfacezone = new InterFaceZone();
            Global::interfacezones.push_back( interfacezone );

            BaseZone baseZone;
            //baseZone.baseId = baseId;
            baseZone.zone_name = zonename;

            int gZoneId = Global::zone_names.FindBaseZone( baseZone );

            for ( int icoord = 0; icoord < ncoords; ++ icoord )
            {
                int coorId = icoord + 1;
                DataType_t dataType;
                char coordname[33];
                cg_coord_info( fileId, baseId, zoneId, coorId, &dataType, coordname );
                std::cout << "coordname = " << coordname << "\n";
                std::cout << "dataType = " << dataType << " DataTypeName = " << DataTypeName[ dataType ] << "\n";
                //std::vector<double> coord( nNodes );
                std::vector<char> coord( nNodes * sizeof(double) );
                grid->x.resize( nNodes );

                cg_coord_read( fileId, baseId, zoneId, coordname, dataType, irmin.data(), irmax.data(), coord.data() );
                double * xd = reinterpret_cast<double *>( const_cast<char *>( coord.data() ) );
                for ( int i = 0; i < nNodes; ++ i )
                {
                    //std::cout << coord[i] << " ";
                    std::cout << xd[i] << " ";
                    grid->x[ i ] = xd[ i ];
                    if ( ( i + 1 ) % 5 == 0 ) std::cout << "\n";
                }
                std::cout << "\n";
            }

            int nbocos = -1;

            cg_nbocos( fileId, baseId, zoneId, &nbocos );
            std::cout << "nbocos = " << nbocos << "\n";
            for ( int iboco = 0; iboco < nbocos; ++ iboco )
            {
                int bccoId = iboco + 1;
                GridLocation_t location;
                cg_boco_gridlocation_read( fileId, baseId, zoneId, bccoId, &location );
                std::cout << "iboco = " << iboco <<  " location = " << location << " GridLocationName = " << GridLocationName[location] << "\n";

                char boconame[ 33 ];
                BCType_t bocotype;
                PointSetType_t ptset_type;
                cgsize_t npnts = 0;
                std::vector<int> normalIndex(index_dim,-1);
                cgsize_t normalListSize = 0;
                DataType_t normalDataType;
                int ndataset = -1;

                cg_boco_info( fileId, baseId, zoneId, bccoId, boconame, &bocotype, &ptset_type,
                    &npnts, normalIndex.data(), &normalListSize, &normalDataType, &ndataset );
                std::cout << "boconame = " << boconame << " bocotype = " << bocotype << " BCTypeName = " << BCTypeName[ bocotype ] << "\n";
                std::cout << "ptset_type = " << ptset_type <<  " PointSetTypeName = " << PointSetTypeName[ptset_type] << "\n";
                std::cout << "npnts = " << npnts << "\n";
                std::cout << "normalIndex = ";
                for ( int i = 0; i < index_dim; ++ i )
                {
                    std::cout << normalIndex[ i ] << " ";
                }
                std::cout << "\n";
                std::cout << "normalListSize = " << normalListSize << "\n";
                std::cout << "normalDataType = " << normalDataType << " DataTypeName = " << DataTypeName[ normalDataType ] << "\n";
                std::cout << "ndataset = " << ndataset << "\n";

                std::vector<char> normalList( nNodes * iphysdim * sizeof( double ) );

                std::vector<cgsize_t> pnts( npnts * index_dim );
                cg_boco_read( fileId, baseId, zoneId, bccoId, pnts.data(), normalList.data() );
                std::cout << "pnts = ";
                for ( int i = 0; i < pnts.size(); ++ i )
                {
                    std::cout << pnts[ i ] << " ";
                }
                std::cout << "\n";
                for ( int i = 0; i < pnts.size(); ++ i )
                {
                    Global::boundary_points.AddPoint( grid->x[ pnts[ i ] - 1 ] );
                }

                if ( index_dim == 1 )
                {
                    int p1 = pnts[ 0 ];
                    double x1 = grid->x[ p1 - 1 ];
                    int id = Global::boundary_points.AddPoint( x1 );
                    Face face;
                    std::vector<int> nodes;
                    nodes.push_back( id );
                    face.AddFace( nodes );
                    Global::boundary_faces.AddFace( face );

                    bc->bctypes.push_back( bocotype );
                    bc->faceids.push_back( p1 - 1 );
                }
                
                double * normal_d = reinterpret_cast<double *>( const_cast<char *>( normalList.data() ) );

                //std::cout << "normalList = ";
                //for ( int i = 0; i < nNodes*iphysdim; ++ i )
                //{
                //    std::cout << normal_d[ i ] << " ";
                //}
                //std::cout << "\n";
            }
            int n1to1 = -1;
            cg_n1to1( fileId, baseId, zoneId, &n1to1 );
            std::cout << "n1to1 = " << n1to1 << "\n";
            for ( int i1to1 = 0; i1to1 < n1to1; ++ i1to1 )
            {
                int i1to1Id = i1to1 + 1;
                char connectname[ 33 ];
                char donorname[ 33 ];
                cgsize_t npnts = 2;
                std::vector<cgsize_t> range( npnts * index_dim );
                std::vector<cgsize_t> donor_range( npnts * index_dim );
                std::vector<int> transform( index_dim );
                cg_1to1_read( fileId, baseId, zoneId, i1to1Id, connectname, donorname, range.data(), donor_range.data(), transform.data() );
                std::cout << "connectname = " << connectname << "\n";
                std::cout << "donorname = " << donorname << "\n";
                for ( int i = 0; i < range.size(); ++ i )
                {
                    Global::boundary_points.AddPoint( grid->x[ range[ i ] - 1 ] );
                }

                BaseZone baseZone;
                baseZone.zone_name = donorname;

                int gDonorZoneId = Global::zone_names.FindBaseZone( baseZone );

                if ( index_dim == 1 )
                {
                    int p1 = range[ 0 ];
                    double x1 = grid->x[ p1 - 1 ];
                    int id = Global::boundary_points.AddPoint( x1 );
                    Face face;
                    std::vector<int> nodes;
                    nodes.push_back( id );
                    face.AddFace( nodes );
                    Global::boundary_faces.AddFace( face );
                    Global::inter_faces.AddFace( face );

                    bc->bctypes.push_back( BCTypeUserDefined );
                    bc->faceids.push_back( p1 - 1 );

                    interfacezone->face_ids.push_back( p1 );

                    interfacezone->left_zones.push_back( gZoneId );
                    interfacezone->left_cells.push_back( p1 - 1 );

                    int q = donor_range[ 0 ];
                    interfacezone->right_zones.push_back( gDonorZoneId );
                    interfacezone->right_cells.push_back( q - 1 );
                }
            }
            int kkk = 1;
        }
    }
    int kkk = 1;
    cg_close(fileId);
}

int main(int argc, char **argv)
{
    Solver solver;
    solver.Run();
    return 0;
}
