#pragma once
#include <vector>
#include <map>
#include <algorithm>
#include <string>

class Grid;

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
    Point( const value_type & x )
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
    bool FindPoint( const point_type & point )
    {
        iterator iter = pointList.find( point );
        if ( iter == pointList.end() )
        {
            return false;
        }
        return true;
    }
private:
    bool FindPoint( const point_type & point, iterator & iter )
    {
        iter = pointList.find( point );
        if ( iter == pointList.end() )
        {
            return false;
        }
        return true;
    }
public:
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

class BC
{
public:
    std::vector<int> bctypes;
    std::vector<int> faceids;
};

class Grid
{
public:
    int zoneIndex;
    std::vector<double> x;
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


class Field;
class BC;
class InterFaceZone;

class Zone;

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
public:
    static std::vector<Zone *> zones;
    static BaseZoneList zone_names;
public:
    static int nt;
    static int cell_dim;
    static int phys_dim;
};

