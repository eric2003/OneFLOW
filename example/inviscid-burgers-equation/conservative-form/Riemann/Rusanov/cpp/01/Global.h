#pragma once
#include "Vec1d.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class Face
{
public:
    int zone = 0;
    int i = 0;
    int j = 1;
    int k = 1;
public:
    bool operator < ( const Face & rhs ) const;
    bool operator == ( const Face & rhs ) const;
public:
    void Print();
};

class FacePair
{
public:
    Face left;
    Face right;
public:
    bool operator < ( const FacePair & rhs ) const;
    void AddPair( const Face &face1, const Face &face2);
public:
    void Print();
};

class Transform;

class InterfaceTopo
{
public:
    std::vector<std::vector<int>> linkmap;
public:
    void InitNeighborInfo();
    void SwapNeighborInfo();
    void SwapNeighborDonorfaces();
};


class Interface
{
public:
    int zoneid;
    std::vector<int> zoneList;
    std::vector<int> global_faceids;
    std::vector<int> mpi_global_faceids;
    std::vector<int> proc_global_faceids;
    std::vector<int> local_faceids;
    std::vector<int> ijk_ghosts;
    std::vector<int> ijk_donors;
    std::vector<double> data_recv;
    std::vector<double> data_send;
    std::unordered_map<int, int> global_local_face_map;
public:
    std::vector<int> neighbor_donor_zones;
    std::vector<std::vector<int>> neighbor_donorfaces;
    std::vector<std::vector<int>> sub_local_faceids;
    std::vector<int> send_to_zones;
    std::vector<std::vector<int>> donorfaces_for_send;
    std::vector<std::vector<int>> donorijk_for_send;
    std::vector<std::vector<double>> donordata_for_send;
public:
    void CalcInterface( Transform * transform, std::vector<int> & start, std::vector<int> & end, int donor_zoneid, int nghost );
    void SendGeom( int zone, std::vector<int> & donorfaces );
};

class Field;
class InterFaceZone;
class Zone;
class Grid;

enum class BasicScheme
{
    FTCS = 0,
    CENTER, //Central Space
    CN, //Crank¨CNicolson
    ICP, //Implicit Compact Pade (ICP) Scheme
    WENO, //Weighted Essentially Non-oscillatory
    CRWENO, //Compact Reconstruction WENO-5 Scheme
    LAX, //Lax-Friedrichs flux splitting
    UpWind, //Lax-Friedrichs flux splitting
    RungeKutta,
    RK1,
    RK2,
    RK3
};

enum class GoverningEquation
{
    Heat = 0,
    Burgers
};

template<typename T>
int to_int( const T & t )
{
    return static_cast<int>( t );
}

template<typename T>
BasicScheme to_BasicScheme( const T & t )
{
    return static_cast<BasicScheme>( t );
}

class Scheme
{
public:
    int inviscid;
    int viscous;
    int time_scheme;
    int reconstruction;
public:
    void read( const json & j );
    void set_inviscid_scheme( const std::string & name );
    void set_reconstruction_scheme( const std::string & name );
    void set_viscous_scheme( const std::string & name );
    void set_time_scheme( const std::string & name );
};

class Global
{
public:
    static std::vector<Grid *> grids;
    static std::vector<Field *> fields;
public:
    static std::vector<Zone *> zones;
    static std::vector<Interface *> interfaces;
public:
    static std::map<Face, int> faceMap;
    static std::map<FacePair, int> facePairMap;
    static std::vector<FacePair> facePairList;
    static std::vector<FacePair> mpi_facePairList;
    static std::vector<std::set<int>> donor_zone_sets;
    static std::vector<std::vector<int>> donor_zones;
    static InterfaceTopo interfaceTopo;
public:
    static Scheme scheme;
    static GoverningEquation governing_equation;
    static double total_time;
    static int iconservation;
    static int iviscous;
    static int nsave;
    static int idump_initial_field;
    static int ifinite_volume;
    static int nt;
    static int iter;
    static int cell_dim;
    static int phys_dim;
    static int nghost;
    static std::string file_string;
public:
    static void InsertFaceMap( const Face & face );
    static int InsertFacePairMap( const FacePair & facePair );
    static void AddFacePairList( std::vector<FacePair> & a, std::vector<FacePair> & b );
};
