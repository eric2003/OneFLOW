#pragma once
#include <vector>

class ZoneState
{
public:
    static int nZones;
    static int zone_id;
    static std::vector<int> pids;
public:
    static bool IsValid( int zoneid );
    static int GetProcID( int zoneid );
};

class LocalZone
{
public:
    static int nZones;
    static std::vector<int> global_zoneids;
};