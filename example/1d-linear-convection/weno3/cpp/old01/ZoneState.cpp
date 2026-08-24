#include "ZoneState.h"
#include "Parallel.h"

int ZoneState::nZones = 0;
std::vector<int> ZoneState::pids;
std::vector<int> ZoneState::g2lzoneids;

bool ZoneState::IsValid( int zoneid )
{
    return ZoneState::pids[ zoneid ] == Parallel::pid;
}

int ZoneState::GetProcID( int zoneid )
{
    return ZoneState::pids[ zoneid ];
}

int LocalZone::nZones = 0;
std::vector<int> LocalZone::global_zoneids;

