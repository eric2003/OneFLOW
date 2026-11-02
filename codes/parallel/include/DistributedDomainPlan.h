/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#pragma once

#include <vector>

namespace ONEFLOW
{

struct BoundaryExchange
{
    int sourceZone = -1;
    int targetZone = -1;

    bool operator==( const BoundaryExchange & ) const = default;
};

class DistributedDomainPlan
{
public:
    static DistributedDomainPlan Build(
        const std::vector< int > & ownerRanks,
        const std::vector< std::vector< int > > & zoneNeighbors,
        int rankCount );

    int ZoneCount() const;
    int OwnerRank( int zoneId ) const;

    const std::vector< int > & OwnerRanks() const;
    const std::vector< BoundaryExchange > & BoundaryExchanges() const;

private:
    std::vector< int > ownerRanks_;
    std::vector< BoundaryExchange > boundaryExchanges_;
};

}
