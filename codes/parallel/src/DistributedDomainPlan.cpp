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

#include "DistributedDomainPlan.h"

#include <set>
#include <stdexcept>
#include <utility>

namespace ONEFLOW
{

DistributedDomainPlan DistributedDomainPlan::Build(
    const std::vector< int > & ownerRanks,
    const std::vector< std::vector< int > > & zoneNeighbors,
    int rankCount )
{
    if ( rankCount <= 0 )
    {
        throw std::invalid_argument(
            "DistributedDomainPlan::Build: rank count must be positive" );
    }
    if ( ownerRanks.size() != zoneNeighbors.size() )
    {
        throw std::invalid_argument(
            "DistributedDomainPlan::Build: ownership and topology sizes differ" );
    }

    const int zoneCount = static_cast< int >( ownerRanks.size() );
    for ( int zoneId = 0; zoneId < zoneCount; ++ zoneId )
    {
        const int ownerRank = ownerRanks[ zoneId ];
        if ( ownerRank < 0 || ownerRank >= rankCount )
        {
            throw std::out_of_range(
                "DistributedDomainPlan::Build: zone owner rank is out of range" );
        }
    }

    std::set< std::pair< int, int > > zonePairs;
    for ( int sourceZone = 0; sourceZone < zoneCount; ++ sourceZone )
    {
        for ( const int targetZone : zoneNeighbors[ sourceZone ] )
        {
            if ( targetZone < 0 || targetZone >= zoneCount )
            {
                throw std::out_of_range(
                    "DistributedDomainPlan::Build: neighbor zone is out of range" );
            }
            // A zone cannot exchange a partition boundary with itself.
            // Ignore duplicate/self topology entries instead of scheduling a
            // second local action pass.
            if ( targetZone == sourceZone ) continue;

            const int first =
                sourceZone < targetZone ? sourceZone : targetZone;
            const int second =
                sourceZone < targetZone ? targetZone : sourceZone;
            zonePairs.emplace( first, second );
        }
    }

    DistributedDomainPlan plan;
    plan.ownerRanks_ = ownerRanks;
    for ( const auto & zonePair : zonePairs )
    {
        plan.boundaryExchanges_.push_back(
            { zonePair.first, zonePair.second } );
        if ( zonePair.first != zonePair.second )
        {
            plan.boundaryExchanges_.push_back(
                { zonePair.second, zonePair.first } );
        }
    }
    return plan;
}

int DistributedDomainPlan::ZoneCount() const
{
    return static_cast< int >( this->ownerRanks_.size() );
}

int DistributedDomainPlan::OwnerRank( int zoneId ) const
{
    if ( zoneId < 0 || zoneId >= this->ZoneCount() )
    {
        throw std::out_of_range(
            "DistributedDomainPlan::OwnerRank: zone id is out of range" );
    }
    return this->ownerRanks_[ zoneId ];
}

const std::vector< int > & DistributedDomainPlan::OwnerRanks() const
{
    return this->ownerRanks_;
}

const std::vector< BoundaryExchange > &
DistributedDomainPlan::BoundaryExchanges() const
{
    return this->boundaryExchanges_;
}

}
