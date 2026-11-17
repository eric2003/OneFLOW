/*---------------------------------------------------------------------------*\\
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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#pragma once

#include <stdexcept>

namespace ONEFLOW
{
namespace SimuProcessLifecycleDetail
{

enum class ProcessRuntimePhase
{
    NotStarted,
    Active,
    Finalized,
    MovedFrom
};

class ProcessRuntimeState
{
public:
    ProcessRuntimeState() = default;
    ProcessRuntimeState( const ProcessRuntimeState & ) = delete;
    ProcessRuntimeState & operator=( const ProcessRuntimeState & ) = delete;

    // Transfer the process-runtime lifecycle marker with its owning context.
    ProcessRuntimeState( ProcessRuntimeState && other ) noexcept
        : phase_( other.phase_ )
    {
        other.phase_ = ProcessRuntimePhase::MovedFrom;
    }

    // Replacing an existing state could silently abandon an active runtime.
    ProcessRuntimeState & operator=( ProcessRuntimeState && ) = delete;

    void Start()
    {
        if ( phase_ != ProcessRuntimePhase::NotStarted )
        {
            throw std::logic_error(
                "Process runtime can only be initialized once" );
        }
        phase_ = ProcessRuntimePhase::Active;
    }

    bool IsActive() const
    {
        return phase_ == ProcessRuntimePhase::Active;
    }

    bool IsFinalized() const
    {
        return phase_ == ProcessRuntimePhase::Finalized;
    }

    void MarkFinalized()
    {
        if ( phase_ == ProcessRuntimePhase::Active )
        {
            phase_ = ProcessRuntimePhase::Finalized;
        }
    }

private:
    ProcessRuntimePhase phase_ = ProcessRuntimePhase::NotStarted;
};

}
}