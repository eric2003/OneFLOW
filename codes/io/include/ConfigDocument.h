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

#include <string>
#include <utility>
#include <vector>

namespace ONEFLOW {

// Format-neutral representation of one named parameter before it is applied
// to a runtime-specific store such as the legacy DataBase.
struct ParameterEntry {
    std::string name;
    std::string typeName;
    std::vector<std::string> values;
};

// Ordered parameter document shared by configuration format adapters.
class ConfigDocument {
public:
    void Add( ParameterEntry entry )
    {
        entries_.push_back( std::move( entry ) );
    }

    void Clear() noexcept { entries_.clear(); }
    bool Empty() const noexcept { return entries_.empty(); }
    const std::vector<ParameterEntry>& Entries() const noexcept { return entries_; }

private:
    std::vector<ParameterEntry> entries_;
};

} // namespace ONEFLOW