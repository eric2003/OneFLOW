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
// FaceLookup.h
#pragma once

#include "HXKey.h"
#include "HXVector.h"
#include "Constant.h"
#include <set>
#include <map>
#include <vector>
#include <cstddef>

BeginNameSpace(ONEFLOW)

template<typename T = int>
class HXLookup
{
public:
    using key_type       = HXKey<T>;
    using container_type = HXVector<T>;   // or std::vector<T>

    //---------------------- Core interface ----------------------
    // Find or insert, return unique id (starts from 0)
    std::pair<int, bool> FindOrAdd(const container_type& nodes)
    {
        return FindOrAddImpl(key_type(nodes));
    }

    std::pair<int, bool> FindOrAdd(const std::vector<T>& nodes)
    {
        return FindOrAddImpl(key_type(nodes));
    }

    // Only find, return INVALID_ID if not exist
    int Find(const container_type& nodes) const
    {
        return FindImpl(key_type(nodes));
    }

    int Find(const std::vector<T>& nodes) const
    {
        return FindImpl(key_type(nodes));
    }

    //---------------------- Helper ----------------------
    std::size_t Size()  const { return map_.size(); }
    void        Clear()       { map_.clear(); }

private:
    std::map<key_type, int> map_;

    std::pair<int, bool> FindOrAddImpl(const key_type& key)
    {
        auto it = map_.find(key);
        if (it != map_.end())
            return {it->second, false};          // already existed

        int newId = static_cast<int>(map_.size());
        map_.emplace(std::move(key), newId);
        return {newId, true};                    // newly inserted
    }

    int FindImpl(const key_type& key) const
    {
        auto it = map_.find(key);
        return (it == map_.end()) ? INVALID_INDEX : it->second;
    }
};

// Convenient aliases
using IntKeyLookup  = HXLookup<int>;
using IntFaceLookup = HXLookup<int>;   // keep old name if needed
using IntCellLookup = HXLookup<int>;
using IntLineLookup = HXLookup<int>;

EndNameSpace