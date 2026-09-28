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
#include <set>
#include <map>
#include <vector>
#include <cstddef>

BeginNameSpace(ONEFLOW)

enum class LookupBackend { Set, Map /*, Hash*/ };

template<typename T = int>
class FaceLookup
{
public:
    using key_type       = HXKey<T>;
    using container_type = HXVector<T>;

    explicit FaceLookup(LookupBackend backend = LookupBackend::Set)
        : backend_(backend)
    {}

    //---------------------- 核心接口 ----------------------
    // 查找或插入，返回 id。支持 HXVector 和 std::vector
    int FindOrAdd(const container_type& nodes)
    {
        key_type key(nodes);          // 自动 Normalize
        return (backend_ == LookupBackend::Set)
                   ? FindOrAdd_Set(key)
                   : FindOrAdd_Map(key);
    }

    int FindOrAdd(const std::vector<T>& nodes)
    {
        key_type key(nodes);
        return (backend_ == LookupBackend::Set)
                   ? FindOrAdd_Set(key)
                   : FindOrAdd_Map(key);
    }

    // 只查找，不插入
    int Find(const container_type& nodes) const
    {
        key_type key(nodes);
        return (backend_ == LookupBackend::Set)
                   ? Find_Set(key)
                   : Find_Map(key);
    }

    int Find(const std::vector<T>& nodes) const
    {
        key_type key(nodes);
        return (backend_ == LookupBackend::Set)
                   ? Find_Set(key)
                   : Find_Map(key);
    }

    //---------------------- 辅助接口 ----------------------
    std::size_t Size() const
    {
        return (backend_ == LookupBackend::Set) ? set_.size() : map_.size();
    }

    void Clear()
    {
        set_.clear();
        map_.clear();
    }

    void SetBackend(LookupBackend b) { backend_ = b; }
    LookupBackend GetBackend() const { return backend_; }

private:
    LookupBackend backend_;

    // -------------------- Set 后端 --------------------
    struct KeyWithId
    {
        key_type key;
        int      id;

        bool operator<(const KeyWithId& rhs) const
        {
            return key < rhs.key;
        }
    };

    std::set<KeyWithId> set_;

    int FindOrAdd_Set(const key_type& key)
    {
        KeyWithId tmp{key, -1};
        auto it = set_.find(tmp);
        if (it != set_.end())
            return it->id;

        int newId = static_cast<int>(set_.size());
        set_.insert(KeyWithId{key, newId});
        return newId;
    }

    int Find_Set(const key_type& key) const
    {
        KeyWithId tmp{key, -1};
        auto it = set_.find(tmp);
        return (it == set_.end()) ? -1 : it->id;   // -1 建议后续换成项目的 INVALID_INDEX
    }

    // -------------------- Map 后端 --------------------
    std::map<key_type, int> map_;

    int FindOrAdd_Map(const key_type& key)
    {
        auto it = map_.find(key);
        if (it != map_.end())
            return it->second;

        int newId = static_cast<int>(map_.size());
        map_.emplace(key, newId);
        return newId;
    }

    int Find_Map(const key_type& key) const
    {
        auto it = map_.find(key);
        return (it == map_.end()) ? -1 : it->second;
    }
};

// 方便使用的别名（与项目风格一致）
using IntFaceLookup = FaceLookup<int>;

EndNameSpace