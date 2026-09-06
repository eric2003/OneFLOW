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
#include "SimuBase.h"

BeginNameSpace(ONEFLOW)

TestRegistry& TestRegistry::Instance()
{
    static TestRegistry inst;
    return inst;
}

void TestRegistry::Register(const std::string& name, TestCreator creator)
{
    m_creators[name] = std::move(creator);
}

std::unique_ptr<SimuBase> TestRegistry::Create(const std::string& name)
{
    auto it = m_creators.find(name);
    if (it == m_creators.end())
        return nullptr;
    return it->second();
}

// ===== 实现列出所有注册名称 =====
std::vector<std::string> TestRegistry::GetAllRegisteredNames() const
{
    std::vector<std::string> names;
    names.reserve(m_creators.size());
    for(const auto& pair : m_creators)
    {
        names.push_back(pair.first);
    }
    return names;
}

EndNameSpace
