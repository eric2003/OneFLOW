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
#include "NamespaceMacros.h"
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <functional>

BeginNameSpace(ONEFLOW)

/**
* @brief unified abstract interface for both full simulation and lightweight test cases
*/
class SimuBase
{
public:
    virtual ~SimuBase() = default;
    virtual void Run() = 0;
};

using TestCreator = std::function<std::unique_ptr<SimuBase>()>;

/**
* @brief global registry for lightweight quick-run test cases
*/
class TestRegistry
{
public:
    static TestRegistry& Instance();

    void Register(const std::string& name, TestCreator creator);
    std::unique_ptr<SimuBase> Create(const std::string& name);
    // ===== 新增：获取全部已注册测试名字列表 =====
    std::vector<std::string> GetAllRegisteredNames() const;

private:
    TestRegistry() = default;
    std::unordered_map<std::string, TestCreator> m_creators;
};

EndNameSpace

/**
* @brief macro for register test case, place in *.cpp, NOT in header
* @param ClassName test case class name derived from ONEFLOW::SimuBase
* @param CaseName unique string identifier
*/
#define REGISTER_TEST_CASE(ClassName, CaseName) \
namespace { \
bool ClassName##_registered = [](){ \
    ONEFLOW::TestRegistry::Instance().Register(CaseName, [](){ \
        return std::make_unique<ONEFLOW::ClassName>(); \
    }); \
    return true; \
}(); \
}
