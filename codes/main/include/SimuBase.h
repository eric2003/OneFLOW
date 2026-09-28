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
#include <memory>
#include <functional>
#include <unordered_map>
#include <vector>
#include <string>

BeginNameSpace(ONEFLOW)

class SimuBase
{
public:
    virtual ~SimuBase() = default;
    virtual void Run() = 0;
};

using TestCreator = std::function<std::unique_ptr<SimuBase>()>;

class TestRegistry
{
public:
    static TestRegistry& Instance();

    void Register(const std::string& caseName, TestCreator creator);
    std::unique_ptr<SimuBase> Create(const std::string& caseName);
    std::vector<std::string> GetAllRegisteredNames() const;

private:
    TestRegistry() = default;
    std::unordered_map<std::string, TestCreator> m_creators;
};

EndNameSpace

// Register test case, ClassName is local class in current translation-unit (anonymous namespace)
#define REGISTER_TEST_CASE(ClassName, CaseName)                          \
namespace {                                                             \
bool ClassName##_registered = [](){                                     \
    ONEFLOW::TestRegistry::Instance().Register(CaseName, [](){          \
        return std::make_unique<ClassName>();                           \
    });                                                                 \
    return true;                                                        \
}();                                                                    \
}

// Macro to generate local wrapper class for existing test type T
// T: original test class name inside ONEFLOW namespace
// case_id: string identifier used for environment variable / registry
// NOTE: invoke this macro at the END of xxxTest.cpp, AFTER EndNameSpace (outside ONEFLOW namespace)
#define WRAP_TEST_CLASS(T, case_id)                                       \
namespace {                                                                \
class T##Wrapper : public ONEFLOW::SimuBase                                \
{                                                                         \
public:                                                                   \
    void Run() override                                                   \
    {                                                                     \
        ONEFLOW::T obj;                                                   \
        obj.Run();                                                        \
    }                                                                     \
};                                                                        \
}                                                                         \
REGISTER_TEST_CASE(T##Wrapper, case_id)
