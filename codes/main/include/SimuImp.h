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
#include "SimuBase.h"
#include "SimuContext.h"
#include <memory>
#include <vector>
#include <string>

BeginNameSpace( ONEFLOW )

// Full simulation path. Owns a SimuContext (phase 2) and dispatches
// work through TaskRegistry (phase 1).
class SimuImp : public SimuBase
{
public:
    explicit SimuImp( std::vector<std::string>& args );
    ~SimuImp() override;

    void Run() override;

    // Exposed for tests that inject a pre-built context path later.
    SimuContext& Context() { return *ctx_; }
    const SimuContext& Context() const { return *ctx_; }

public:
    void PreProcess();
    void MainProcess();
    void PostProcess();

protected:
    void InitSimu();
    void RunSimu();

private:
    std::unique_ptr<SimuContext> ctx_;

public:
    // Kept for source compatibility with any code reading simu.args.
    // Prefer Context().Args().
    std::vector<std::string> args;
};

EndNameSpace
