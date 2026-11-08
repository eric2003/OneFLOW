/*
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.

    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
*/
#pragma once

#include <string>

namespace ONEFLOW
{

struct VersionCompatibilityAssessment
{
    std::string status;
    std::string basis;
    std::string detail;
};

// An empty string or "unknown" means that the corresponding version is not
// independent evidence.
VersionCompatibilityAssessment AssessVersionCompatibility(
    const std::string & headerVersion,
    const std::string & runtimeVersion );

}
