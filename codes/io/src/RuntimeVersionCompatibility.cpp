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

#include "RuntimeVersionCompatibility.h"

namespace ONEFLOW
{
namespace
{

bool HasVersionEvidence( const std::string & version )
{
    return !version.empty() && version != "unknown";
}

}

VersionCompatibilityAssessment AssessVersionCompatibility(
    const std::string & headerVersion,
    const std::string & runtimeVersion )
{
    const bool hasHeader = HasVersionEvidence( headerVersion );
    const bool hasRuntime = HasVersionEvidence( runtimeVersion );

    if ( hasHeader && hasRuntime )
    {
        const bool matches = headerVersion == runtimeVersion;
        return {
            matches ? "verified_match" : "mismatch",
            "header_and_runtime_versions",
            matches
                ? "Header and runtime versions were independently obtained and match."
                : "Header and runtime versions were independently obtained and differ."
        };
    }

    if ( hasHeader )
    {
        return {
            "assumed_match",
            "header_version",
            "Runtime version was not independently obtained; the header version is used as the default assumption. This does not verify ABI compatibility."
        };
    }

    return {
        "unknown",
        "insufficient_version_evidence",
        "There is insufficient independent version information to assess header/runtime version consistency."
    };
}

}
