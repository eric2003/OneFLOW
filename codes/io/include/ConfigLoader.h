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
// ConfigLoader.h
#pragma once
#include <string>
#include <vector>

namespace ONEFLOW {

    class TextFileParser;

    // Represents a single parsed parameter entry before committing to DataBase
    struct ParameterEntry {
        std::string name;
        int type; // Maps to HX_INT, HX_REAL, HX_STRING, etc.
        std::vector<std::string> values;
    };

    // Responsible for parsing script files and loading parameters
    class ConfigLoader {
    public:
        ConfigLoader() = default;
        ~ConfigLoader() = default;

        // Parse directly from a file path
        void ParseFile(const std::string& fileName);

        // Parse from an existing TextFileParser (Used for legacy interface integration)
        void ParseFromParser(TextFileParser& parser);

        // Commit all parsed entries to the global DataBase
        void CommitToDataBase() const;

    private:
        std::vector<ParameterEntry> entries_;

        void ParseScalarParameter(TextFileParser& parser, ParameterEntry& entry);
        void ParseArrayParameter(TextFileParser& parser, ParameterEntry& entry);
    };

} // namespace ONEFLOW
