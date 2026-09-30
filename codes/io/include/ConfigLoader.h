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
#include "ConfigDocument.h"
#include <functional>
#include <string>
#include <utility>

namespace ONEFLOW {

    class TextFileParser;

    // Reads the legacy C-like parameter syntax into a ConfigDocument.
    class ConfigLoader {
    public:
        using ArraySizeResolver = std::function<int(const std::string&)>;

        // The caller supplies runtime lookup for array sizes named by variables.
        explicit ConfigLoader(ArraySizeResolver arraySizeResolver = {})
            : arraySizeResolver_(std::move(arraySizeResolver)) {}
        ~ConfigLoader() = default;

        // Parse directly from a file path
        void ParseFile(const std::string& fileName);

        // Parse from an existing TextFileParser (Used for legacy interface integration)
        void ParseFromParser(TextFileParser& parser);

        const ConfigDocument& Document() const noexcept { return document_; }

    private:
        ConfigDocument document_;
        ArraySizeResolver arraySizeResolver_;

        void ParseScalarParameter(TextFileParser& parser, ParameterEntry& entry);
        void ParseArrayParameter(TextFileParser& parser, ParameterEntry& entry);
    };

} // namespace ONEFLOW
