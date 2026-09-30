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
#include "ConfigLoader.h"
#include "ParaFile.h"       // Reuse IsArrayParameter, GetParameterArraySize
#include "TextFileParser.h"
#include "Word.h"
#include "Fatal.h"
#include <utility>

namespace ONEFLOW {

    void ConfigLoader::ParseFile(const std::string& fileName) {
        TextFileParser parser;
        parser.OpenFile(fileName, std::ios_base::in);
        ParseFromParser(parser);
        parser.CloseFile();
    }

    void ConfigLoader::ParseFromParser(TextFileParser& parser) {
        document_.Clear();
        std::string keyWordSeparator = " =\r\n\t#$,;\"";
        parser.SetDefaultSeparator(keyWordSeparator);

        while (!parser.ReachTheEndOfFile()) {
            if (!parser.ReadNextMeaningfulLine()) break;

            std::string keyWord = parser.ReadNextWord();
            if (keyWord.empty()) continue;

            std::string currentLine = parser.GetCurrentLine();
            ParameterEntry entry;
            entry.typeName = keyWord;

            if (IsArrayParameter(currentLine)) {
                ParseArrayParameter(parser, entry);
            } else {
                ParseScalarParameter(parser, entry);
            }

            if (!entry.name.empty()) {
                document_.Add( std::move( entry ) );
            }
        }
    }

    void ConfigLoader::ParseScalarParameter(TextFileParser& parser, ParameterEntry& entry) {
        std::string separator = " =\r\n\t#$,;\"";
        entry.name = parser.ReadNextWord(separator);
        std::string value = parser.ReadNextWord(separator);
        entry.values.push_back(value);
    }

    void ConfigLoader::ParseArrayParameter(TextFileParser& parser, ParameterEntry& entry) {
        std::string commSeparator = "=\r\n\t#$,;\"";
        std::string arraySeparator = " =\r\n\t#$,;\"[]";

        std::string arrayInfo = parser.ReadNextWord(commSeparator);
        entry.name = Word::FindNextWord(arrayInfo, arraySeparator);
        std::string arraySizeName = Word::FindNextWord(arrayInfo, arraySeparator);

        // Reuse legacy logic: supports literal digits or variable names from DataBase
        int arraySize = GetParameterArraySize(arraySizeName);

        for (int i = 0; i < arraySize; ++i) {
            std::string val = parser.ReadNextWord(arraySeparator);
            // Legacy cross-line reading compatibility
            if (val.empty()) {
                parser.ReadNextNonEmptyLine();
                val = parser.ReadNextWord(arraySeparator);
                if (val.empty()) {
                    Fatal("error in parameter file: array element missing");
                }
            }
            entry.values.push_back(val);
        }
    }

} // namespace ONEFLOW

