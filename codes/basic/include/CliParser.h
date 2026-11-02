/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/
#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ONEFLOW
{

class CliParser;

class ParsedCommandLine
{
public:
    bool Has( const std::string & name ) const
    {
        return values_.find( name ) != values_.end();
    }

    const std::string & GetOne( const std::string & name ) const
    {
        const std::vector<std::string> & values = GetAll( name );
        if ( values.size() != 1 )
        {
            throw std::logic_error(
                "Command line argument is not singular: " + name );
        }
        return values.front();
    }

    const std::vector<std::string> & GetAll(
        const std::string & name ) const
    {
        const auto iter = values_.find( name );
        if ( iter == values_.end() )
        {
            throw std::out_of_range(
                "Command line argument was not parsed: " + name );
        }
        return iter->second;
    }

private:
    friend class CliParser;
    std::map<std::string, std::vector<std::string>> values_;
};

class CliParser
{
public:
    static constexpr std::size_t Unlimited =
        std::numeric_limits<std::size_t>::max();

    explicit CliParser( std::string programName )
        : programName_( std::move( programName ) )
    {
    }

    void AddPositional(
        std::string name,
        std::string description,
        std::size_t minimumCount = 1,
        std::size_t maximumCount = 1 )
    {
        if ( name.empty() )
        {
            throw std::invalid_argument(
                "Command line positional name must not be empty" );
        }
        if ( minimumCount > maximumCount || maximumCount == 0 )
        {
            throw std::invalid_argument(
                "Invalid cardinality for command line positional: " + name );
        }
        if ( std::any_of(
            positionals_.begin(), positionals_.end(),
            [&]( const Positional & positional )
            {
                return positional.name == name;
            } ) )
        {
            throw std::invalid_argument(
                "Duplicate command line positional: " + name );
        }
        if ( ! positionals_.empty()
            && positionals_.back().maximumCount == Unlimited )
        {
            throw std::invalid_argument(
                "A repeatable positional must be the final argument" );
        }

        positionals_.push_back( {
            std::move( name ), std::move( description ),
            minimumCount, maximumCount } );
    }

    ParsedCommandLine Parse(
        const std::vector<std::string> & arguments ) const
    {
        ParsedCommandLine parsed;
        std::size_t argumentIndex = 0;

        for ( std::size_t specIndex = 0;
            specIndex < positionals_.size(); ++specIndex )
        {
            const Positional & positional = positionals_[ specIndex ];
            std::size_t requiredAfter = 0;
            for ( std::size_t later = specIndex + 1;
                later < positionals_.size(); ++later )
            {
                requiredAfter += positionals_[ later ].minimumCount;
            }

            const std::size_t available = arguments.size() - argumentIndex;
            if ( available < requiredAfter
                || available - requiredAfter < positional.minimumCount )
            {
                throw std::invalid_argument(
                    "Missing required command line argument: "
                    + positional.name + "\n" + Usage() );
            }

            const std::size_t availableHere = available - requiredAfter;
            const std::size_t count = std::min(
                positional.maximumCount, availableHere );
            if ( count > 0 )
            {
                std::vector<std::string> & values =
                    parsed.values_[ positional.name ];
                values.insert(
                    values.end(),
                    arguments.begin() + argumentIndex,
                    arguments.begin() + argumentIndex + count );
            }
            argumentIndex += count;
        }

        if ( argumentIndex < arguments.size() )
        {
            throw std::invalid_argument(
                "Unexpected command line argument: "
                + arguments[ argumentIndex ] + "\n" + Usage() );
        }

        return parsed;
    }

    std::string Usage() const
    {
        std::ostringstream usage;
        usage << "Usage: " << programName_;
        for ( const Positional & positional : positionals_ )
        {
            const bool optional = positional.minimumCount == 0;
            const std::string token = "<" + positional.name + ">";
            if ( optional )
            {
                usage << " [" << token << "]";
            }
            else
            {
                usage << " " << token;
            }
            if ( positional.maximumCount == Unlimited
                || positional.maximumCount > 1 )
            {
                usage << "...";
            }
        }

        if ( ! positionals_.empty() )
        {
            usage << "\nPositional arguments:";
            for ( const Positional & positional : positionals_ )
            {
                usage << "\n  " << positional.name;
                if ( ! positional.description.empty() )
                {
                    usage << "  " << positional.description;
                }
            }
        }
        return usage.str();
    }

private:
    struct Positional
    {
        std::string name;
        std::string description;
        std::size_t minimumCount;
        std::size_t maximumCount;
    };

    std::string programName_;
    std::vector<Positional> positionals_;
};

}