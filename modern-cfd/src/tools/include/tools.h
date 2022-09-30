/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include <string>
#include <sstream>

//using uint = std::size_t;

template <class... Ts>
std::string add_string( Ts const&... args )
{
    std::ostringstream oss;
    ((oss << args), ... );
    return oss.str();
}

template <class... Ts>
void print_all(std::ostream& os, Ts const&... args) {
    ((os << args << "\n" ), ... );
}

inline std::string alignl(const std::size_t n, const std::string& x="")
{
	// converts x to string with spaces behind such that length is n if x is not longer than n
	std::string s = x;
	for ( std::size_t i = x.length(); i < n; ++ i )
	{
		s += " ";
	}
	return s;
}

inline std::string alignr(const std::size_t n, const std::string& x="")
{
	// converts x to string with spaces behind such that length is n if x is not longer than n
	std::string s = "";
	for ( std::size_t i = x.length(); i < n; ++ i )
	{
		s += " ";
	}
	s += x;
	return s;
}


