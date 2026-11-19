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

#include "HXTypeBasic.h"

#include <cmath>
#include <cstddef>
#include <vector>

BeginNameSpace( ONEFLOW )

// One face reported by the device-versus-host viscous flux comparison.
struct ViscousFluxDifference
{
    int face = -1;
    // Equation index, or -1 when the face is reported for a non-finite value
    // rather than for exceeding the relative tolerance.
    int equation = -1;
    Real relative = 0.0;
    bool nonFinite = false;
};

// Result of comparing the device viscous face flux against the host reference
// on the same state.
//
// The comparison is deliberately explicit about non-finite values instead of
// relying on the relative-difference test: every comparison against NaN is
// false, so a NaN face silently satisfies "rel <= tolerance".  That is how a
// diverging device arm can report zero differing faces.  Non-finite values are
// counted separately and always fail the check.
struct ViscousFluxSelfCheckSummary
{
    // False when the caller passed a shape or pointer the comparison cannot
    // interpret; the summary then carries no measurement.
    bool inputValid = false;

    // Faces where at least one value is non-finite, per side.  A face can be
    // counted on both sides.
    long nonFiniteDeviceFaces = 0;
    long nonFiniteHostFaces = 0;

    // Faces with a non-zero scale that were compared numerically.
    long comparedFaces = 0;

    // Faces that either contain a non-finite value or exceed the tolerance.
    long differingFaces = 0;

    Real worstRelative = 0.0;
    int worstFace = -1;
    int worstEquation = -1;
    int firstFace = -1;
    int firstEquation = -1;

    // Caller-requested leading differences, in face order.
    std::vector< ViscousFluxDifference > differences;

    bool Finite() const
    {
        return nonFiniteDeviceFaces == 0 && nonFiniteHostFaces == 0;
    }

    bool Passed() const
    {
        return inputValid && Finite() && differingFaces == 0;
    }
};

// True when the value can take part in a relative comparison.
inline bool IsViscousFluxValueFinite( Real value )
{
    return std::isfinite( static_cast< double >( value ) ) != 0;
}

// Compare two equation-major flux sets on identical inputs.
//
//   device[eq][face], host[eq][face]   caller-owned planes, both nEquations
//   nEquations                          must be positive
//   nFaces                              number of faces per plane
//   tolerance                           relative tolerance, must be finite
//                                       and non-negative
//   maxDifferences                      how many leading differing faces to
//                                       record in summary.differences
//
// The scale for a face is the largest |value| over both sides and all
// equations, the same normalisation the previous in-place check used.  Faces
// where that scale is zero carry no information and are skipped.
inline ViscousFluxSelfCheckSummary CompareViscousFaceFlux(
    const Real * const * device,
    const Real * const * host,
    int nEquations,
    int nFaces,
    Real tolerance,
    std::size_t maxDifferences = 4 )
{
    ViscousFluxSelfCheckSummary summary;

    if ( device == nullptr || host == nullptr ) return summary;
    if ( nEquations <= 0 || nFaces < 0 ) return summary;
    if ( ! IsViscousFluxValueFinite( tolerance ) || tolerance < 0.0 )
    {
        return summary;
    }
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( device[ equation ] == nullptr || host[ equation ] == nullptr )
        {
            return summary;
        }
    }

    summary.inputValid = true;

    for ( int face = 0; face < nFaces; ++ face )
    {
        bool deviceNonFinite = false;
        bool hostNonFinite = false;
        Real scale = 0.0;

        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            const Real a = device[ equation ][ face ];
            const Real b = host[ equation ][ face ];
            if ( ! IsViscousFluxValueFinite( a ) ) deviceNonFinite = true;
            if ( ! IsViscousFluxValueFinite( b ) ) hostNonFinite = true;
            if ( deviceNonFinite || hostNonFinite ) continue;
            const Real magnitudeA = std::fabs( a );
            const Real magnitudeB = std::fabs( b );
            const Real magnitude = magnitudeA > magnitudeB ? magnitudeA : magnitudeB;
            if ( magnitude > scale ) scale = magnitude;
        }

        if ( deviceNonFinite ) ++ summary.nonFiniteDeviceFaces;
        if ( hostNonFinite ) ++ summary.nonFiniteHostFaces;

        if ( deviceNonFinite || hostNonFinite )
        {
            if ( summary.firstFace < 0 )
            {
                summary.firstFace = face;
                summary.firstEquation = -1;
            }
            ++ summary.differingFaces;
            if ( summary.differences.size() < maxDifferences )
            {
                ViscousFluxDifference difference;
                difference.face = face;
                difference.equation = -1;
                difference.nonFinite = true;
                summary.differences.push_back( difference );
            }
            continue;
        }

        // Both sides are zero for every equation: nothing to compare.
        if ( scale == 0.0 ) continue;

        ++ summary.comparedFaces;

        int firstDifferingEquation = -1;
        Real firstDifferingRelative = 0.0;
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            const Real difference =
                std::fabs( device[ equation ][ face ] - host[ equation ][ face ] );
            const Real relative = difference / scale;
            // relative is finite here: both inputs and the scale are finite and
            // the scale is non-zero.
            if ( relative <= tolerance ) continue;

            if ( firstDifferingEquation < 0 )
            {
                firstDifferingEquation = equation;
                firstDifferingRelative = relative;
            }
            if ( relative > summary.worstRelative )
            {
                summary.worstRelative = relative;
                summary.worstFace = face;
                summary.worstEquation = equation;
            }
        }

        if ( firstDifferingEquation < 0 ) continue;

        if ( summary.firstFace < 0 )
        {
            summary.firstFace = face;
            summary.firstEquation = firstDifferingEquation;
        }
        ++ summary.differingFaces;
        if ( summary.differences.size() < maxDifferences )
        {
            ViscousFluxDifference difference;
            difference.face = face;
            difference.equation = firstDifferingEquation;
            difference.relative = firstDifferingRelative;
            summary.differences.push_back( difference );
        }
    }

    return summary;
}

EndNameSpace
