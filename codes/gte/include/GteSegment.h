// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

//#include <Mathematics/GteVector.h>
#include <GteVector.h>

// The segment is represented by (1-t)*P0 + t*P1, where P0 and P1 are the
// endpoints of the segment and 0 <= t <= 1.  Some algorithms prefer a
// centered representation that is similar to how oriented bounding boxes are
// defined.  This representation is C + s*D, where C = (P0 + P1)/2 is the
// center of the segment, D = (P1 - P0)/|P1 - P0| is a unit-length direction
// vector for the segment, and |t| <= e.  The value e = |P1 - P0|/2 is the
// extent (or radius or half-length) of the segment.

namespace gte
{

template <int N, typename Real>
class Segment
{
public:
    // Construction and destruction.  The default constructor sets p0 to
    // (-1,0,...,0) and p1 to (1,0,...,0).  NOTE:  If you set p0 and p1;
    // compute C, D, and e; and then recompute q0 = C-e*D and q1 = C+e*D,
    // numerical round-off errors can lead to q0 not exactly equal to p0
    // and q1 not exactly equal to p1.
    Segment();
    Segment(Vector<N, Real> const& p0, Vector<N, Real> const& p1);
    //Segment(std::array<Vector<N, Real>, 2> const& inP);
    Segment(testarray<Vector<N, Real>, 2> const& inP);
    Segment(Vector<N, Real> const& center, Vector<N, Real> const& direction,
        Real extent);

    // Manipulation via the centered form.
    void SetCenteredForm(Vector<N, Real> const& center,
        Vector<N, Real> const& direction, Real extent);

    void GetCenteredForm(Vector<N, Real>& center, Vector<N, Real>& direction,
        Real& extent) const;

    // Public member access.
    //std::array<Vector<N, Real>, 2> p;
    testarray<Vector<N, Real>, 2> p;
    

public:
    // Comparisons to support sorted containers.
    bool operator==(Segment const& segment) const;
    bool operator!=(Segment const& segment) const;
    bool operator< (Segment const& segment) const;
    bool operator<=(Segment const& segment) const;
    bool operator> (Segment const& segment) const;
    bool operator>=(Segment const& segment) const;
};

// Template aliases for convenience.
//template <typename Real>
//using Segment2 = Segment<2, Real>;
//
//template <typename Real>
//using Segment3 = Segment<3, Real>;


template <int N, typename Real>
Segment<N, Real>::Segment()
{
    p[1].MakeUnit(0);
    p[0] = -p[1];
}

template <int N, typename Real>
Segment<N, Real>::Segment(Vector<N, Real> const& p0,
    Vector<N, Real> const& p1)
{
    p[0] = p0;
    p[1] = p1;
}

//template <int N, typename Real>
//Segment<N, Real>::Segment(std::array<Vector<N, Real>, 2> const& inP)
//{
//    p = inP;
//}

template <int N, typename Real>
Segment<N, Real>::Segment(testarray<Vector<N, Real>, 2> const& inP)
{
    p = inP;
}

template <int N, typename Real>
Segment<N, Real>::Segment(Vector<N, Real> const& center,
    Vector<N, Real> const& direction, Real extent)
{
    SetCenteredForm(center, direction, extent);
}

template <int N, typename Real>
void Segment<N, Real>::SetCenteredForm(Vector<N, Real> const& center,
    Vector<N, Real> const& direction, Real extent)
{
    p[0] = center - extent * direction;
    p[1] = center + extent * direction;
}

template <int N, typename Real>
void Segment<N, Real>::GetCenteredForm(Vector<N, Real>& center,
    Vector<N, Real>& direction, Real& extent) const
{
    center = ((Real)0.5)*(p[0] + p[1]);
    direction = p[1] - p[0];
    extent = ((Real)0.5)*Normalize(direction);
}

template <int N, typename Real>
bool Segment<N, Real>::operator==(Segment const& segment) const
{
    return p == segment.p;
}

template <int N, typename Real>
bool Segment<N, Real>::operator!=(Segment const& segment) const
{
    return !operator==(segment);
}

template <int N, typename Real>
bool Segment<N, Real>::operator<(Segment const& segment) const
{
    return p < segment.p;
}

template <int N, typename Real>
bool Segment<N, Real>::operator<=(Segment const& segment) const
{
    return operator<(segment) || operator==(segment);
}

template <int N, typename Real>
bool Segment<N, Real>::operator>(Segment const& segment) const
{
    return !operator<=(segment);
}

template <int N, typename Real>
bool Segment<N, Real>::operator>=(Segment const& segment) const
{
    return !operator<(segment);
}


}
