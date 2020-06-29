/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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

#include "UnitElement.h"
#include "HXCgns.h"
#include "Dimension.h"
using namespace std;

BeginNameSpace( ONEFLOW )

UnitElement::UnitElement()
{
    ;
}

UnitElement::~UnitElement()
{
    ;
}

int UnitElement::GetElementNodeNumbers( int elementType )
{
    int nodeNumberOfElement = 0;
    if ( elementType == ONEFLOW::NODE )
    {
        nodeNumberOfElement = 1;
    }
    else if ( elementType == ONEFLOW::BAR_2 )
    {
        nodeNumberOfElement = 2;
    }
    else if ( elementType == ONEFLOW::BAR_3 )
    {
        nodeNumberOfElement = 3;
    }
    else if ( elementType == ONEFLOW::TRI_3 )
    {
        nodeNumberOfElement = 3;
    }
    else if ( elementType == ONEFLOW::TRI_6 )
    {
        nodeNumberOfElement = 6;
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        nodeNumberOfElement = 4;
    }
    else if ( elementType == ONEFLOW::QUAD_8 )
    {
        nodeNumberOfElement = 8;
    }
    else if ( elementType == ONEFLOW::QUAD_9 )
    {
        nodeNumberOfElement = 9;
    }
    else if ( elementType == ONEFLOW::TETRA_4 )
    {
        nodeNumberOfElement = 4;
    }
    else if ( elementType == ONEFLOW::TETRA_10 )
    {
        nodeNumberOfElement = 10;
    }
    else if ( elementType == ONEFLOW::PYRA_5 )
    {
        nodeNumberOfElement = 5;
    }
    else if ( elementType == ONEFLOW::PYRA_14 )
    {
        nodeNumberOfElement = 14;
    }
    else if ( elementType == ONEFLOW::PENTA_6 )
    {
        nodeNumberOfElement = 6;
    }
    else if ( elementType == ONEFLOW::PENTA_15 )
    {
        nodeNumberOfElement = 15;
    }
    else if ( elementType == ONEFLOW::PENTA_18 )
    {
        nodeNumberOfElement = 18;
    }
    else if ( elementType == ONEFLOW::HEXA_8 )
    {
        nodeNumberOfElement = 8;
    }
    else if ( elementType == ONEFLOW::HEXA_20 )
    {
        nodeNumberOfElement = 20;
    }
    else if ( elementType == ONEFLOW::HEXA_27 )
    {
        nodeNumberOfElement = 27;
    }
    return nodeNumberOfElement;
}

void UnitElement::Initialize( int elementType )
{
    this->elementType = elementType;

    int elementNodeNumbers = this->GetElementNodeNumbers( elementType );
    middlePointStruct.resize( elementNodeNumbers );

    nodeId.resize( elementNodeNumbers );

    for ( int iNode = 0; iNode < elementNodeNumbers; ++ iNode )
    {
        nodeId[ iNode ] = iNode;
    }

    if ( elementType == ONEFLOW::NODE )
    {
        //Face 1
        // 1
        this->PushElementFace( ONEFLOW::NODE, 1 );
    }
    else if ( elementType == ONEFLOW::BAR_2 )
    {
        //Face 1
        // 1 2
        this->PushElementFace( ONEFLOW::NODE, 1 );
        this->PushElementFace( ONEFLOW::NODE, 2 );
    }
    else if ( elementType == ONEFLOW::BAR_3 )
    {
        //Face 2
        // 1 3
        // 3 2

        this->PushElementFace( ONEFLOW::NODE, 1 );
        this->PushElementFace( ONEFLOW::NODE, 2 );

        //BAR_3
        // 1.........3..........2
        // 对于中点3来说相邻点为1、2
        this->PushMiddlePoint( 3, 1, 2 );

        //subcell
        // BAR_2 1 3
        // BAR_2 3 2

        this->PushChildElement( ONEFLOW::BAR_2, 1, 3 );
        this->PushChildElement( ONEFLOW::BAR_2, 3, 2 );
    }
    else if ( elementType == ONEFLOW::TRI_3 )
    {
        //TRI_3
        // 1.........2.........3

        // Face 3
        // 1 2
        // 2 3
        // 3 1
        this->PushElementFace( ONEFLOW::BAR_2, 1, 2 );
        this->PushElementFace( ONEFLOW::BAR_2, 2, 3 );
        this->PushElementFace( ONEFLOW::BAR_2, 3, 1 );
    }
    else if ( elementType == ONEFLOW::TRI_6 )
    {
        //TRI_6
        // 1....4....2 ...5....3....6....1

        // Composite Face 3
        // BAR_3 1 2 4
        // BAR_3 2 3 5
        // BAR_3 3 1 6

        this->PushCompositeFace( ONEFLOW::BAR_3, 1, 2, 4 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 2, 3, 5 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 3, 1, 6 );
        // Face 6
        // 1 4
        // 4 2
        // 2 5
        // 5 3
        // 3 6
        // 6 1

        this->PushElementFace( ONEFLOW::BAR_2, 1, 4 );
        this->PushElementFace( ONEFLOW::BAR_2, 4, 2 );
        this->PushElementFace( ONEFLOW::BAR_2, 2, 5 );
        this->PushElementFace( ONEFLOW::BAR_2, 5, 3 );
        this->PushElementFace( ONEFLOW::BAR_2, 3, 6 );
        this->PushElementFace( ONEFLOW::BAR_2, 6, 1 );

        //TRI_6
        // 1.........4..........2
        // 2.........5..........3
        // 3.........6..........1

        this->PushMiddlePoint( 4, 1, 2 );
        this->PushMiddlePoint( 5, 2, 3 );
        this->PushMiddlePoint( 6, 3, 1 );

        //subcell
        // TRI_3 1 4 6
        // TRI_3 4 2 5
        // TRI_3 4 5 6
        // TRI_3 6 5 3 

        this->PushChildElement( ONEFLOW::TRI_3, 1, 4, 6 );
        this->PushChildElement( ONEFLOW::TRI_3, 4, 2, 5 );
        this->PushChildElement( ONEFLOW::TRI_3, 4, 5, 6 );
        this->PushChildElement( ONEFLOW::TRI_3, 6, 5, 3 );
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        // QUAD_4
        // 1....2....3 ...4....1

        // Face 4
        // 1 2
        // 2 3
        // 3 4
        // 4 1
        this->PushElementFace( ONEFLOW::BAR_2, 1, 2 );
        this->PushElementFace( ONEFLOW::BAR_2, 2, 3 );
        this->PushElementFace( ONEFLOW::BAR_2, 3, 4 );
        this->PushElementFace( ONEFLOW::BAR_2, 4, 1 );
    }
    else if ( elementType == ONEFLOW::QUAD_8 )
    {
        // QUAD_8
        // 1..5..2..6..3..7..4..8..1

        // Composite Face 4
        // BAR_3 1 2 5
        // BAR_3 2 3 6
        // BAR_3 3 4 7
        // BAR_3 4 1 8

        this->PushCompositeFace( ONEFLOW::BAR_3, 1, 2, 5 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 2, 3, 6 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 3, 4, 7 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 4, 1, 8 );

        // Face 8
        // 1 5
        // 5 2
        // 2 6
        // 6 3

        // 3 7
        // 7 4
        // 4 8
        // 8 1

        this->PushElementFace( ONEFLOW::BAR_2, 1, 5 );
        this->PushElementFace( ONEFLOW::BAR_2, 5, 2 );
        this->PushElementFace( ONEFLOW::BAR_2, 2, 6 );
        this->PushElementFace( ONEFLOW::BAR_2, 6, 3 );

        this->PushElementFace( ONEFLOW::BAR_2, 3, 7 );
        this->PushElementFace( ONEFLOW::BAR_2, 7, 4 );
        this->PushElementFace( ONEFLOW::BAR_2, 4, 8 );
        this->PushElementFace( ONEFLOW::BAR_2, 8, 1 );

        //QUAD_8
        // 1.........5..........2
        // 2.........6..........3
        // 3.........7..........4
        // 4.........8..........1

        this->PushMiddlePoint( 5, 1, 2 );
        this->PushMiddlePoint( 6, 2, 3 );
        this->PushMiddlePoint( 7, 3, 4 );
        this->PushMiddlePoint( 8, 4, 1 );
    }
    else if ( elementType == ONEFLOW::QUAD_9 )
    {
        // QUAD_9
        // 1..5.2..6..3..7..4..8..1..9

        // Composite Face 4
        // BAR_3 1 2 5
        // BAR_3 2 3 6
        // BAR_3 3 4 7
        // BAR_3 4 1 8
        this->PushCompositeFace( ONEFLOW::BAR_3, 1, 2, 5 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 2, 3, 6 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 3, 4, 7 );
        this->PushCompositeFace( ONEFLOW::BAR_3, 4, 1, 8 );

        // Face 8
        // 1 5
        // 5 2
        // 2 6
        // 6 3

        // 3 7
        // 7 4
        // 4 8
        // 8 1
        this->PushElementFace( ONEFLOW::BAR_2, 1, 5 );
        this->PushElementFace( ONEFLOW::BAR_2, 5, 2 );
        this->PushElementFace( ONEFLOW::BAR_2, 2, 6 );
        this->PushElementFace( ONEFLOW::BAR_2, 6, 3 );

        this->PushElementFace( ONEFLOW::BAR_2, 3, 7 );
        this->PushElementFace( ONEFLOW::BAR_2, 7, 4 );
        this->PushElementFace( ONEFLOW::BAR_2, 4, 8 );
        this->PushElementFace( ONEFLOW::BAR_2, 8, 1 );
        //QUAD_9
        // 1.........5..........2
        // 2.........6..........3
        // 3.........7..........4
        // 4.........8..........1
        // 1 2.......9........3 4

        this->PushMiddlePoint( 5, 1, 2 );
        this->PushMiddlePoint( 6, 2, 3 );
        this->PushMiddlePoint( 7, 3, 4 );
        this->PushMiddlePoint( 8, 4, 1 );
        this->PushMiddlePoint( 9, 1, 2, 3, 4 );
        //subcell
        // QUAD_4 : 1 5 9 8
        // QUAD_4 : 5 2 6 9
        // QUAD_4 : 9 6 3 7
        // QUAD_4 : 8 9 7 4

        this->PushChildElement( ONEFLOW::QUAD_4, 1, 5, 9, 8 );
        this->PushChildElement( ONEFLOW::QUAD_4, 5, 2, 6, 9 );
        this->PushChildElement( ONEFLOW::QUAD_4, 9, 6, 3, 7 );
        this->PushChildElement( ONEFLOW::QUAD_4, 8, 9, 7, 4 );
    }
    else if ( elementType == ONEFLOW::TETRA_4 )
    {
        // TETRA_4
        // 1..2...3...4

        // Face 4
        // 1 3 2
        // 1 2 4
        // 2 3 4
        // 3 1 4
        this->PushElementFace( ONEFLOW::TRI_3, 1, 3, 2 );
        this->PushElementFace( ONEFLOW::TRI_3, 1, 2, 4 );
        this->PushElementFace( ONEFLOW::TRI_3, 2, 3, 4 );
        this->PushElementFace( ONEFLOW::TRI_3, 3, 1, 4 );
    }
    else if ( elementType == ONEFLOW::TETRA_10 )
    {
        // TETRA_10
        // 1..7..3..6..2..5..8..10..9..4

        // Composite Face 4
        // TRI_6 1  3  2  7  6  5
        // TRI_6 2  3  4  6  10 9
        // TRI_6 3  1  4  7  8  10
        // TRI_6 1  2  4  5  9  8

        this->PushCompositeFace( ONEFLOW::TRI_6, 1, 3, 2, 7, 6, 5 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 2, 3, 4, 6, 10, 9 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 3, 1, 4, 7, 8, 10 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 1, 2, 4, 5, 9, 8 );
        // Face 16
        // 1  7  5
        // 7  6  5
        // 5  6  2
        // 6  7  3

        this->PushElementFace( ONEFLOW::TRI_3, 1, 7, 5 );
        this->PushElementFace( ONEFLOW::TRI_3, 7, 6, 5 );
        this->PushElementFace( ONEFLOW::TRI_3, 5, 6, 2 );
        this->PushElementFace( ONEFLOW::TRI_3, 6, 7, 3 );


        // 2..6..3..10..4..9..2
        // 2  6  9
        // 9  6  10
        // 10 6  3
        // 10 4  9

        this->PushElementFace( ONEFLOW::TRI_3, 2, 6, 9 );
        this->PushElementFace( ONEFLOW::TRI_3, 9, 6, 10 );
        this->PushElementFace( ONEFLOW::TRI_3, 10, 6, 3 );
        this->PushElementFace( ONEFLOW::TRI_3, 10, 4, 9 );

        // 1..8..4..10..3..7..1
        // 1  8  7
        // 8  10 7
        // 10 3  7
        // 10 8  4

        this->PushElementFace( ONEFLOW::TRI_3, 1, 8, 7 );
        this->PushElementFace( ONEFLOW::TRI_3, 8, 10, 7 );
        this->PushElementFace( ONEFLOW::TRI_3, 10, 3, 7 );
        this->PushElementFace( ONEFLOW::TRI_3, 10, 8, 4 );

        // 1..5..2..9..4..8..1
        // 1  5  8
        // 8  5  9
        // 9  5  2
        // 9  4  8

        this->PushElementFace( ONEFLOW::TRI_3, 1, 5, 8 );
        this->PushElementFace( ONEFLOW::TRI_3, 8, 5, 9 );
        this->PushElementFace( ONEFLOW::TRI_3, 9, 5, 2 );
        this->PushElementFace( ONEFLOW::TRI_3, 9, 4, 8 );

        //TETRA_10
        // 1.........5..........2
        // 2.........6..........3
        // 3.........7..........1

        // 1.........8..........4
        // 2.........9..........4
        // 3.........10.........4

        this->PushMiddlePoint( 5, 1, 2 );
        this->PushMiddlePoint( 6, 2, 3 );
        this->PushMiddlePoint( 7, 3, 4 );
        this->PushMiddlePoint( 8, 1, 4 );
        this->PushMiddlePoint( 9, 2, 4 );
        this->PushMiddlePoint( 10, 3, 4 );
        //subcell
        // TETRA_4 : 8  5  1  7
        // TETRA_4 : 8  9  5  7
        // TETRA_4 : 9  10 6  7
        // TETRA_4 : 10 3  6  7
        // TETRA_4 : 2  6  5  9
        // TETRA_4 : 5  6  7  9
        // TETRA_4 : 4  8  10 9
        // TETRA_4 : 10 8  7  9

        this->PushChildElement( ONEFLOW::TETRA_4, 8, 5, 1, 7 );
        this->PushChildElement( ONEFLOW::TETRA_4, 8, 9, 5, 7 );
        this->PushChildElement( ONEFLOW::TETRA_4, 9, 10, 6, 7 );
        this->PushChildElement( ONEFLOW::TETRA_4, 10, 3, 6, 7 );
        this->PushChildElement( ONEFLOW::TETRA_4, 2, 6, 5, 9 );
        this->PushChildElement( ONEFLOW::TETRA_4, 5, 6, 7, 9 );
        this->PushChildElement( ONEFLOW::TETRA_4, 4, 8, 10, 9 );
        this->PushChildElement( ONEFLOW::TETRA_4, 10, 8, 7, 9 );
    }
    else if ( elementType == ONEFLOW::PYRA_5 )
    {
        // PYRA_5
        // 1..2..3..4..5

        // Face 5
        // 1  4  3  2
        // 1  5  4
        // 1  2  5
        // 2  3  5
        // 3  4  5

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 4, 3, 2 );
        this->PushElementFace( ONEFLOW::TRI_3, 1, 5, 4 );
        this->PushElementFace( ONEFLOW::TRI_3, 1, 2, 5 );
        this->PushElementFace( ONEFLOW::TRI_3, 2, 3, 5 );
        this->PushElementFace( ONEFLOW::TRI_3, 3, 4, 5 );

    }
    else if ( elementType == ONEFLOW::PYRA_14 )
    {
        // PYRA_14
        // Face 20

        // Composite Face 5
        // QUAD_9 : 1  4  3  2  9   8   7  6  14
        // TRI_6  : 1  2  5  6  11  10
        // TRI_6  : 2  3  5  7  12  11
        // TRI_6  : 3  4  5  8  13  12
        // TRI_6  : 4  1  5  9  10  13

        this->PushCompositeFace( ONEFLOW::QUAD_9, 1, 4, 3, 2, 9, 8, 7, 6, 14 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 1, 2, 5, 6, 11, 10 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 2, 3, 5, 7, 12, 11 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 3, 4, 5, 8, 13, 12 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 4, 1, 5, 9, 10, 13 );

        // 1  9  4  8  3  7  2  6  1   middle 14
        // Face 4
        // 1  9  14 6
        // 9  4  8  14
        // 6  14 7  2
        // 14 8  3  7

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 9, 14, 6 );
        this->PushElementFace( ONEFLOW::QUAD_4, 9, 4, 8, 14 );
        this->PushElementFace( ONEFLOW::QUAD_4, 6, 14, 7, 2 );
        this->PushElementFace( ONEFLOW::QUAD_4, 14, 8, 3, 7 );

        // 1  10  5  13  4  9  1
        // Face 4
        // 1  10  9
        // 9  10  13
        // 4  9   13
        // 10 5   13

        this->PushElementFace( ONEFLOW::TRI_3, 1, 10, 9 );
        this->PushElementFace( ONEFLOW::TRI_3, 9, 10, 13 );
        this->PushElementFace( ONEFLOW::TRI_3, 4, 9, 13 );
        this->PushElementFace( ONEFLOW::TRI_3, 10, 5, 13 );

        // 1  6  2  11  5  10  1
        // Face 4
        // 1  6  10
        // 6  11 10
        // 6  2  11
        // 10 11 5

        this->PushElementFace( ONEFLOW::TRI_3, 1, 6, 10 );
        this->PushElementFace( ONEFLOW::TRI_3, 6, 11, 10 );
        this->PushElementFace( ONEFLOW::TRI_3, 6, 2, 11 );
        this->PushElementFace( ONEFLOW::TRI_3, 10, 11, 5 );

        // 2  7   3  12  5  11  2
        // Face 4
        // 2  7  11
        // 7  12 11
        // 7  3  12
        // 11 12 5

        this->PushElementFace( ONEFLOW::TRI_3, 2, 7, 11 );
        this->PushElementFace( ONEFLOW::TRI_3, 7, 12, 11 );
        this->PushElementFace( ONEFLOW::TRI_3, 7, 3, 12 );
        this->PushElementFace( ONEFLOW::TRI_3, 11, 12, 5 );

        // 3  8  4  13  5  12  3
        // Face 4
        // 3  8  12
        // 8  13 12
        // 8  4  13
        // 12 13 5

        this->PushElementFace( ONEFLOW::TRI_3, 3, 8, 12 );
        this->PushElementFace( ONEFLOW::TRI_3, 8, 13, 12 );
        this->PushElementFace( ONEFLOW::TRI_3, 8, 4, 13 );
        this->PushElementFace( ONEFLOW::TRI_3, 12, 13, 5 );

        //PYRA_14
        // 1.........6..........2
        // 2.........7..........3
        // 3.........8..........4
        // 4.........9..........1

        // 1.........10.........5
        // 2.........11.........5
        // 3.........12.........5
        // 4.........13.........5

        // 1 4.......14.......3 2//

        this->PushMiddlePoint( 6, 1, 2 );
        this->PushMiddlePoint( 7, 2, 3 );
        this->PushMiddlePoint( 8, 3, 4 );
        this->PushMiddlePoint( 9, 4, 1 );
        this->PushMiddlePoint( 10, 1, 5 );
        this->PushMiddlePoint( 11, 2, 5 );
        this->PushMiddlePoint( 12, 3, 5 );
        this->PushMiddlePoint( 13, 4, 5 );
        this->PushMiddlePoint( 14, 1, 4, 3, 2 );

        //subcell
        // PYRA_5  : 10 11 12 13 5
        // PYRA_5  : 1  6  14 9  10
        // PYRA_5  : 6  2  7  14 11
        // PYRA_5  : 4  9  14 8  13
        // PYRA_5  : 14 7  3  8  12
        // PYRA_5  : 10 13 12 11 14
        // TETRA_4 : 6  10 11 14
        // TETRA_4 : 8  12 13 14
        // TETRA_4 : 9  13 10 14
        // TETRA_4 : 7  11 12 14

        this->PushChildElement( ONEFLOW::PYRA_5, 10, 11, 12, 13, 5 );
        this->PushChildElement( ONEFLOW::PYRA_5, 1, 6, 14, 9, 10 );
        this->PushChildElement( ONEFLOW::PYRA_5, 6, 2, 7, 14, 11 );
        this->PushChildElement( ONEFLOW::PYRA_5, 4, 9, 14, 8, 13 );
        this->PushChildElement( ONEFLOW::PYRA_5, 14, 7, 3, 8, 12 );
        this->PushChildElement( ONEFLOW::PYRA_5, 10, 13, 12, 11, 14 );
        this->PushChildElement( ONEFLOW::TETRA_4, 6, 10, 11, 14 );
        this->PushChildElement( ONEFLOW::TETRA_4, 8, 12, 13, 14 );
        this->PushChildElement( ONEFLOW::TETRA_4, 9, 13, 10, 14 );
        this->PushChildElement( ONEFLOW::TETRA_4, 7, 11, 12, 14 );
    }
    else if ( elementType == ONEFLOW::PENTA_6 )
    {
        // PENTA_6

        // 1  2  3  4  5  6
        // Face 5

        // 1  3  2
        // 4  5  6
        // 1  2  5  4
        // 2  3  6  5
        // 3  1  4  6

        this->PushElementFace( ONEFLOW::TRI_3, 1, 3, 2 );
        this->PushElementFace( ONEFLOW::TRI_3, 4, 5, 6 );

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 2, 5, 4 );
        this->PushElementFace( ONEFLOW::QUAD_4, 2, 3, 6, 5 );
        this->PushElementFace( ONEFLOW::QUAD_4, 3, 1, 4, 6 );
    }
    else if ( elementType == ONEFLOW::PENTA_15 )
    {
        //PENTA_15
        // 1.........7..........2
        // 2.........8..........3
        // 3.........9..........1

        // 1.........10.........4
        // 2.........11.........5
        // 3.........12.........6

        // 4.........13..........5
        // 5.........14..........6
        // 6.........15..........4

        this->PushMiddlePoint( 7, 1, 2 );
        this->PushMiddlePoint( 8, 2, 3 );
        this->PushMiddlePoint( 9, 3, 1 );

        this->PushMiddlePoint( 10, 1, 4 );
        this->PushMiddlePoint( 11, 2, 5 );
        this->PushMiddlePoint( 12, 3, 6 );

        this->PushMiddlePoint( 13, 4, 5 );
        this->PushMiddlePoint( 14, 5, 6 );
        this->PushMiddlePoint( 15, 6, 4 );

    }
    else if ( elementType == ONEFLOW::PENTA_18 )
    {
        // PENTA_18

        // Composite Face 5
        // QUAD_9 : 1  2  5  4  7  11  13  10  16
        // QUAD_9 : 2  3  6  5  8  12  14  11  17
        // QUAD_9 : 3  1  4  6  9  10  15  12  18
        // TRI_6  : 1  3  2  9  8  7
        // TRI_6  : 4  5  6  13 14 15

        this->PushCompositeFace( ONEFLOW::QUAD_9, 1, 2, 5, 4, 7, 11, 13, 10, 16 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 2, 3, 6, 5, 8, 12, 14, 11, 17 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 3, 1, 4, 6, 9, 10, 15, 12, 18 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 1, 3, 2, 9, 8, 7 );
        this->PushCompositeFace( ONEFLOW::TRI_6, 4, 5, 6, 13, 14, 15 );

        // Face 20

        // 1  9  3  8  2  7  1
        // Face 4
        // 1  9  7
        // 9  8  7
        // 7  8  2
        // 9  3  8

        this->PushElementFace( ONEFLOW::TRI_3, 1, 9, 7 );
        this->PushElementFace( ONEFLOW::TRI_3, 9, 8, 7 );
        this->PushElementFace( ONEFLOW::TRI_3, 7, 8, 2 );
        this->PushElementFace( ONEFLOW::TRI_3, 9, 3, 8 );

        // 4  13  5  14  6  15  4
        // Face 4
        // 4  13 15
        // 13 14 15
        // 13 5  14
        // 15 14 6

        this->PushElementFace( ONEFLOW::TRI_3, 4, 13, 15 );
        this->PushElementFace( ONEFLOW::TRI_3, 13, 14, 15 );
        this->PushElementFace( ONEFLOW::TRI_3, 13, 5, 14 );
        this->PushElementFace( ONEFLOW::TRI_3, 15, 14, 6 );

        // 1  7  2  11  5  13  4  10  1  middle16
        // Face 4
        // 1  7  16 10
        // 7  2  11 16
        // 10 16 13 4
        // 16 11 5  13

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 7, 16, 10 );
        this->PushElementFace( ONEFLOW::QUAD_4, 7, 2, 11, 16 );
        this->PushElementFace( ONEFLOW::QUAD_4, 10, 16, 13, 4 );
        this->PushElementFace( ONEFLOW::QUAD_4, 16, 11, 5, 13 );
        // 2  8  3  12  6  14  5  11  2  middle17
        // Face 4
        // 2  8  17 11
        // 8  3  12 17
        // 11 17 14 5
        // 17 12 6  14

        this->PushElementFace( ONEFLOW::QUAD_4, 2, 8, 17, 11 );
        this->PushElementFace( ONEFLOW::QUAD_4, 8, 3, 12, 17 );
        this->PushElementFace( ONEFLOW::QUAD_4, 11, 17, 14, 5 );
        this->PushElementFace( ONEFLOW::QUAD_4, 17, 12, 6, 14 );

        // 1  10  4  15  6  12  3  9  1  middle18
        // Face 4
        // 1  10 18 9
        // 9  18 12 3
        // 10 4  15 18
        // 18 15 6  12

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 10, 18, 9 );
        this->PushElementFace( ONEFLOW::QUAD_4, 9, 18, 12, 3 );
        this->PushElementFace( ONEFLOW::QUAD_4, 10, 4, 15, 18 );
        this->PushElementFace( ONEFLOW::QUAD_4, 18, 15, 6, 12 );

        //PENTA_18
        // 1.........7..........2
        // 2.........8..........3
        // 3.........9..........1

        // 1.........10.........4
        // 2.........11.........5
        // 3.........12.........6

        // 4.........13.........5
        // 5.........14.........6
        // 6.........15.........4

        // 1.2.......16.......5.4//
        // 2.3.......17.......6.5//
        // 3.1.......18.......4.6//

        this->PushMiddlePoint( 7, 1, 2 );
        this->PushMiddlePoint( 8, 2, 3 );
        this->PushMiddlePoint( 9, 3, 1 );

        this->PushMiddlePoint( 10, 1, 4 );
        this->PushMiddlePoint( 11, 2, 5 );
        this->PushMiddlePoint( 12, 3, 6 );

        this->PushMiddlePoint( 13, 4, 5 );
        this->PushMiddlePoint( 14, 5, 6 );
        this->PushMiddlePoint( 15, 6, 4 );

        this->PushMiddlePoint( 16, 1, 2, 5, 4 );
        this->PushMiddlePoint( 17, 2, 3, 6, 5 );
        this->PushMiddlePoint( 18, 3, 1, 4, 6 );
        //subcell
        // PENTA_6  : 1  7  9  10 16 18
        // PENTA_6  : 7  8  9  16 17 18
        // PENTA_6  : 7  2  8  16 11 17
        // PENTA_6  : 9  8  3  18 17 12
        // PENTA_6  : 10 16 18 4  13 15
        // PENTA_6  : 16 17 18 13 14 15
        // PENTA_6  : 16 11 17 13 5  14
        // PENTA_6  : 18 17 12 15 14 6

        this->PushChildElement( ONEFLOW::PENTA_6, 1, 7, 9, 10, 16, 18 );
        this->PushChildElement( ONEFLOW::PENTA_6, 7, 8, 9, 16, 17, 18 );
        this->PushChildElement( ONEFLOW::PENTA_6, 7, 2, 8, 16, 11, 17 );
        this->PushChildElement( ONEFLOW::PENTA_6, 9, 8, 3, 18, 17, 12 );
        this->PushChildElement( ONEFLOW::PENTA_6, 10, 16, 18, 4, 13, 15 );
        this->PushChildElement( ONEFLOW::PENTA_6, 16, 17, 18, 13, 14, 15 );
        this->PushChildElement( ONEFLOW::PENTA_6, 16, 11, 17, 13, 5, 14 );
        this->PushChildElement( ONEFLOW::PENTA_6, 18, 17, 12, 15, 14, 6 );

    }
    else if ( elementType == ONEFLOW::HEXA_8 )
    {
        // HEXA_8

        // 1 2 3 4 5 6 7 8
        // Face 6
        // 1, 4, 3, 2
        // 1, 2, 6, 5
        // 2, 3, 7, 6
        // 3, 4, 8, 7
        // 1, 5, 8, 4
        // 5, 6, 7, 8
        this->PushElementFace( ONEFLOW::QUAD_4, 1, 4, 3, 2 );
        this->PushElementFace( ONEFLOW::QUAD_4, 1, 2, 6, 5 );
        this->PushElementFace( ONEFLOW::QUAD_4, 2, 3, 7, 6 );
        this->PushElementFace( ONEFLOW::QUAD_4, 3, 4, 8, 7 );
        this->PushElementFace( ONEFLOW::QUAD_4, 1, 5, 8, 4 );
        this->PushElementFace( ONEFLOW::QUAD_4, 5, 6, 7, 8 );
    }
    else if ( elementType == ONEFLOW::HEXA_20 )
    {
        //HEXA_20
        // 1.........9 .........2
        // 2.........10.........3
        // 3.........11.........4
        // 4.........12.........1

        // 1.........13.........5
        // 2.........14.........6
        // 3.........15.........7
        // 4.........16.........8

        // 5.........17.........6
        // 6.........18.........7
        // 7.........19.........8
        // 8.........20.........5

        this->PushMiddlePoint( 9, 1, 2 );
        this->PushMiddlePoint( 10, 2, 3 );
        this->PushMiddlePoint( 11, 3, 4 );
        this->PushMiddlePoint( 12, 4, 1 );

        this->PushMiddlePoint( 13, 1, 5 );
        this->PushMiddlePoint( 14, 2, 6 );
        this->PushMiddlePoint( 15, 3, 7 );
        this->PushMiddlePoint( 16, 4, 8 );

        this->PushMiddlePoint( 17, 5, 6 );
        this->PushMiddlePoint( 18, 6, 7 );
        this->PushMiddlePoint( 19, 7, 8 );
        this->PushMiddlePoint( 20, 8, 5 );
    }
    else if ( elementType == ONEFLOW::HEXA_27 )
    {
        // HEXA_27

        // Composite Face 6
        // QUAD_9 : 1  4  3  2  12  11  10  9   21
        // QUAD_9 : 5  6  7  8  17  18  19  20  26
        // QUAD_9 : 1  2  6  5  9   14  17  13  22
        // QUAD_9 : 2  3  7  6  10  15  18  14  23
        // QUAD_9 : 3  4  8  7  11  16  19  15  24
        // QUAD_9 : 1  5  8  4  13  20  16  12  25

        this->PushCompositeFace( ONEFLOW::QUAD_9, 1, 4, 3, 2, 12, 11, 10, 9, 21 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 5, 6, 7, 8, 17, 18, 19, 20, 26 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 1, 2, 6, 5, 9, 14, 17, 13, 22 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 2, 3, 7, 6, 10, 15, 18, 14, 23 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 3, 4, 8, 7, 11, 16, 19, 15, 24 );
        this->PushCompositeFace( ONEFLOW::QUAD_9, 1, 5, 8, 4, 13, 20, 16, 12, 25 );

        // Face 24

        // 1 12 4 11 3  10  2  9  1  middle 21
        // Face 4
        // 1  12 21 9
        // 12 4  11 21
        // 9  21 10 2
        // 21 11 3  10

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 12, 21, 9 );
        this->PushElementFace( ONEFLOW::QUAD_4, 12, 4, 11, 21 );
        this->PushElementFace( ONEFLOW::QUAD_4, 9, 21, 10, 2 );
        this->PushElementFace( ONEFLOW::QUAD_4, 21, 11, 3, 10 );

        // 5 17 6  18 7 19 8 20 5  middle 26
        // Face 4
        // 5  17 26 20
        // 20 26 19 8
        // 17 6  18 26
        // 26 18 7  19

        this->PushElementFace( ONEFLOW::QUAD_4, 5, 17, 26, 20 );
        this->PushElementFace( ONEFLOW::QUAD_4, 20, 26, 19, 8 );
        this->PushElementFace( ONEFLOW::QUAD_4, 17, 6, 18, 26 );
        this->PushElementFace( ONEFLOW::QUAD_4, 26, 18, 7, 19 );

        // 1  13  5  20  8  16 4 12 1  middle 25
        // Face 4
        // 1  13 25 12
        // 12 25 16  4
        // 13 5  20 25
        // 25 20 8  16

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 13, 25, 12 );
        this->PushElementFace( ONEFLOW::QUAD_4, 12, 25, 16, 4 );
        this->PushElementFace( ONEFLOW::QUAD_4, 13, 5, 20, 25 );
        this->PushElementFace( ONEFLOW::QUAD_4, 25, 20, 8, 16 );

        // 2  10  3  15  7 18 6 14 2  middle 23
        // Face 4
        // 2  10 23 14
        // 10 3  15 23
        // 14 23 18 6
        // 23 15 7  18

        this->PushElementFace( ONEFLOW::QUAD_4, 2, 10, 23, 14 );
        this->PushElementFace( ONEFLOW::QUAD_4, 10, 3, 15, 23 );
        this->PushElementFace( ONEFLOW::QUAD_4, 14, 23, 18, 6 );
        this->PushElementFace( ONEFLOW::QUAD_4, 23, 15, 7, 18 );
        // 4 16 8 19  7 15 3 11 4  middle 24
        // Face 4
        // 4  16 24 11
        // 11 24 15 3
        // 16 8  19 24
        // 24 19 7  15

        this->PushElementFace( ONEFLOW::QUAD_4, 4, 16, 24, 11 );
        this->PushElementFace( ONEFLOW::QUAD_4, 11, 24, 15, 3 );
        this->PushElementFace( ONEFLOW::QUAD_4, 16, 8, 19, 24 );
        this->PushElementFace( ONEFLOW::QUAD_4, 24, 19, 7, 15 );

        // 1 9 2 14 6 17 5 13 1  middle 22
        // Face 4
        // 1  9  22 13
        // 9  2  14 22
        // 13 22 17 5
        // 22 14 6  17

        this->PushElementFace( ONEFLOW::QUAD_4, 1, 9, 22, 13 );
        this->PushElementFace( ONEFLOW::QUAD_4, 9, 2, 14, 22 );
        this->PushElementFace( ONEFLOW::QUAD_4, 13, 22, 17, 5 );
        this->PushElementFace( ONEFLOW::QUAD_4, 22, 14, 6, 17 );
        //HEXA_27
        // 1.........9. ........2
        // 2.........10.........3
        // 3.........11.........4
        // 4.........12.........1

        // 1.........13.........5
        // 2.........14.........6
        // 3.........15.........7
        // 4.........16.........8

        // 5.........17.........6
        // 6.........18.........7
        // 7.........19.........8
        // 8.........20.........5

        this->PushMiddlePoint( 9, 1, 2 );
        this->PushMiddlePoint( 10, 2, 3 );
        this->PushMiddlePoint( 11, 3, 4 );
        this->PushMiddlePoint( 12, 4, 1 );

        this->PushMiddlePoint( 13, 1, 5 );
        this->PushMiddlePoint( 14, 2, 6 );
        this->PushMiddlePoint( 15, 3, 7 );
        this->PushMiddlePoint( 16, 4, 8 );

        this->PushMiddlePoint( 17, 5, 6 );
        this->PushMiddlePoint( 18, 6, 7 );
        this->PushMiddlePoint( 19, 7, 8 );
        this->PushMiddlePoint( 20, 8, 5 );

        // 1 4.......21.........3 2
        // 1 2.......22.........6 5
        // 2 3.......23.........7 6
        // 3 4.......24.........8 7
        // 1 5.......25.........8 4
        // 5 6.......26.........7 8
        // 1 2 3 4...27.....5 6 7 8

        this->PushMiddlePoint( 21, 1, 4, 3, 2 );
        this->PushMiddlePoint( 22, 1, 2, 6, 5 );
        this->PushMiddlePoint( 23, 2, 3, 7, 6 );
        this->PushMiddlePoint( 24, 3, 4, 8, 7 );
        this->PushMiddlePoint( 25, 1, 5, 8, 4 );
        this->PushMiddlePoint( 26, 5, 6, 7, 8 );

        this->PushMiddlePoint( 27, 1, 2, 3, 4, 5, 6, 7, 8 );
        //subcell
        // HEXA_8  : 4  12 21 11 16 25 27 24
        // HEXA_8  : 12 1  9  21 25 13 22 27
        // HEXA_8  : 11 21 10 3  24 27 23 15
        // HEXA_8  : 21 9  2  10 27 22 14 23
        // HEXA_8  : 16 25 27 24 8  20 26 19
        // HEXA_8  : 25 13 22 27 20 5  17 26
        // HEXA_8  : 24 27 23 15 19 26 18 7
        // HEXA_8  : 27 22 14 23 26 17 6  18

        this->PushChildElement( ONEFLOW::HEXA_8, 4, 12, 21, 11, 16, 25, 27, 24 );
        this->PushChildElement( ONEFLOW::HEXA_8, 12, 1, 9, 21, 25, 13, 22, 27 );
        this->PushChildElement( ONEFLOW::HEXA_8, 11, 21, 10, 3, 24, 27, 23, 15 );
        this->PushChildElement( ONEFLOW::HEXA_8, 21, 9, 2, 10, 27, 22, 14, 23 );
        this->PushChildElement( ONEFLOW::HEXA_8, 16, 25, 27, 24, 8, 20, 26, 19 );
        this->PushChildElement( ONEFLOW::HEXA_8, 25, 13, 22, 27, 20, 5, 17, 26 );
        this->PushChildElement( ONEFLOW::HEXA_8, 24, 27, 23, 15, 19, 26, 18, 7 );
        this->PushChildElement( ONEFLOW::HEXA_8, 27, 22, 14, 23, 26, 17, 6, 18 );
    }
}

int UnitElement::GetRefinedElementType( int elementType )
{
    int nextElementType = -1;
    if ( elementType == ONEFLOW::NODE )
    {
    }
    else if ( elementType == ONEFLOW::BAR_2 )
    {
        nextElementType = ONEFLOW::BAR_3;
    }
    else if ( elementType == ONEFLOW::BAR_3 )
    {
    }
    else if ( elementType == ONEFLOW::TRI_3 )
    {
        nextElementType = ONEFLOW::TRI_6;
    }
    else if ( elementType == ONEFLOW::TRI_6 )
    {
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        nextElementType = ONEFLOW::QUAD_9;
    }
    else if ( elementType == ONEFLOW::QUAD_8 )
    {
    }
    else if ( elementType == ONEFLOW::QUAD_9 )
    {
    }
    else if ( elementType == ONEFLOW::TETRA_4 )
    {
        nextElementType = ONEFLOW::TETRA_10;
    }
    else if ( elementType == ONEFLOW::TETRA_10 )
    {
    }
    else if ( elementType == ONEFLOW::PYRA_5 )
    {
        nextElementType = ONEFLOW::PYRA_14;
    }
    else if ( elementType == ONEFLOW::PYRA_14 )
    {
    }
    else if ( elementType == ONEFLOW::PENTA_6 )
    {
        nextElementType = ONEFLOW::PENTA_18;
    }
    else if ( elementType == ONEFLOW::PENTA_15 )
    {
    }
    else if ( elementType == ONEFLOW::PENTA_18 )
    {
    }
    else if ( elementType == ONEFLOW::HEXA_8 )
    {
        nextElementType = ONEFLOW::HEXA_27;
    }
    else if ( elementType == ONEFLOW::HEXA_20 )
    {
    }
    else if ( elementType == ONEFLOW::HEXA_27 )
    {
    }
    return nextElementType;
}

int UnitElement::GetSimpleElementType( int elementType )
{
    int previousElementType = -1;
    if ( elementType == ONEFLOW::NODE )
    {
    }
    else if ( elementType == ONEFLOW::BAR_2 )
    {
    }
    else if ( elementType == ONEFLOW::BAR_3 )
    {
        previousElementType = ONEFLOW::BAR_2;
    }
    else if ( elementType == ONEFLOW::TRI_3 )
    {
    }
    else if ( elementType == ONEFLOW::TRI_6 )
    {
        previousElementType = ONEFLOW::TRI_3;
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
    }
    else if ( elementType == ONEFLOW::QUAD_8 )
    {
    }
    else if ( elementType == ONEFLOW::QUAD_9 )
    {
        previousElementType = ONEFLOW::QUAD_4;
    }
    else if ( elementType == ONEFLOW::TETRA_4 )
    {
    }
    else if ( elementType == ONEFLOW::TETRA_10 )
    {
        previousElementType = ONEFLOW::TETRA_4;
    }
    else if ( elementType == ONEFLOW::PYRA_5 )
    {
    }
    else if ( elementType == ONEFLOW::PYRA_14 )
    {
        previousElementType = ONEFLOW::PYRA_5;
    }
    else if ( elementType == ONEFLOW::PENTA_6 )
    {
    }
    else if ( elementType == ONEFLOW::PENTA_15 )
    {
    }
    else if ( elementType == ONEFLOW::PENTA_18 )
    {
        previousElementType = ONEFLOW::PENTA_6;
    }
    else if ( elementType == ONEFLOW::HEXA_8 )
    {
    }
    else if ( elementType == ONEFLOW::HEXA_20 )
    {
    }
    else if ( elementType == ONEFLOW::HEXA_27 )
    {
        previousElementType = ONEFLOW::HEXA_8;
    }

    return previousElementType;
}

bool UnitElement::IsUnitElementType( int elementType )
{
    bool results = false;
    if ( elementType == ONEFLOW::NODE )
    {
    }
    else if ( elementType == ONEFLOW::BAR_2 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::BAR_3 )
    {
    }
    else if ( elementType == ONEFLOW::TRI_3 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::TRI_6 )
    {
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::QUAD_8 )
    {
    }
    else if ( elementType == ONEFLOW::QUAD_9 )
    {
    }
    else if ( elementType == ONEFLOW::TETRA_4 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::TETRA_10 )
    {
    }
    else if ( elementType == ONEFLOW::PYRA_5 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::PYRA_14 )
    {
    }
    else if ( elementType == ONEFLOW::PENTA_6 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::PENTA_15 )
    {
    }
    else if ( elementType == ONEFLOW::PENTA_18 )
    {
    }
    else if ( elementType == ONEFLOW::HEXA_8 )
    {
        results = true;
    }
    else if ( elementType == ONEFLOW::HEXA_20 )
    {
    }
    else if ( elementType == ONEFLOW::HEXA_27 )
    {
    }
    return results;
}

bool UnitElement::IsFaceElementType( int elementType )
{
    bool results = false;

    if ( elementType == ONEFLOW::BAR_2 )
    {
        if ( ONEFLOW::IsTwoD() )
        {
            results = true;
        }
        else
        {
            results = false;
        }
    }
    else if ( elementType == ONEFLOW::TRI_3 )
    {
        if ( ONEFLOW::IsTwoD() )
        {
            results = false;
        }
        else
        {
            results = true;
        }
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        if ( ONEFLOW::IsTwoD() )
        {
            results = false;
        }
        else
        {
            results = true;
        }
    }
    return results;
}

bool UnitElement::IsFaceElementTypeAtLeast( int elementType )
{
    bool results = IsUnitElementType( elementType );

    if ( elementType == ONEFLOW::BAR_2 )
    {
        if ( ONEFLOW::IsOneD() )
        {
            results = true;
        }
        else
        {
            results = false;
        }
    }
    return results;
}

bool UnitElement::IsBasicVolumeElementType( int elementType )
{
    if ( ONEFLOW::IsOneD() )
    {
        bool results = false;
        if ( elementType == ONEFLOW::CGNS_ENUMV( BAR_2 ) )
        {
            results = true;
        }
        return results;
    }
    else if ( ONEFLOW::IsTwoD() )
    {
        bool results = false;
        if ( elementType == ONEFLOW::CGNS_ENUMV( TRI_3 ) )
        {
            results = true;
        }
        else if ( elementType == ONEFLOW::CGNS_ENUMV( QUAD_4 ) )
        {
            results = true;
        }
        return results;
    }
    else
    {
        bool results = false;
        if ( elementType == ONEFLOW::TETRA_4 )
        {
            results = true;
        }
        else if ( elementType == ONEFLOW::PYRA_5 )
        {
            results = true;
        }
        else if ( elementType == ONEFLOW::PENTA_6 )
        {
            results = true;
        }
        else if ( elementType == ONEFLOW::HEXA_8 )
        {
            results = true;
        }
        return results;
    }
}

void UnitElement::PushElementFace( int faceType, int p1 )
{
    faceTypeContainer.push_back( faceType );

    IntField face;
    face.push_back( p1 );
    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }
    faceList.push_back( face );
}

void UnitElement::PushElementFace( int faceType, int p1, int p2 )
{
    faceTypeContainer.push_back( faceType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }
    faceList.push_back( face );
}

void UnitElement::PushElementFace( int faceType, int p1, int p2, int p3 )
{
    faceTypeContainer.push_back( faceType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );
    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }
    faceList.push_back( face );
}

void UnitElement::PushElementFace( int faceType, int p1, int p2, int p3, int p4 )
{
    faceTypeContainer.push_back( faceType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );
    face.push_back( p4 );
    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }
    faceList.push_back( face );
}

void UnitElement::PushChildElement( int elementType, int p1, int p2 )
{
    childElementType.push_back( elementType );

    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );

    for ( int iPoint = 0; iPoint < element.size(); ++ iPoint )
    {
        element[ iPoint ] -= 1;
    }

    childElementIndex.push_back( element );
}

void UnitElement::PushChildElement( int elementType, int p1, int p2, int p3 )
{
    childElementType.push_back( elementType );

    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );

    for ( int iPoint = 0; iPoint < element.size(); ++ iPoint )
    {
        element[ iPoint ] -= 1;
    }

    childElementIndex.push_back( element );
}

void UnitElement::PushChildElement( int elementType, int p1, int p2, int p3, int p4 )
{
    childElementType.push_back( elementType );

    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );
    element.push_back( p4 );

    for ( int iPoint = 0; iPoint < element.size(); ++ iPoint )
    {
        element[ iPoint ] -= 1;
    }

    childElementIndex.push_back( element );
}

void UnitElement::PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5 )
{
    childElementType.push_back( elementType );

    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );
    element.push_back( p4 );
    element.push_back( p5 );

    for ( int iPoint = 0; iPoint < element.size(); ++ iPoint )
    {
        element[ iPoint ] -= 1;
    }

    childElementIndex.push_back( element );
}

void UnitElement::PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5, int p6 )
{
    childElementType.push_back( elementType );

    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );
    element.push_back( p4 );
    element.push_back( p5 );
    element.push_back( p6 );

    for ( int iPoint = 0; iPoint < element.size(); ++ iPoint )
    {
        element[ iPoint ] -= 1;
    }

    childElementIndex.push_back( element );
}


IntField & UnitElement::GetElementPhysicsFace(int iFace) { return compositeFaceList[iFace]; };
int UnitElement::GetElementPhysicsFaceType(int iFace) const { return compositeFaceElementType[iFace]; }
int UnitElement::GetElementFaceNumber() const { return faceList.size(); }
int UnitElement::GetElementPhysicsFaceNumber() const { return compositeFaceList.size(); }
int UnitElement::GetElementType() const { return elementType; };


IntField & UnitElement::GetElementFace(int iFace) { return faceList[iFace]; };
int UnitElement::GetFaceType(int iFace) const { return faceTypeContainer[iFace]; };

IntField & UnitElement::GetRelatedPointListForMiddlePointCalcutation(int iMiddlePoint) { return middlePointStruct[iMiddlePoint]; }
IntField & UnitElement::GetChildElementRelativeNodeIndex(int iChildElement) { return childElementIndex[iChildElement]; };

int UnitElement::GetChildElementNumbers() const { return childElementType.size(); };
int UnitElement::GetChildElementType(int iChild) { return childElementType[iChild]; }

void UnitElement::PushChildElement( int elementType, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8 )
{
    childElementType.push_back( elementType );

    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );
    element.push_back( p4 );
    element.push_back( p5 );
    element.push_back( p6 );
    element.push_back( p7 );
    element.push_back( p8 );

    for ( int iPoint = 0; iPoint < element.size(); ++ iPoint )
    {
        element[ iPoint ] -= 1;
    }

    childElementIndex.push_back( element );
}

void UnitElement::PushMiddlePoint( int pm, int p1, int p2 )
{
    middlePointStruct[ pm - 1 ].push_back( p1 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p2 - 1 );
}

void UnitElement::PushMiddlePoint( int pm, int p1, int p2, int p3 )
{
    middlePointStruct[ pm - 1 ].push_back( p1 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p2 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p3 - 1 );
}

void UnitElement::PushMiddlePoint( int pm, int p1, int p2, int p3, int p4 )
{
    middlePointStruct[ pm - 1 ].push_back( p1 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p2 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p3 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p4 - 1 );
}

void UnitElement::PushMiddlePoint( int pm, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8 )
{
    middlePointStruct[ pm - 1 ].push_back( p1 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p2 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p3 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p4 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p5 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p6 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p7 - 1 );
    middlePointStruct[ pm - 1 ].push_back( p8 - 1 );
}

void UnitElement::PushCompositeFace( int faceElementType, int p1, int p2, int p3 )
{
    compositeFaceElementType.push_back( faceElementType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );

    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }

    compositeFaceList.push_back( face );
}

void UnitElement::PushCompositeFace( int faceElementType, int p1, int p2, int p3, int p4, int p5, int p6 )
{
    compositeFaceElementType.push_back( faceElementType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );
    face.push_back( p4 );
    face.push_back( p5 );
    face.push_back( p6 );

    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }

    compositeFaceList.push_back( face );
}

void UnitElement::PushCompositeFace( int faceElementType, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8, int p9 )
{
    compositeFaceElementType.push_back( faceElementType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );
    face.push_back( p4 );
    face.push_back( p5 );
    face.push_back( p6 );
    face.push_back( p7 );
    face.push_back( p8 );
    face.push_back( p9 );

    for ( int iPoint = 0; iPoint < face.size(); ++ iPoint )
    {
        face[ iPoint ] -= 1;
    }

    compositeFaceList.push_back( face );
}

EndNameSpace