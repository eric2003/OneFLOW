#include "Vec1d.h"

void VecWrap::Allocate( int nEqu, int ist, int ied, double value )
{
    this->data.resize( nEqu );
    for ( int m = 0; m < nEqu; ++ m )
    {
        this->data[ m ].Allocate( ist, ied, value );
    }
}
