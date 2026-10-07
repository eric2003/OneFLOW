#include "DataObject.h"

BeginNameSpace( ONEFLOW )

template < typename T >
void TDataObjectDump( std::fstream & file, const std::vector< T > & data )
{
    if ( data.empty() ) return;
    file << data[ 0 ];
    for ( int i = 1; i < static_cast<int>(data.size()); ++ i )
    {
        file << " , ";
        file << data[ i ];
    }
}

EndNameSpace
