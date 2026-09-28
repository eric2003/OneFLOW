#include "EulerDomainMrFieldAdapter.h"

#include <algorithm>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

EulerDomainMrFieldSnapshot::EulerDomainMrFieldSnapshot(
    const MRField & field,
    int nCells )
    : nCells_( nCells )
    , nEquations_( static_cast< int >( field.GetNEqu() ) )
{
    if ( nCells_ <= 0 || ( nEquations_ != 3 && nEquations_ != 5 ) )
    {
        throw std::invalid_argument( "invalid MRField Euler snapshot shape" );
    }

    values_.resize(
        static_cast< std::size_t >( nCells_ )
        * static_cast< std::size_t >( nEquations_ ) );
    for ( int equation = 0; equation < nEquations_; ++ equation )
    {
        const auto & source = field[ equation ];
        if ( static_cast< int >( source.size() ) < nCells_ )
        {
            throw std::invalid_argument(
                "MRField does not contain the requested internal cells" );
        }
        std::copy(
            source.begin(),
            source.begin() + nCells_,
            values_.begin() + equation * nCells_ );
    }
}

EulerDomainConstFieldView EulerDomainMrFieldSnapshot::View() const
{
    return { nCells_, nEquations_, values_.data() };
}

EndNameSpace
