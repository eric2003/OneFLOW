#pragma once

#include "EulerDomain.h"
#include "HXArray.h"

#include <vector>

BeginNameSpace( ONEFLOW )

// Packs the internal-cell portion of the OneFLOW equation-major MRField
// into the contiguous view expected by EulerDomainBackend. Boundary-cell
// storage remains metadata on EulerDomainProblem::nGhostCells and is not
// uploaded by this snapshot.
class EulerDomainMrFieldSnapshot
{
public:
    EulerDomainMrFieldSnapshot( const MRField & field, int nCells );

    int NCells() const { return nCells_; }
    int NEquations() const { return nEquations_; }
    const std::vector< Real > & Values() const { return values_; }
    EulerDomainConstFieldView View() const;

private:
    int nCells_ = 0;
    int nEquations_ = 0;
    std::vector< Real > values_;
};

EndNameSpace
