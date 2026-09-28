#include "reconstructor/basic/FirstOrderReconstructor.h"

namespace cfd {

void FirstOrderReconstructor::reconstruct(const Vector& q, 
                                        Vector& q_face_left,
                                        Vector& q_face_right,
                                        const ComputationalDomain& domain) {
    const int ist = domain.ist();
    const int ied = domain.ied();
    
    for (int i = ist; i <= ied; ++i) {
        int j = i - ist;
        q_face_left[j] = q[i-1];
        q_face_right[j] = q[i];
    }
}

} // namespace cfd