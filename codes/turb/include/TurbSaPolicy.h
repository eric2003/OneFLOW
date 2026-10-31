#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

Real SaChi( Real nuTilde, Real invLaminarKinematicViscosity );
Real SaFv1( Real chi, Real cv1Cubed );
Real SaFv2( Real chi, Real fv1 );
Real SaSBar( Real rs, Real fv2 );
Real SaModifiedStrainRate( Real omega, Real sBar, Real c2, Real c3 );
Real SaRRatio( Real rs, Real modifiedStrain );
Real SaG( Real r, Real cw2 );
Real SaFw( Real g, Real cw3ToSix );
Real SaFt2( Real chi, Real ct3, Real ct4, bool enabled );

// @doc-algorithm sa-source
// @doc-stage production
// @doc-pure
// @doc-formula P = cb1 (1-ft2) S_tilde nu_t
Real SaProduction( Real cb1, Real ft2, Real modifiedStrain, Real nuTilde );

// @doc-algorithm sa-source
// @doc-stage destruction
// @doc-pure
// @doc-formula D = (cw1 fw - cb1 ft2) nu_t rs
Real SaDestruction( Real cw1k, Real fw, Real cb1, Real ft2, Real nuTildeRs );

// @doc-algorithm sa-source
// @doc-stage diffusion
// @doc-pure
// @doc-formula Diff = cb2/sigma * |grad nu_t|^2 / Re
Real SaDiffusion( Real cb2OverSigma, Real gradNuSquared, Real invReynolds );

// @doc-algorithm sa-source
// @doc-stage limit-diffusion
// @doc-pure
Real SaLimitedDiffusion( Real diffusion, Real rprod, Real production );

// @doc-algorithm sa-source
// @doc-stage residual
// @doc-pure
// @doc-formula S = P - D + Diff
Real SaSourceResidual( Real production, Real destruction, Real diffusion );

// @doc-algorithm sa-source
// @doc-stage negative-part
// @doc-pure
// @doc-formula -1/2 (x - |x|)
Real SaNegativePart( Real x );

EndNameSpace
