#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// --- closure ---
Real SaChi( Real nuTilde, Real invLaminarKinematicViscosity );
Real SaFv1( Real chi, Real cv1Cubed );
Real SaFv2( Real chi, Real fv1 );
Real SaSBar( Real rs, Real fv2 );
Real SaModifiedStrainRate( Real omega, Real sBar, Real c2, Real c3 );
Real SaRRatio( Real rs, Real modifiedStrain );
Real SaG( Real r, Real cw2 );
Real SaFw( Real g, Real cw3ToSix );
Real SaFt2( Real chi, Real ct3, Real ct4, bool enabled );

// --- source residual ---
Real SaProduction( Real cb1, Real ft2, Real modifiedStrain, Real nuTilde );
Real SaDestruction( Real cw1k, Real fw, Real cb1, Real ft2, Real nuTildeRs );
Real SaDiffusion( Real cb2OverSigma, Real gradNuSquared, Real invReynolds );
Real SaLimitedDiffusion( Real diffusion, Real rprod, Real production );
Real SaSourceResidual( Real production, Real destruction, Real diffusion );
Real SaNegativePart( Real x );

// --- eddy viscosity ---
Real SaEddyViscosity( Real density, Real nuTilde, Real fv1 );
Real SaCappedEddyViscosity(
    Real eddyViscosity, Real laminarViscosity, Real maxViscosityRatio );


// --- SA production strain measure ---
// @doc-algorithm sa-prod-strain
// @doc-stage classical
// @doc-pure
// @doc-formula S = sqrt(|2 work2 + work3 - 2/3 work1^2|)
Real SaStrainProductionClassical( Real work1, Real work2, Real work3 );

// @doc-algorithm sa-prod-strain
// @doc-stage rotational-v4
// @doc-pure
// @doc-formula S = Omega - 2 max(0, Omega - sqrt(|stp|))
Real SaStrainProductionRotationalV4( Real omega, Real strainProductionTensor );

// @doc-algorithm sa-prod-strain
// @doc-stage rotational-v5
// @doc-pure
Real SaStrainProductionRotationalV5( Real omega, Real absStrainProductionTensor );

EndNameSpace
