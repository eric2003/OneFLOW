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

EndNameSpace
