!* ------------------------------------------------------------------------- *
!* CGNS - CFD General Notation System (http://www.cgns.org)                  *
!* CGNS/MLL - Mid-Level Library header file                                  *
!* Please see cgnsconfig.h file for this local installation configuration    *
!* ------------------------------------------------------------------------- *
!
!* ------------------------------------------------------------------------- *
!
! This software is provided 'as-is', without any express or implied warranty.
! In no event will the authors be held liable for any damages arising from
! the use of this software.
!
! Permission is granted to anyone to use this software for any purpose,
! including commercial applications, and to alter it and redistribute it
! freely, subject to the following restrictions:
!
! 1. The origin of this software must not be misrepresented; you must not
!    claim that you wrote the original software. If you use this software
!    in a product, an acknowledgment in the product documentation would be
!    appreciated but is not required.
!
! 2. Altered source versions must be plainly marked as such, and must not
!    be misrepresented as being the original software.
!
! 3. This notice may not be removed or altered from any source distribution.
!
!* ------------------------------------------------------------------------- *

#ifndef CGNSTYPES_F03_H
#define CGNSTYPES_F03_H

#define CG_BUILD_64BIT_F 0
#define HAVE_FORTRAN_95 0
#define HAVE_FORTRAN_2003 0
#define HAVE_FORTRAN_2008TS 0
#define HDF5_HAVE_MULTI_DATASETS 0

#define CG_BUILD_SCOPE  1

#define CGNS_BASESCOPE  0

#if CG_BUILD_SCOPE
# ifndef CGNS_SCOPE_ENUMS
#  define CGNS_SCOPE_ENUMS
# endif
#else
# ifdef CGNS_SCOPE_ENUMS
#  undef CGNS_SCOPE_ENUMS
# endif
#endif

#ifdef NO_CONCATENATION
# define IDENTITY(x)      x
# define CONCATENATE(a,b) IDENTITY(a)b
#else
# define CONCATENATE(a,b) a##b
#endif

!#if  __STDC__
!#define CONCATENATE(a,b) a##b
!#else
!#define CONCATENATE(a,b) a/**/b
!#endif

#ifdef CGNS_SCOPE_ENUMS
! set scope prefix (case does not matter in Fortran)
#define CGNS_ENUMV(x) CONCATENATE(CG_,x)
#else
#define CGNS_ENUMV(x) x
#endif

#endif
