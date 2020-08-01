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

#ifndef CGNSTYPES_F_H
#define CGNSTYPES_F_H

#define CG_BUILD_64BIT_F 1
!#define HAVE_FORTRAN_2003 0
#define HAVE_FORTRAN_2003 0

#if HAVE_FORTRAN_2003
#  if CG_BUILD_64BIT_F
#   define cgsize_t INTEGER(KIND=C_LONG_LONG)
#   define CGSIZE_T INTEGER(KIND=C_LONG_LONG)
#  else
#   define cgsize_t INTEGER(KIND=C_INT)
#   define CGSIZE_T INTEGER(KIND=C_INT)
#  endif

#  define cglong_t INTEGER(KIND=C_LONG_LONG)
#  define CGLONG_T INTEGER(KIND=C_LONG_LONG)
#  define cgid_t   REAL(KIND=C_DOUBLE)
#  define CGID_T   REAL(KIND=C_DOUBLE)
#else
#  if CG_BUILD_64BIT_F
#   define cgsize_t integer*8
#   define CGSIZE_T integer*8
#  else
#   define cgsize_t integer*4
#   define CGSIZE_T integer*4
#  endif

#  define cglong_t integer*8
#  define CGLONG_T integer*8
#  define cgid_t   real*8
#  define CGID_T   real*8
#endif
#endif

