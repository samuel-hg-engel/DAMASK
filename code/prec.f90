! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief setting precision for real and int type, using double precision for real
!--------------------------------------------------------------------------------------------------
#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER<1200
#define LEGACY_COMPILER
#endif
#endif

!--------------------------------------------------------------------------------------------------
module prec


 implicit none
 private 

 integer,     parameter, public :: pReal     = selected_real_kind(15,300)                           !< floating point number with 15 significant digits, up to 1e+-300 (double precision)
 integer,     parameter, public :: pInt      = selected_int_kind(9)                                 !< integer representation with at least up to +- 1e9 (32 bit)
 integer,     parameter, public :: pLongInt  = selected_int_kind(12)                                !< integer representation with at least up to +- 1e12 (64 bit)
 real(pReal), parameter, public :: tol_math_check = 1.0e-8_pReal
 real(pReal), parameter, public :: tol_gravityNodePos = 1.0e-100_pReal
 
! NaN is precision dependent 
! from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html
! copy can be found in documentation/Code/Fortran
#ifdef LEGACY_COMPILER
 real(pReal), parameter, public :: DAMASK_NaN = Z'7FF8000000000000'                                 !< when using old compiler without standard check
#else
 real(pReal), parameter, public :: DAMASK_NaN = real(Z'7FF8000000000000', pReal)                    !< quiet NaN for double precision
#endif

 type, public :: p_vec
   real(pReal), dimension(:), pointer :: p
 end type p_vec

 public :: prec_init
 
contains
!--------------------------------------------------------------------------------------------------
!> @brief reporting precision and checking if DAMASK_NaN is set correctly
!--------------------------------------------------------------------------------------------------
subroutine prec_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 
 implicit none

!$OMP CRITICAL (write2out)
#ifndef LEGACY_COMPILER
 open (6, encoding='UTF-8')  
#endif

 write(6,*)
 write(6,*) '<<<+-  prec init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 write(6,'(a,i3)')    ' Bytes for pReal:    ',pReal
 write(6,'(a,i3)')    ' Bytes for pInt:     ',pInt
 write(6,'(a,i3)')    ' Bytes for pLongInt: ',pLongInt
 write(6,'(a,e10.3)') ' NaN:           ',     DAMASK_NaN
 write(6,'(a,l3)')    ' NaN /= NaN:         ',DAMASK_NaN/=DAMASK_NaN
 if (DAMASK_NaN == DAMASK_NaN) call quit(9000)
 write(6,*)
!$OMP END CRITICAL (write2out)

end subroutine prec_init

end module prec
