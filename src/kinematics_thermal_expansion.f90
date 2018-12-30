!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from thermal expansion
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_thermal_expansion
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private

 public :: &
   kinematics_thermal_expansion_init, &
   kinematics_thermal_expansion_initialStrain, &
   kinematics_thermal_expansion_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use IO, only: &
   IO_timeStamp
 use material, only: &
   phase_kinematics, &
   KINEMATICS_thermal_expansion_label, &
   KINEMATICS_thermal_expansion_ID
 use config, only: &
   config_phase

 implicit none
 integer(pInt) maxNinstance
 
 write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_thermal_expansion_LABEL//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_kinematics == KINEMATICS_thermal_expansion_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

! ToDo: this subroutine should read in lattice_thermal_expansion. No need to make it a global array

end subroutine kinematics_thermal_expansion_init

!--------------------------------------------------------------------------------------------------
!> @brief  report initial thermal strain based on current temperature deviation from reference
!--------------------------------------------------------------------------------------------------
pure function kinematics_thermal_expansion_initialStrain(ipc, ip, el)
 use material, only: &
   material_phase, &
   material_homog, &
   temperature, &
   thermalMapping
 use lattice, only: &
   lattice_thermalExpansion33, &
   lattice_referenceTemperature
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   kinematics_thermal_expansion_initialStrain                                                       !< initial thermal strain (should be small strain, though)
 integer(pInt) :: &
   phase, &
   homog, offset
   
 phase = material_phase(ipc,ip,el)
 homog = material_homog(ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 
 kinematics_thermal_expansion_initialStrain = &
   (temperature(homog)%p(offset) - lattice_referenceTemperature(phase))**1 / 1. * &
   lattice_thermalExpansion33(1:3,1:3,1,phase) + &                                                  ! constant  coefficient
   (temperature(homog)%p(offset) - lattice_referenceTemperature(phase))**2 / 2. * &
   lattice_thermalExpansion33(1:3,1:3,2,phase) + &                                                  ! linear    coefficient
   (temperature(homog)%p(offset) - lattice_referenceTemperature(phase))**3 / 3. * &
   lattice_thermalExpansion33(1:3,1:3,3,phase)                                                      ! quadratic coefficient
  
end function kinematics_thermal_expansion_initialStrain

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_LiAndItsTangent(Li, dLi_dTstar, ipc, ip, el)
 use material, only: &
   material_phase, &
   material_homog, &
   temperature, &
   temperatureRate, &
   thermalMapping
 use lattice, only: &
   lattice_thermalExpansion33, &
   lattice_referenceTemperature
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< thermal velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLi_dTstar                                                                                       !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
 integer(pInt) :: &
   phase, &
   homog, offset
 real(pReal) :: &
   T, TRef, TDot  
   
 phase = material_phase(ipc,ip,el)
 homog = material_homog(ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 T = temperature(homog)%p(offset)
 TDot = temperatureRate(homog)%p(offset)
 TRef = lattice_referenceTemperature(phase)
 
 Li = TDot * ( &
               lattice_thermalExpansion33(1:3,1:3,1,phase)*(T - TRef)**0 &                           ! constant  coefficient
             + lattice_thermalExpansion33(1:3,1:3,2,phase)*(T - TRef)**1 &                           ! linear    coefficient
             + lattice_thermalExpansion33(1:3,1:3,3,phase)*(T - TRef)**2 &                           ! quadratic coefficient
             ) / &
      (1.0_pReal &
            + lattice_thermalExpansion33(1:3,1:3,1,phase)*(T - TRef)**1 / 1. &
            + lattice_thermalExpansion33(1:3,1:3,2,phase)*(T - TRef)**2 / 2. &
            + lattice_thermalExpansion33(1:3,1:3,3,phase)*(T - TRef)**3 / 3. &
      )
 dLi_dTstar = 0.0_pReal 
  
end subroutine kinematics_thermal_expansion_LiAndItsTangent

end module kinematics_thermal_expansion
