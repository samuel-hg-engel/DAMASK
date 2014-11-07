!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module constitutive_dislokmc
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   constitutive_dislokmc_sizePostResults                                                            !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   constitutive_dislokmc_sizePostResult                                                             !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   constitutive_dislokmc_output                                                                     !< name of each post result output
   
 character(len=12),                   dimension(3),           parameter,           private :: &
   CONSTITUTIVE_DISLOKMC_listBasicSlipStates = &
   ['rhoEdge     ',      'rhoEdgeDip  ',      'accshearslip']

 character(len=12),                   dimension(2),           parameter,           private :: &
   CONSTITUTIVE_DISLOKMC_listBasicTwinStates = & 
   ['twinFraction',      'accsheartwin']

 character(len=17),                   dimension(4),           parameter,           private :: &
   CONSTITUTIVE_DISLOKMC_listDependentSlipStates = &
   ['invLambdaSlip    ', 'invLambdaSlipTwin', 'meanFreePathSlip ', 'tauSlipThreshold ']

 character(len=16),                   dimension(4),           parameter,           private :: &
   CONSTITUTIVE_DISLOKMC_listDependentTwinStates = & 
   ['invLambdaTwin   ',  'meanFreePathTwin',  'tauTwinThreshold',  'twinVolume      ']

 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin

 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   constitutive_dislokmc_Noutput                                                                    !< number of outputs per instance of this plasticity 

 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   constitutive_dislokmc_totalNslip, &                                                              !< total number of active slip systems for each instance
   constitutive_dislokmc_totalNtwin                                                                 !< total number of active twin systems for each instance

 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   constitutive_dislokmc_Nslip, &                                                                   !< number of active slip systems for each family and instance
   constitutive_dislokmc_Ntwin                                                                      !< number of active twin systems for each family and instance

 real(pReal),                         dimension(:),           allocatable,         private :: &
   constitutive_dislokmc_CAtomicVolume, &                                                           !< atomic volume in Bugers vector unit
   constitutive_dislokmc_D0, &                                                                      !< prefactor for self-diffusion coefficient
   constitutive_dislokmc_Qsd, &                                                                     !< activation energy for dislocation climb
   constitutive_dislokmc_GrainSize, &                                                               !< grain size
   constitutive_dislokmc_MaxTwinFraction, &                                                         !< maximum allowed total twin volume fraction
   constitutive_dislokmc_CEdgeDipMinDistance, &                                                     !<
   constitutive_dislokmc_Cmfptwin, &                                                                !<
   constitutive_dislokmc_Cthresholdtwin, &                                                          !<
   constitutive_dislokmc_SolidSolutionStrength, &                                                   !< Strength due to elements in solid solution
   constitutive_dislokmc_L0, &                                                                      !< Length of twin nuclei in Burgers vectors
   constitutive_dislokmc_xc, &                                                                      !< critical distance for formation of twin nucleus
   constitutive_dislokmc_VcrossSlip, &                                                              !< cross slip volume
   constitutive_dislokmc_SFE_0K, &                                                                  !< stacking fault energy at zero K
   constitutive_dislokmc_dSFE_dT, &                                                                 !< temperature dependance of stacking fault energy
   constitutive_dislokmc_dipoleFormationFactor, &                                                   !< scaling factor for dipole formation: 0: off, 1: on. other values not useful
   constitutive_dislokmc_aTolRho, &                                                                 !< absolute tolerance for integration of dislocation density
   constitutive_dislokmc_aTolTwinFrac                                                               !< absolute tolerance for integration of twin volume fraction

 real(pReal),                         dimension(:,:,:,:),     allocatable,         private :: &
   constitutive_dislokmc_Ctwin66                                                                    !< twin elasticity matrix in Mandel notation for each instance
 real(pReal),                         dimension(:,:,:,:,:,:), allocatable,         private :: &
   constitutive_dislokmc_Ctwin3333                                                                  !< twin elasticity matrix for each instance
 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   constitutive_dislokmc_rhoEdge0, &                                                                !< initial edge dislocation density per slip system for each family and instance
   constitutive_dislokmc_rhoEdgeDip0, &                                                             !< initial edge dipole density per slip system for each family and instance
   constitutive_dislokmc_burgersPerSlipFamily, &                                                    !< absolute length of burgers vector [m] for each slip family and instance
   constitutive_dislokmc_burgersPerSlipSystem, &                                                    !< absolute length of burgers vector [m] for each slip system and instance
   constitutive_dislokmc_burgersPerTwinFamily, &                                                    !< absolute length of burgers vector [m] for each twin family and instance
   constitutive_dislokmc_burgersPerTwinSystem, &                                                    !< absolute length of burgers vector [m] for each twin system and instance
   constitutive_dislokmc_QedgePerSlipFamily, &                                                      !< activation energy for glide [J] for each slip family and instance
   constitutive_dislokmc_QedgePerSlipSystem, &                                                      !< activation energy for glide [J] for each slip system and instance
   constitutive_dislokmc_v0PerSlipFamily, &                                                         !< dislocation velocity prefactor [m/s] for each family and instance
   constitutive_dislokmc_v0PerSlipSystem, &                                                         !< dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_dislokmc_tau_peierlsPerSlipFamily, &                                                !< Peierls stress [Pa] for each family and instance
   constitutive_dislokmc_Ndot0PerTwinFamily, &                                                      !< twin nucleation rate [1/m³s] for each twin family and instance
   constitutive_dislokmc_Ndot0PerTwinSystem, &                                                      !< twin nucleation rate [1/m³s] for each twin system and instance
   constitutive_dislokmc_tau_r, &                                                                   !< stress to bring partial close together for each twin system and instance
   constitutive_dislokmc_twinsizePerTwinFamily, &                                                   !< twin thickness [m] for each twin family and instance
   constitutive_dislokmc_twinsizePerTwinSystem, &                                                   !< twin thickness [m] for each twin system and instance
   constitutive_dislokmc_CLambdaSlipPerSlipFamily, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_dislokmc_CLambdaSlipPerSlipSystem, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_dislokmc_interaction_SlipSlip, &                                                    !< coefficients for slip-slip interaction for each interaction type and instance
   constitutive_dislokmc_interaction_SlipTwin, &                                                    !< coefficients for slip-twin interaction for each interaction type and instance
   constitutive_dislokmc_interaction_TwinSlip, &                                                    !< coefficients for twin-slip interaction for each interaction type and instance
   constitutive_dislokmc_interaction_TwinTwin, &                                                    !< coefficients for twin-twin interaction for each interaction type and instance
   constitutive_dislokmc_pPerSlipFamily, &                                                          !< p-exponent in glide velocity
   constitutive_dislokmc_qPerSlipFamily, &                                                          !< q-exponent in glide velocity
   constitutive_dislokmc_uPerSlipFamily, &                                                          !< u-exponent in glide velocity           
   constitutive_dislokmc_sPerSlipFamily, &                                                          !< self-hardening in glide velocity
   constitutive_dislokmc_rPerTwinFamily, &                                                          !< r-exponent in twin nucleation rate
   constitutive_dislokmc_nonSchmidCoeff                                                             !< non-Schmid coefficients (bcc)
 real(pReal),                         dimension(:,:,:),       allocatable,         private :: &
   constitutive_dislokmc_interactionMatrix_SlipSlip, &                                              !< interaction matrix of the different slip systems for each instance
   constitutive_dislokmc_interactionMatrix_SlipTwin, &                                              !< interaction matrix of slip systems with twin systems for each instance
   constitutive_dislokmc_interactionMatrix_TwinSlip, &                                              !< interaction matrix of twin systems with slip systems for each instance
   constitutive_dislokmc_interactionMatrix_TwinTwin, &                                              !< interaction matrix of the different twin systems for each instance
   constitutive_dislokmc_forestProjectionEdge                                                       !< matrix of forest projections of edge dislocations for each instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 edge_density_ID, &
                 dipole_density_ID, &
                 shear_rate_slip_ID, &
                 accumulated_shear_slip_ID, &
                 mfp_slip_ID, &
                 resolved_stress_slip_ID, &
                 threshold_stress_slip_ID, &
                 edge_dipole_distance_ID, &
                 stress_exponent_ID, &
                 twin_fraction_ID, &
                 shear_rate_twin_ID, &
                 accumulated_shear_twin_ID, &
                 mfp_twin_ID, &
                 resolved_stress_twin_ID, &
                 threshold_stress_twin_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   constitutive_dislokmc_outputID                                                                  !< ID of each post result output


 public :: &
   constitutive_dislokmc_init, &
   constitutive_dislokmc_homogenizedC, &
   constitutive_dislokmc_microstructure, &
   constitutive_dislokmc_LpAndItsTangent, &
   constitutive_dislokmc_dotState, &
   constitutive_dislokmc_getAccumulatedSlip, &
   constitutive_dislokmc_getSlipRate, &
   constitutive_dislokmc_postResults
 private :: &
   constitutive_dislokmc_stateInit, &
   constitutive_dislokmc_aTolState

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
   math_mul3x3
 use mesh, only: &
   mesh_NcpElems
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_DISLOKMC_label, &
   PLASTICITY_DISLOKMC_ID, &
   material_phase, &  
   plasticState, & 
   MATERIAL_partPhase
 use lattice
 use numerics,only: &
   worldrank, &
   numerics_integrator
 
 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = LATTICE_maxNinteraction + 1_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,maxTotalNslip,maxTotalNtwin,&
                  f,instance,j,k,l,m,n,o,p,q,r,s,ns,nt, &
                  Nchunks_SlipSlip, Nchunks_SlipTwin, Nchunks_TwinSlip, Nchunks_TwinTwin, &
                  Nchunks_SlipFamilies, Nchunks_TwinFamilies, Nchunks_nonSchmid, &
                  index_myFamily, index_otherFamily
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase
 character(len=65536) :: &
   tag  = '', &
   line = ''
 real(pReal), dimension(:), allocatable :: tempPerSlip, tempPerTwin
  
 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOKMC_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_DISLOKMC_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(constitutive_dislokmc_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(constitutive_dislokmc_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(constitutive_dislokmc_output(maxval(phase_Noutput),maxNinstance))
          constitutive_dislokmc_output = ''
 allocate(constitutive_dislokmc_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(constitutive_dislokmc_Noutput(maxNinstance),                             source=0_pInt)
 allocate(constitutive_dislokmc_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(constitutive_dislokmc_Ntwin(lattice_maxNtwinFamily,maxNinstance),        source=0_pInt)
 allocate(constitutive_dislokmc_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(constitutive_dislokmc_totalNtwin(maxNinstance),                          source=0_pInt)
 allocate(constitutive_dislokmc_CAtomicVolume(maxNinstance),                       source=0.0_pReal)
 allocate(constitutive_dislokmc_D0(maxNinstance),                                  source=0.0_pReal)
 allocate(constitutive_dislokmc_Qsd(maxNinstance),                                 source=0.0_pReal)
 allocate(constitutive_dislokmc_GrainSize(maxNinstance),                           source=0.0_pReal)
 allocate(constitutive_dislokmc_MaxTwinFraction(maxNinstance),                     source=0.0_pReal)
 allocate(constitutive_dislokmc_CEdgeDipMinDistance(maxNinstance),                 source=0.0_pReal)
 allocate(constitutive_dislokmc_Cmfptwin(maxNinstance),                            source=0.0_pReal)
 allocate(constitutive_dislokmc_Cthresholdtwin(maxNinstance),                      source=0.0_pReal)
 allocate(constitutive_dislokmc_SolidSolutionStrength(maxNinstance),               source=0.0_pReal)
 allocate(constitutive_dislokmc_L0(maxNinstance),                                  source=0.0_pReal)
 allocate(constitutive_dislokmc_xc(maxNinstance),                                  source=0.0_pReal)
 allocate(constitutive_dislokmc_VcrossSlip(maxNinstance),                          source=0.0_pReal)
 allocate(constitutive_dislokmc_aTolRho(maxNinstance),                             source=0.0_pReal)
 allocate(constitutive_dislokmc_aTolTwinFrac(maxNinstance),                        source=0.0_pReal)
 allocate(constitutive_dislokmc_SFE_0K(maxNinstance),                              source=0.0_pReal)
 allocate(constitutive_dislokmc_dSFE_dT(maxNinstance),                             source=0.0_pReal)
 allocate(constitutive_dislokmc_dipoleFormationFactor(maxNinstance),               source=1.0_pReal) !should be on by default
 allocate(constitutive_dislokmc_rhoEdge0(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)
 allocate(constitutive_dislokmc_rhoEdgeDip0(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal)
 allocate(constitutive_dislokmc_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_burgersPerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_QedgePerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_v0PerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_tau_peierlsPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_pPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(constitutive_dislokmc_qPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(constitutive_dislokmc_uPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(constitutive_dislokmc_sPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(constitutive_dislokmc_Ndot0PerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_twinsizePerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_CLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_rPerTwinFamily(lattice_maxNtwinFamily,maxNinstance),source=0.0_pReal)
 allocate(constitutive_dislokmc_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_interaction_SlipTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_interaction_TwinSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_interaction_TwinTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislokmc_nonSchmidCoeff(lattice_maxNnonSchmid,maxNinstance), &
                                                                                  source=0.0_pReal)
 

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     if (phase_plasticity(phase) == PLASTICITY_DISLOKMC_ID) then
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)
       Nchunks_TwinFamilies = count(lattice_NtwinSystem(:,phase) > 0_pInt)
       Nchunks_SlipSlip =     maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_SlipTwin =     maxval(lattice_interactionSlipTwin(:,:,phase))
       Nchunks_TwinSlip =     maxval(lattice_interactionTwinSlip(:,:,phase))
       Nchunks_TwinTwin =     maxval(lattice_interactionTwinTwin(:,:,phase))
       Nchunks_nonSchmid =    lattice_NnonSchmid(phase)
       if(allocated(tempPerSlip)) deallocate(tempPerSlip)
       if(allocated(tempPerTwin)) deallocate(tempPerTwin)
       allocate(tempPerSlip(Nchunks_SlipFamilies))
       allocate(tempPerTwin(Nchunks_TwinFamilies))
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_DISLOKMC_ID) then           ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
      select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('edge_density')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = edge_density_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('dipole_density')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = dipole_density_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('shear_rate_slip')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = shear_rate_slip_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('accumulated_shear_slip')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = accumulated_shear_slip_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('mfp_slip')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = mfp_slip_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('resolved_stress_slip')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = resolved_stress_slip_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('edge_dipole_distance')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = edge_dipole_distance_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('stress_exponent')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = stress_exponent_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('twin_fraction')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = twin_fraction_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('shear_rate_twin')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = shear_rate_twin_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('accumulated_shear_twin')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = accumulated_shear_twin_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('mfp_twin')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = mfp_twin_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('resolved_stress_twin')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = resolved_stress_twin_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('threshold_stress_twin')
             constitutive_dislokmc_Noutput(instance) = constitutive_dislokmc_Noutput(instance) + 1_pInt
             constitutive_dislokmc_outputID(constitutive_dislokmc_Noutput(instance),instance) = threshold_stress_twin_ID
             constitutive_dislokmc_output(constitutive_dislokmc_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of slip system families
       case ('nslip')
         if (positions(1) < Nchunks_SlipFamilies + 1_pInt) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         if (positions(1) > Nchunks_SlipFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_dislokmc_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('rhoedge0','rhoedgedip0','slipburgers','qedge','v0','clambdaslip','tau_peierls','p_slip','q_slip',&
            'u_slip','v_slip','s_slip')
         do j = 1_pInt, Nchunks_SlipFamilies
           tempPerSlip(j) = IO_floatValue(line,positions,1_pInt+j)
         enddo
         select case(tag)
           case ('rhoedge0')
             constitutive_dislokmc_rhoEdge0(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('rhoedgedip0')
             constitutive_dislokmc_rhoEdgeDip0(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('slipburgers')
             constitutive_dislokmc_burgersPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('qedge')
             constitutive_dislokmc_QedgePerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('v0')
             constitutive_dislokmc_v0PerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('clambdaslip')
             constitutive_dislokmc_CLambdaSlipPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('tau_peierls')
             if (lattice_structure(phase) /= LATTICE_bcc_ID) &
               call IO_warning(42_pInt,ext_msg=trim(tag)//' for non-bcc ('//PLASTICITY_DISLOKMC_label//')')
             constitutive_dislokmc_tau_peierlsPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('p_slip')
             constitutive_dislokmc_pPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('q_slip')
             constitutive_dislokmc_qPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('u_slip')
             constitutive_dislokmc_uPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('s_slip')
             constitutive_dislokmc_sPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)   
         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on slip number of twin families
       case ('ntwin')
         if (positions(1) < Nchunks_TwinFamilies + 1_pInt) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         if (positions(1) > Nchunks_TwinFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         Nchunks_TwinFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_dislokmc_Ntwin(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('ndot0','twinsize','twinburgers','r_twin')
         do j = 1_pInt, Nchunks_TwinFamilies
           tempPerTwin(j) = IO_floatValue(line,positions,1_pInt+j)
         enddo
         select case(tag)
           case ('ndot0')
             if (lattice_structure(phase) == LATTICE_fcc_ID) &
               call IO_warning(42_pInt,ext_msg=trim(tag)//' for fcc ('//PLASTICITY_DISLOKMC_label//')')
             constitutive_dislokmc_Ndot0PerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
           case ('twinsize')
             constitutive_dislokmc_twinsizePerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
           case ('twinburgers')
             constitutive_dislokmc_burgersPerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
           case ('r_twin')
             constitutive_dislokmc_rPerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of interactions
       case ('interaction_slipslip','interactionslipslip')
         if (positions(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           constitutive_dislokmc_interaction_SlipSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_sliptwin','interactionsliptwin')
         if (positions(1) < 1_pInt + Nchunks_SlipTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         do j = 1_pInt, Nchunks_SlipTwin
           constitutive_dislokmc_interaction_SlipTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twinslip','interactiontwinslip')
         if (positions(1) < 1_pInt + Nchunks_TwinSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         do j = 1_pInt, Nchunks_TwinSlip
           constitutive_dislokmc_interaction_TwinSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twintwin','interactiontwintwin')
         if (positions(1) < 1_pInt + Nchunks_TwinTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         do j = 1_pInt, Nchunks_TwinTwin
           constitutive_dislokmc_interaction_TwinTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('nonschmid_coefficients')
         if (positions(1) < 1_pInt + Nchunks_nonSchmid) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOKMC_label//')')
         do j = 1_pInt,Nchunks_nonSchmid
           constitutive_dislokmc_nonSchmidCoeff(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
!--------------------------------------------------------------------------------------------------
! parameters independent of number of slip/twin systems
       case ('grainsize')
         constitutive_dislokmc_GrainSize(instance) = IO_floatValue(line,positions,2_pInt)
       case ('maxtwinfraction')
         constitutive_dislokmc_MaxTwinFraction(instance) = IO_floatValue(line,positions,2_pInt)
       case ('d0')
         constitutive_dislokmc_D0(instance) = IO_floatValue(line,positions,2_pInt)
       case ('qsd')
         constitutive_dislokmc_Qsd(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_rho')
         constitutive_dislokmc_aTolRho(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_twinfrac')
         constitutive_dislokmc_aTolTwinFrac(instance) = IO_floatValue(line,positions,2_pInt)
       case ('cmfptwin')
         constitutive_dislokmc_Cmfptwin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('cthresholdtwin')
         constitutive_dislokmc_Cthresholdtwin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('solidsolutionstrength')
         constitutive_dislokmc_SolidSolutionStrength(instance) = IO_floatValue(line,positions,2_pInt)
       case ('l0')
         constitutive_dislokmc_L0(instance) = IO_floatValue(line,positions,2_pInt)
       case ('xc')
              constitutive_dislokmc_xc(instance) = IO_floatValue(line,positions,2_pInt)
       case ('vcrossslip')
              constitutive_dislokmc_VcrossSlip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('cedgedipmindistance')
         constitutive_dislokmc_CEdgeDipMinDistance(instance) = IO_floatValue(line,positions,2_pInt)
       case ('catomicvolume')
         constitutive_dislokmc_CAtomicVolume(instance) = IO_floatValue(line,positions,2_pInt)
       case ('sfe_0k')
         constitutive_dislokmc_SFE_0K(instance) = IO_floatValue(line,positions,2_pInt)
       case ('dsfe_dt')
         constitutive_dislokmc_dSFE_dT(instance) = IO_floatValue(line,positions,2_pInt)
       case ('dipoleformationfactor')
         constitutive_dislokmc_dipoleFormationFactor(instance) = IO_floatValue(line,positions,2_pInt)
     end select
   endif; endif
 enddo parsingFile
 
 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
    myPhase: if (phase_plasticity(phase) == PLASTICITY_dislokmc_ID) then
      instance = phase_plasticityInstance(phase)
      if (sum(constitutive_dislokmc_Nslip(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Nslip ('//PLASTICITY_DISLOKMC_label//')')
      if (sum(constitutive_dislokmc_Ntwin(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Ntwin ('//PLASTICITY_DISLOKMC_label//')')
      do f = 1_pInt,lattice_maxNslipFamily
        if (constitutive_dislokmc_Nslip(f,instance) > 0_pInt) then
          if (constitutive_dislokmc_rhoEdge0(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOKMC_label//')')
          if (constitutive_dislokmc_rhoEdgeDip0(f,instance) < 0.0_pReal) & 
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOKMC_label//')')
          if (constitutive_dislokmc_burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOKMC_label//')')
          if (constitutive_dislokmc_v0PerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOKMC_label//')')
          if (constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='tau_peierls ('//PLASTICITY_DISLOKMC_label//')')
        endif
      enddo
      do f = 1_pInt,lattice_maxNtwinFamily
        if (constitutive_dislokmc_Ntwin(f,instance) > 0_pInt) then
          if (constitutive_dislokmc_burgersPerTwinFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='twinburgers ('//PLASTICITY_DISLOKMC_label//')')
          if (constitutive_dislokmc_Ndot0PerTwinFamily(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='ndot0 ('//PLASTICITY_DISLOKMC_label//')')
        endif
      enddo
      if (constitutive_dislokmc_CAtomicVolume(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOKMC_label//')')
      if (constitutive_dislokmc_D0(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='D0 ('//PLASTICITY_DISLOKMC_label//')')
      if (constitutive_dislokmc_Qsd(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='Qsd ('//PLASTICITY_DISLOKMC_label//')')
      if (sum(constitutive_dislokmc_Ntwin(:,instance)) > 0_pInt) then
        if (constitutive_dislokmc_SFE_0K(instance) == 0.0_pReal .and. &
            constitutive_dislokmc_dSFE_dT(instance) == 0.0_pReal .and. &
            lattice_structure(phase) == LATTICE_fcc_ID) &
          call IO_error(211_pInt,el=instance,ext_msg='SFE0K ('//PLASTICITY_DISLOKMC_label//')')
        if (constitutive_dislokmc_aTolRho(instance) <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_DISLOKMC_label//')')   
        if (constitutive_dislokmc_aTolTwinFrac(instance) <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolTwinFrac ('//PLASTICITY_DISLOKMC_label//')')
      endif
      if (constitutive_dislokmc_dipoleFormationFactor(instance) /= 0.0_pReal .and. &
          constitutive_dislokmc_dipoleFormationFactor(instance) /= 1.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='dipoleFormationFactor ('//PLASTICITY_DISLOKMC_label//')')

!--------------------------------------------------------------------------------------------------
! Determine total number of active slip or twin systems
      constitutive_dislokmc_Nslip(:,instance) = min(lattice_NslipSystem(:,phase),constitutive_dislokmc_Nslip(:,instance))
      constitutive_dislokmc_Ntwin(:,instance) = min(lattice_NtwinSystem(:,phase),constitutive_dislokmc_Ntwin(:,instance))
      constitutive_dislokmc_totalNslip(instance) = sum(constitutive_dislokmc_Nslip(:,instance))
      constitutive_dislokmc_totalNtwin(instance) = sum(constitutive_dislokmc_Ntwin(:,instance))
   endif myPhase
 enddo sanityChecks
 
!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 maxTotalNslip = maxval(constitutive_dislokmc_totalNslip)
 maxTotalNtwin = maxval(constitutive_dislokmc_totalNtwin)
 
 allocate(constitutive_dislokmc_burgersPerSlipSystem(maxTotalNslip, maxNinstance),    source=0.0_pReal)
 allocate(constitutive_dislokmc_burgersPerTwinSystem(maxTotalNtwin, maxNinstance),    source=0.0_pReal)
 allocate(constitutive_dislokmc_QedgePerSlipSystem(maxTotalNslip, maxNinstance),      source=0.0_pReal)
 allocate(constitutive_dislokmc_v0PerSlipSystem(maxTotalNslip, maxNinstance),         source=0.0_pReal)
 allocate(constitutive_dislokmc_Ndot0PerTwinSystem(maxTotalNtwin, maxNinstance),      source=0.0_pReal)
 allocate(constitutive_dislokmc_tau_r(maxTotalNtwin, maxNinstance),                   source=0.0_pReal)
 allocate(constitutive_dislokmc_twinsizePerTwinSystem(maxTotalNtwin, maxNinstance),   source=0.0_pReal)
 allocate(constitutive_dislokmc_CLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance),source=0.0_pReal)
 
 allocate(constitutive_dislokmc_interactionMatrix_SlipSlip(maxval(constitutive_dislokmc_totalNslip),&  ! slip resistance from slip activity
                                                            maxval(constitutive_dislokmc_totalNslip),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(constitutive_dislokmc_interactionMatrix_SlipTwin(maxval(constitutive_dislokmc_totalNslip),&  ! slip resistance from twin activity
                                                            maxval(constitutive_dislokmc_totalNtwin),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(constitutive_dislokmc_interactionMatrix_TwinSlip(maxval(constitutive_dislokmc_totalNtwin),&  ! twin resistance from slip activity
                                                            maxval(constitutive_dislokmc_totalNslip),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(constitutive_dislokmc_interactionMatrix_TwinTwin(maxval(constitutive_dislokmc_totalNtwin),&  ! twin resistance from twin activity
                                                            maxval(constitutive_dislokmc_totalNtwin),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(constitutive_dislokmc_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(constitutive_dislokmc_Ctwin66(6,6,maxTotalNtwin,maxNinstance),              source=0.0_pReal)
 allocate(constitutive_dislokmc_Ctwin3333(3,3,3,3,maxTotalNtwin,maxNinstance),        source=0.0_pReal)

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
    myPhase2: if (phase_plasticity(phase) == PLASTICITY_dislokmc_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_plasticityInstance(phase)
 
     ns = constitutive_dislokmc_totalNslip(instance)
     nt = constitutive_dislokmc_totalNtwin(instance)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputs: do o = 1_pInt,constitutive_dislokmc_Noutput(instance)
       select case(constitutive_dislokmc_outputID(o,instance))
         case(edge_density_ID, &
              dipole_density_ID, &
              shear_rate_slip_ID, &
              accumulated_shear_slip_ID, &
              mfp_slip_ID, &
              resolved_stress_slip_ID, &
              threshold_stress_slip_ID, &
              edge_dipole_distance_ID, &
              stress_exponent_ID &
              )
           mySize = ns
         case(twin_fraction_ID, &
              shear_rate_twin_ID, &
              accumulated_shear_twin_ID, &
              mfp_twin_ID, &
              resolved_stress_twin_ID, &
              threshold_stress_twin_ID &
              )
           mySize = nt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          constitutive_dislokmc_sizePostResult(o,instance) = mySize
          constitutive_dislokmc_sizePostResults(instance)  = constitutive_dislokmc_sizePostResults(instance) + mySize
       endif
     enddo outputs
  
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeDotState              =   int(size(CONSTITUTIVE_DISLOKMC_listBasicSlipStates),pInt) * ns &
                                 + int(size(CONSTITUTIVE_DISLOKMC_listBasicTwinStates),pInt) * nt
     sizeState                 =   sizeDotState &
                                 + int(size(CONSTITUTIVE_DISLOKMC_listDependentSlipStates),pInt) * ns &
                                 + int(size(CONSTITUTIVE_DISLOKMC_listDependentTwinStates),pInt) * nt
                
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     allocate(plasticState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(plasticState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(phase)%deltaState          (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(phase)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)
    !* Process slip related parameters ------------------------------------------------ 
 
     mySlipFamilies: do f = 1_pInt,lattice_maxNslipFamily
       index_myFamily = sum(constitutive_dislokmc_Nslip(1:f-1_pInt,instance))                      ! index in truncated slip system list
       mySlipSystems: do j = 1_pInt,constitutive_dislokmc_Nslip(f,instance)

      !* Burgers vector, 
      !  dislocation velocity prefactor,
      !  mean free path prefactor,
      !  and minimum dipole distance
 
         constitutive_dislokmc_burgersPerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislokmc_burgersPerSlipFamily(f,instance)
 
         constitutive_dislokmc_QedgePerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislokmc_QedgePerSlipFamily(f,instance)
 
         constitutive_dislokmc_v0PerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislokmc_v0PerSlipFamily(f,instance)
 
         constitutive_dislokmc_CLambdaSlipPerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislokmc_CLambdaSlipPerSlipFamily(f,instance)
  
       !* Calculation of forest projections for edge dislocations
       !* Interaction matrices
  
         otherSlipFamilies: do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(constitutive_dislokmc_Nslip(1:o-1_pInt,instance))
           otherSlipSystems: do k = 1_pInt,constitutive_dislokmc_Nslip(o,instance)
             constitutive_dislokmc_forestProjectionEdge(index_myFamily+j,index_otherFamily+k,instance) = &
               abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,phase))+j,phase), &
                               lattice_st(:,sum(lattice_NslipSystem(1:o-1,phase))+k,phase)))
             constitutive_dislokmc_interactionMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislokmc_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                 phase), instance )
         enddo otherSlipSystems; enddo otherSlipFamilies
  
         otherTwinFamilies: do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(constitutive_dislokmc_Ntwin(1:o-1_pInt,instance))
           otherTwinSystems: do k = 1_pInt,constitutive_dislokmc_Ntwin(o,instance)
             constitutive_dislokmc_interactionMatrix_SlipTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislokmc_interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo otherTwinSystems; enddo otherTwinFamilies
  
       enddo mySlipSystems
     enddo mySlipFamilies
  
    !* Process twin related parameters ------------------------------------------------
    
     myTwinFamilies: do f = 1_pInt,lattice_maxNtwinFamily
       index_myFamily = sum(constitutive_dislokmc_Ntwin(1:f-1_pInt,instance))                      ! index in truncated twin system list
       myTwinSystems: do j = 1_pInt,constitutive_dislokmc_Ntwin(f,instance)
 
       !* Burgers vector,
       !  nucleation rate prefactor,
       !  and twin size
 
         constitutive_dislokmc_burgersPerTwinSystem(index_myFamily+j,instance)  = &
         constitutive_dislokmc_burgersPerTwinFamily(f,instance)

         constitutive_dislokmc_Ndot0PerTwinSystem(index_myFamily+j,instance)  = &
         constitutive_dislokmc_Ndot0PerTwinFamily(f,instance)

         constitutive_dislokmc_twinsizePerTwinSystem(index_myFamily+j,instance) = &
         constitutive_dislokmc_twinsizePerTwinFamily(f,instance)
 
       !* Rotate twin elasticity matrices
 
         index_otherFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase))                             ! index in full lattice twin list
         do l = 1_pInt,3_pInt; do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt
           do p = 1_pInt,3_pInt; do q = 1_pInt,3_pInt; do r = 1_pInt,3_pInt; do s = 1_pInt,3_pInt
             constitutive_dislokmc_Ctwin3333(l,m,n,o,index_myFamily+j,instance) = &
             constitutive_dislokmc_Ctwin3333(l,m,n,o,index_myFamily+j,instance) + &
               lattice_C3333(p,q,r,s,instance) * &
               lattice_Qtwin(l,p,index_otherFamily+j,phase) * &
               lattice_Qtwin(m,q,index_otherFamily+j,phase) * &
               lattice_Qtwin(n,r,index_otherFamily+j,phase) * &
               lattice_Qtwin(o,s,index_otherFamily+j,phase)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo; enddo
         constitutive_dislokmc_Ctwin66(1:6,1:6,index_myFamily+j,instance) = &
           math_Mandel3333to66(constitutive_dislokmc_Ctwin3333(1:3,1:3,1:3,1:3,index_myFamily+j,instance))
 
      !* Interaction matrices
         otherSlipFamilies2: do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(constitutive_dislokmc_Nslip(1:o-1_pInt,instance))
           otherSlipSystems2: do k = 1_pInt,constitutive_dislokmc_Nslip(o,instance)
             constitutive_dislokmc_interactionMatrix_TwinSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislokmc_interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo otherSlipSystems2; enddo otherSlipFamilies2
 
         otherTwinFamilies2: do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(constitutive_dislokmc_Ntwin(1:o-1_pInt,instance))
           otherTwinSystems2:  do k = 1_pInt,constitutive_dislokmc_Ntwin(o,instance)
             constitutive_dislokmc_interactionMatrix_TwinTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislokmc_interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo otherTwinSystems2; enddo otherTwinFamilies2
 
       enddo myTwinSystems
     enddo myTwinFamilies
     call constitutive_dislokmc_stateInit(phase,instance)
     call constitutive_dislokmc_aTolState(phase,instance)
   endif myPhase2
 
 enddo initializeInstances
 
end subroutine constitutive_dislokmc_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_stateInit(ph,instance)
 use math, only: &
   pi
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_mu
 use material, only: &
   plasticState

 implicit none
 integer(pInt), intent(in) :: &
   instance, &                                                                                      !< number specifying the instance of the plasticity
   ph 

  real(pReal), dimension(plasticState(ph)%sizeState) :: tempState

 integer(pInt) :: i,j,f,ns,nt, index_myFamily
 real(pReal), dimension(constitutive_dislokmc_totalNslip(instance)) :: &
   rhoEdge0, &
   rhoEdgeDip0, &
   invLambdaSlip0, &
   MeanFreePathSlip0, &
   tauSlipThreshold0
 real(pReal), dimension(constitutive_dislokmc_totalNtwin(instance)) :: &
   MeanFreePathTwin0,TwinVolume0
 tempState = 0.0_pReal
 ns = constitutive_dislokmc_totalNslip(instance)
 nt = constitutive_dislokmc_totalNtwin(instance)

!--------------------------------------------------------------------------------------------------
! initialize basic slip state variables
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily   = sum(constitutive_dislokmc_Nslip(1:f-1_pInt,instance))                        ! index in truncated slip system list
   rhoEdge0(index_myFamily+1_pInt: &
            index_myFamily+constitutive_dislokmc_Nslip(f,instance)) = &
     constitutive_dislokmc_rhoEdge0(f,instance)
   rhoEdgeDip0(index_myFamily+1_pInt: &
               index_myFamily+constitutive_dislokmc_Nslip(f,instance)) = &
     constitutive_dislokmc_rhoEdgeDip0(f,instance)
 enddo
 
 tempState(1_pInt:ns)           = rhoEdge0
 tempState(ns+1_pInt:2_pInt*ns) = rhoEdgeDip0
 
!--------------------------------------------------------------------------------------------------
! initialize dependent slip microstructural variables
 forall (i = 1_pInt:ns) &
   invLambdaSlip0(i) = sqrt(dot_product((rhoEdge0+rhoEdgeDip0),constitutive_dislokmc_forestProjectionEdge(1:ns,i,instance)))/ &
                       constitutive_dislokmc_CLambdaSlipPerSlipSystem(i,instance)
 tempState(3_pInt*ns+2_pInt*nt+1:4_pInt*ns+2_pInt*nt) = invLambdaSlip0
 
 forall (i = 1_pInt:ns) &
   MeanFreePathSlip0(i) = &
     constitutive_dislokmc_GrainSize(instance)/(1.0_pReal+invLambdaSlip0(i)*constitutive_dislokmc_GrainSize(instance))
 tempState(5_pInt*ns+3_pInt*nt+1:6_pInt*ns+3_pInt*nt) = MeanFreePathSlip0
 
 forall (i = 1_pInt:ns) &
   tauSlipThreshold0(i) = &
     lattice_mu(ph)*constitutive_dislokmc_burgersPerSlipSystem(i,instance) * &
     sqrt(dot_product((rhoEdge0+rhoEdgeDip0),constitutive_dislokmc_interactionMatrix_SlipSlip(i,1:ns,instance)))

 tempState(6_pInt*ns+4_pInt*nt+1:7_pInt*ns+4_pInt*nt) = tauSlipThreshold0


 
!--------------------------------------------------------------------------------------------------
! initialize dependent twin microstructural variables
 forall (j = 1_pInt:nt) &
   MeanFreePathTwin0(j) = constitutive_dislokmc_GrainSize(instance)
 tempState(6_pInt*ns+3_pInt*nt+1_pInt:6_pInt*ns+4_pInt*nt) = MeanFreePathTwin0
 
 forall (j = 1_pInt:nt) &
   TwinVolume0(j) = &
     (pi/4.0_pReal)*constitutive_dislokmc_twinsizePerTwinSystem(j,instance)*MeanFreePathTwin0(j)**(2.0_pReal)
 tempState(7_pInt*ns+5_pInt*nt+1_pInt:7_pInt*ns+6_pInt*nt) = TwinVolume0
 
plasticState(ph)%state0 = spread(tempState,2,size(plasticState(ph)%state(1,:)))

end subroutine constitutive_dislokmc_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_aTolState(ph,instance)
 use material, only: &
  plasticState

 implicit none
 integer(pInt), intent(in) ::  &
   ph, &
   instance                                                                                         ! number specifying the current instance of the plasticity
 
 ! Tolerance state for dislocation densities
 plasticState(ph)%aTolState(1_pInt:2_pInt*constitutive_dislokmc_totalNslip(instance)) = &
   constitutive_dislokmc_aTolRho(instance)

 ! Tolerance state for accumulated shear due to slip 
 plasticState(ph)%aTolState(2_pInt*constitutive_dislokmc_totalNslip(instance)+1_pInt: &
                                  3_pInt*constitutive_dislokmc_totalNslip(instance))=1e6_pReal
   
 
 ! Tolerance state for twin volume fraction
 plasticState(ph)%aTolState(3_pInt*constitutive_dislokmc_totalNslip(instance)+1_pInt: &
                                  3_pInt*constitutive_dislokmc_totalNslip(instance)+&
                                   constitutive_dislokmc_totalNtwin(instance)) = &
   constitutive_dislokmc_aTolTwinFrac(instance)

! Tolerance state for accumulated shear due to twin
 plasticState(ph)%aTolState(3_pInt*constitutive_dislokmc_totalNslip(instance)+ &
                                  constitutive_dislokmc_totalNtwin(instance)+1_pInt: &
                                  3_pInt*constitutive_dislokmc_totalNslip(instance)+ &
                                  2_pInt*constitutive_dislokmc_totalNtwin(instance)) = 1e6_pReal

end subroutine constitutive_dislokmc_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
function constitutive_dislokmc_homogenizedC(ipc,ip,el)
 use material, only: &
  homogenization_maxNgrains, &
  phase_plasticityInstance, &
  plasticState, &
  mappingConstitutive
 use lattice, only: &
  lattice_C66
 
  implicit none
  real(pReal), dimension(6,6) :: &
    constitutive_dislokmc_homogenizedC
  integer(pInt), intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

 integer(pInt) :: instance,ns,nt,i, &
                  ph, &
                  of
 real(pReal) :: sumf

 !* Shortened notation
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = constitutive_dislokmc_totalNslip(instance)
 nt = constitutive_dislokmc_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(plasticState(ph)%state((3_pInt*ns+1_pInt):(3_pInt*ns+nt),of))           ! safe for nt == 0
 !* Homogenized elasticity matrix
 constitutive_dislokmc_homogenizedC = (1.0_pReal-sumf)*lattice_C66(1:6,1:6,ph)
 do i=1_pInt,nt
    constitutive_dislokmc_homogenizedC = constitutive_dislokmc_homogenizedC &
                   + plasticState(ph)%state(3_pInt*ns+i, of)*constitutive_dislokmc_Ctwin66(1:6,1:6,i,instance)
 enddo 
 
 end function constitutive_dislokmc_homogenizedC
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_microstructure(temperature,ipc,ip,el)
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &
   mappingConstitutive
 use lattice, only: &
   lattice_mu, &
   lattice_nu

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in) :: &
   temperature                                                                                      !< temperature at IP 

 integer(pInt) :: &
   instance, &
   ns,nt,s,t, &
   ph, &
   of
 real(pReal) :: &
   sumf,sfe,x0
 real(pReal), dimension(constitutive_dislokmc_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: fOverStacksize

 !* Shortened notation
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = constitutive_dislokmc_totalNslip(instance)
 nt = constitutive_dislokmc_totalNtwin(instance)
 !* State: 1           :  ns         rho_edge
 !* State: ns+1        :  2*ns       rho_dipole
 !* State: 2*ns+1      :  3*ns       accumulated shear due to slip
 !* State: 3*ns+1      :  3*ns+nt    f
 !* State: 3*ns+nt+1   :  3*ns+2*nt  accumulated shear due to twin
 !* State: 3*ns+2*nt+1 :  4*ns+2*nt  1/lambda_slip
 !* State: 4*ns+2*nt+1 :  5*ns+2*nt  1/lambda_sliptwin
 !* State: 5*ns+2*nt+1 :  5*ns+3*nt  1/lambda_twin
 !* State: 5*ns+3*nt+1 :  6*ns+3*nt  mfp_slip
 !* State: 6*ns+3*nt+1 :  6*ns+4*nt  mfp_twin
 !* State: 6*ns+4*nt+1 :  7*ns+4*nt  threshold_stress_slip
 !* State: 7*ns+4*nt+1 :  7*ns+5*nt  threshold_stress_twin
 !* State: 7*ns+5*nt+1 :  7*ns+6*nt  twin volume
 
 !* Total twin volume fraction
 sumf = sum(plasticState(ph)%state((3*ns+1):(3*ns+nt), of)) ! safe for nt == 0
 
 !* Stacking fault energy
 sfe = constitutive_dislokmc_SFE_0K(instance) + & 
       constitutive_dislokmc_dSFE_dT(instance) * Temperature
 
 !* rescaled twin volume fraction for topology
 forall (t = 1_pInt:nt) &
   fOverStacksize(t) = &
     plasticState(ph)%state(3_pInt*ns+t, of)/constitutive_dislokmc_twinsizePerTwinSystem(t,instance)
 
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (s = 1_pInt:ns) &
   plasticState(ph)%state(3_pInt*ns+2_pInt*nt+s, of) = &
     sqrt(dot_product((plasticState(ph)%state(1:ns,of)+plasticState(ph)%state(ns+1_pInt:2_pInt*ns,of)),&
                      constitutive_dislokmc_forestProjectionEdge(1:ns,s,instance)))/ &
     constitutive_dislokmc_CLambdaSlipPerSlipSystem(s,instance)
 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 plasticState(ph)%state((4_pInt*ns+2_pInt*nt+1_pInt):(5_pInt*ns+2_pInt*nt), of) = 0.0_pReal
 if (nt > 0_pInt .and. ns > 0_pInt) &
   plasticState(ph)%state((4_pInt*ns+2_pInt*nt+1):(5_pInt*ns+2_pInt*nt), of) = &
     matmul(constitutive_dislokmc_interactionMatrix_SlipTwin(1:ns,1:nt,instance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
 !$OMP END CRITICAL (evilmatmul)
 
 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
 !$OMP CRITICAL (evilmatmul)
 if (nt > 0_pInt) &
   plasticState(ph)%state((5_pInt*ns+2_pInt*nt+1_pInt):(5_pInt*ns+3_pInt*nt), of) = &
     matmul(constitutive_dislokmc_interactionMatrix_TwinTwin(1:nt,1:nt,instance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
 !$OMP END CRITICAL (evilmatmul)
 
 !* mean free path between 2 obstacles seen by a moving dislocation
 do s = 1_pInt,ns
    if (nt > 0_pInt) then
       plasticState(ph)%state(5_pInt*ns+3_pInt*nt+s, of) = &
         constitutive_dislokmc_GrainSize(instance)/(1.0_pReal+constitutive_dislokmc_GrainSize(instance)*&
         (plasticState(ph)%state(3_pInt*ns+2_pInt*nt+s, of)+plasticState(ph)%state(4_pInt*ns+2_pInt*nt+s, of)))
    else
       plasticState(ph)%state(5_pInt*ns+s, of) = &
         constitutive_dislokmc_GrainSize(instance)/&
         (1.0_pReal+constitutive_dislokmc_GrainSize(instance)*(plasticState(ph)%state(3_pInt*ns+s, of)))
    endif
 enddo

 !* mean free path between 2 obstacles seen by a growing twin
 forall (t = 1_pInt:nt) &
   plasticState(ph)%state(6_pInt*ns+3_pInt*nt+t, of) = &
     (constitutive_dislokmc_Cmfptwin(instance)*constitutive_dislokmc_GrainSize(instance))/&
     (1.0_pReal+constitutive_dislokmc_GrainSize(instance)*plasticState(ph)%state(5_pInt*ns+2_pInt*nt+t, of))
 
 !* threshold stress for dislocation motion
 forall (s = 1_pInt:ns) &
   plasticState(ph)%state(6_pInt*ns+4_pInt*nt+s, of) = &
     lattice_mu(ph)*constitutive_dislokmc_burgersPerSlipSystem(s,instance)*&
     sqrt(dot_product((plasticState(ph)%state(1:ns, of)+plasticState(ph)%state(ns+1_pInt:2_pInt*ns, of)),&
                      constitutive_dislokmc_interactionMatrix_SlipSlip(s,1:ns,instance)))
 
 !* threshold stress for growing twin
 forall (t = 1_pInt:nt) &
   plasticState(ph)%state(7_pInt*ns+4_pInt*nt+t, of) = &
     constitutive_dislokmc_Cthresholdtwin(instance)*&
     (sfe/(3.0_pReal*constitutive_dislokmc_burgersPerTwinSystem(t,instance))+&
     3.0_pReal*constitutive_dislokmc_burgersPerTwinSystem(t,instance)*lattice_mu(ph)/&
     (constitutive_dislokmc_L0(instance)*constitutive_dislokmc_burgersPerSlipSystem(t,instance)))
 
 !* final twin volume after growth
 forall (t = 1_pInt:nt) &
   plasticState(ph)%state(7_pInt*ns+5_pInt*nt+t, of) = &
     (pi/4.0_pReal)*constitutive_dislokmc_twinsizePerTwinSystem(t,instance)*plasticState(ph)%state(6*ns+3*nt+t, of)**(2.0_pReal)

 !* equilibrium seperation of partial dislocations
 do t = 1_pInt,nt
   x0 = lattice_mu(ph)*constitutive_dislokmc_burgersPerTwinSystem(t,instance)**(2.0_pReal)/&
     (sfe*8.0_pReal*pi)*(2.0_pReal+lattice_nu(ph))/(1.0_pReal-lattice_nu(ph))
   constitutive_dislokmc_tau_r(t,instance)= &
        lattice_mu(ph)*constitutive_dislokmc_burgersPerTwinSystem(t,instance)/(2.0_pReal*pi)*&     
        (1/(x0+constitutive_dislokmc_xc(instance))+cos(pi/3.0_pReal)/x0)                              !!! used where??
 enddo

end subroutine constitutive_dislokmc_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,Temperature,slipDamage,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_spectralDecompositionSym33, &
   math_tensorproduct, &
   math_symmetric33, &
   math_mul33x3
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &
   mappingConstitutive
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily,&
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid, &
   lattice_shearTwin, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID
 
 implicit none
 integer(pInt), intent(in)                  :: ipc,ip,el
 real(pReal), intent(in)                    :: Temperature
 real(pReal), dimension(6),   intent(in)    :: Tstar_v
 real(pReal), &
 dimension(constitutive_dislokmc_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))), &
 intent(in) :: &
   slipDamage
 real(pReal), dimension(3,3), intent(out)   :: Lp
 real(pReal), dimension(9,9), intent(out)   :: dLp_dTstar99

 integer(pInt) :: instance,ph,of,ns,nt,f,i,j,k,l,m,n,index_myFamily,s1,s2
 real(pReal) :: sumf,StressRatio_p,StressRatio_pminus1,StressRatio_r,BoltzmannRatio,DotGamma0,Ndot0, &
                StressRatio_u,StressRatio_uminus1,tau_slip_pos,tau_slip_neg,vel_slip,dvel_slip,&
                dgdot_dtauslip_pos,dgdot_dtauslip_neg,dgdot_dtautwin,tau_twin,gdot_twin,stressRatio
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333
 real(pReal), dimension(constitutive_dislokmc_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,gdot_slip_neg
   
 !* Shortened notation
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = constitutive_dislokmc_totalNslip(instance)
 nt = constitutive_dislokmc_totalNtwin(instance)
 
 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 
!--------------------------------------------------------------------------------------------------
! Dislocation glide part
 gdot_slip_pos = 0.0_pReal
 gdot_slip_neg = 0.0_pReal
 dgdot_dtauslip_pos = 0.0_pReal
 dgdot_dtauslip_neg = 0.0_pReal

 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystems: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
     j = j+1_pInt
     !* Boltzmann ratio
     BoltzmannRatio = constitutive_dislokmc_QedgePerSlipSystem(j,instance)/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = &
        plasticState(ph)%state(j, of)*constitutive_dislokmc_burgersPerSlipSystem(j,instance)*&
        constitutive_dislokmc_v0PerSlipSystem(j,instance)
     !* Resolved shear stress on slip system
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)

     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
       tau_slip_pos = tau_slip_pos + constitutive_dislokmc_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + constitutive_dislokmc_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + constitutive_dislokmc_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k,index_myFamily+i,ph)
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + constitutive_dislokmc_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph)
     enddo nonSchmidSystems
     !* Applying damage to slip system
     tau_slip_pos = tau_slip_pos/slipDamage(j)
     tau_slip_neg = tau_slip_neg/slipDamage(j)

     significantPostitiveStress: if((abs(tau_slip_pos)-plasticState(ph)%state(6*ns+4*nt+j, of)) > tol_math_check) then
       !* Stress ratios
       stressRatio = ((abs(tau_slip_pos)-plasticState(ph)%state(6*ns+4*nt+j, of))/&
                      (constitutive_dislokmc_SolidSolutionStrength(instance)+&
                       constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p       = stressRatio** constitutive_dislokmc_pPerSlipFamily(f,instance)
       stressRatio_pminus1 = stressRatio**(constitutive_dislokmc_pPerSlipFamily(f,instance)-1.0_pReal)
       stressRatio_u       = stressRatio** constitutive_dislokmc_uPerSlipFamily(f,instance)
       stressRatio_uminus1 = stressRatio**(constitutive_dislokmc_uPerSlipFamily(f,instance)-1.0_pReal)
       !* Shear rates due to slip                                                                                                                                              
       vel_slip = exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)) &
                  * (1.0_pReal-constitutive_dislokmc_sPerSlipFamily(f,instance) &
                  * exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))
                       
       gdot_slip_pos(j) = DotGamma0 &
                       * StressRatio_u * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos) 

       !* Derivatives of shear rates                                                                                                              
       dvel_slip = &
         (abs(exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))&
         *BoltzmannRatio*constitutive_dislokmc_pPerSlipFamily(f,instance)&
         *constitutive_dislokmc_qPerSlipFamily(f,instance)/&
         (constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))*&
         StressRatio_pminus1*(1.0_pReal-StressRatio_p)**(constitutive_dislokmc_qPerSlipFamily(f,instance)-1.0_pReal)  )&
         *(1.0_pReal - 2.0_pReal*constitutive_dislokmc_sPerSlipFamily(f,instance)&
         *abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance))))   
                       
       dgdot_dtauslip_pos = DotGamma0 * &
         ( constitutive_dislokmc_uPerSlipFamily(f,instance)*StressRatio_uminus1 &
         /(constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))&
         * vel_slip &
         + StressRatio_u * dvel_slip)
     endif significantPostitiveStress
     significantNegativeStress: if((abs(tau_slip_neg)-plasticState(ph)%state(6*ns+4*nt+j, of)) > tol_math_check) then
       !* Stress ratios
       stressRatio = ((abs(tau_slip_neg)-plasticState(ph)%state(6*ns+4*nt+j, of))/&
                      (constitutive_dislokmc_SolidSolutionStrength(instance)+&
                       constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p       = stressRatio** constitutive_dislokmc_pPerSlipFamily(f,instance)
       stressRatio_pminus1 = stressRatio**(constitutive_dislokmc_pPerSlipFamily(f,instance)-1.0_pReal)
       stressRatio_u       = stressRatio** constitutive_dislokmc_uPerSlipFamily(f,instance)
       stressRatio_uminus1 = stressRatio**(constitutive_dislokmc_uPerSlipFamily(f,instance)-1.0_pReal)
       !* Shear rates due to slip                                                                                                                                              
       vel_slip = exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)) &
                  * (1.0_pReal-constitutive_dislokmc_sPerSlipFamily(f,instance) &
                  * exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))
                       
       gdot_slip_neg(j) = DotGamma0 &
                       * StressRatio_u * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg) 

       !* Derivatives of shear rates                                                                                                              
       dvel_slip = &
         (abs(exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))&
         *BoltzmannRatio*constitutive_dislokmc_pPerSlipFamily(f,instance)&
         *constitutive_dislokmc_qPerSlipFamily(f,instance)/&
         (constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))*&
         StressRatio_pminus1*(1.0_pReal-StressRatio_p)**(constitutive_dislokmc_qPerSlipFamily(f,instance)-1.0_pReal)  )&
         *(1.0_pReal - 2.0_pReal*constitutive_dislokmc_sPerSlipFamily(f,instance)&
         *abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance))))   
                       
       dgdot_dtauslip_neg = DotGamma0 * &
         ( constitutive_dislokmc_uPerSlipFamily(f,instance)*StressRatio_uminus1 &
         /(constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))&
         * vel_slip &
         + StressRatio_u * dvel_slip)
     endif significantNegativeStress
     !* Plastic velocity gradient for dislocation glide
     Lp = Lp + (gdot_slip_pos(j)+gdot_slip_neg(j))*0.5_pReal*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
 
     !* Calculation of the tangent of Lp
     forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + (dgdot_dtauslip_pos*nonSchmid_tensor(m,n,1)+&
                                   dgdot_dtauslip_neg*nonSchmid_tensor(m,n,2))*0.5_pReal*&
                                   lattice_Sslip(k,l,1,index_myFamily+i,ph)
   enddo slipSystems
 enddo slipFamilies

!--------------------------------------------------------------------------------------------------
! correct Lp and dLp_dTstar3333 for twinned fraction 
 !* Total twin volume fraction
 sumf = sum(plasticState(ph)%state((3_pInt*ns+1_pInt):(3_pInt*ns+nt), of)) ! safe for nt == 0
 Lp = Lp * (1.0_pReal - sumf)
 dLp_dTstar3333 = dLp_dTstar3333 * (1.0_pReal - sumf)
 
!--------------------------------------------------------------------------------------------------
! Mechanical twinning part
 gdot_twin = 0.0_pReal
 dgdot_dtautwin = 0.0_pReal
 j = 0_pInt
 twinFamilies: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph)) ! at which index starts my family
   twinSystems: do i = 1_pInt,constitutive_dislokmc_Ntwin(f,instance)
      j = j+1_pInt
      !* Resolved shear stress on twin system
      tau_twin = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_twin > tol_math_check) then
        StressRatio_r = (plasticState(ph)%state(7*ns+4*nt+j, of)/tau_twin)**constitutive_dislokmc_rPerTwinFamily(f,instance) 
      !* Shear rates and their derivatives due to twin
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin < constitutive_dislokmc_tau_r(j,instance)) then
              Ndot0=(abs(gdot_slip_pos(s1))*(plasticState(ph)%state(s2,of)+plasticState(ph)%state(ns+s2, of))+& !no non-Schmid behavior for fcc, just take the not influenced positive gdot_slip_pos (= gdot_slip_neg)
                     abs(gdot_slip_pos(s2))*(plasticState(ph)%state(s1,of)+plasticState(ph)%state(ns+s1, of)))/&
                    (constitutive_dislokmc_L0(instance)*constitutive_dislokmc_burgersPerSlipSystem(j,instance))*&
                    (1.0_pReal-exp(-constitutive_dislokmc_VcrossSlip(instance)/(kB*Temperature)*&
                    (constitutive_dislokmc_tau_r(j,instance)-tau_twin)))
            else
              Ndot0=0.0_pReal
            end if
          case default
            Ndot0=constitutive_dislokmc_Ndot0PerTwinSystem(j,instance)
        end select
        gdot_twin = &
          (constitutive_dislokmc_MaxTwinFraction(instance)-sumf)*lattice_shearTwin(index_myFamily+i,ph)*&
          plasticState(ph)%state(7*ns+5*nt+j, of)*Ndot0*exp(-StressRatio_r)
        dgdot_dtautwin = ((gdot_twin*constitutive_dislokmc_rPerTwinFamily(f,instance))/tau_twin)*StressRatio_r
      endif
 
      !* Plastic velocity gradient for mechanical twinning
      Lp = Lp + gdot_twin*lattice_Stwin(1:3,1:3,index_myFamily+i,ph)
 
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin*&
                                  lattice_Stwin(k,l,index_myFamily+i,ph)*&
                                  lattice_Stwin(m,n,index_myFamily+i,ph)
   enddo twinSystems
 enddo twinFamilies
 
 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)
 
end subroutine constitutive_dislokmc_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_dotState(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &
   mappingConstitutive
 use lattice,  only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_Sslip, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid, &
   lattice_sheartwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID

 implicit none
 real(pReal), dimension(6),  intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: instance,ns,nt,f,i,j,k,index_myFamily,s1,s2, &
                  ph, &
                  of
 real(pReal) :: &
   sumf, &
   stressRatio_p,&
   BoltzmannRatio,&
   DotGamma0,&
   stressRatio_u,&
   stressRatio, &
   EdgeDipMinDistance,&
   AtomicVolume,&
   VacancyDiffusion,&
   StressRatio_r,&
   Ndot0,&
   tau_slip_pos,&
   tau_slip_neg,&
   DotRhoMultiplication,&
   EdgeDipDistance, &
   DotRhoEdgeDipAnnihilation, &
   DotRhoEdgeEdgeAnnihilation, &
   ClimbVelocity, &
   DotRhoEdgeDipClimb, &
   DotRhoDipFormation, &
   tau_twin, &
   vel_slip, &
   gdot_slip
 real(pReal), dimension(constitutive_dislokmc_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos, gdot_slip_neg

 !* Shortened notation
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = constitutive_dislokmc_totalNslip(instance)
 nt = constitutive_dislokmc_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(plasticState(ph)%state((3_pInt*ns+1_pInt):(3_pInt*ns+nt), of)) ! safe for nt == 0
 plasticState(ph)%dotState(:,of) = 0.0_pReal
 
 !* Dislocation density evolution
 gdot_slip_pos = 0.0_pReal
 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystems: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
     j = j+1_pInt
     !* Boltzmann ratio
     BoltzmannRatio = constitutive_dislokmc_QedgePerSlipSystem(j,instance)/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = &
        plasticState(ph)%state(j, of)*constitutive_dislokmc_burgersPerSlipSystem(j,instance)*&
        constitutive_dislokmc_v0PerSlipSystem(j,instance)
     !* Resolved shear stress on slip system
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos

     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
       tau_slip_pos = tau_slip_pos + constitutive_dislokmc_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,  index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + constitutive_dislokmc_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems

     significantPositiveStress: if((abs(tau_slip_pos)-plasticState(ph)%state(6*ns+4*nt+j, of)) > tol_math_check) then
       !* Stress ratios
       stressRatio = ((abs(tau_slip_pos)-plasticState(ph)%state(6*ns+4*nt+j, of))/&
                      (constitutive_dislokmc_SolidSolutionStrength(instance)+&
                       constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p = stressRatio** constitutive_dislokmc_pPerSlipFamily(f,instance)
       stressRatio_u = stressRatio** constitutive_dislokmc_uPerSlipFamily(f,instance)
       !* Shear rates due to slip                                                                                                                                                                                                                                                                           
       vel_slip = exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)) &
                  * (1.0_pReal-constitutive_dislokmc_sPerSlipFamily(f,instance) &
                  * exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))
                       
       gdot_slip_pos(j) = DotGamma0 &
                       * StressRatio_u * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos) 
     endif significantPositiveStress
     significantNegativeStress: if((abs(tau_slip_neg)-plasticState(ph)%state(6*ns+4*nt+j, of)) > tol_math_check) then
       !* Stress ratios
       stressRatio = ((abs(tau_slip_neg)-plasticState(ph)%state(6*ns+4*nt+j, of))/&
                      (constitutive_dislokmc_SolidSolutionStrength(instance)+&
                       constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p = stressRatio** constitutive_dislokmc_pPerSlipFamily(f,instance)
       stressRatio_u = stressRatio** constitutive_dislokmc_uPerSlipFamily(f,instance)
       !* Shear rates due to slip                                                                                                                                                                                                                                                                           
       vel_slip = exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)) &
                  * (1.0_pReal-constitutive_dislokmc_sPerSlipFamily(f,instance) &
                  * exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))
                       
       gdot_slip_neg(j) = DotGamma0 &
                       * StressRatio_u * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg) 
     endif significantNegativeStress
     gdot_slip = (gdot_slip_pos(j)+gdot_slip_neg(j))*0.5_pReal
     !* Multiplication
     DotRhoMultiplication = abs(gdot_slip)/&
                               (constitutive_dislokmc_burgersPerSlipSystem(j,instance)* &
                                plasticState(ph)%state(5*ns+3*nt+j, of))
 
     !* Dipole formation
     EdgeDipMinDistance = &
       constitutive_dislokmc_CEdgeDipMinDistance(instance)*constitutive_dislokmc_burgersPerSlipSystem(j,instance)
     if (tau_slip_pos == 0.0_pReal) then
       DotRhoDipFormation = 0.0_pReal
     else
       EdgeDipDistance = &
         (3.0_pReal*lattice_mu(ph)*constitutive_dislokmc_burgersPerSlipSystem(j,instance))/&
         (16.0_pReal*pi*abs(tau_slip_pos))
       if (EdgeDipDistance>plasticState(ph)%state(5*ns+3*nt+j, of)) EdgeDipDistance=plasticState(ph)%state(5*ns+3*nt+j, of)
       if (EdgeDipDistance<EdgeDipMinDistance) EdgeDipDistance=EdgeDipMinDistance
       DotRhoDipFormation = &
         ((2.0_pReal*EdgeDipDistance)/constitutive_dislokmc_burgersPerSlipSystem(j,instance))*&
         plasticState(ph)%state(j, of)*abs(gdot_slip)*constitutive_dislokmc_dipoleFormationFactor(instance)
     endif
 
    !* Spontaneous annihilation of 2 single edge dislocations
    DotRhoEdgeEdgeAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/constitutive_dislokmc_burgersPerSlipSystem(j,instance))*&
        plasticState(ph)%state(j, of)*abs(gdot_slip)
 
    !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
    DotRhoEdgeDipAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/constitutive_dislokmc_burgersPerSlipSystem(j,instance))*&
        plasticState(ph)%state(ns+j, of)*abs(gdot_slip)
 
      !* Dislocation dipole climb
     AtomicVolume = &
        constitutive_dislokmc_CAtomicVolume(instance)*constitutive_dislokmc_burgersPerSlipSystem(j,instance)**(3.0_pReal)
     VacancyDiffusion = &
        constitutive_dislokmc_D0(instance)*exp(-constitutive_dislokmc_Qsd(instance)/(kB*Temperature))
     if (tau_slip_pos == 0.0_pReal) then
       DotRhoEdgeDipClimb = 0.0_pReal
     else
       ClimbVelocity = &
          ((3.0_pReal*lattice_mu(ph)*VacancyDiffusion*AtomicVolume)/(2.0_pReal*pi*kB*Temperature))*&
          (1/(EdgeDipDistance+EdgeDipMinDistance))
       DotRhoEdgeDipClimb = &
          (4.0_pReal*ClimbVelocity*plasticState(ph)%state(ns+j, of))/(EdgeDipDistance-EdgeDipMinDistance)
     endif
 
     !* Edge dislocation density rate of change
     plasticState(ph)%dotState(j, of) = &
        DotRhoMultiplication-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
 
     !* Edge dislocation dipole density rate of change
     plasticState(ph)%dotState(ns+j, of) = &
        DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb
 
     !* Dotstate for accumulated shear due to slip
     plasticState(ph)%dotState(2_pInt*ns+j, of) = gdot_slip
 
   enddo slipSystems
 enddo slipFamilies
 
 !* Twin volume fraction evolution
 j = 0_pInt
 twinFamilies: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph)) ! at which index starts my family
   twinSystems: do i = 1_pInt,constitutive_dislokmc_Ntwin(f,instance)
      j = j+1_pInt
      !* Resolved shear stress on twin system
      tau_twin = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
      !* Stress ratios
      if (tau_twin > tol_math_check) then
        StressRatio_r = (plasticState(ph)%state(7*ns+4*nt+j, of)/tau_twin)**constitutive_dislokmc_rPerTwinFamily(f,instance)
      !* Shear rates and their derivatives due to twin

        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin < constitutive_dislokmc_tau_r(j,instance)) then
              Ndot0=(abs(gdot_slip_pos(s1))*(plasticState(ph)%state(s2, of)+plasticState(ph)%state(ns+s2, of))+&  !no non-Schmid behavior for fcc, just take the not influenced positive slip (gdot_slip_pos = gdot_slip_neg)
                     abs(gdot_slip_pos(s2))*(plasticState(ph)%state(s1, of)+plasticState(ph)%state(ns+s1, of)))/&
                    (constitutive_dislokmc_L0(instance)*constitutive_dislokmc_burgersPerSlipSystem(j,instance))*&
                    (1.0_pReal-exp(-constitutive_dislokmc_VcrossSlip(instance)/(kB*Temperature)*&
                    (constitutive_dislokmc_tau_r(j,instance)-tau_twin)))
            else
              Ndot0=0.0_pReal
            end if
          case default
            Ndot0=constitutive_dislokmc_Ndot0PerTwinSystem(j,instance)
        end select

        plasticState(ph)%dotState(3_pInt*ns+j, of) = &
          (constitutive_dislokmc_MaxTwinFraction(instance)-sumf)*&
          plasticState(ph)%state(7_pInt*ns+5_pInt*nt+j, of)*Ndot0*exp(-StressRatio_r)
        !* Dotstate for accumulated shear due to twin
        plasticState(ph)%dotState(3_pInt*ns+nt+j, of) = plasticState(ph)%dotState(3_pInt*ns+j, of) * &
                                                          lattice_sheartwin(index_myfamily+i,ph)
      endif
   enddo twinSystems
 enddo twinFamilies
 
end subroutine constitutive_dislokmc_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns accumulated slip
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_getAccumulatedSlip(nSlip,accumulatedSlip,ipc, ip, el)
 use lattice, only: &
   lattice_maxNslipFamily
 use material, only: &
   mappingConstitutive, &
   plasticState, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(:), allocatable :: &
   accumulatedSlip
 integer(pInt) :: &
   nSlip
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   offset, &
   phase, &
   instance, &
   offset_accshear_slip, &
   f, j, i

 offset = mappingConstitutive(1,ipc,ip,el)
 phase = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(phase)
 nSlip = constitutive_dislokmc_totalNslip(instance)
 allocate(accumulatedSlip(nSlip))
 offset_accshear_slip = 2_pInt*nSlip

 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   slipSystems: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
      j = j+1_pInt
      accumulatedSlip(j) = plasticState(phase)%state(offset_accshear_slip+j,offset)
   enddo slipSystems
 enddo slipFamilies
   
end subroutine constitutive_dislokmc_getAccumulatedSlip


!--------------------------------------------------------------------------------------------------
!> @brief returns accumulated slip
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislokmc_getSlipRate(nSlip,slipRate,ipc, ip, el)
 use lattice, only: &
   lattice_maxNslipFamily
 use material, only: &
   mappingConstitutive, &
   plasticState, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(:), allocatable :: &
   slipRate
 integer(pInt) :: &
   nSlip
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   offset, &
   phase, &
   instance, &
   offset_accshear_slip, &
   f, j, i

 offset = mappingConstitutive(1,ipc,ip,el)
 phase = mappingConstitutive(2,ipc,ip,el)
 instance = phase_plasticityInstance(phase)
 nSlip = constitutive_dislokmc_totalNslip(instance)
 allocate(slipRate(nSlip))
 offset_accshear_slip = 2_pInt*nSlip

 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   slipSystems: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
      j = j+1_pInt
      slipRate(j) = plasticState(phase)%dotState(offset_accshear_slip+j,offset)
   enddo slipSystems
 enddo slipFamilies
   
end subroutine constitutive_dislokmc_getSlipRate

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_dislokmc_postResults(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance,& 
   plasticState, &
   mappingConstitutive
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_Sslip, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid, &
   lattice_shearTwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(constitutive_dislokmc_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           constitutive_dislokmc_postResults

 integer(pInt) :: &
   instance,&
   ns,nt,&
   f,o,i,c,j,k,index_myFamily,&
   s1,s2, &
   ph, &
   of
 real(pReal) :: sumf,tau_twin,StressRatio_p,StressRatio_pminus1,&
                BoltzmannRatio,DotGamma0,StressRatio_r,Ndot0,stressRatio
 real(pReal) :: dvel_slip, vel_slip
 real(pReal) :: StressRatio_u,StressRatio_uminus1
 real(pReal), dimension(constitutive_dislokmc_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg
 
 !* Shortened notation
 of = mappingConstitutive(1,ipc,ip,el)
 ph = mappingConstitutive(2,ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = constitutive_dislokmc_totalNslip(instance)
 nt = constitutive_dislokmc_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(plasticState(ph)%state((3_pInt*ns+1_pInt):(3_pInt*ns+nt), of))                          ! safe for nt == 0
 
 !* Required output
 c = 0_pInt
 constitutive_dislokmc_postResults = 0.0_pReal

 do o = 1_pInt,constitutive_dislokmc_Noutput(instance)
    select case(constitutive_dislokmc_outputID(o,instance))
 
      case (edge_density_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+ns) = plasticState(ph)%state(1_pInt:ns, of)
        c = c + ns
      case (dipole_density_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+ns) = plasticState(ph)%state(ns+1_pInt:2_pInt*ns, of)
        c = c + ns
      case (shear_rate_slip_ID,shear_rate_twin_ID,stress_exponent_ID)
        gdot_slip_pos = 0.0_pReal
        gdot_slip_neg = 0.0_pReal
        dgdot_dtauslip_pos = 0.0_pReal
        dgdot_dtauslip_neg = 0.0_pReal
        j = 0_pInt
        slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
          index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
          slipSystems: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
            j = j + 1_pInt
            !* Boltzmann ratio
            BoltzmannRatio = constitutive_dislokmc_QedgePerSlipSystem(j,instance)/(kB*Temperature)
            !* Initial shear rates
            DotGamma0 = &
              plasticState(ph)%state(j, of)*constitutive_dislokmc_burgersPerSlipSystem(j,instance)*&
              constitutive_dislokmc_v0PerSlipSystem(j,instance)
            !* Resolved shear stress on slip system
            tau_slip_pos(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
            tau_slip_neg(j)  = tau_slip_pos(j)

            nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
              tau_slip_pos = tau_slip_pos + constitutive_dislokmc_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
              tau_slip_neg = tau_slip_neg + constitutive_dislokmc_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
            enddo nonSchmidSystems

            significantPostitiveStress: if((abs(tau_slip_pos(j))-plasticState(ph)%state(6*ns+4*nt+j, of)) > tol_math_check) then
              !* Stress ratios
              stressRatio = ((abs(tau_slip_pos(j))-plasticState(ph)%state(6*ns+4*nt+j, of))/&
                      (constitutive_dislokmc_SolidSolutionStrength(instance)+&
                       constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance)))
              stressRatio_p       = stressRatio** constitutive_dislokmc_pPerSlipFamily(f,instance)
              stressRatio_pminus1 = stressRatio**(constitutive_dislokmc_pPerSlipFamily(f,instance)-1.0_pReal)
              stressRatio_u       = stressRatio** constitutive_dislokmc_uPerSlipFamily(f,instance)
              stressRatio_uminus1 = stressRatio**(constitutive_dislokmc_uPerSlipFamily(f,instance)-1.0_pReal)
              !* Shear rates due to slip                                                                                                                                                                                                                                                                           
              vel_slip = exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)) &
                     * (1.0_pReal-constitutive_dislokmc_sPerSlipFamily(f,instance) &
                     * exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))
                       
              gdot_slip_pos(j) = DotGamma0 &
                       * StressRatio_u * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos(j))
              !* Derivatives of shear rates 
              dvel_slip = &
                (abs(exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))&
                *BoltzmannRatio*constitutive_dislokmc_pPerSlipFamily(f,instance)&
                *constitutive_dislokmc_qPerSlipFamily(f,instance)/&
                (constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))*&
                StressRatio_pminus1*(1.0_pReal-StressRatio_p)**(constitutive_dislokmc_qPerSlipFamily(f,instance)-1.0_pReal)  )&
                *(1.0_pReal - 2.0_pReal*constitutive_dislokmc_sPerSlipFamily(f,instance)&
                *abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance))))   
                              
              dgdot_dtauslip_pos(j) = DotGamma0 * &
                ( constitutive_dislokmc_uPerSlipFamily(f,instance)*StressRatio_uminus1 &
                /(constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))&
                * vel_slip &
                + StressRatio_u * dvel_slip)
            endif significantPostitiveStress
            significantNegativeStress: if((abs(tau_slip_neg(j))-plasticState(ph)%state(6*ns+4*nt+j, of)) > tol_math_check) then
              !* Stress ratios
              stressRatio = ((abs(tau_slip_neg(j))-plasticState(ph)%state(6*ns+4*nt+j, of))/&
                      (constitutive_dislokmc_SolidSolutionStrength(instance)+&
                       constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance)))
              stressRatio_p       = stressRatio** constitutive_dislokmc_pPerSlipFamily(f,instance)
              stressRatio_pminus1 = stressRatio**(constitutive_dislokmc_pPerSlipFamily(f,instance)-1.0_pReal)
              stressRatio_u       = stressRatio** constitutive_dislokmc_uPerSlipFamily(f,instance)
              stressRatio_uminus1 = stressRatio**(constitutive_dislokmc_uPerSlipFamily(f,instance)-1.0_pReal)
              !* Shear rates due to slip                                                                                                                                                                                                                                                                           
              vel_slip = exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)) &
                     * (1.0_pReal-constitutive_dislokmc_sPerSlipFamily(f,instance) &
                     * exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))
                       
              gdot_slip_neg(j) = DotGamma0 &
                       * StressRatio_u * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg(j))
              !* Derivatives of shear rates 
              dvel_slip = &
                (abs(exp(-BoltzmannRatio*(1.0_pReal-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance)))&
                *BoltzmannRatio*constitutive_dislokmc_pPerSlipFamily(f,instance)&
                *constitutive_dislokmc_qPerSlipFamily(f,instance)/&
                (constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))*&
                StressRatio_pminus1*(1.0_pReal-StressRatio_p)**(constitutive_dislokmc_qPerSlipFamily(f,instance)-1.0_pReal)  )&
                *(1.0_pReal - 2.0_pReal*constitutive_dislokmc_sPerSlipFamily(f,instance)&
                *abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislokmc_qPerSlipFamily(f,instance))))   
                              
              dgdot_dtauslip_neg(j) = DotGamma0 * &
                ( constitutive_dislokmc_uPerSlipFamily(f,instance)*StressRatio_uminus1 &
                /(constitutive_dislokmc_SolidSolutionStrength(instance)+constitutive_dislokmc_tau_peierlsPerSlipFamily(f,instance))&
                * vel_slip &
                + StressRatio_u * dvel_slip)
            endif significantNegativeStress
          enddo slipSystems
        enddo slipFamilies

        if     (constitutive_dislokmc_outputID(o,instance) == shear_rate_slip_ID) then
          constitutive_dislokmc_postResults(c+1:c+ns) = (gdot_slip_pos + gdot_slip_neg)*0.5_pReal
          c = c + ns
        elseif (constitutive_dislokmc_outputID(o,instance) == shear_rate_twin_ID) then
          if (nt > 0_pInt) then
            j = 0_pInt
            twinFamilies1: do f = 1_pInt,lattice_maxNtwinFamily
              index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
              twinSystems1: do i = 1,constitutive_dislokmc_Ntwin(f,instance)
                j = j + 1_pInt
 
                !* Resolved shear stress on twin system
                tau_twin = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
                !* Stress ratios
                StressRatio_r = (plasticState(ph)%state(7_pInt*ns+4_pInt*nt+j, of)/ &
                                            tau_twin)**constitutive_dislokmc_rPerTwinFamily(f,instance)
 
                !* Shear rates due to twin
                if ( tau_twin > 0.0_pReal ) then
                  select case(lattice_structure(ph))
                    case (LATTICE_fcc_ID)
                      s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
                      s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
                      if (tau_twin < constitutive_dislokmc_tau_r(j,instance)) then
                        Ndot0=(abs(gdot_slip_pos(s1))*(plasticState(ph)%state(s2, of)+plasticState(ph)%state(ns+s2, of))+& !no non-Schmid behavior for fcc, just take the not influenced positive slip (gdot_slip_pos = gdot_slip_neg)
                               abs(gdot_slip_pos(s2))*(plasticState(ph)%state(s1, of)+plasticState(ph)%state(ns+s1, of)))/&
                              (constitutive_dislokmc_L0(instance)*&
                               constitutive_dislokmc_burgersPerSlipSystem(j,instance))*&
                              (1.0_pReal-exp(-constitutive_dislokmc_VcrossSlip(instance)/(kB*Temperature)*&
                              (constitutive_dislokmc_tau_r(j,instance)-tau_twin)))
                      else
                        Ndot0=0.0_pReal
                      end if

                    case default
                      Ndot0=constitutive_dislokmc_Ndot0PerTwinSystem(j,instance)
                  end select
                  constitutive_dislokmc_postResults(c+j) = &
                    (constitutive_dislokmc_MaxTwinFraction(instance)-sumf)*lattice_shearTwin(index_myFamily+i,ph)*&
                    plasticState(ph)%state(7_pInt*ns+5_pInt*nt+j, of)*Ndot0*exp(-StressRatio_r)
                endif
              enddo twinSystems1
            enddo twinFamilies1
          endif
          c = c + nt
        elseif(constitutive_dislokmc_outputID(o,instance) == stress_exponent_ID) then
          do j = 1_pInt, ns
            if ((gdot_slip_pos(j)+gdot_slip_neg(j))*0.5_pReal==0.0_pReal) then
              constitutive_dislokmc_postResults(c+j) = 0.0_pReal
            else
              constitutive_dislokmc_postResults(c+j) = (tau_slip_pos(j)+tau_slip_neg(j))/&
                                                       (gdot_slip_pos(j)+gdot_slip_neg(j))*&
                                                       (dgdot_dtauslip_pos(j)+dgdot_dtauslip_neg(j))* 0.5_pReal
            endif
          enddo
           c = c + ns
        endif

      case (accumulated_shear_slip_ID)
       constitutive_dislokmc_postResults(c+1_pInt:c+ns) = &
                      plasticState(ph)%state((2_pInt*ns+1_pInt):(3_pInt*ns), of)
        c = c + ns
      case (mfp_slip_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+ns) =&
                      plasticState(ph)%state((5_pInt*ns+3_pInt*nt+1_pInt):(6_pInt*ns+3_pInt*nt), of)
        c = c + ns
      case (resolved_stress_slip_ID)
        j = 0_pInt
        slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           slipSystems1: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
              j = j + 1_pInt
              constitutive_dislokmc_postResults(c+j) =&
                                dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
        enddo slipSystems1; enddo slipFamilies1
        c = c + ns
      case (threshold_stress_slip_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+ns) = &
                                plasticState(ph)%state((6_pInt*ns+4_pInt*nt+1_pInt):(7_pInt*ns+4_pInt*nt), of)
        c = c + ns
      case (edge_dipole_distance_ID)
        j = 0_pInt
        slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           slipSystems2: do i = 1_pInt,constitutive_dislokmc_Nslip(f,instance)
              j = j + 1_pInt
              constitutive_dislokmc_postResults(c+j) = &
                (3.0_pReal*lattice_mu(ph)*constitutive_dislokmc_burgersPerSlipSystem(j,instance))/&
                (16.0_pReal*pi*abs(dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))))
              constitutive_dislokmc_postResults(c+j)=min(constitutive_dislokmc_postResults(c+j),&
                                                            plasticState(ph)%state(5*ns+3*nt+j, of))
        enddo slipSystems2; enddo slipFamilies2
        c = c + ns
      case (twin_fraction_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+nt) = plasticState(ph)%state((3_pInt*ns+1_pInt):(3_pInt*ns+nt), of)
        c = c + nt

      case (accumulated_shear_twin_ID)
       constitutive_dislokmc_postResults(c+1_pInt:c+nt) = plasticState(ph)% &
                                        state((3_pInt*ns+nt+1_pInt)       :(3_pInt*ns+2_pInt*nt), of)
        c = c + nt     
      case (mfp_twin_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+nt) = plasticState(ph)% &
                                        state((6_pInt*ns+3_pInt*nt+1_pInt):(6_pInt*ns+4_pInt*nt), of)
        c = c + nt
      case (resolved_stress_twin_ID)
        if (nt > 0_pInt) then
          j = 0_pInt
          twinFamilies2: do f = 1_pInt,lattice_maxNtwinFamily
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
            twinSystems2: do i = 1_pInt,constitutive_dislokmc_Ntwin(f,instance)
              j = j + 1_pInt
              constitutive_dislokmc_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
          enddo twinSystems2; enddo twinFamilies2
        endif
        c = c + nt
      case (threshold_stress_twin_ID)
        constitutive_dislokmc_postResults(c+1_pInt:c+nt) = plasticState(ph)% &
                                       state((7_pInt*ns+4_pInt*nt+1_pInt):(7_pInt*ns+5_pInt*nt), of)
        c = c + nt
    end select
 enddo
end function constitutive_dislokmc_postResults

end module constitutive_dislokmc
