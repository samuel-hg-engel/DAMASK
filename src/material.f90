!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parses material config file, either solverJobName.materialConfig or material.config
!> @details reads the material configuration file, where solverJobName.materialConfig takes
!! precedence over material.config and parses the sections 'homogenization', 'crystallite',
!! 'phase', 'texture', and 'microstucture'
!--------------------------------------------------------------------------------------------------
module material
 use prec
 use math
 use config

 implicit none
 private
 character(len=*),                         parameter,            public :: &
   ELASTICITY_hooke_label               = 'hooke', &
   PLASTICITY_none_label                = 'none', &
   PLASTICITY_isotropic_label           = 'isotropic', &
   PLASTICITY_phenopowerlaw_label       = 'phenopowerlaw', &
   PLASTICITY_kinehardening_label       = 'kinehardening', &
   PLASTICITY_dislotwin_label           = 'dislotwin', &
   PLASTICITY_disloucla_label           = 'disloucla', &
   PLASTICITY_nonlocal_label            = 'nonlocal', &
   SOURCE_thermal_dissipation_label     = 'thermal_dissipation', &
   SOURCE_thermal_externalheat_label    = 'thermal_externalheat', &
   SOURCE_damage_isoBrittle_label       = 'damage_isobrittle', &
   SOURCE_damage_isoDuctile_label       = 'damage_isoductile', &
   SOURCE_damage_anisoBrittle_label     = 'damage_anisobrittle', &
   SOURCE_damage_anisoDuctile_label     = 'damage_anisoductile', &
   KINEMATICS_thermal_expansion_label   = 'thermal_expansion', &
   KINEMATICS_cleavage_opening_label    = 'cleavage_opening', &
   KINEMATICS_slipplane_opening_label   = 'slipplane_opening', &
   STIFFNESS_DEGRADATION_damage_label   = 'damage', &
   THERMAL_isothermal_label             = 'isothermal', &
   THERMAL_adiabatic_label              = 'adiabatic', &
   THERMAL_conduction_label             = 'conduction', &
   DAMAGE_none_label                    = 'none', &
   DAMAGE_local_label                   = 'local', &
   DAMAGE_nonlocal_label                = 'nonlocal', &
   HOMOGENIZATION_none_label            = 'none', &
   HOMOGENIZATION_isostrain_label       = 'isostrain', &
   HOMOGENIZATION_rgc_label             = 'rgc'



 enum, bind(c)
   enumerator :: ELASTICITY_undefined_ID, &
                 ELASTICITY_hooke_ID
 end enum
 enum, bind(c)
   enumerator :: PLASTICITY_undefined_ID, &
                 PLASTICITY_none_ID, &
                 PLASTICITY_isotropic_ID, &
                 PLASTICITY_phenopowerlaw_ID, &
                 PLASTICITY_kinehardening_ID, &
                 PLASTICITY_dislotwin_ID, &
                 PLASTICITY_disloucla_ID, &
                 PLASTICITY_nonlocal_ID
 end enum

 enum, bind(c)
   enumerator :: SOURCE_undefined_ID, &
                 SOURCE_thermal_dissipation_ID, &
                 SOURCE_thermal_externalheat_ID, &
                 SOURCE_damage_isoBrittle_ID, &
                 SOURCE_damage_isoDuctile_ID, &
                 SOURCE_damage_anisoBrittle_ID, &
                 SOURCE_damage_anisoDuctile_ID
 end enum

 enum, bind(c)
   enumerator :: KINEMATICS_undefined_ID, &
                 KINEMATICS_cleavage_opening_ID, &
                 KINEMATICS_slipplane_opening_ID, &
                 KINEMATICS_thermal_expansion_ID
 end enum

 enum, bind(c)
   enumerator :: STIFFNESS_DEGRADATION_undefined_ID, &
                 STIFFNESS_DEGRADATION_damage_ID
 end enum

 enum, bind(c)
   enumerator :: THERMAL_isothermal_ID, &
                 THERMAL_adiabatic_ID, &
                 THERMAL_conduction_ID
 end enum

 enum, bind(c)
   enumerator :: DAMAGE_none_ID, &
                 DAMAGE_local_ID, &
                 DAMAGE_nonlocal_ID
 end enum

 enum, bind(c)
   enumerator :: HOMOGENIZATION_undefined_ID, &
                 HOMOGENIZATION_none_ID, &
                 HOMOGENIZATION_isostrain_ID, &
                 HOMOGENIZATION_rgc_ID
 end enum

 integer(kind(ELASTICITY_undefined_ID)),     dimension(:),   allocatable, public, protected :: &
   phase_elasticity                                                                                 !< elasticity of each phase
 integer(kind(PLASTICITY_undefined_ID)),     dimension(:),   allocatable, public, protected :: &
   phase_plasticity                                                                                 !< plasticity of each phase
 integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
   thermal_type                                                                                     !< thermal transport model
 integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
   damage_type                                                                                      !< nonlocal damage model

 integer(kind(SOURCE_undefined_ID)),         dimension(:,:), allocatable, public, protected :: &
   phase_source, &                                                                                  !< active sources mechanisms of each phase
   phase_kinematics, &                                                                              !< active kinematic mechanisms of each phase
   phase_stiffnessDegradation                                                                       !< active stiffness degradation mechanisms of each phase

 integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
   homogenization_type                                                                              !< type of each homogenization

 integer(pInt), public, protected :: &
   homogenization_maxNgrains                                                                        !< max number of grains in any USED homogenization

 integer(pInt), dimension(:), allocatable, public, protected :: &
   phase_Nsources, &                                                                                !< number of source mechanisms active in each phase
   phase_Nkinematics, &                                                                             !< number of kinematic mechanisms active in each phase
   phase_NstiffnessDegradations, &                                                                  !< number of stiffness degradation mechanisms active in each phase
   phase_Noutput, &                                                                                 !< number of '(output)' items per phase
   phase_elasticityInstance, &                                                                      !< instance of particular elasticity of each phase
   phase_plasticityInstance, &                                                                      !< instance of particular plasticity of each phase
   crystallite_Noutput, &                                                                           !< number of '(output)' items per crystallite setting
   homogenization_Ngrains, &                                                                        !< number of grains in each homogenization
   homogenization_Noutput, &                                                                        !< number of '(output)' items per homogenization
   homogenization_typeInstance, &                                                                   !< instance of particular type of each homogenization
   thermal_typeInstance, &                                                                          !< instance of particular type of each thermal transport
   damage_typeInstance, &                                                                           !< instance of particular type of each nonlocal damage
   microstructure_crystallite                                                                       !< crystallite setting ID of each microstructure ! DEPRECATED !!!!

 real(pReal), dimension(:), allocatable, public, protected :: &
   thermal_initialT, &                                                                              !< initial temperature per each homogenization
   damage_initialPhi                                                                                !< initial damage per each homogenization

! NEW MAPPINGS 
 integer, dimension(:),     allocatable, public, protected :: &                                     ! (elem)
   material_homogenizationAt                                                                        !< homogenization ID of each element (copy of mesh_homogenizationAt)
 integer, dimension(:,:),   allocatable, public, protected :: &                                     ! (ip,elem)
   material_homogenizationMemberAt                                                                  !< position of the element within its homogenization instance
 integer, dimension(:,:), allocatable, public, protected :: &                                       ! (constituent,elem)
   material_phaseAt                                                                                 !< phase ID of each element
 integer, dimension(:,:,:), allocatable, public, protected :: &                                     ! (constituent,ip,elem)
   material_phaseMemberAt                                                                           !< position of the element within its phase instance
! END NEW MAPPINGS
 
! DEPRECATED: use material_phaseAt
 integer(pInt), dimension(:,:,:), allocatable, public :: &
   material_phase                                                                                   !< phase (index) of each grain,IP,element

 type(tPlasticState), allocatable, dimension(:), public :: &
   plasticState
 type(tSourceState),  allocatable, dimension(:), public :: &
   sourceState
 type(tState),        allocatable, dimension(:), public :: &
   homogState, &
   thermalState, &
   damageState

 integer(pInt), dimension(:,:,:), allocatable, public, protected :: &
   material_texture                                                                                 !< texture (index) of each grain,IP,element

 real(pReal), dimension(:,:,:,:), allocatable, public, protected :: &
   material_EulerAngles                                                                             !< initial orientation of each grain,IP,element

 logical, dimension(:), allocatable, public, protected :: &
   microstructure_active, &
   microstructure_elemhomo, &                                                                       !< flag to indicate homogeneous microstructure distribution over element's IPs
   phase_localPlasticity                                                                            !< flags phases with local constitutive law

 integer(pInt), private :: &
   microstructure_maxNconstituents, &                                                               !< max number of constituents in any phase
   texture_maxNgauss                                                                                !< max number of Gauss components in any texture

 integer(pInt), dimension(:), allocatable, private :: &
   microstructure_Nconstituents, &                                                                  !< number of constituents in each microstructure
   texture_Ngauss                                                                                   !< number of Gauss components per texture

 integer(pInt), dimension(:,:), allocatable, private :: &
   microstructure_phase, &                                                                          !< phase IDs of each microstructure
   microstructure_texture                                                                           !< texture IDs of each microstructure

 real(pReal), dimension(:,:), allocatable, private :: &
   microstructure_fraction                                                                          !< vol fraction of each constituent in microstructure

 real(pReal), dimension(:,:,:), allocatable, private :: &
   texture_Gauss, &                                                                                 !< data of each Gauss component
   texture_transformation                                                                           !< transformation for each texture

 logical, dimension(:), allocatable, private :: &
   homogenization_active

! BEGIN DEPRECATED
 integer(pInt), dimension(:,:,:), allocatable, public :: phaseAt                                    !< phase ID of every material point (ipc,ip,el)
 integer(pInt), dimension(:,:,:), allocatable, public :: phasememberAt                              !< memberID of given phase at every material point (ipc,ip,el)

 integer(pInt), dimension(:,:,:), allocatable, public, target :: mappingHomogenization              !< mapping from material points to offset in heterogenous state/field
 integer(pInt), dimension(:,:),   allocatable, private, target :: mappingHomogenizationConst         !< mapping from material points to offset in constant state/field
! END DEPRECATED

 type(tHomogMapping), allocatable, dimension(:), public :: &
   thermalMapping, &                                                                                !< mapping for thermal state/fields
   damageMapping                                                                                    !< mapping for damage state/fields

 type(group_float),  allocatable, dimension(:), public :: &
   temperature, &                                                                                   !< temperature field
   damage, &                                                                                        !< damage field
   temperatureRate                                                                                  !< temperature change rate field

 public :: &
   material_init, &
   material_allocatePlasticState, &
   material_allocateSourceState, &
   ELASTICITY_hooke_ID ,&
   PLASTICITY_none_ID, &
   PLASTICITY_isotropic_ID, &
   PLASTICITY_phenopowerlaw_ID, &
   PLASTICITY_kinehardening_ID, &
   PLASTICITY_dislotwin_ID, &
   PLASTICITY_disloucla_ID, &
   PLASTICITY_nonlocal_ID, &
   SOURCE_thermal_dissipation_ID, &
   SOURCE_thermal_externalheat_ID, &
   SOURCE_damage_isoBrittle_ID, &
   SOURCE_damage_isoDuctile_ID, &
   SOURCE_damage_anisoBrittle_ID, &
   SOURCE_damage_anisoDuctile_ID, &
   KINEMATICS_cleavage_opening_ID, &
   KINEMATICS_slipplane_opening_ID, &
   KINEMATICS_thermal_expansion_ID, &
   STIFFNESS_DEGRADATION_damage_ID, &
   THERMAL_isothermal_ID, &
   THERMAL_adiabatic_ID, &
   THERMAL_conduction_ID, &
   DAMAGE_none_ID, &
   DAMAGE_local_ID, &
   DAMAGE_nonlocal_ID, &
   HOMOGENIZATION_none_ID, &
   HOMOGENIZATION_isostrain_ID, &
   HOMOGENIZATION_RGC_ID

 private :: &
   material_parseHomogenization, &
   material_parseMicrostructure, &
   material_parseCrystallite, &
   material_parsePhase, &
   material_parseTexture, &
   material_populateGrains

contains


!--------------------------------------------------------------------------------------------------
!> @brief parses material configuration file
!> @details figures out if solverJobName.materialConfig is present, if not looks for
!> material.config
!--------------------------------------------------------------------------------------------------
subroutine material_init
#if defined(PETSc) || defined(DAMASK_HDF5)
 use results
#endif
 use IO, only: &
   IO_error
 use debug, only: &
   debug_level, &
   debug_material, &
   debug_levelBasic, &
   debug_levelExtensive
 use mesh, only: &
   theMesh

 integer(pInt), parameter :: FILEUNIT = 210_pInt
 integer(pInt)            :: m,c,h, myDebug, myPhase, myHomog
 integer(pInt) :: &
  g, &                                                                                              !< grain number
  i, &                                                                                              !< integration point number
  e                                                                                                 !< element number
 integer(pInt), dimension(:), allocatable :: &
  CounterPhase, &
  CounterHomogenization

 myDebug = debug_level(debug_material)

 write(6,'(/,a)') ' <<<+-  material init  -+>>>'

 call material_parsePhase()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Phase          parsed'; flush(6)
 
 call material_parseMicrostructure()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Microstructure parsed'; flush(6)
 
 call material_parseCrystallite()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Crystallite    parsed'; flush(6)
 
 call material_parseHomogenization()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Homogenization parsed'; flush(6)
 
 call material_parseTexture()
 if (iand(myDebug,debug_levelBasic) /= 0_pInt) write(6,'(a)') ' Texture        parsed'; flush(6)

 allocate(plasticState       (size(config_phase)))
 allocate(sourceState        (size(config_phase)))
 do myPhase = 1,size(config_phase)
   allocate(sourceState(myPhase)%p(phase_Nsources(myPhase)))
 enddo

 allocate(homogState         (size(config_homogenization)))
 allocate(thermalState       (size(config_homogenization)))
 allocate(damageState        (size(config_homogenization)))

 allocate(thermalMapping     (size(config_homogenization)))
 allocate(damageMapping      (size(config_homogenization)))

 allocate(temperature        (size(config_homogenization)))
 allocate(damage             (size(config_homogenization)))

 allocate(temperatureRate    (size(config_homogenization)))

 do m = 1_pInt,size(config_microstructure)
   if(microstructure_crystallite(m) < 1_pInt .or. &
      microstructure_crystallite(m) > size(config_crystallite)) &
        call IO_error(150_pInt,m,ext_msg='crystallite')
   if(minval(microstructure_phase(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
      maxval(microstructure_phase(1:microstructure_Nconstituents(m),m)) > size(config_phase)) &
        call IO_error(150_pInt,m,ext_msg='phase')
   if(minval(microstructure_texture(1:microstructure_Nconstituents(m),m)) < 1_pInt .or. &
      maxval(microstructure_texture(1:microstructure_Nconstituents(m),m)) > size(config_texture)) &
        call IO_error(150_pInt,m,ext_msg='texture')
   if(microstructure_Nconstituents(m) < 1_pInt) &
        call IO_error(151_pInt,m)
 enddo

 debugOut: if (iand(myDebug,debug_levelExtensive) /= 0_pInt) then
   write(6,'(/,a,/)') ' MATERIAL configuration'
   write(6,'(a32,1x,a16,1x,a6)') 'homogenization                  ','type            ','grains'
   do h = 1_pInt,size(config_homogenization)
     write(6,'(1x,a32,1x,a16,1x,i6)') homogenization_name(h),homogenization_type(h),homogenization_Ngrains(h)
   enddo
   write(6,'(/,a14,18x,1x,a11,1x,a12,1x,a13)') 'microstructure','crystallite','constituents','homogeneous'
   do m = 1_pInt,size(config_microstructure)
     write(6,'(1x,a32,1x,i11,1x,i12,1x,l13)') microstructure_name(m), &
                                        microstructure_crystallite(m), &
                                        microstructure_Nconstituents(m), &
                                        microstructure_elemhomo(m)
     if (microstructure_Nconstituents(m) > 0_pInt) then
       do c = 1_pInt,microstructure_Nconstituents(m)
         write(6,'(a1,1x,a32,1x,a32,1x,f7.4)') '>',phase_name(microstructure_phase(c,m)),&
                                                   texture_name(microstructure_texture(c,m)),&
                                                   microstructure_fraction(c,m)
       enddo
       write(6,*)
     endif
   enddo
 endif debugOut

 call material_populateGrains
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! new mappings
 allocate(material_homogenizationAt,source=theMesh%homogenizationAt)
 allocate(material_homogenizationMemberAt(theMesh%elem%nIPs,theMesh%Nelems),source=0)

 allocate(CounterHomogenization(size(config_homogenization)),source=0)
 do e = 1, theMesh%Nelems
   do i = 1, theMesh%elem%nIPs
     CounterHomogenization(material_homogenizationAt(e)) = &
     CounterHomogenization(material_homogenizationAt(e)) + 1
     material_homogenizationMemberAt(i,e) = CounterHomogenization(material_homogenizationAt(e))
   enddo
 enddo


 allocate(material_phaseAt(homogenization_maxNgrains,theMesh%Nelems), source=material_phase(:,1,:))
 allocate(material_phaseMemberAt(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0)
 
 allocate(CounterPhase(size(config_phase)),source=0)
 do e = 1, theMesh%Nelems
   do i = 1, theMesh%elem%nIPs
     do c = 1, homogenization_maxNgrains
       CounterPhase(material_phaseAt(c,e)) = &
       CounterPhase(material_phaseAt(c,e)) + 1
       material_phaseMemberAt(c,i,e) = CounterPhase(material_phaseAt(c,e))
     enddo
   enddo
 enddo
 
#if defined(PETSc) || defined(DAMASK_HDF5)
 call results_openJobFile
 call results_mapping_constituent(material_phaseAt,material_phaseMemberAt,phase_name)
 call results_mapping_materialpoint(material_homogenizationAt,material_homogenizationMemberAt,homogenization_name)
 call results_closeJobFile
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN DEPRECATED
 allocate(phaseAt                   (  homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0_pInt)
 allocate(phasememberAt             (  homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0_pInt)
 allocate(mappingHomogenization     (2,                          theMesh%elem%nIPs,theMesh%Nelems),source=0_pInt)
 allocate(mappingHomogenizationConst(                            theMesh%elem%nIPs,theMesh%Nelems),source=1_pInt)
 
 CounterHomogenization=0
 CounterPhase         =0


 do e = 1_pInt,theMesh%Nelems
 myHomog = theMesh%homogenizationAt(e)
   do i = 1_pInt, theMesh%elem%nIPs
     CounterHomogenization(myHomog) = CounterHomogenization(myHomog) + 1_pInt
     mappingHomogenization(1:2,i,e) = [CounterHomogenization(myHomog),huge(1)]
     do g = 1_pInt,homogenization_Ngrains(myHomog)
       myPhase = material_phase(g,i,e)
       CounterPhase(myPhase) = CounterPhase(myPhase)+1_pInt                             ! not distinguishing between instances of same phase
       phaseAt(g,i,e)              = myPhase
       phasememberAt(g,i,e)        = CounterPhase(myPhase)
     enddo
   enddo
 enddo
! END DEPRECATED

! REMOVE !!!!!
! hack needed to initialize field values used during constitutive and crystallite initializations
 do myHomog = 1,size(config_homogenization)
   thermalMapping     (myHomog)%p => mappingHomogenizationConst
   damageMapping      (myHomog)%p => mappingHomogenizationConst
   allocate(temperature     (myHomog)%p(1), source=thermal_initialT(myHomog))
   allocate(damage          (myHomog)%p(1), source=damage_initialPhi(myHomog))
   allocate(temperatureRate (myHomog)%p(1), source=0.0_pReal)
 enddo

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization
 use mesh, only: &
   theMesh
 use IO, only: &
   IO_error

 integer(pInt)        :: h
 character(len=65536) :: tag

 allocate(homogenization_type(size(config_homogenization)),           source=HOMOGENIZATION_undefined_ID)
 allocate(thermal_type(size(config_homogenization)),                  source=THERMAL_isothermal_ID)
 allocate(damage_type (size(config_homogenization)),                  source=DAMAGE_none_ID)
 allocate(homogenization_typeInstance(size(config_homogenization)),   source=0_pInt)
 allocate(thermal_typeInstance(size(config_homogenization)),          source=0_pInt)
 allocate(damage_typeInstance(size(config_homogenization)),           source=0_pInt)
 allocate(homogenization_Ngrains(size(config_homogenization)),        source=0_pInt)
 allocate(homogenization_Noutput(size(config_homogenization)),        source=0_pInt)
 allocate(homogenization_active(size(config_homogenization)),         source=.false.)  !!!!!!!!!!!!!!!
 allocate(thermal_initialT(size(config_homogenization)),              source=300.0_pReal)
 allocate(damage_initialPhi(size(config_homogenization)),             source=1.0_pReal)

 forall (h = 1_pInt:size(config_homogenization)) &
   homogenization_active(h) = any(theMesh%homogenizationAt == h)


 do h=1_pInt, size(config_homogenization)
   homogenization_Noutput(h) = config_homogenization(h)%countKeys('(output)')

   tag = config_homogenization(h)%getString('mech')
   select case (trim(tag))
     case(HOMOGENIZATION_NONE_label)
       homogenization_type(h) = HOMOGENIZATION_NONE_ID
       homogenization_Ngrains(h) = 1_pInt
     case(HOMOGENIZATION_ISOSTRAIN_label)
       homogenization_type(h) = HOMOGENIZATION_ISOSTRAIN_ID
       homogenization_Ngrains(h) = config_homogenization(h)%getInt('nconstituents')
     case(HOMOGENIZATION_RGC_label)
       homogenization_type(h) = HOMOGENIZATION_RGC_ID
       homogenization_Ngrains(h) = config_homogenization(h)%getInt('nconstituents')
     case default
       call IO_error(500_pInt,ext_msg=trim(tag))
   end select
   
   homogenization_typeInstance(h) = count(homogenization_type==homogenization_type(h))

   if (config_homogenization(h)%keyExists('thermal')) then
     thermal_initialT(h) =  config_homogenization(h)%getFloat('t0',defaultVal=300.0_pReal)

     tag = config_homogenization(h)%getString('thermal')
     select case (trim(tag))
       case(THERMAL_isothermal_label)
         thermal_type(h) = THERMAL_isothermal_ID
       case(THERMAL_adiabatic_label)
         thermal_type(h) = THERMAL_adiabatic_ID
       case(THERMAL_conduction_label)
         thermal_type(h) = THERMAL_conduction_ID
       case default
         call IO_error(500_pInt,ext_msg=trim(tag))
     end select

   endif

   if (config_homogenization(h)%keyExists('damage')) then
     damage_initialPhi(h) =  config_homogenization(h)%getFloat('initialdamage',defaultVal=1.0_pReal)

     tag = config_homogenization(h)%getString('damage')
     select case (trim(tag))
       case(DAMAGE_NONE_label)
         damage_type(h) = DAMAGE_none_ID
       case(DAMAGE_LOCAL_label)
         damage_type(h) = DAMAGE_local_ID
       case(DAMAGE_NONLOCAL_label)
         damage_type(h) = DAMAGE_nonlocal_ID
       case default
         call IO_error(500_pInt,ext_msg=trim(tag))
     end select

   endif

 enddo

 do h=1_pInt, size(config_homogenization)
   homogenization_typeInstance(h)  = count(homogenization_type(1:h)  == homogenization_type(h))
   thermal_typeInstance(h)         = count(thermal_type       (1:h)  == thermal_type       (h))
   damage_typeInstance(h)          = count(damage_type        (1:h)  == damage_type        (h))
 enddo

 homogenization_maxNgrains = maxval(homogenization_Ngrains,homogenization_active)

end subroutine material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the microstructure part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseMicrostructure
 use IO, only: &
   IO_floatValue, &
   IO_intValue, &
   IO_stringValue, &
   IO_stringPos, &
   IO_error
 use mesh, only: &
   theMesh

 character(len=65536), dimension(:), allocatable :: &
   strings
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: e, m, c, i
 character(len=65536) :: &
   tag

 allocate(microstructure_crystallite(size(config_microstructure)),          source=0_pInt)
 allocate(microstructure_Nconstituents(size(config_microstructure)),        source=0_pInt)
 allocate(microstructure_active(size(config_microstructure)),               source=.false.)
 allocate(microstructure_elemhomo(size(config_microstructure)),             source=.false.)

 if(any(theMesh%microstructureAt > size(config_microstructure))) &
  call IO_error(155_pInt,ext_msg='More microstructures in geometry than sections in material.config')

 forall (e = 1_pInt:theMesh%Nelems) &
   microstructure_active(theMesh%microstructureAt(e)) = .true.                                         ! current microstructure used in model? Elementwise view, maximum N operations for N elements

 do m=1_pInt, size(config_microstructure)
   microstructure_Nconstituents(m) =  config_microstructure(m)%countKeys('(constituent)')
   microstructure_crystallite(m)   =  config_microstructure(m)%getInt('crystallite')
   microstructure_elemhomo(m)      =  config_microstructure(m)%keyExists('/elementhomogeneous/')
 enddo

 microstructure_maxNconstituents = maxval(microstructure_Nconstituents)
 allocate(microstructure_phase   (microstructure_maxNconstituents,size(config_microstructure)),source=0_pInt)
 allocate(microstructure_texture (microstructure_maxNconstituents,size(config_microstructure)),source=0_pInt)
 allocate(microstructure_fraction(microstructure_maxNconstituents,size(config_microstructure)),source=0.0_pReal)

 allocate(strings(1))                                                                               ! Intel 16.0 Bug
 do m=1_pInt, size(config_microstructure)
   strings = config_microstructure(m)%getStrings('(constituent)',raw=.true.)
   do c = 1_pInt, size(strings)
     chunkPos = IO_stringPos(strings(c))

     do i = 1_pInt,5_pInt,2_pInt
        tag = IO_stringValue(strings(c),chunkPos,i)

        select case (tag)
          case('phase')
            microstructure_phase(c,m) =    IO_intValue(strings(c),chunkPos,i+1_pInt)
          case('texture')
            microstructure_texture(c,m) =  IO_intValue(strings(c),chunkPos,i+1_pInt)
          case('fraction')
            microstructure_fraction(c,m) =  IO_floatValue(strings(c),chunkPos,i+1_pInt)
        end select
     
     enddo
   enddo
 enddo

 do m = 1_pInt, size(config_microstructure)
   if (dNeq(sum(microstructure_fraction(:,m)),1.0_pReal)) &
     call IO_error(153_pInt,ext_msg=microstructure_name(m))
 enddo
 
end subroutine material_parseMicrostructure


!--------------------------------------------------------------------------------------------------
!> @brief parses the crystallite part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseCrystallite

 integer(pInt)        :: c

 allocate(crystallite_Noutput(size(config_crystallite)),source=0_pInt)
 do c=1_pInt, size(config_crystallite)
   crystallite_Noutput(c) =  config_crystallite(c)%countKeys('(output)')
 enddo

end subroutine material_parseCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief parses the phase part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parsePhase
 use IO, only: &
   IO_error, &
   IO_getTag, &
   IO_stringValue

 integer(pInt) :: sourceCtr, kinematicsCtr, stiffDegradationCtr, p
 character(len=65536), dimension(:), allocatable ::  str 


 allocate(phase_elasticity(size(config_phase)),source=ELASTICITY_undefined_ID)
 allocate(phase_plasticity(size(config_phase)),source=PLASTICITY_undefined_ID)
 allocate(phase_Nsources(size(config_phase)),              source=0_pInt)
 allocate(phase_Nkinematics(size(config_phase)),           source=0_pInt)
 allocate(phase_NstiffnessDegradations(size(config_phase)),source=0_pInt)
 allocate(phase_Noutput(size(config_phase)),               source=0_pInt)
 allocate(phase_localPlasticity(size(config_phase)),       source=.false.)

 do p=1_pInt, size(config_phase)
   phase_Noutput(p) =                 config_phase(p)%countKeys('(output)')
   phase_Nsources(p) =                config_phase(p)%countKeys('(source)')
   phase_Nkinematics(p) =             config_phase(p)%countKeys('(kinematics)')
   phase_NstiffnessDegradations(p) =  config_phase(p)%countKeys('(stiffness_degradation)')
   phase_localPlasticity(p) = .not.   config_phase(p)%KeyExists('/nonlocal/')

   select case (config_phase(p)%getString('elasticity'))
     case (ELASTICITY_HOOKE_label)
       phase_elasticity(p) = ELASTICITY_HOOKE_ID
     case default
       call IO_error(200_pInt,ext_msg=trim(config_phase(p)%getString('elasticity')))
   end select

   select case (config_phase(p)%getString('plasticity'))
     case (PLASTICITY_NONE_label)
       phase_plasticity(p) = PLASTICITY_NONE_ID
     case (PLASTICITY_ISOTROPIC_label)
       phase_plasticity(p) = PLASTICITY_ISOTROPIC_ID
     case (PLASTICITY_PHENOPOWERLAW_label)
       phase_plasticity(p) = PLASTICITY_PHENOPOWERLAW_ID
     case (PLASTICITY_KINEHARDENING_label)
       phase_plasticity(p) = PLASTICITY_KINEHARDENING_ID
     case (PLASTICITY_DISLOTWIN_label)
       phase_plasticity(p) = PLASTICITY_DISLOTWIN_ID
     case (PLASTICITY_DISLOUCLA_label)
       phase_plasticity(p) = PLASTICITY_DISLOUCLA_ID
     case (PLASTICITY_NONLOCAL_label)
       phase_plasticity(p) = PLASTICITY_NONLOCAL_ID
     case default
       call IO_error(201_pInt,ext_msg=trim(config_phase(p)%getString('plasticity')))
   end select

 enddo

 allocate(phase_source(maxval(phase_Nsources),size(config_phase)), source=SOURCE_undefined_ID)
 allocate(phase_kinematics(maxval(phase_Nkinematics),size(config_phase)), source=KINEMATICS_undefined_ID)
 allocate(phase_stiffnessDegradation(maxval(phase_NstiffnessDegradations),size(config_phase)), &
          source=STIFFNESS_DEGRADATION_undefined_ID)
 do p=1_pInt, size(config_phase)
#if defined(__GFORTRAN__) || defined(__PGI)
   str = ['GfortranBug86277']
   str = config_phase(p)%getStrings('(source)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_phase(p)%getStrings('(source)',defaultVal=[character(len=65536)::])
#endif
   do sourceCtr = 1_pInt, size(str)
     select case (trim(str(sourceCtr)))
       case (SOURCE_thermal_dissipation_label)
         phase_source(sourceCtr,p) = SOURCE_thermal_dissipation_ID
       case (SOURCE_thermal_externalheat_label)
         phase_source(sourceCtr,p) = SOURCE_thermal_externalheat_ID
       case (SOURCE_damage_isoBrittle_label)
         phase_source(sourceCtr,p) = SOURCE_damage_isoBrittle_ID
       case (SOURCE_damage_isoDuctile_label)
         phase_source(sourceCtr,p) = SOURCE_damage_isoDuctile_ID
       case (SOURCE_damage_anisoBrittle_label)
         phase_source(sourceCtr,p) = SOURCE_damage_anisoBrittle_ID
       case (SOURCE_damage_anisoDuctile_label)
         phase_source(sourceCtr,p) = SOURCE_damage_anisoDuctile_ID
     end select
   enddo

#if defined(__GFORTRAN__) || defined(__PGI)
   str = ['GfortranBug86277']
   str = config_phase(p)%getStrings('(kinematics)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_phase(p)%getStrings('(kinematics)',defaultVal=[character(len=65536)::])
#endif
   do kinematicsCtr = 1_pInt, size(str)
     select case (trim(str(kinematicsCtr)))
       case (KINEMATICS_cleavage_opening_label)
         phase_kinematics(kinematicsCtr,p) = KINEMATICS_cleavage_opening_ID
       case (KINEMATICS_slipplane_opening_label)
         phase_kinematics(kinematicsCtr,p) = KINEMATICS_slipplane_opening_ID
       case (KINEMATICS_thermal_expansion_label)
         phase_kinematics(kinematicsCtr,p) = KINEMATICS_thermal_expansion_ID
     end select
   enddo
#if defined(__GFORTRAN__) || defined(__PGI)
   str = ['GfortranBug86277']
   str = config_phase(p)%getStrings('(stiffness_degradation)',defaultVal=str)
   if (str(1) == 'GfortranBug86277') str = [character(len=65536)::]
#else
   str = config_phase(p)%getStrings('(stiffness_degradation)',defaultVal=[character(len=65536)::])
#endif
   do stiffDegradationCtr = 1_pInt, size(str)
     select case (trim(str(stiffDegradationCtr)))
       case (STIFFNESS_DEGRADATION_damage_label)
         phase_stiffnessDegradation(stiffDegradationCtr,p) = STIFFNESS_DEGRADATION_damage_ID
    end select
   enddo
 enddo

 allocate(phase_plasticityInstance(size(config_phase)),   source=0_pInt)
 allocate(phase_elasticityInstance(size(config_phase)),   source=0_pInt)

 do p=1_pInt, size(config_phase)
   phase_elasticityInstance(p)  = count(phase_elasticity(1:p)  == phase_elasticity(p))
   phase_plasticityInstance(p)  = count(phase_plasticity(1:p)  == phase_plasticity(p))
 enddo

end subroutine material_parsePhase

!--------------------------------------------------------------------------------------------------
!> @brief parses the texture part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseTexture
 use IO, only: &
   IO_error, &
   IO_stringPos, &
   IO_floatValue, &
   IO_stringValue

 integer(pInt) :: section, gauss, j, t, i
 character(len=65536), dimension(:), allocatable ::  strings                                     ! Values for given key in material config 
 integer(pInt), dimension(:), allocatable :: chunkPos

 allocate(texture_Ngauss(size(config_texture)),   source=0_pInt)

 do t=1_pInt, size(config_texture)
   texture_Ngauss(t) =  config_texture(t)%countKeys('(gauss)')
   if (config_texture(t)%keyExists('symmetry')) call IO_error(147,ext_msg='symmetry')
   if (config_texture(t)%keyExists('(random)')) call IO_error(147,ext_msg='(random)')
   if (config_texture(t)%keyExists('(fiber)'))  call IO_error(147,ext_msg='(fiber)')
 enddo

 texture_maxNgauss = maxval(texture_Ngauss)
 allocate(texture_Gauss (5,texture_maxNgauss,size(config_texture)), source=0.0_pReal)
 allocate(texture_transformation(3,3,size(config_texture)),         source=0.0_pReal)
          texture_transformation = spread(math_I3,3,size(config_texture))

 do t=1_pInt, size(config_texture)
   section = t
   gauss = 0_pInt
   
   if (config_texture(t)%keyExists('axes')) then
     strings = config_texture(t)%getStrings('axes')
     do j = 1_pInt, 3_pInt                                                                          ! look for "x", "y", and "z" entries
       select case (strings(j))
         case('x', '+x')
           texture_transformation(j,1:3,t) = [ 1.0_pReal, 0.0_pReal, 0.0_pReal]                     ! original axis is now +x-axis
         case('-x')
           texture_transformation(j,1:3,t) = [-1.0_pReal, 0.0_pReal, 0.0_pReal]                     ! original axis is now -x-axis
         case('y', '+y')
           texture_transformation(j,1:3,t) = [ 0.0_pReal, 1.0_pReal, 0.0_pReal]                     ! original axis is now +y-axis
         case('-y')
           texture_transformation(j,1:3,t) = [ 0.0_pReal,-1.0_pReal, 0.0_pReal]                     ! original axis is now -y-axis
         case('z', '+z')
           texture_transformation(j,1:3,t) = [ 0.0_pReal, 0.0_pReal, 1.0_pReal]                     ! original axis is now +z-axis
         case('-z')
           texture_transformation(j,1:3,t) = [ 0.0_pReal, 0.0_pReal,-1.0_pReal]                     ! original axis is now -z-axis
         case default
           call IO_error(157_pInt,t)
       end select
     enddo
     if(dNeq(math_det33(texture_transformation(1:3,1:3,t)),1.0_pReal)) call IO_error(157_pInt,t)
   endif

   if (config_texture(t)%keyExists('(gauss)')) then
     gauss = gauss + 1_pInt
     strings = config_texture(t)%getStrings('(gauss)',raw= .true.)
     do i = 1_pInt , size(strings)
       chunkPos = IO_stringPos(strings(i))
       do j = 1_pInt,9_pInt,2_pInt
         select case (IO_stringValue(strings(i),chunkPos,j))
             case('phi1')
                 texture_Gauss(1,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
             case('phi')
                 texture_Gauss(2,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
             case('phi2')
                 texture_Gauss(3,gauss,t) = IO_floatValue(strings(i),chunkPos,j+1_pInt)*inRad
          end select
      enddo
     enddo
   endif
 enddo    
 
 call config_deallocate('material.config/texture')

end subroutine material_parseTexture


!--------------------------------------------------------------------------------------------------
!> @brief allocates the plastic state of a phase
!--------------------------------------------------------------------------------------------------
subroutine material_allocatePlasticState(phase,NofMyPhase,&
                                         sizeState,sizeDotState,sizeDeltaState,&
                                         Nslip,Ntwin,Ntrans)
 use numerics, only: &
   numerics_integrator

 integer(pInt), intent(in) :: &
   phase, &
   NofMyPhase, &
   sizeState, &
   sizeDotState, &
   sizeDeltaState, &
   Nslip, &
   Ntwin, &
   Ntrans

 plasticState(phase)%sizeState        = sizeState
 plasticState(phase)%sizeDotState     = sizeDotState
 plasticState(phase)%sizeDeltaState   = sizeDeltaState
 plasticState(phase)%offsetDeltaState = sizeState-sizeDeltaState                                    ! deltaState occupies latter part of state by definition
 plasticState(phase)%Nslip = Nslip
 plasticState(phase)%Ntwin = Ntwin
 plasticState(phase)%Ntrans= Ntrans

 allocate(plasticState(phase)%aTolState           (sizeState),                source=0.0_pReal)
 allocate(plasticState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(plasticState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(plasticState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(plasticState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

 allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 1_pInt) then
   allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
   allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
 endif
 if (numerics_integrator == 4_pInt) &
   allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 5_pInt) &
   allocate(plasticState(phase)%RKCK45dotState  (6,sizeDotState,NofMyPhase),  source=0.0_pReal)

 allocate(plasticState(phase)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)

end subroutine material_allocatePlasticState


!--------------------------------------------------------------------------------------------------
!> @brief allocates the source state of a phase
!--------------------------------------------------------------------------------------------------
subroutine material_allocateSourceState(phase,of,NofMyPhase,&
                                        sizeState,sizeDotState,sizeDeltaState)
 use numerics, only: &
   numerics_integrator

 integer(pInt), intent(in) :: &
   phase, &
   of, &
   NofMyPhase, &
   sizeState, sizeDotState,sizeDeltaState

 sourceState(phase)%p(of)%sizeState        = sizeState
 sourceState(phase)%p(of)%sizeDotState     = sizeDotState
 sourceState(phase)%p(of)%sizeDeltaState   = sizeDeltaState
 sourceState(phase)%p(of)%offsetDeltaState = sizeState-sizeDeltaState                               ! deltaState occupies latter part of state by definition

 allocate(sourceState(phase)%p(of)%aTolState           (sizeState),                source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
 allocate(sourceState(phase)%p(of)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

 allocate(sourceState(phase)%p(of)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 1_pInt) then
   allocate(sourceState(phase)%p(of)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
   allocate(sourceState(phase)%p(of)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
 endif
 if (numerics_integrator == 4_pInt) &
   allocate(sourceState(phase)%p(of)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
 if (numerics_integrator == 5_pInt) &
   allocate(sourceState(phase)%p(of)%RKCK45dotState  (6,sizeDotState,NofMyPhase),  source=0.0_pReal)

 allocate(sourceState(phase)%p(of)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)

end subroutine material_allocateSourceState


!--------------------------------------------------------------------------------------------------
!> @brief populates the grains
!> @details populates the grains by identifying active microstructure/homogenization pairs,
!! calculates the volume of the grains and deals with texture components
!--------------------------------------------------------------------------------------------------
subroutine material_populateGrains
 use mesh, only: &
   theMesh

 integer(pInt) :: e,i,c,homog,micro

 allocate(material_phase(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),        source=0_pInt)
 allocate(material_texture(homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),      source=0_pInt)
 allocate(material_EulerAngles(3,homogenization_maxNgrains,theMesh%elem%nIPs,theMesh%Nelems),source=0.0_pReal)

  do e = 1, theMesh%Nelems
    do i = 1, theMesh%elem%nIPs
      homog = theMesh%homogenizationAt(e)
      micro = theMesh%microstructureAt(e)
      do c = 1, homogenization_Ngrains(homog)
        material_phase(c,i,e)   = microstructure_phase(c,micro)
        material_texture(c,i,e) = microstructure_texture(c,micro)
        material_EulerAngles(1:3,c,i,e) = texture_Gauss(1:3,1,material_texture(c,i,e))
        material_EulerAngles(1:3,c,i,e) = math_RtoEuler( &                                       ! translate back to Euler angles
                                             math_mul33x33( &                                       ! pre-multiply
                                               math_EulertoR(material_EulerAngles(1:3,c,i,e)), &    ! face-value orientation
                                               texture_transformation(1:3,1:3,material_texture(c,i,e)) &          ! and transformation matrix
                                             ) &
                                             )
      enddo
    enddo
  enddo

 deallocate(texture_transformation)

 call config_deallocate('material.config/microstructure')

end subroutine material_populateGrains

end module material
