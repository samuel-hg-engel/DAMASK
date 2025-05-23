!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Su Leen Wong, Max-Planck-Institut für Eisenforschung GmbH
!> @author Nan Jia, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) dislobasic

  type :: tParameters
    real(pREAL),               allocatable, dimension(:) :: &
      b_sl, &                                                                                       !< magnitude of Burgers vector (m)
      delta_F, &                                                                                     !< activation energy for glide (J)
      k_1, &                                                                                        !< Dislocation multiplication
      alpha_n, &                                                                                    !< Slip-system interaction strength
      tau_0, &                                                                                      !< Intrinsic strength
      k_2, &                                                                                        !< Dislocation annhilation
      nu_g, &                                                                                       !< Dislocation jump frequency
      delta_V, &                                                                                    !< Dislocation activation volume
      rho_mob_0                                                                                     !< Mobile Dislocation Density
    real(pREAL),               allocatable, dimension(:,:) :: &
      forestProjection, &
      h_sl_sl                                                                                       !< components of slip-slip interaction matrix
    real(pREAL),               allocatable, dimension(:,:,:) :: &
      P_sl
    integer :: &
      sum_N_sl                                                                                      !< total number of active slip systems
    character(len=:),          allocatable :: &
      isotropic_bound
    character(len=pSTRLEN),    allocatable, dimension(:) :: &
      output
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl
  end type tParameters                                                                              !< container type for internal constitutive parameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      rho_ssd, &
      gamma_sl
  end type tIndexDotState

  type :: tDislobasicState
    real(pREAL),                  dimension(:,:),   pointer :: &
      rho_ssd, &
      gamma_sl
  end type tDislobasicState

  type :: tDislobasicDependentState
    real(pREAL),                  dimension(:,:),   allocatable :: &
      tau_pass                                                                                      !< threshold stress for slip
  end type tDislobasicDependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),              allocatable, dimension(:) :: param
  type(tIndexDotState),           allocatable, dimension(:) :: indexDotState
  type(tDislobasicState),          allocatable, dimension(:) :: state
  type(tDislobasicDependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_dislobasic_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl
  real(pREAL), allocatable, dimension(:) :: &
    f_edge, &
    rho_ssd_0                                                                                       !< initial SSD dislocation density per slip system
  character(len=:), allocatable :: &
    refs, &
    extmsg
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('dislobasic')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:dislobasic init  -+>>>'

  print'(/,1x,a)', 'A. Ma and F. Roters, Acta Materialia 52(12):3603–3612, 2004'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2004.04.012'

  print'(/,1x,a)', 'F. Roters et al., Computational Materials Science 39:91–95, 2007'
  print'(  1x,a)', 'https://doi.org/10.1016/j.commatsci.2006.04.014'

  print'(/,1x,a)', 'S.L. Wong et al., Acta Materialia 118:140–151, 2016'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2016.07.032'

  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  allocate(dependentState(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), &
              stt => state(ph), dst => dependentState(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get_dict(ph)
    mech  => phase%get_dict('mechanical')
    pl    => mech%get_dict('plastic')

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(pl,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs

#if defined (__GFORTRAN__)
    prm%output = output_as1dStr(pl)
#else
    prm%output = pl%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

   prm%isotropic_bound = pl%get_asStr('isotropic_bound',defaultVal='isostrain')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = crystal_labels_slip(N_sl,phase_lattice(ph))
      prm%P_sl       = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))
      f_edge       = math_expand(pl%get_as1dReal('f_edge',    requiredSize=size(N_sl), &
                                                 defaultVal=[(0.5_pREAL,i=1,size(N_sl))]),N_sl)
#ifdef __GFORTRAN__
      rho_ssd_0    = pl%get_as1dReal('rho_ssd_0', requiredChunks=N_sl)
#else
      rho_ssd_0    = math_expand(pl%get_as1dReal('rho_ssd_0', requiredSize=size(N_sl)),N_sl)
#endif

      prm%b_sl      = math_expand(pl%get_as1dReal('b_sl',      requiredSize=size(N_sl)),N_sl)
      prm%delta_F   = math_expand(pl%get_as1dReal('delta_F',   requiredSize=size(N_sl)),N_sl)
      prm%k_1       = math_expand(pl%get_as1dReal('k_1',       requiredSize=size(N_sl)),N_sl)
      prm%tau_0     = math_expand(pl%get_as1dReal('tau_0',     requiredSize=size(N_sl)),N_sl)
      prm%k_2       = math_expand(pl%get_as1dReal('k_2',       requiredSize=size(N_sl)),N_sl)
      prm%alpha_n   = math_expand(pl%get_as1dReal('alpha_n',   requiredSize=size(N_sl)),N_sl)
      prm%nu_g      = math_expand(pl%get_as1dReal('nu_g',      requiredSize=size(N_sl)),N_sl)
      prm%delta_V   = math_expand(pl%get_as1dReal('delta_V',   requiredSize=size(N_sl)),N_sl)
      prm%rho_mob_0 = math_expand(pl%get_as1dReal('rho_mob_0', requiredSize=size(N_sl)),N_sl)

      prm%h_sl_sl = crystal_interaction_SlipBySlip(N_sl,pl%get_as1dReal('h_sl-sl'),phase_lattice(ph))

      prm%forestProjection = spread(          f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_edge (N_sl,phase_lattice(ph),phase_cOverA(ph)) &
                           + spread(1.0_pREAL-f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_screw(N_sl,phase_lattice(ph),phase_cOverA(ph))

      ! sanity checks
      if (any(rho_ssd_0         <  0.0_pREAL))         extmsg = trim(extmsg)//' rho_ssd_0'
      if (any(prm%rho_mob_0     <  0.0_pREAL))         extmsg = trim(extmsg)//' rho_mob_0'
      if (any(prm%b_sl          <= 0.0_pREAL))         extmsg = trim(extmsg)//' b_sl'
      if (any(prm%delta_F       <= 0.0_pREAL))         extmsg = trim(extmsg)//' delta_F'
      if (any(prm%k_1           <= 0.0_pREAL))         extmsg = trim(extmsg)//' k_1'
      if (any(prm%alpha_n       <= 0.0_pREAL))         extmsg = trim(extmsg)//' alpha_n'
      if (any(prm%k_2           <  0.0_pREAL))         extmsg = trim(extmsg)//' k_2'
      if (any(prm%nu_g          <  0.0_pREAL))         extmsg = trim(extmsg)//' n_g'
      if (any(prm%delta_V       <  0.0_pREAL))         extmsg = trim(extmsg)//' delta_V'

    else slipActive
      rho_ssd_0 = emptyRealArray
      allocate(prm%b_sl, &
               prm%delta_F, &
               prm%delta_V, &
               prm%k_1, &
               prm%tau_0, &
               prm%k_2, &
               prm%alpha_n, &
               prm%nu_g, &
               prm%rho_mob_0, &
               source=emptyRealArray)
      allocate(prm%forestProjection(0,0), &
               prm%h_sl_sl(0,0))

    end if slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers  = count(material_ID_phase == ph)
    sizeDotState = size(['rho_ssd ','gamma_sl']) * prm%sum_N_sl

    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%rho_ssd = [startIndex,endIndex]
    stt%rho_ssd => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_ssd = spread(rho_ssd_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_gamma'

    allocate(dst%tau_pass (prm%sum_N_sl,Nmembers),source=0.0_pREAL)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_dislobasic_init

!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
module subroutine dislobasic_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pREAL), dimension(3,3),     intent(out) :: Lp
  real(pREAL), dimension(3,3,3,3), intent(out) :: dLp_dMp
  real(pREAL), dimension(3,3),     intent(in)  :: Mp
  integer,                         intent(in)  :: ph,en

  integer :: i,k,l,m,n
  real(pREAL) :: &
    T
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl,ddot_gamma_dtau_sl

  T = thermal_T(ph,en)
  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph), stt => state(ph))

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl,ddot_gamma_dtau_sl)
    slipContribution: do i = 1, prm%sum_N_sl
      Lp = Lp + dot_gamma_sl(i)*prm%P_sl(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_sl(i) * prm%P_sl(k,l,i) * prm%P_sl(m,n,i)
    end do slipContribution

    end associate

end subroutine dislobasic_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function dislobasic_dotState(Mp,ph,en) result(dotState)

  real(pREAL), dimension(3,3),  intent(in):: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl
  real(pREAL) :: &
    mu, nu, &
    T

  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph), &
            dot_rho_ssd => dotState(indexDotState(ph)%rho_ssd(1):indexDotState(ph)%rho_ssd(2)), &
            abs_dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)))

    mu = elastic_mu(ph,en,prm%isotropic_bound)
    nu = elastic_nu(ph,en,prm%isotropic_bound)
    T = thermal_T(ph,en)

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl)
    abs_dot_gamma_sl = abs(dot_gamma_sl)

    !dot_rho_ssd = abs_dot_gamma_sl * (prm%k_1 / prm%b_sl * prm%alpha_n * sqrt(matmul(prm%forestProjection,stt%rho_ssd(:,en)))) &
    !            - abs_dot_gamma_sl * (prm%k_2 * stt%rho_ssd(:,en))

    dot_rho_ssd = abs_dot_gamma_sl * (prm%k_1 / prm%b_sl * sqrt(matmul(prm%forestProjection,stt%rho_ssd(:,en)))) &
                - abs_dot_gamma_sl * (prm%k_2 * stt%rho_ssd(:,en))

  end associate

end function dislobasic_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine dislobasic_dependentState(ph,en)

  integer,       intent(in) :: &
    ph, &
    en
  real(pREAL) :: &
    mu


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    mu = elastic_mu(ph,en,prm%isotropic_bound)

    !* threshold stress for dislocation motion
    dst%tau_pass(:,en) = prm%tau_0 + mu * prm%b_sl * prm%alpha_n * sqrt(matmul(prm%forestProjection,stt%rho_ssd(:,en)))
    !dst%tau_pass(:,en) = prm%tau_0 + mu * prm%b_sl * prm%alpha_n * sqrt(matmul(prm%h_sl_sl,stt%rho_ssd(:,en)))
    !dst%tau_pass(:,en) = prm%tau_0 + mu * prm%b_sl * prm%alpha_n * sqrt(stt%rho_ssd(:,en))

  end associate

end subroutine dislobasic_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dislobasic_result(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('rho_ssd')
          call result_writeDataset(stt%rho_ssd,group,trim(prm%output(ou)), &
                                   'SSD density','1/m²',prm%systems_sl)
        case('gamma_sl')
          call result_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                   'plastic shear','1',prm%systems_sl)
        case('tau_pass')
          call result_writeDataset(dst%tau_pass,group,trim(prm%output(ou)), &
                                   'passing stress for slip','Pa',prm%systems_sl)
      end select

    end do

  end associate

end subroutine plastic_dislobasic_result

!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems, their derivatives with respect to resolved
!         stress, and the resolved stress.
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_sl(Mp,T,ph,en, &
                            dot_gamma_sl,ddot_gamma_dtau_sl,tau_sl)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pREAL),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(param(ph)%sum_N_sl), intent(out) :: &
    dot_gamma_sl
  real(pREAL), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
    ddot_gamma_dtau_sl, &
    tau_sl
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    ddot_gamma_dtau
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau, &
    v_g, &
    tau_eff, &                                                                                      !< effective resolved stress
    dv_g_dtau
  integer :: i

  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    tau = [(math_tensordot(Mp,prm%P_sl(1:3,1:3,i)),i = 1, prm%sum_N_sl)]

    tau_eff = abs(tau)-dst%tau_pass(:,en)

    significantStress: where(tau_eff > tol_math_check)

      v_g = prm%nu_g * prm%b_sl * exp(-1.0_pREAL*prm%delta_F/(K_B*T)) * sinh((tau_eff * prm%delta_V) /(K_B*T))

      dot_gamma_sl = sign(prm%rho_mob_0 * prm%b_sl * v_g, tau)

      dv_g_dtau = prm%nu_g * prm%b_sl * exp(-1.0_pREAL*prm%delta_F/(K_B*T)) * cosh((tau_eff * prm%delta_V) /(K_B*T)) * prm%delta_V / (K_B*T)

      ddot_gamma_dtau = prm%rho_mob_0 * prm%b_sl * dv_g_dtau

    else where significantStress
      dot_gamma_sl    = 0.0_pREAL
      ddot_gamma_dtau = 0.0_pREAL
    end where significantStress

  end associate

  if (present(ddot_gamma_dtau_sl)) ddot_gamma_dtau_sl = ddot_gamma_dtau
  if (present(tau_sl))             tau_sl             = tau

end subroutine kinetics_sl

end submodule dislobasic