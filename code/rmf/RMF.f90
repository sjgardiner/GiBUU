!******************************************************************************
!****m* /RMF
! NAME
! module RMF
! PURPOSE
! Includes all information about relativistic mean-field potential for baryons
! and mesons.
! NOTES
! * When hyperon coupling is scaled by the well known factor of 2/3, the kaons
!   are scaled by the factor 1/3 in order to compensate the missing self energy
!   between incoming and outgoing channel.
!   This is because the threshold condition, e.g. for pi N->YK,
!   sqrt(s*)>m*_y+m*_k, in the medium assumes no changes in the self energy
!   between initial and final states
!   (cf. https://inspirehep.net/literature/748586)
! * This prescription should better not be used at energies
!   near the kaon-production threshold, since in this method the kaon potential
!   is not consistent within Chiral Perturbation Theory or
!   One-Boson-Exchange models.
! * The same non-trivial feature appears if the baryon self energies depend on
!   isospin.
!   Presently no isospin-dependent part is included in the baryon fields.
! * Going beyond this simple approximation means to explicitly include different
!   threshold conditions for all channels considered in the collision term.
!******************************************************************************
module RMF

  use CallStack, only: TRACEBACK


  implicit none
  private

  !****************************************************************************
  !****g* RMF/RMF_flag
  ! SOURCE
  !
  logical, save :: RMF_flag = .false.
  ! PURPOSE
  ! If .true. then use relativistic mean fields.
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/N_set
  ! SOURCE
  !
  integer, save :: N_set = 1
  ! PURPOSE
  ! Select parameter set to use:
  ! *  1 --- NL1 [Lalazissis] (K=211.29 MeV, m*/m=0.57)
  ! *  2 --- NL3 [Lalazissis] (K=271.76 MeV, m*/m=0.60)
  ! *  3 --- NL2 [Lang]       (K=210 MeV,    m*/m=0.83)
  ! *  4 --- NLZ2 [Bender]    (K=172 MeV,    m*/m=0.583)
  ! *  5 --- NL3* [Lalazissis, priv. comm.] (K=258.28 MeV, m*/m=0.594)
  ! *  6 --- Same as N_set=3, but including the rho meson.
  ! *  7 --- NL1 [Lee]        (K=212 MeV,    m*/m=0.57)
  ! *  8 --- NL2 [Lee]        (K=399 MeV,    m*/m=0.67)
  ! *  9 --- Set I [ Liu]     (K=240 MeV,    m*/m=0.75)
  ! * 10 --- NL1 [Lang]       (K=380 MeV,    m*/m=0.83)
  ! * 11 --- NL3 [Lang]       (K=380 MeV,    m*/m=0.70)
  ! * 31 --- Parity doublet model Set P3 [Zschiesche] (K=374 MeV)
  ! * 32 --- Parity doublet model Set P2 [Zschiesche] (K=374 MeV)
  ! * 33 --- Parity doublet model Set 1 [Shin] (K=240 MeV)
  ! * 34 --- Parity doublet model Set 2 [Shin] (K=215 MeV)
  !
  ! References:
  ! * Bender et al., PRC 60, 34304 (1999)
  ! * Lalazissis et al., PRC 55, 540 (1997),
  ! * Lang et al., NPA 541, 507 (1992)
  ! * Lee et al., PRL 57, 2916 (1986)
  ! * Liu et al., PRC 65, 045201 (2002)
  ! * Shin et al., arXiv:1805.03402
  ! * Zschiesche et al., PRC 75, 055202 (2007)
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/grad_flag
  ! SOURCE
  !
  logical, save, public :: grad_flag = .false.
  ! PURPOSE
  ! If .true. then include space derivatives of the fields.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/lorentz_flag
  ! SOURCE
  !
  logical, save, public :: lorentz_flag = .true.
  ! PURPOSE
  ! If .false. then the space components of the omega and rho fields are put
  ! to zero.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/Tens_flag
  ! SOURCE
  !
  logical, save, public :: Tens_flag = .false.
  ! PURPOSE
  ! If .true. then compute the energy-momentum tensor and four-momentum density field
  ! (not used in propagation).
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/flagCorThr
  ! SOURCE
  !
  logical, save, public :: flagCorThr=.false.
  ! PURPOSE
  ! If .true. then the srtfree of colliding particles is corrected
  ! to ensure in-medium thresholds of BB -> BB and MB -> B
  !****************************************************************************

  !*** Modification factors for the coupling constants with mean meson fields:

  !****************************************************************************
  !****g* RMF/fact_pbar
  ! SOURCE
  real, save :: fact_pbar    = 1.
  ! PURPOSE
  ! Modification factor for the antiproton coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_Delta
  ! SOURCE
  real, save :: fact_Delta    = 1.
  ! PURPOSE
  ! Modification factor for the Delta(1232) coupling constants.
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/fact_hyp
  ! SOURCE
  real, save :: fact_hyp     = 1.
  ! PURPOSE
  ! Modification factor for the hyperon coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_antihyp
  ! SOURCE
  real, save :: fact_antihyp = 1.
  ! PURPOSE
  ! Modification factor for the antihyperon coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_Xi
  ! SOURCE
  real, save :: fact_Xi      = 1.
  ! PURPOSE
  ! Modification factor for the Xi and XiStar coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_antiXi
  ! SOURCE
  real, save :: fact_antiXi  = 1.
  ! PURPOSE
  ! Modification factor for the antiXi and antiXiStar coupling constants.
  !****************************************************************************

  !****************************************************************************
  !****g* RMF/fact_kaon
  ! SOURCE
  real, save :: fact_kaon    = 0.
  ! PURPOSE
  ! Modification factor for the Kaon and antikaon coupling constants.
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/kaonpot_flag
  ! SOURCE
  !
  logical, save, public :: kaonpot_flag = .false.
  ! PURPOSE
  ! This switch turns on the Kaon potential in RMF mode.
  !****************************************************************************


  !****************************************************************************
  !****g* RMF/flagVectMod
  ! SOURCE
  !
  logical, save, public :: flagVectMod=.true.
  ! PURPOSE
  ! This switch turns on the modification factors for vector couplings.
  !****************************************************************************

  !**** variables related to the kaon-nucleon potential in RMF mode:
  real, parameter, public :: gs_kaon = 52.0291   ! Sigma_KN/f_pi^2 in GeV^-1 (with Sigma_KN=0.450 GeV and f_pi=0.093 GeV)
  real, parameter, public :: gv_kaon = 72.2627   ! 3/(8f_pi*^2) in GeV^-2 (with f_pi=sqrt(0.6)*f_pi)

  real, parameter, public :: Sigma_N=0.045     ! nucleon sigma term (GeV) needed for Feynman-Hellmann theorem
  ! used to initialize sigma field in the PDM

  !**** Mean field parameters in the notations of Ref.
  !**** G.A. Lalazissis et al., PRC 55, 540 (1997):

  real, public, save :: m_nucleon   ! -- nucleon mass (GeV),
  real, public, save :: m_sigma     ! -- sigma-meson mass (GeV),
  real, public, save :: m_omega     ! -- omega-meson mass (GeV),
  real, public, save :: m_rho       ! -- rho-meson mass (GeV),
  real, public, save :: g_sigma     ! -- sigma-nucleon coupling constant,
  real, public, save :: g_omega     ! -- omega-nucleon coupling constant,
  real, public, save :: g_rho=0.    ! -- rho-nucleon coupling constant,
  real, public, save :: g_2 ! -- coefficient at sigma^3 in the Lagrangian (GeV),
  real, public, save :: g_3 ! -- coefficient at sigma^4 in the Lagrangian

  real, public, save :: a_0, a_1, a_2, a_3, a_4, a_5, a_6, a_7    ! auxiliary parameters

  !***  Additional parameters needed for the parity doublet (mirror) model in the notations of Ref.
  !***  D. Jido et al., Prog. Theor. Phys. 106, 873 (2001)

  real, public, save :: m_minus     ! -- mass of the partner state with negative parity (GeV)
  real, public, save :: m0          ! -- chiral invariant mass (GeV)
  real, public, save :: g1          ! -- sigma N1 coupling constant
  real, public, save :: g2          ! -- sigma N2 coupling constant


  !***  Parameters of the meson lagrangian of the parity doublet model
  !*** (see concrete parameterizations below)
  real, public, save :: mu      ! -- 'mass parameter' of the phi field (GeV)
  real, public, save :: lambda  ! -- coefficient at (-phi^4/4) in the Lagrangian
  real, public, save :: epsilon ! -- coefficient at (sigma) in the Lagrangian
  real, public, save :: lambda6 ! -- coefficient at (phi^6/6) in the Lagrangian
  real, public, save :: g4      ! -- coefficient at (omega_mu*omega^mu)^2 in the Lagrangian

  real, public, save :: endens_vac, pressure_vac     ! energy density and pressure in vacuum (GeV/fm^3)


  logical, public, save :: flagWalecka=.true.   ! if false - the parity doublet model is used

  logical, save :: initFlag=.true.


  public :: mDiracNucleon_Approx
  public :: mDirac1535_Approx
  public :: waleckaShift
  public :: PD
  public :: mPD
  public :: dmPD
  public :: d2mPD
  public :: f
  public :: fprime
  public :: fshift
  public :: getRMF_flag
  public :: getRMF_parSet
  public :: ModificationFactor

  public :: PD_rhoB
  public :: PD_2

contains


  !****************************************************************************
  !****f* RMF/ModificationFactor
  ! NAME
  ! real function ModificationFactor(Id,antiFlag)
  ! PURPOSE
  ! Returns the modification factor of the RMF coupling constants for a given
  ! particle.
  ! INPUTS
  ! * integer :: Id         -- Id of particle
  ! * logical :: antiFlag   -- if .true. the particle is an antiparticle
  ! * logical, optional :: vectFlag  -- if .true. the modification for vector couplings
  !****************************************************************************
  pure real function ModificationFactor(Id,antiFlag,vectFlag)
    use IdTable

    integer, intent(in) :: Id
    logical, intent(in) :: antiFlag
    logical, optional, intent(in) :: vectFlag

    if (present(vectFlag)) then
       if (vectFlag .and. .not.flagVectMod) then
          ModificationFactor=1.
          return
       end if
    end if

    if ( nucleon <= Id .and. Id <= F37_1950 .and. antiFlag ) then ! Nonstrange antibaryon
       ModificationFactor=fact_pbar
    else if (Id==delta .and. .not.antiFlag) then ! Delta(1232)
       ModificationFactor=fact_Delta
    else if ( Lambda  <= Id .and. Id <= sigma_1915 ) then  ! Baryon with s=-1
       if (.not.antiFlag) then
          ModificationFactor=fact_hyp
       else
          ModificationFactor=fact_antihyp
       end if
    else if ( Xi  <= Id .and. Id <= XiStar ) then  ! Baryon with s=-2
       if (.not.antiFlag) then
          ModificationFactor=fact_Xi
       else
          ModificationFactor=fact_antiXi
       end if
    else if (id==Kaon .or. id==kaonBar) then
       ModificationFactor=fact_kaon
    else
       ModificationFactor=1.
    end if
  end function ModificationFactor

  !****************************************************************************
  !****f* RMF/getRMF_flag()
  ! NAME
  ! logical function getRMF_flag()
  ! PURPOSE
  ! return the value of the variable RMF_flag. Ensures that jobcard is read.
  !****************************************************************************
  logical function getRMF_flag()
    if (initFlag) call init
    getRMF_flag = RMF_flag
  end function getRMF_flag

  !****************************************************************************
  !****f* RMF/getRMF_parSet()
  ! NAME
  ! integer function getRMF_parSet()
  ! PURPOSE
  ! return the value of the variable N_set. Ensures that jobcard is read.
  !****************************************************************************
  integer function getRMF_parSet()
    if (initFlag) call init
    getRMF_parSet = N_set
  end function getRMF_parSet


  !****************************************************************************
  !****s* RMF/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads input switches. Initializes the mean field parameters.
  !****************************************************************************
  subroutine init

    use constants, only: hbarc,f_pi,mpi
    use output, only: Write_ReadingInput

    use RMFflag, only: setRMF

    integer :: ios
    real :: tmp1

    !**************************************************************************
    !****n* RMF/RMF_input
    ! NAME
    ! NAMELIST /RMF_input/
    ! PURPOSE
    ! Includes the following input switches:
    ! * RMF_flag
    ! * N_set
    ! * grad_flag
    ! * lorentz_flag
    ! * Tens_flag
    ! * flagCorThr
    ! * kaonpot_flag
    ! * fact_pbar
    ! * fact_Delta
    ! * fact_hyp
    ! * fact_antihyp
    ! * fact_Xi
    ! * fact_antiXi
    ! * fact_kaon
    ! * flagVectMod
    !**************************************************************************
    NAMELIST /RMF_input/ RMF_flag, N_set, grad_flag, lorentz_flag, &
         Tens_flag, flagCorThr, kaonpot_flag, &
         fact_pbar, fact_Delta, fact_hyp, fact_antihyp, &
         fact_Xi, fact_antiXi, fact_kaon, flagVectMod

    call Write_ReadingInput('RMF_input',0)
    rewind(5)
    read(5,nml=RMF_input,iostat=ios)
    call Write_ReadingInput('RMF_input',0,ios)
    write(*,*) 'RMF_flag =', RMF_flag,'.'
    if (RMF_flag) then
       write(*,*) '  N_set =', N_set
       write(*,*) '  grad_flag    =', grad_flag
       write(*,*) '  lorentz_flag =', lorentz_flag
       write(*,*) '  kaonpot_flag =', kaonpot_flag
       write(*,*) '  Tens_flag    =', Tens_flag
       write(*,*) '  flagCorThr   =', flagCorThr
       write(*,*) '  fact_pbar    =', fact_pbar
       write(*,*) '  fact_Delta   =', fact_Delta
       write(*,*) '  fact_hyp     =', fact_hyp
       write(*,*) '  fact_antihyp =', fact_antihyp
       write(*,*) '  fact_Xi      =', fact_Xi
       write(*,*) '  fact_antiXi  =', fact_antiXi
       write(*,*) '  fact_kaon    =', fact_kaon
       write(*,*) '  flagVectMod  =', flagVectMod
    end if
    call Write_ReadingInput('RMF_input',1)

    initFlag = .false.

    call setRMF(RMF_flag)

    if (RMF_flag) then

       select case (N_set)

       case (1)
          ! NL1 set from G.A. Lalazissis et al., PRC 55, 540 (1997)
          ! (K=211.29 MeV, m*/m=0.57):
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.492250  ! -- sigma-meson mass (GeV)
          m_omega=0.795359  ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=10.138    ! -- sigma-nucleon coupling constant
          g_omega=13.285    ! -- omega-nucleon coupling constant
          g_rho=4.976       ! -- rho-nucleon coupling constant
          g_2=-12.172*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-36.265       ! -- coefficient at sigma^4 in the Lagrangian

       case (2)
          ! NL3 set from G.A. Lalazissis et al., PRC 55, 540 (1997)
          ! (K=271.76 MeV, m*/m=0.60) :
          m_nucleon=0.939   ! -- nucleon mass (GeV)
          m_sigma=0.508194  ! -- sigma-meson mass (GeV)
          m_omega=0.782501  ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=10.217    ! -- sigma-nucleon coupling constant
          g_omega=12.868    ! -- omega-nucleon coupling constant
          g_rho=4.474       ! -- rho-nucleon coupling constant
          g_2=-10.431*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-28.885       ! -- coefficient at sigma^4 in the Lagrangian

       case (3)
          ! NL2 set from A. Lang et al., NPA 541, 507 (1992)
          ! (K=210 MeV, m*/m=0.83) :
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.5505    ! -- sigma-meson mass (GeV)
          m_omega=0.7833    ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=8.50      ! -- sigma-nucleon coupling constant
          g_omega=7.54      ! -- omega-nucleon coupling constant
          g_rho=0.          ! -- rho-nucleon coupling constant
          g_2=-50.37*hbarc  ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-6.26         ! -- coefficient at sigma^4 in the Lagrangian

       case (4)
          ! NLZ2 set from M. Bender et al., PRC 60, 34304 (1999)
          ! (K=172 MeV, m*/m=0.583) :
          m_nucleon=0.9389  ! -- nucleon mass (GeV)
          m_sigma=0.493150  ! -- sigma-meson mass (GeV)
          m_omega=0.7800    ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=10.1369   ! -- sigma-nucleon coupling constant
          g_omega=12.9084   ! -- omega-nucleon coupling constant
          g_rho=4.55627     ! -- rho-nucleon coupling constant
          g_2=-13.7561*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-41.4013      ! -- coefficient at sigma^4 in the Lagrangian

       case (5)
          ! NL3* set from G.A. Lalazissis, private communication.
          ! (K=258.28 MeV, m*/m=0.594) :
          m_nucleon=0.939   ! -- nucleon mass (GeV)
          m_sigma=0.5026    ! -- sigma-meson mass (GeV)
          m_omega=0.7826    ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=10.0944   ! -- sigma-nucleon coupling constant
          g_omega=12.8065   ! -- omega-nucleon coupling constant
          g_rho=4.5748      ! -- rho-nucleon coupling constant
          g_2=-10.8093*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-30.1486      ! -- coefficient at sigma^4 in the Lagrangian

       case (6)
          ! Same as N_set=3, but including the rho meson.
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.5505    ! -- sigma-meson mass (GeV)
          m_omega=0.7833    ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=8.50      ! -- sigma-nucleon coupling constant
          g_omega=7.54      ! -- omega-nucleon coupling constant
          g_rho=4.271       ! -- rho-nucleon coupling constant
          g_2=-50.37*hbarc  ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-6.26         ! -- coefficient at sigma^4 in the Lagrangian

       case (7)
          ! NL1 set from S.J. Lee et al., PRL 57, 2916 (1986)
          ! (K=212 MeV, m*/m=0.57) :
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.49225   ! -- sigma-meson mass (GeV)
          m_omega=0.795359  ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=10.138    ! -- sigma-nucleon coupling constant
          g_omega=13.285    ! -- omega-nucleon coupling constant
          g_rho=4.976       ! -- rho-nucleon coupling constant
          g_2=-12.172*hbarc ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=-36.265       ! -- coefficient at sigma^4 in the Lagrangian

       case (8)
          ! NL2 set from S.J. Lee et al., PRL 57, 2916 (1986)
          ! (K=399 MeV, m*/m=0.67) :
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.50489   ! -- sigma-meson mass (GeV)
          m_omega=0.78000   ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=9.111     ! -- sigma-nucleon coupling constant
          g_omega=11.493    ! -- omega-nucleon coupling constant
          g_rho=5.507       ! -- rho-nucleon coupling constant
          g_2=-2.304*hbarc  ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=13.783        ! -- coefficient at sigma^4 in the Lagrangian

       case (9)
          ! Set I from B. Liu et al., PRC 65, 045201 (2002)
          ! (K=240 MeV, m*/m=0.75)
          m_nucleon=0.939   ! -- nucleon mass (GeV)
          m_sigma=0.550     ! -- sigma-meson mass (GeV)
          m_omega=0.783     ! -- omega-meson mass (GeV)
          m_rho=0.763       ! -- rho-meson mass (GeV)
          g_sigma=m_sigma*sqrt(10.33)/hbarc   ! -- sigma-nucleon coupling
          g_omega=m_omega*sqrt(5.42)/hbarc    ! -- omega-nucleon coupling
          g_rho=m_rho*sqrt(0.95)/ hbarc       ! -- rho-nucleon coupling
          g_2=-0.033*g_sigma**3*hbarc         ! -- coefficient at sigma^3
          g_3=-0.0048*g_sigma**4              ! -- coefficient at sigma^4

       case (10)
          ! NL1 set from A. Lang et al., NPA 541, 507 (1992)
          ! (K=380 MeV, m*/m=0.83) :
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.5505    ! -- sigma-meson mass (GeV)
          m_omega=0.7833    ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=6.91      ! -- sigma-nucleon coupling constant
          g_omega=7.54      ! -- omega-nucleon coupling constant
          g_rho=0.          ! -- rho-nucleon coupling constant
          g_2=40.6*hbarc    ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=384.4         ! -- coefficient at sigma^4 in the Lagrangian

       case (11)
          ! NL3 set from A. Lang et al., NPA 541, 507 (1992)
          ! (K=380 MeV, m*/m=0.70) :
          m_nucleon=0.938   ! -- nucleon mass (GeV)
          m_sigma=0.5505    ! -- sigma-meson mass (GeV)
          m_omega=0.7833    ! -- omega-meson mass (GeV)
          m_rho=0.763000    ! -- rho-meson mass (GeV)
          g_sigma=9.50      ! -- sigma-nucleon coupling constant
          g_omega=10.95     ! -- omega-nucleon coupling constant
          g_rho=0.          ! -- rho-nucleon coupling constant
          g_2=-1.589*hbarc  ! -- coefficient at sigma^3 in the Lagrangian (GeV)
          g_3=34.23         ! -- coefficient at sigma^4 in the Lagrangian

       case (31)
          ! Parity doublet model Set P3 from D. Zschiesche et al., PRC 75, 055202 (2007)
          ! (K=510 MeV)
          flagWalecka=.false.
          m0=0.790            ! -- chiral invariant mass (GeV)
          m_sigma=0.37063     ! -- sigma-meson mass (GeV)
          m_omega=0.783       ! -- omega-meson mass (GeV)
          g_omega=6.79        ! -- omega-nucleon coupling constant
          g1=13.00            ! -- sigma N1 coupling constant (g1=a)
          g2=6.97             ! -- sigma N2 coupling constant (g2=b)
          lambda6=0. ! -- coefficient at (phi^6/6) in the Lagrangian (GeV^-2)
          g4=0.      ! -- coefficient at (omega_mu*omega^mu)^2 in the Lagrangian
          ! m_nucleon=0.939   ! -- nucleon mass (GeV)
          ! m_minus=1.500     ! -- mass of the partner state with negative parity (GeV)

       case (32)
          ! Parity doublet model Set P2 from D. Zschiesche et al., PRC 75, 055202 (2007)
          ! (K=374 MeV)
          flagWalecka=.false.
          m0=0.790            ! -- chiral invariant mass (GeV)
          m_sigma=0.302       ! -- sigma-meson mass (GeV)
          m_omega=0.783       ! -- omega-meson mass (GeV)
          g_omega=6.77        ! -- omega-nucleon coupling constant
          g1=9.16             ! -- sigma N1 coupling constant (g1=a)
          g2=6.35             ! -- sigma N2 coupling constant (g2=b)
          lambda6=0. ! -- coefficient at (phi^6/6) in the Lagrangian (GeV^-2)
          g4=3.8     ! -- coefficient at (omega_mu*omega^mu)^2 in the Lagrangian
          ! m_nucleon=0.939   ! -- nucleon mass (GeV)
          ! m_minus=1.200     ! -- mass of the partner state with negative parity (GeV)

       case (33)
          ! Parity doublet model Set 1 from I.J. Shin et al., arXiv:1805.03402
          ! (K=240 MeV)
          flagWalecka=.false.
          m0=0.700            ! -- chiral invariant mass (GeV)
          m_sigma=0.385805    ! -- sigma-meson mass (GeV)
          m_omega=0.783       ! -- omega-meson mass (GeV)
          m_rho=0.776         ! -- rho-meson mass (GeV)
          g_omega=7.30465     ! -- omega-nucleon coupling constant
          g_rho=4.06502       ! -- rho-nucleon coupling constant
          g1=14.1708          ! -- sigma N1 coupling constant (g1=a)
          g2=7.76222          ! -- sigma N2 coupling constant (g2=b)
          lambda6=13.5401/f_pi**2  ! -- coefficient at (phi^6/6)  (GeV^-2)
          g4=0.     ! -- coefficient at (omega_mu*omega^mu)^2 in the Lagrangian
          ! m_nucleon=0.939   ! -- nucleon mass (GeV)
          ! m_minus=1.535     ! -- mass of the partner state with negative parity (GeV)

       case (34)
          ! Parity doublet model Set 2 from I.J. Shin et al., arXiv:1805.03402
          ! (K=215 MeV)
          flagWalecka=.false.
          m0=0.700            ! -- chiral invariant mass (GeV)
          m_sigma=0.384428    ! -- sigma-meson mass (GeV)
          m_omega=0.783       ! -- omega-meson mass (GeV)
          m_rho=0.776         ! -- rho-meson mass (GeV)
          g_omega=7.05508     ! -- omega-nucleon coupling constant
          g_rho=4.07986       ! -- rho-nucleon coupling constant
          g1=14.1708          ! -- sigma N1 coupling constant (g1=a)
          g2=7.76222          ! -- sigma N2 coupling constant (g2=b)
          lambda6=15.7393/f_pi**2  ! -- coefficient at (phi^6/6)  (GeV^-2)
          g4=0.     ! -- coefficient at (omega_mu*omega^mu)^2 in the Lagrangian
          ! m_nucleon=0.939   ! -- nucleon mass (GeV)
          ! m_minus=1.535     ! -- mass of the partner state with negative parity (GeV)

       case (35)
          ! Parity doublet model Set P4 from D. Zschiesche et al., PRC 75, 055202 (2007)
          ! (K=440 MeV)
          flagWalecka=.false.
          m0=0.790            ! -- chiral invariant mass (GeV)
          m_sigma=0.34659     ! -- sigma-meson mass (GeV)
          m_omega=0.783       ! -- omega-meson mass (GeV)
          g_omega=7.75        ! -- omega-nucleon coupling constant
          g1=13.00            ! -- sigma N1 coupling constant (g1=a)
          g2=6.97             ! -- sigma N2 coupling constant (g2=b)
          lambda6=0.          ! -- coefficient at (phi^6/6)  (GeV^-2)
          g4=3.8    ! -- coefficient at (omega_mu*omega^mu)^2 in the Lagrangian
          ! m_nucleon=0.939   ! -- nucleon mass (GeV)
          ! m_minus=1.500     ! -- mass of the partner state with negative parity (GeV)

       case default
          write(*,*) 'invalid value for N_set: ', N_set
          call TRACEBACK()

       end select


       if (flagWalecka) then
          a_0=0.
          a_1=(g_sigma/m_sigma)**2*hbarc**3
          a_2=g_2*g_sigma/m_sigma**2
          a_3=g_3*g_sigma/m_sigma**2
          a_4=2.*g_2/m_sigma**2
          a_5=3.*g_3/m_sigma**2
       else
          mu=sqrt((m_sigma**2-3.*mpi**2)/2.+lambda6*f_pi**4)         ! -- 'mass parameter' of the phi field (GeV)
          lambda=(m_sigma**2-mpi**2)/2./f_pi**2+2.*lambda6*f_pi**2   ! -- coefficient at (-phi^4/4.) in the Lagrangian
          epsilon=mpi**2*f_pi                                        ! -- coefficient at (sigma) in the Lagrangian (GeV^3)
          tmp1=sqrt((f_pi*(g1+g2))**2+4.*m0**2)
          m_nucleon=0.5*(tmp1-f_pi*(g1-g2))
          m_minus=0.5*(tmp1+f_pi*(g1-g2))
          endens_vac = (-mu**2/2.*f_pi**2+lambda/4.*f_pi**4-lambda6/6.*f_pi**6-epsilon*f_pi)/hbarc**3
          pressure_vac=-endens_vac
          g_sigma=0.5*f_pi*(g1+g2)**2/tmp1          ! g_sigma = (dmPlus/dsigma + dmMinus/dsigma)/2
          a_0=epsilon/mu**2*g_sigma
          a_1=(g_sigma/mu)**2*hbarc**3
          a_2=0.
          a_3=lambda/mu**2*g_sigma
          a_4=0.
          a_5=lambda6/mu**2*g_sigma

          write(*,*) 'Parity doublet model is used with parameters: '
          write(*,*) '  m0, GeV :      ', m0
          write(*,*) '  m_sigma, GeV : ', m_sigma
          write(*,*) '  m_omega, GeV : ', m_omega
          write(*,*) '  g_omega : ', g_omega
          write(*,*) '  g1 : ', g1
          write(*,*) '  g2 : ', g2
          write(*,*) '  mu, GeV : ', mu
          write(*,*) '  lambda : ', lambda
          write(*,*) '  epsilon, GeV^3 : ', epsilon
          write(*,*) '  lambda6, GeV^-2 : ', lambda6
          write(*,*) '  g4 : ', g4
          write(*,*) '  effective coupling constant g_sigma : ',  g_sigma
          write(*,*) '  m_nucleon, m_minus, GeV : ', m_nucleon, m_minus
          write(*,*) '  vacuum energy density and pressure,  GeV/fm^3 : ', &
               endens_vac, pressure_vac
       end if

       a_6=(g_omega/m_omega)**2*hbarc**3
       a_7=(g_rho/m_rho)**2*hbarc**3

    end if

  end subroutine init


  !****************************************************************************
  !****f* RMF/waleckaShift
  ! NAME
  ! function waleckaShift(rhobar,em0,rhoscalar,endens,S,V,potential)
  ! return(shift)
  ! PURPOSE
  ! Determine the mass shift of the nucleon in equilibrated
  ! isospin symmetric nuclear matter at zero temperature
  ! within Walecka model with nonlinear sigma-coupling.
  ! INPUTS
  ! * real :: rhobar          ! -- baryon density (fm^-3)
  ! * real, optional :: em0   ! -- starting value of mass for iterations (GeV)
  ! OUTPUT
  ! * real  :: shift              ! = m - m^* -- mass shift (GeV)
  ! * real, optional :: rhoscalar ! -- scalar density (fm^-3)
  ! * real, optional :: endens    ! -- energy density (GeV/fm^3)
  ! * real, optional :: pressure  ! -- pressure (GeV/fm^3)
  ! * real, optional :: S         ! -- scalar potential (GeV)
  ! * real, optional :: V         ! -- vector potential (GeV)
  ! * real, optional :: potential ! -- Schroedinger equivalent potential (GeV)
  !
  ! NOTES:
  ! * The potential is for a particle at rest, so there is no gamma factor
  !   E/m_N in front of U_v for the SE potential
  !****************************************************************************
  function waleckaShift(rhobar,em0,rhoscalar,endens,pressure,S,V,potential) result(shift)

    use constants, only: hbarc, pi, mN

    real :: shift

    real, intent(in) :: rhobar
    real, optional, intent(in) :: em0
    real, optional, intent(out) :: rhoscalar
    real, optional, intent(out) :: endens
    real, optional, intent(out) :: pressure
    real, optional, intent(out) :: S
    real, optional, intent(out) :: V
    real, optional, intent(out) :: potential

    real :: pf, dmstm, sigma, a, rhos, drhos, fun, dfun, U_s, U_v
    integer :: niter
    logical, parameter :: debug=.false.

    if (initFlag) call init

    if (rhobar < 0.001) then
       shift = 0.
       if (present(rhoscalar)) rhoscalar = 0.
       if (present(endens)) endens = 0.
       if (present(pressure)) pressure = 0.
       if (present(S)) S=0.
       if (present(V)) V=0.
       if (present(potential)) potential=0.
       return
    end if

    pf = (1.5*pi**2*rhobar)**0.333333*hbarc

    ! Here we reset m_nucleon using the default value for the
    ! nucleon mass:
    m_nucleon = mN

    if (present(em0)) then
       dmstm = em0 - m_nucleon
    else
       dmstm = -fshift(rhobar)
    end if

    ! Test: *****************
    !shift = fshift(rhobar)
    !return
    !************************

    niter = 0
    do
       niter = niter + 1

       sigma = dmstm/g_sigma    ! mean value of the sigma-meson field
       a = (dmstm + m_nucleon)/pf
       rhos = rhobar*f(a)        ! scalar density
       drhos = rhobar/pf*fprime(a) ! d rhos / d mst
       fun = dmstm + a_1*rhos + a_2*sigma**2 + a_3*sigma**3   ! we want to have fun = 0.
       dfun = 1. + a_1*drhos + a_4*sigma + a_5*sigma**2  ! d fun / d mst

       if (dfun.ne.0.) dmstm = dmstm - fun/dfun

       if ( abs(fun) <= 1.e-04 ) exit ! ==> success

       if ( niter == 10 ) then
          write(*,*) 'bad convergence after ',niter,' iterations:',&
               & rhobar,dmstm,abs(fun)
          call TRACEBACK()
       end if

       if (debug) then
          write(*,'(a26,1x,f4.2,1x,i2,1x,e13.6,1x,e13.6)') &
               'rhobar, niter, dmstm, fun:', rhobar, niter, dmstm, fun
       end if

    end do

    shift = -dmstm
    a = (dmstm + m_nucleon)/pf

    if (present(rhoscalar)) rhoscalar=rhobar*f(a)

    if (present(endens)) then
       ! Compute also the energy density and pressure:
       sigma = dmstm/g_sigma
       endens = ( 2.*pf**4/pi**2*g(a) &
            &+   0.5*(m_sigma*sigma)**2 + g_2*sigma**3/3. &
            &+   g_3*sigma**4/4. )/hbarc**3 &
            &+ 0.5*a_6*rhobar**2
       pressure = rhobar*sqrt(pf**2+(dmstm + m_nucleon)**2) &
            &-( 2.*pf**4/pi**2*g(a) &
            &+  0.5*(m_sigma*sigma)**2 + g_2*sigma**3/3. &
            &+  g_3*sigma**4/4. )/hbarc**3 &
            &+ 0.5*a_6*rhobar**2
    end if

    U_s=g_sigma*sigma
    U_v=a_6*rhobar

    if (present(S)) S=U_s
    if (present(V)) V=U_v
    if (present(potential)) potential = U_s + U_v + (U_s**2-U_v**2)/(2.*m_nucleon)

  end function waleckaShift


  !****************************************************************************
  !****f* RMF/fshift
  ! NAME
  ! real function fshift(rho)
  ! PURPOSE
  ! Fit of the nucleon mass shift m - m* for the various RMF parameter sets.
  ! INPUTS
  ! * real :: rho  -- baryon density (fm**-3)
  ! OUTPUT
  ! * real :: fshift  -- m - m* (GeV)
  ! NOTES
  ! This is a very rough fit which is only good to provide
  ! the starting value for iterations in walecka.
  ! The density rho must be in the interval from 0 up to 12*rhoNull.
  !****************************************************************************
  real function fshift(rho)

    use constants, only: rhoNull

    real, intent(in) :: rho

    ! constants for each parameter set (index = N_set)
    real, parameter :: A(1:11) = (/ &
         1.03099, 0.783659, 0.201221, 1.01769, 0.75293, &
         0.201221, 1.00994, 0.599025, 0.359, 0.22723, &
         0.539378 /)
    real, parameter :: B(1:11) = (/ &
         1.52439, 1.41926,  0.888537, 1.57845, 1.63993, &
         0.888537, 1.67593, 1.16527,  1.105, 0.536761, &
         0.987796 /)

    if (initFlag) call init

    if (N_set.gt.11) &
         call TRACEBACK('no fitted m*(rho) for this rmf set')

    fshift = m_nucleon * ( 1.-1./(1.+A(N_set)*(rho/rhoNull)**B(N_set)) )

  end function fshift



  !****************************************************************************
  !****s* RMF/PD
  ! NAME
  ! subroutine PD(mubStar,sigma,shift,flagPlot,sigma_inp,mub,rhoPlus,rhoMinus,rhoscalar,endens,pressure,S,V,potential)
  ! PURPOSE
  ! Determine the scalar field and mass shifts of the nucleon and its negative
  ! parity partner in equilibrated isospin symmetric nuclear matter at zero
  ! temperature within parity doublet model.
  !
  ! INPUTS
  ! * real :: mubStar               ! = sqrt(pf_pm**2+m_pm**2)
  !   -- kinetic part of baryon chemical potential (GeV)
  ! * logical, optional :: flagPlot ! if true: equation for sigma field
  !   f(sigma)=0 is not solved, only function f vs sigma is plotted
  ! * real, optional :: sigma_inp   ! -- starting value of sigma field (GeV)
  !
  ! OUTPUT
  ! * real :: sigma               ! -- scalar field (GeV)
  ! * real :: shift(1:2)          ! = m - m^* -- mass shift (GeV),
  !   1 - nucleon, 2 - negative parity partner
  ! * real, optional :: mub       ! -- baryon chemical potential (GeV)
  ! * real, optional :: rhoPlus   ! -- nucleon density (fm^-3)
  ! * real, optional :: rhoMinus  ! -- partner density (fm^-3)
  ! * real, optional :: rhoscalar ! -- scalar density (fm^-3)
  ! * real, optional :: endens    ! -- energy density (GeV/fm^3)
  ! * real, optional :: pressure  ! -- pressure (GeV/fm^3)
  ! * real, optional :: S(1:2)         ! -- scalar potential (GeV)
  ! * real, optional :: V(1:2)         ! -- vector potential (GeV)
  ! * real, optional :: potential(1:2) ! -- Schroedinger equivalent pot. (GeV)
  !
  ! NOTES:
  ! * The potential is for a particle at rest, so there is no gamma factor
  !   E/m in front of U_v for the SE potential
  !****************************************************************************
  subroutine PD(mubStar,sigma,shift,flagPlot,sigma_inp,mub,rhoPlus,rhoMinus,rhoscalar,endens,pressure,S,V,potential)

    use constants, only: hbarc, pi, f_pi

    real, intent(in) :: mubStar
    real, intent(out) :: sigma
    real, intent(out) :: shift(1:2)
    logical, optional, intent(in) :: flagPlot
    real, optional, intent(in) :: sigma_inp
    real, optional, intent(out) :: mub
    real, optional, intent(out) :: rhoPlus
    real, optional, intent(out) :: rhoMinus
    real, optional, intent(out) :: rhoscalar
    real, optional, intent(out) :: endens
    real, optional, intent(out) :: pressure
    real, optional, intent(out) :: S(1:2)
    real, optional, intent(out) :: V(1:2)
    real, optional, intent(out) :: potential(1:2)

    real :: pfPlus,pfMinus,rhob,mPlus,mMinus,rhosPlus,rhosMinus,dmPlus,dmMinus,fun,drhosPlus,drhosMinus,d2m,dfun,U_s(1:2),U_v(1:2),omega
    integer :: niter
    logical, parameter :: debug=.true.
    logical :: flagPlt

    if (initFlag) call init

!!$    if (abs(mubStar-m_nucleon) < 1.e-05) then ! ?????
!!$       sigma=f_pi
!!$       shift = 0.
!!$       if (present(mub)) mub = mubStar
!!$       if (present(rhoPlus)) rhoPlus = 0.
!!$       if (present(rhoMinus)) rhoMinus = 0.
!!$       if (present(rhoscalar)) rhoscalar = 0.
!!$       if (present(endens)) endens = 0.
!!$       if (present(pressure)) pressure = 0.
!!$       if (present(S)) S=0.
!!$       if (present(V)) V=0.
!!$       if (present(potential)) potential=0.
!!$       return
!!$    end if

    if (present(sigma_inp)) then
       sigma=sigma_inp
    else
       sigma=f_pi     ! vacuum value
    end if

    flagPlt=.false.
    if (present(flagPlot)) flagPlt = flagPlot

    if (flagPlt) open(2,file='test.dat',status='unknown')

    niter = 0
    do
       niter = niter + 1

       if (flagPlt) sigma = sigma - f_pi/1000.

       mPlus=mPD(sigma,1)     ! nucleon
       mMinus=mPD(sigma,-1)   ! partner

       dmPlus = dmPD(sigma,1)  ! d mPlus / d sigma
       dmMinus = dmPD(sigma,-1)  ! d mMinus / d sigma

       d2m = d2mPD(sigma)   ! d^2 m / d sigma^2

       pfPlus=sqrt(max(0.,mubStar**2-mPlus**2))
       pfMinus=sqrt(max(0.,mubStar**2-mMinus**2))

       rhosPlus=0.
       drhosPlus=0.
       rhosMinus=0.
       drhosMinus=0.

       if (pfPlus.gt.0.) then
          rhosPlus = 2./(3.*pi**2)*pfPlus**3*f(mPlus/pfPlus)      ! nucleon scalar density
          drhosPlus = -2.*dmPlus/pi**2*( mPlus*pfPlus*f(mPlus/pfPlus) - mubStar**2/3.*fprime(mPlus/pfPlus) ) ! d rhosPlus / d sigma
       end if

       if (pfMinus.gt.0.) then
          rhosMinus = 2./(3.*pi**2)*pfMinus**3*f(mMinus/pfMinus)  ! partner scalar density
          drhosMinus = -2.*dmMinus/pi**2*( mMinus*pfMinus*f(mMinus/pfMinus) - mubStar**2/3.*fprime(mMinus/pfMinus) ) ! d rhosMinus / d sigma
       end if

       fun = dmPlus*rhosPlus + dmMinus*rhosMinus  - mu**2*sigma + lambda*sigma**3 - lambda6*sigma**5 - epsilon    ! we want to have fun = 0.

       dfun = d2m*rhosPlus + dmPlus*drhosPlus + d2m*rhosMinus + dmMinus*drhosMinus - mu**2 + 3.*lambda*sigma**2 - 5.*lambda6*sigma**4     ! d fun / d sigma

       if (debug) then
          write(*,'(a26,1x,e13.6,1x,i3,1x,e13.6,1x,e13.6)') &
               'mubStar, niter, sigma, fun:', mubStar, niter, sigma, fun
       end if


       if (flagPlt) then
          write(2,*) sigma,fun,dfun,rhosPlus,drhosPlus,mPlus,dmPlus,d2m
          if (niter.eq.1000) then
             exit
          else
             cycle
          end if
       end if

       if ( abs(fun) <= 1.e-06 ) exit

       if ( niter == 11 ) then
          write(*,*) 'bad convergence after', niter-1,' iterations:',&
               & mubStar,sigma,abs(fun)
          call Traceback()
       end if

       if (dfun.ne.0.) sigma = sigma - fun/dfun

    end do

    shift(1) = m_nucleon - mPlus
    shift(2) = m_minus - mMinus

    rhob=2./(3.*pi**2)*(pfPlus**3+pfMinus**3)   ! baryon density (GeV**3)

    if (present(rhoPlus)) rhoPlus =2./(3.*pi**2)*pfPlus**3/hbarc**3     ! fm^-3
    if (present(rhoMinus)) rhoMinus =2./(3.*pi**2)*pfMinus**3/hbarc**3  ! fm^-3

    if (present(rhoscalar)) rhoscalar=(rhosPlus+rhosMinus)/hbarc**3     ! fm^-3

    omega=g_omega/m_omega**2*rhob   ! GeV

    if (abs(g4).gt.1.e-03) then

       niter = 0
       do
          niter = niter + 1

          fun = 4.*g4**4*omega**3 + m_omega**2*omega - g_omega*rhob

          if (debug) then
             write(*,'(a26,1x,e13.6,1x,i3,1x,e13.6,1x,e13.6)') &
                  'rhob, niter, omega, fun:', rhob, niter, omega, fun
          end if

          if ( abs(fun) <= 1.e-06 ) exit

          if ( niter == 11 ) then
             write(*,*) 'bad convergence after', niter-1,' iterations:',&
                  & rhob,omega,abs(fun)
             call Traceback()
          end if

          dfun = 12.*g4**4*omega**2 + m_omega**2   ! d fun/ d omega

          if (dfun.ne.0.) omega = omega - fun/dfun

       end do

    end if


    if (present(endens)) then
       endens = 0.
       if (pfPlus.gt.0.) endens = endens + 2.*pfPlus**4/pi**2*g(mPlus/pfPlus)
       if (pfMinus.gt.0.) endens = endens + 2.*pfMinus**4/pi**2*g(mMinus/pfMinus)
       endens = ( endens  &
            &   - mu**2/2.*sigma**2 + lambda/4.*sigma**4 - lambda6/6.*sigma**6 - epsilon*sigma &
            &   + g_omega*omega*rhob - 0.5*m_omega**2*omega**2 - g4**4*omega**4)/hbarc**3 - endens_vac
    end if

    if (present(pressure)) then
       pressure = 0.
       if (pfPlus.gt.0.) pressure = pressure - 2.*pfPlus**4/pi**2*g(mPlus/pfPlus)
       if (pfMinus.gt.0.) pressure = pressure - 2.*pfMinus**4/pi**2*g(mMinus/pfMinus)
       pressure = ( pressure  + rhob*mubStar &
            &   + mu**2/2.*sigma**2 - lambda/4.*sigma**4 + lambda6/6.*sigma**6 + epsilon*sigma &
            &   + 0.5*m_omega**2*omega**2 + g4**4*omega**4)/hbarc**3 - pressure_vac
    end if

    U_s(1:2) = -shift(1:2)
    U_v(1:2) = g_omega*omega

    if (present(mub)) mub = mubStar + U_v(1)

    if (present(S)) S=U_s
    if (present(V)) V=U_v
    if (present(potential)) then
       potential(1) = U_s(1) + U_v(1) + (U_s(1)**2-U_v(1)**2)/(2.*m_nucleon)
       potential(2) = U_s(2) + U_v(2) + (U_s(2)**2-U_v(2)**2)/(2.*m_minus)
    end if

  end subroutine PD


  !****************************************************************************
  !****f* RMF/mPD
  ! NAME
  ! real function mPD
  ! PURPOSE
  ! Computes effective mass (in GeV) in the parity doublet model
  !
  ! INPUTS
  ! * real :: sigma     -- scalar field (GeV)
  ! * integer :: parity -- +1 for nucleon, -1 for S11_1535
  !****************************************************************************
  real function mPD(sigma,parity)

    real, intent(in) :: sigma
    integer, intent(in) :: parity

    if (abs(parity).ne.1) then
       call TRACEBACK('wrong parity')
    end if

    mPD=0.5*(sqrt((sigma*(g1+g2))**2+4.*m0**2) - parity*sigma*(g1-g2))

  end function mPD


  !****************************************************************************
  !****f* RMF/dmPD
  ! NAME
  ! real function dmPD
  ! PURPOSE
  ! Computes derivative of effective mass over sigma field d m/d sigma in the
  ! parity doublet model
  !
  ! INPUTS
  ! * real :: sigma     -- scalar field (GeV)
  ! * integer :: parity -- +1 for nucleon, -1 for S11_1535
  !****************************************************************************
  real function dmPD(sigma,parity)

    real, intent(in) :: sigma
    integer, intent(in) :: parity

    if (abs(parity).ne.1) then
       call TRACEBACK('wrong parity')
    end if

    dmPD = 0.5*(sigma*(g1+g2)**2/sqrt((sigma*(g1+g2))**2+4.*m0**2)-parity*(g1-g2))

  end function dmPD


  !****************************************************************************
  !****f* RMF/d2mPD
  ! NAME
  ! real function d2mPD
  ! PURPOSE
  ! Computes second derivative of effective mass over sigma field
  ! d^2 m/d sigma^2 (in GeV^-1) in the parity doublet model
  !
  ! INPUTS
  ! * real :: sigma  -- scalar field (GeV)
  !****************************************************************************
  real function d2mPD(sigma)

    real, intent(in) :: sigma

    d2mPD=2.*(m0*(g1+g2))**2/(sqrt((sigma*(g1+g2))**2+4.*m0**2))**3

  end function d2mPD


  !****************************************************************************
  !****if* RMF/f
  ! NAME
  ! real function f(a)
  ! PURPOSE
  ! Computes analytically the expression
  ! 3*a*\int_0^1 dx x^2/\sqrt(x^2+a^2)
  ! INPUTS
  ! * real, intent(in) :: a  -- dimensionless parameter equal to m^*/p_F
  !****************************************************************************
  real function f(a)

    real, intent(in) :: a
    real :: tmp

    tmp = sqrt(1.+a**2)
    f = 1.5*a * ( tmp - 0.5*a**2*log((tmp + 1.)/(tmp - 1.)) )

  end function f


  !****************************************************************************
  !****if* RMF/fprime
  ! NAME
  ! real function fprime(a)
  ! PURPOSE
  ! Computes analytically the derivative
  ! of function f(a) with respect to a.
  ! INPUTS
  ! * real :: a  -- dimensionless parameter equal to m^*/p_F
  !****************************************************************************
  real function fprime(a)

    real, intent(in) :: a
    real :: tmp

    tmp = sqrt(1.+a**2)
    fprime = f(a)/a + 3.*a**2*( 1./tmp - log((tmp+1.)/a) )

  end function fprime


  !****************************************************************************
  !****if* RMF/g
  ! NAME
  ! real function g(a)
  ! PURPOSE
  ! Computes analytically the expression
  ! \int_0^1 dx x^2*\sqrt(x^2+a^2)
  ! INPUTS
  ! * real :: a  -- dimensionless parameter equal to m^*/p_F
  !****************************************************************************
  real function g(a)

    real, intent(in) :: a
    real :: tmp

    tmp = sqrt(1.+a**2)
    g = ( tmp**3 + tmp - 0.5*a**4*log((tmp + 1.)/(tmp - 1.)) )/8.

  end function g

  !****************************************************************************
  !****f* RMF/mDiracNucleon_Approx
  ! NAME
  ! real function mDiracNucleon_Approx(rhoBar)
  ! PURPOSE
  ! shortcut to calculate the Dirac mass in Walecka and Parity Doublet Model
  ! for the nucleon mass used in 'initNucPhaseSpace'
  ! INPUTS
  ! * real :: rhobar          -- baryon density (fm^-3)
  ! OUTPUT
  ! * the Dirac mass (in GeV)
  !
  ! NOTES
  ! * a shortcut to
  !      mDirac = mN - waleckaShift(...)
  !   resp.
  !      mDirac = mPD(sigma,+1)
  ! * For PDM, the sigma field is only approximated
  !****************************************************************************
  function mDiracNucleon_Approx(rhoBar) result(mDirac)
    use constants, only: mN, hbarc, mpi, f_pi

    real :: mDirac
    real, intent(in) :: rhobar

    real :: sigma

    if (flagWalecka) then
       mDirac = mN - waleckaShift(rhoBar)
    else ! Parity doublet model
       sigma = f_pi*(1.-Sigma_N/(mpi*f_pi)**2*rhoBar*hbarc**3)
       mDirac = mPD(sigma,+1)
    end if

  end function mDiracNucleon_Approx

  !****************************************************************************
  !****f* RMF/mDirac1535_Approx
  ! NAME
  ! real function mDirac1535_Approx(rhoBar)
  ! PURPOSE
  ! shortcut to calculate the Dirac mass in Walecka and Parity Doublet Model
  ! for the mass of the N*(1535) resonance
  ! INPUTS
  ! * real :: rhobar          -- baryon density (fm^-3)
  ! OUTPUT
  ! * the Dirac mass (in GeV)
  ! NOTES
  ! * For PDM, the sigma field is only approximated
  !****************************************************************************
  function mDirac1535_Approx(rhoBar) result(mDirac)
    use constants, only: hbarc, mpi, f_pi
    use particleProperties, only: hadron
    use IdTable, only: S11_1535

    real :: mDirac
    real, intent(in) :: rhobar

    real :: sigma

    if (flagWalecka) then
       mDirac = hadron(S11_1535)%mass - waleckaShift(rhoBar)
    else ! Parity doublet model
       sigma = f_pi*(1.-Sigma_N/(mpi*f_pi)**2*rhoBar*hbarc**3)
       mDirac = mPD(sigma,-1)
    end if

  end function mDirac1535_Approx


  subroutine PD_rhoB(rhoB_in)
    use constants, only: hbarc, pi, mpi, f_pi

    real, intent(in) :: rhoB_in

    real :: rhoB, sigma, omega, mubStar, sigmaApprox
    real :: mPlus, mMinus, pFplus, pFminus
    real :: dmPlus, dmMinus, d2m
    real :: rhosPlus, rhosMinus, drhosPlus, drhosMinus
    real :: fun, dfun
    real :: endens

    integer :: niter




    rhoB = rhoB_in*hbarc**3
    sigmaApprox = f_pi*(1.-Sigma_N/(mpi*f_pi)**2*rhoB)

    !    sigma = f_pi
    sigma = sigmaApprox

    pFplus = (rhoB * 3*pi**2/2)**(1./3.)
    pFminus = 0.

    niter = 0
    do
       exit ! do not iterate at all

       niter = niter + 1

       mPlus=mPD(sigma,1)     ! nucleon
       mMinus=mPD(sigma,-1)   ! partner

       dmPlus = dmPD(sigma,1)  ! d mPlus / d sigma
       dmMinus = dmPD(sigma,-1)  ! d mMinus / d sigma

       d2m = d2mPD(sigma)   ! d^2 m / d sigma^2

       mubStar = sqrt(pFplus**2+mPlus**2)

       rhosPlus=0.
       drhosPlus=0.
       rhosMinus=0.
       drhosMinus=0.

       if (pfPlus.gt.0.) then
          rhosPlus = 2./(3.*pi**2)*pfPlus**3*f(mPlus/pfPlus)      ! nucleon scalar density
          drhosPlus = -2.*dmPlus/pi**2*( mPlus*pfPlus*f(mPlus/pfPlus) - mubStar**2/3.*fprime(mPlus/pfPlus) ) ! d rhosPlus / d sigma
       end if

       fun = dmPlus*rhosPlus + dmMinus*rhosMinus  - mu**2*sigma + lambda*sigma**3 - lambda6*sigma**5 - epsilon    ! we want to have fun = 0.

       dfun = d2m*rhosPlus + dmPlus*drhosPlus + d2m*rhosMinus + dmMinus*drhosMinus - mu**2 + 3.*lambda*sigma**2 - 5.*lambda6*sigma**4     ! d fun / d sigma

       write(*,*) niter,sigmaApprox,sigma,fun,dfun

       if ( abs(fun) <= 1.e-06 ) exit ! ==> success

       if ( niter == 11 ) then
          write(*,*) 'bad convergence after', niter-1,' iterations:',&
               & rhoB,sigma,abs(fun)
          call Traceback()
       end if

       if (dfun.ne.0.) sigma = sigma - fun/dfun

    end do

    omega=g_omega/m_omega**2*rhob   ! GeV

    mPlus=mPD(sigma,1)     ! nucleon
    mMinus=mPD(sigma,-1)   ! partner
    mubStar = sqrt(pFplus**2+mPlus**2)

    endens = 0.
    if (pfPlus.gt.0.) endens = endens + 2.*pfPlus**4/pi**2*g(mPlus/pfPlus)
    if (pfMinus.gt.0.) endens = endens + 2.*pfMinus**4/pi**2*g(mMinus/pfMinus)
    endens = ( endens  &
         &   - mu**2/2.*sigma**2 + lambda/4.*sigma**4 - lambda6/6.*sigma**6 - epsilon*sigma &
         &   + g_omega*omega*rhob - 0.5*m_omega**2*omega**2 - g4**4*omega**4)/hbarc**3

    write(*,*) rhoB/hbarc**3,niter, sigma, mPlus,mMinus, endens,endens_vac, mubStar
    write(679,*) rhoB/hbarc**3,niter, sigma, mPlus,mMinus, endens,endens_vac, mubStar

  end subroutine PD_rhoB

  subroutine PD_2(mubStar)

    use constants, only: hbarc, pi, f_pi

    real, intent(in) :: mubStar

    real :: rhoB, sigma, omega
    real :: mPlus, mMinus, pFplus, pFminus
    real :: dmPlus, dmMinus
    real :: rhosPlus, rhosMinus
    real :: fun
    real :: fun0
    real :: endens

    integer :: i, iS
    integer :: nStep


    if (initFlag) call init

!    sigma=f_pi     ! vacuum value
    sigma=f_pi +0.010    ! vacuum value + dummy

    if (mubStar < 0.091) then
       nStep = 1000
    else
       nStep = 10000
    end if

    iS = 0

    do i=0,nStep

       sigma = sigma - f_pi/nStep

       mPlus=mPD(sigma,1)     ! nucleon
       mMinus=mPD(sigma,-1)   ! partner

       dmPlus = dmPD(sigma,1)  ! d mPlus / d sigma
       dmMinus = dmPD(sigma,-1)  ! d mMinus / d sigma
!!$
!!$       d2m = d2mPD(sigma)   ! d^2 m / d sigma^2

       pfPlus=sqrt(max(0.,mubStar**2-mPlus**2))
       pfMinus=sqrt(max(0.,mubStar**2-mMinus**2))

       rhosPlus=0.
!!$       drhosPlus=0.
       rhosMinus=0.
!!$       drhosMinus=0.

       if (pfPlus.gt.0.) then
          rhosPlus = 2./(3.*pi**2)*pfPlus**3*f(mPlus/pfPlus)      ! nucleon scalar density
!!$          drhosPlus = -2.*dmPlus/pi**2*( mPlus*pfPlus*f(mPlus/pfPlus) - mubStar**2/3.*fprime(mPlus/pfPlus) ) ! d rhosPlus / d sigma
       end if

       if (pfMinus.gt.0.) then
          rhosMinus = 2./(3.*pi**2)*pfMinus**3*f(mMinus/pfMinus)  ! partner scalar density
!!$          drhosMinus = -2.*dmMinus/pi**2*( mMinus*pfMinus*f(mMinus/pfMinus) - mubStar**2/3.*fprime(mMinus/pfMinus) ) ! d rhosMinus / d sigma
       end if

       fun = dmPlus*rhosPlus + dmMinus*rhosMinus  - mu**2*sigma + lambda*sigma**3 - lambda6*sigma**5 - epsilon    ! we want to have fun = 0.

!!$       dfun = d2m*rhosPlus + dmPlus*drhosPlus + d2m*rhosMinus + dmMinus*drhosMinus - mu**2 + 3.*lambda*sigma**2 - 5.*lambda6*sigma**4     ! d fun / d sigma

       if (i.eq.0) then
          fun0 = fun
!!$          dfun0 = dfun
          cycle
       end if

       write(717,*) i,sigma,fun

       if (fun*fun0 < 0) then
          ! we bracket the zero

          write(*,*) mubStar,sigma,fun,mPlus,mMinus
          write(718,*) mubStar,sigma,fun,mPlus,mMinus

          iS = iS+1

          if (iS==2) exit ! this is the solution we want

       end if

       fun0 = fun

    end do

    write(717,*) " "
    write(717,*) " "
    write(718,*)

    rhob=2./(3.*pi**2)*(pfPlus**3+pfMinus**3)   ! baryon density (GeV**3)
    omega=g_omega/m_omega**2*rhob   ! GeV

        endens = 0.
    if (pfPlus.gt.0.) endens = endens + 2.*pfPlus**4/pi**2*g(mPlus/pfPlus)
    if (pfMinus.gt.0.) endens = endens + 2.*pfMinus**4/pi**2*g(mMinus/pfMinus)
    endens = ( endens  &
         &   - mu**2/2.*sigma**2 + lambda/4.*sigma**4 - lambda6/6.*sigma**6 - epsilon*sigma &
         &   + g_omega*omega*rhob - 0.5*m_omega**2*omega**2 - g4**4*omega**4)/hbarc**3

    write(*,*) rhoB/hbarc**3, 0, sigma, mPlus,mMinus, endens,endens_vac, mubStar
    write(719,*) rhoB/hbarc**3, 0, sigma, mPlus,mMinus, endens,endens_vac, mubStar,mubStar+g_omega*omega

  end subroutine PD_2

end module RMF
