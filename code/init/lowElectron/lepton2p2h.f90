!******************************************************************************
!****m* /lepton2p2h
! NAME
! module lepton2p2h
!
! PURPOSE
! Do all the internals for 2p2h scattering:
!
! EM:
! * ell N1 N2 --> ell' N1' N2'  == gamma* N1 N2 --> N1' N2'
! * ell N1 N2 --> ell' N Delta  == gamma* N1 N2 --> N Delta
! NC:
! * nu  N1 N2 --> nu'  N1' N2'
! * nu  N1 N2 --> nu'  N Delta
! CC:
! * nu  N1 N2 --> ell- N1' N2' (sum of hadronic charges increases by +1)
! * nu  N1 N2 --> ell- N Delta (  -- " --                              )
!
! antiEM, antiNC and antiCC are the same as EM, NC, CC.
!
! cases 1 - 3 give parametrizations for 2p2h part of structure function W1
! in terms of Q^2, no distinction for neutrinos and antineutrinos
!
! Cases 1 and 2 are those in: Lalakulich Gallmeister Mosel PRC86(2012)014614

! Case 3 gives a reasonable description of MiniBooNE dd neutrino data
!
! cases 4 - 5 give parametrization for MEC part of W1 from Christy and Bosted
! In this case also W3 is related to W1 (acc. Martini and Ericsson)
!
! Case 4 describes double-differential data from MiniBooNE for neutrino
! and antineutrino scattering. It also describes the dd inclusive Xsection for
! neutrinos from T2K.
!******************************************************************************
module lepton2p2h

  use particleDefinition
  use eN_eventDefinition
  use leptonicID
  use CALLSTACK, only: TRACEBACK

  implicit none

  private
  public :: lepton2p2h_DoQE
  public :: lepton2p2h_DoDelta

  logical, save:: initflag = .true.

contains

  !****************************************************************************
  !****s* lepton2p2h/lepton2p2h_DoQE
  ! NAME
  ! subroutine lepton2p2h_DoQE(eN,outPart,XS)
  !
  ! PURPOSE
  ! Do all the electron induced 2p2h-QE scattering gamma* N1 N2 -> N1' N2'
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: OutPart -- the two produced nucleons
  ! * real :: XS -- the cross section
  !****************************************************************************
  subroutine lepton2p2h_DoQE(eN,outPart,XS)

    use neutrinoParms, only: readInput_2p2h

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    logical :: flagOK


    if (initFlag) then
       call readInput_2p2h
       initFlag = .false.
    end if

    XS = 0.0

    call lepton2p2h_SelectN2(eN,flagOK)
    if (.not.flagOK)  return ! ==> failure

    call lepton2p2h_FinalState(eN,outPart,.true.,flagOK)
    if (.not.flagOK) return ! ==> failure

    XS = lepton2p2h_XS(eN,outPart,.true.)
    outPart%perWeight=XS

  end subroutine lepton2p2h_DoQE

  !****************************************************************************
  !****s* lepton2p2h/lepton2p2h_DoDelta
  ! NAME
  ! subroutine lepton2p2h_DoDelta(eN,outPart,XS)
  !
  ! PURPOSE
  ! Do all the electron induced 2p2h-QE scattering gamma* N1 N2 -> N Delta
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: OutPart -- the two produced hadrons
  ! * real :: XS -- the cross section
  !****************************************************************************
  subroutine lepton2p2h_DoDelta(eN,outPart,XS)

    use neutrinoParms, only: readInput_2p2h

    type(electronNucleon_event), intent(inout) :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    real, intent(out) :: XS

    logical :: flagOK

    if (initFlag) then
       call readInput_2p2h
       initFlag = .false.
    end if

    XS = 0.0

    call lepton2p2h_SelectN2(eN,flagOK)
    if (.not.flagOK) return ! ==> failure

    call lepton2p2h_FinalState(eN,outPart,.false.,flagOK)
    if (.not.flagOK) return ! ==> failure

    XS = lepton2p2h_XS(eN,outPart,.false.)
    outPart%perWeight=XS

  end subroutine lepton2p2h_DoDelta

  !**************************************************************************
  !****s* lepton2p2h/lepton2p2h_SelectN2
  ! NAME
  ! subroutine lepton2p2h_SelectN2(eN)
  !
  ! PURPOSE
  ! Finds the second nucleon for the 2p2h collision
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  !
  ! OUTPUT
  ! * type(electronNucleon_event) :: eN -- a second nucleon is added
  !
  ! NOTES
  ! * The seond particle is generated analytically, not by selecting
  !   a testparticle from the real particle vector.
  ! * This is at a very basic level. You may add more sophisticated features
  !   as eq. two-particle correlatione etc.
  ! * A threshold check  Wfree>(2*mN+1MeV) is performed
  !****************************************************************************
  subroutine lepton2p2h_SelectN2(eN,flagOK)

    use mediumDefinition
    use mediumModule, only: mediumAt
    use densitymodule, only: FermiMomAt
    use random, only: rn, rnOmega
    use constants, only: mN
    use energyCalc, only: energyDetermination
    use minkowski, only: abs4
    use lorentzTrafo, only: lorentz

    type(electronNucleon_event), intent(inout) :: eN
    logical, intent(out) :: flagOK

    type(medium)   :: media
    type(particle) :: partN2
    real :: p,pF
    type(particle), dimension(2) :: nucleon
    real, dimension(0:3) :: momentum
    integer :: i

    ! 0) Set some defaults:
    call setToDefault(partN2)
    flagOK = .false.

    partN2%ID = 1
    partN2%mass=mN

    ! 1) select charge:
    media=mediumAt(eN%nucleon%pos)
    if (rn()*media%density.gt.media%densityProton) then
       partN2%charge = 0
    else
       partN2%charge = 1
    end if

    ! 2) select position:
    partN2%pos = eN%nucleon%pos

    ! 3) select 3-momentum:
    pF = FermiMomAt(partN2%pos)
    p = pF * rn()**(1./3.)
    partN2%mom(1:3) = p * rnOmega()
    partN2%mom(0) = sqrt(mN**2+p**2)

    call energyDetermination(partN2)

    ! 4) change the eN information:

    eN%nucleon2 = partN2

    nucleon(1) = eN%nucleon
    nucleon(2) = eN%nucleon2
    momentum = eN%boson%mom+nucleon(1)%mom+nucleon(2)%mom
    eN%betacm = momentum(1:3)/momentum(0)
    eN%W = abs4(momentum)

    ! we calculate Wfree in the CM system:
    do i=1,2
       nucleon(i)%pos = 999999999
       call lorentz(eN%betacm,nucleon(i)%mom)
       nucleon(i)%mom(0) = FreeEnergy(nucleon(i))
       call lorentz(-eN%betacm,nucleon(i)%mom)
    end do

    momentum = eN%boson%mom+nucleon(1)%mom+nucleon(2)%mom
    eN%W_free = abs4(momentum)
    ! invariant mass of two-nucleon + boson system

    if (eN%W_free.le.2*mN+0.001) return ! ===> failure

    ! 5) abuse of 'offshellParameter' for storage of density,
    !    needed for the 'cross section' calculation in cases 1 - 3

    eN%nucleon2%offshellPar = media%density
    flagOK = .true.

  end subroutine lepton2p2h_SelectN2


  !****************************************************************************
  !****s* lepton2p2h/lepton2p2h_FinalState
  ! NAME
  ! subroutine lepton2p2h_FinalState(eN,outPart,DoQE,flagOK)
  ! PURPOSE
  ! Generate the final state of the electron 2p2h event
  !****************************************************************************
  subroutine lepton2p2h_FinalState(eN,outPart,DoQE,flagOK)
    use mediumDefinition
    use mediumModule, only: mediumAt
    use collisionNumbering, only: pert_numbering
    use master_2Body, only: setKinematics
    use propagation, only: updateVelocity
    use IDtable, only: nucleon, delta
    use random, only: rn
    use particleProperties, only: hadron
    use baryonWidthMedium, only: get_MediumSwitch_coll

    type(electronNucleon_event), intent(in)    :: eN
    type(particle),dimension(:), intent(inout) :: OutPart
    logical,                     intent(in)    :: DoQE
    logical,                     intent(out)   :: flagOK

    type(particle), dimension(2) :: pairIN
    type(medium)         :: media
    real, dimension(1:3) :: betaToLRF
    integer :: i, ChargeIn

    flagOK = .false.

    if (size(OutPart).lt.2) call TRACEBACK('OutPart array too small.')

    pairIN = (/ eN%boson, eN%nucleon /)
    media=mediumAt(eN%nucleon%pos)

    ChargeIn = eN%nucleon%Charge+eN%nucleon2%Charge

    call setToDefault(OutPart)

    if (DoQE) then !=== N N final state ===
       OutPart%ID =     (/ nucleon, nucleon /)

       select case (eN%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          OutPart%Charge = (/ eN%nucleon%Charge, eN%nucleon2%Charge /)

       case (2)      ! == CC
          select case (ChargeIn)
          case (0)
             OutPart%Charge = (/ 0, 1 /)
          case (1)
             OutPart%Charge = (/ 1, 1 /)
          case (2)
             return ! ==> failure
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select

       case (-2)     ! == antiCC
          select case (ChargeIn)
          case (0)
             return ! ==> failure
          case (1)
             OutPart%Charge = (/ 0, 0 /)
          case (2)
             OutPart%Charge = (/ 0, 1 /)
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select

       end select

    else           !=== N Delta final state ===
       OutPart%ID =     (/ nucleon, delta /)

       select case (eN%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          OutPart(1)%Charge = nint(rn())
          OutPart(2)%Charge = eN%nucleon%Charge+eN%nucleon2%Charge &
               &  - OutPart(1)%Charge
       case (2)      ! == CC
          select case (ChargeIn)
          case (0)
             OutPart(1)%Charge = nint(rn())
             OutPart(2)%Charge = 1 - OutPart(1)%Charge
          case (1)
             OutPart(1)%Charge = 1
             if (rn()<1./3.) OutPart(1)%Charge = 0  ! Delta++ n : Delta+ p = 1:2
             !based on counting diagrams, but better ideas are needed
             OutPart(2)%Charge = 2 - OutPart(1)%Charge
          case (2)
             OutPart%Charge = (/ 1, 2 /)
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select
       case (-2)     ! == antiCC
          select case (ChargeIn)
          case (2)
             OutPart(1)%Charge = nint(rn())
             OutPart(2)%Charge = 1 - OutPart(1)%Charge
          case (1)
             OutPart(1)%Charge = 0
             if (rn()<3./4.) OutPart(1)%Charge = 1  ! Delta0 n : Delta- p = 4:3
             !based on counting diagrams, but better ideas are needed
             OutPart(2)%Charge =  - OutPart(1)%Charge
          case (0)
             OutPart%Charge = (/ 0, -1 /)
          case default
             call TRACEBACK('ChargeIn not allowed')
          end select
       end select

       ! The following is in order to avoid problems in massass:
       if (.not.get_MediumSwitch_coll()) then
          ! minimal value: 0.938 + Delta-MinMass + epsilon
          if (eN%W_free .lt. hadron(1)%mass+hadron(2)%minmass+0.005) then
             return ! ==> failure
          end if
       end if


    end if

    OutPart%anti=.false.
    OutPart%pert=.true.

    OutPart%perWeight=0. ! perturbative weight = XS (here only dummy)

    do i=1,2
       OutPart(i)%pos=eN%nucleon%pos
       OutPart(i)%event=pert_numbering(eN%nucleon)
    end do

    ! Now the final-state two nucleons get started with proper 4-momentum

    betaToLRF=0.
    call setKinematics (eN%W, eN%W_free, betaToLRF, eN%betacm, media, pairIn, &
         & OutPart(1:2), flagOK, .false.)
    if (.not.flagOK) return ! ==> failure

    call updateVelocity(OutPart)

  end subroutine lepton2p2h_FinalState

  !****************************************************************************
  !****f* lepton2p2h/lepton2p2h_XS
  ! NAME
  ! real function lepton2p2h_XS(eN,outPart,DoQE)
  ! PURPOSE
  ! calculate the electron induced 2p2h-QE cross section
  ! INPUTS
  ! * type(electronNucleon_event) :: eN -- electron-Nucleon event info
  ! * type(particle),dimension(:) :: OutPart -- the outgoing particles
  ! * logical :: DoQE -- .true. for NN final state, .false. for N Delta
  ! OUTPUT
  ! * the function value
  ! NOTES
  ! * One has to give a realistic parametrization of the matrix element
  ! * If one randomly selects the position of the second particle, one
  !   has to respect this in the XS calculation (and  not
  !   multiply it with the density at the position)
  !****************************************************************************
  real function lepton2p2h_XS(eN,outPart,DoQE)

    use minkowski, only: abs4Sq
    use constants, only: pi,twopi,mN,hbarc,alphaQED,GF,coscab,mW,mZ
    use twoBodyTools, only: pCM_sqr
    use AZN
    use neutrinoParms, only: ME_Version, ME_Norm_QE, ME_Norm_Delta, &
         ME_Mass_QE, ME_Mass_Delta, ME_Transversity, ME_LONG, &
         ME_ODW, inmedW, T, Adep

    type(electronNucleon_event), intent(in) :: eN
    type(particle),dimension(:), intent(in) :: OutPart
    logical,                     intent(in) :: DoQE

    real :: mf1_2,mf2_2  !squares of final state nucleon masses in 2p2h process
    real :: k1 !,k absolute values of the 3-momentum of the in/outgoing leptons
    real :: sqpp,pcm2 ! sqpp = (q + p + p2)^2, pcm2 = (p_cm)^2
    !    real :: d2PS      !d2PS = two-body phase space, includes factor (2pi)^4
    real :: ME        ! contraction of lepton and hadron tensor


    ! NOTE: lepton tensor contains coupling constants and propagator**2
    real :: couplProp    !coupling constant times propagator^2
    real :: Q2 !,omega
    integer :: iP

    call AZNsub(Atarget,Ntarget,Ztarget)

    if (Atarget < 2) then
       lepton2p2h_XS = 0.
       return
    end if


    Q2=eN%Q2
    !    omega=eN%boson%mom(0)              !omega = energy transfer

    mf1_2 = abs4Sq(outPart(1)%mom)
    mf2_2 = abs4Sq(outPart(2)%mom)
    k1=absMom(en%lepton_out)
    !    k =absMom(en%lepton_in)

    if (ME_version < 4) then

       lepton2p2h_XS = en%lepton_out%mom(0)/en%lepton_in%mom(0)
       ! correct for both e and nu

       sqpp = abs4Sq(eN%boson%mom+eN%nucleon%mom+eN%nucleon2%mom)
       !   initial state Mandelstam s of boson + 2 nucleons
       pcm2=pCM_sqr(sqpp, mf1_2, mf2_2)
       !   final state nucleon momentum squared in 2N cm system

       !   dOmega = 4. * pi   ! for isotropic phase space of 2 outgoing nucleons
       !   d2PS = 1.0/(4.*pi)**2  * sqrt(pcm2/sqpp) * dOmega
       !   includes factor (2pi)^4, different from PDG

       !   d2PS = 1.0/(4.*pi)  * sqrt(pcm2/sqpp)
       !   2-body final state phase space in cm system, integrated over NN angle

       lepton2p2h_XS=lepton2p2h_XS &
            * sqrt(pcm2/(16.0*sqpp)) & ! <-- the deltas
            * 2*twopi               ! <-- the angular integration

       ! Now XS times probability for 2nd nucleon to be at same position
       lepton2p2h_XS=lepton2p2h_XS &
            * eN%nucleon2%offshellPar ! <-- abuse !! = media%density


       ! Now we have to calculate the Matrixelement:
       select case (ME_Version)
       case (1)
          ME = ME_const(eN)
       case (2)
          ME = ME_transverse(eN)*0.635/( eN%lepton_in%mom(0) )**2
       case (3)
          ME = ME_Dipole_transverse(eN)
       end select

       lepton2p2h_XS=lepton2p2h_XS* ME / (2. *eN%nucleon2%mom(0))  &
            & * 1/(twopi**5 *8. *eN%nucleon%mom(0))

       ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2
       !      = (197/1000)**2 * 10 mb
       ! Now the cross section is given in units of mb/GeV:

       lepton2p2h_XS = lepton2p2h_XS*hbarc**2*10.

       ! Symmetry-Factor:
       if (IsSamePart(OutPart(1),OutPart(2))) lepton2p2h_XS = lepton2p2h_XS *0.5
       if (IsSamePart(eN%nucleon,eN%nucleon2)) lepton2p2h_XS = lepton2p2h_XS *0.5


    else       ! for Bosted parametrization of MEC in W1

       iP = abs(eN%idProcess)
       select case (iP)

       case (1)
          couplProp = (4*alphaQED**2)/Q2**2 * en%lepton_out%mom(0)**2
       case (2)
          couplProp = (GF*coscab)**2/(2*pi**2) * (mW**2/(mW**2 + Q2))**2 * en%lepton_out%mom(0)**2
       case (3)
          couplProp = GF**2/(2*pi**2) * (mZ**2/(mZ**2 + Q2))**2 * en%lepton_out%mom(0)**2 *1./2.
          ! factor 1/2 for NC because of outgoing neutrino with only one handedness
       case default
          write(*,*) 'reaction type must be EM, CC or NC'
          stop
       end select

       ! Now we have to calculate the Matrixelement:
       ME = ME_W1W2W3(eN)
       !   In this case the parametrisation of the structure function is taken from Bosted and Mamyan
       !   Fit to electron data. The fit does not contain the Pauli-blocking!


       ! use parametrization of W1, W2, W3


       lepton2p2h_XS=couplProp * ME/Atarget
       !    The scaling of ME with 1/Atarget is necessary
       !    since all other cross sections are given in 1/Atarget
       !    Parametrizations are given for nucleus with Atarget


       ! 1/GeV**2=1/1000**2/MeV**2=1/1000**2/(1/197 fm**2)=(197/1000)**2 fm**2
       ! = (197/1000)**2 * 10 mb
       ! Now the cross section is given in units of mb/GeV:

       lepton2p2h_XS = lepton2p2h_XS*hbarc**2*10.

       !   No symmetry factors since the Bosted-Christy
       !   X-sections are fitted to exp.

    end if


  contains
    !**************************************************************************
    !****if* lepton2p2h_XS/ME_const
    ! NAME
    ! real function ME_const(eN)
    ! PURPOSE
    !**************************************************************************
    real function ME_const(eN)

      type(electronNucleon_event), intent(in) :: eN
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME_const=1.0e-5*ME_Norm_QE(iP)
      else           !=== N Delta final state ===
         ME_const=1.0e-5*ME_Norm_Delta(iP)
      end if

    end function ME_const

    !**************************************************************************
    !****if* lepton2p2h_XS/ME_transverse
    ! NAME
    ! real function ME_transverse(eN)
    ! PURPOSE
    !**************************************************************************
    real function ME_transverse(eN)
      use minkowski, only: metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=4.8e4*ME_Norm_QE(iP)
      else           !=== N Delta final state ===
         ME=4.8e4*ME_Norm_Delta(iP)
      end if

      do mu=0,3
         do nu=0,3
            dummy(mu,nu)=eN%boson%mom(mu)*eN%boson%mom(nu)/eN%Q2
         end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%mom,  &
           & eN%lepton_out%mom)

      ME_transverse = Contract(hadronTens,leptonTens)

    end function ME_transverse


    !**************************************************************************
    !****if* lepton2p2h_XS/ME_Dipole_transverse
    ! NAME
    ! real function ME_Dipole_transverse(eN)
    ! PURPOSE
    ! calculate the 2p2h matrix element according to
    ! W_1(g_munu -q_um q_nu /Q2) * L^munu
    ! so that the contribution is only to the transverse part
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !**************************************************************************
    real function ME_Dipole_transverse(eN)
      use minkowski, only: metricTensor, Contract
      use leptonTensor

      type(electronNucleon_event), intent(in) :: eN
      complex, dimension(0:3,0:3) :: leptonTens, hadronTens, dummy
      integer :: mu,nu
      real :: ME
      integer :: iP

      iP = abs(eN%idProcess) ! anti-... same as EM, NC and CC

      if (DoQE) then !=== N N final state ===
         ME=8.0e4*ME_Norm_QE(iP)*(1. + eN%Q2/ME_Mass_QE(iP)**2)**(-4)
      else           !=== N Delta final state ===
         ME=8.0e4*ME_Norm_Delta(iP)*(1+eN%Q2/ME_Mass_Delta(iP)**2)**(-4)
      end if


      do mu=0,3
         do nu=0,3
            dummy(mu,nu)=eN%boson%mom(mu)*eN%boson%mom(nu)/eN%Q2
         end do
      end do

      hadronTens = ME *( - metricTensor - dummy )
      leptonTens = leptonicTensor(eN%idProcess,eN%lepton_in%mom,  &
           & eN%lepton_out%mom)

      ME_Dipole_transverse = Contract(hadronTens,leptonTens)

    end function ME_Dipole_transverse

    !**************************************************************************
    !****if* lepton2p2h_XS/ME_W1W2W3
    ! NAME
    ! real function ME_W1W2W3(eN)
    ! PURPOSE
    ! to calculate the 2p2h contribution to the inclusive cross sections for
    ! electrons and neutrinos, cross section depends on all 3 structure functs
    !
    ! NOTES
    !
    ! You have full access to all incoming and outgoing particles:
    ! * eN%lepton_in  -- incoming lepton
    ! * eN%nucleon    -- incoming nucleon 1
    ! * eN%nucleon2   -- incoming nucleon 2
    !
    ! exchanged boson:
    ! * eN%boson      -- exchanged boson
    !
    ! even without considering the final state particles, you know the kind
    ! of process via 'eN%idProcess', which may take the values EM,NC,CC and
    ! also antiEM,antiNC,antiCC
    !**************************************************************************
    real function ME_W1W2W3(eN)

      use particleDefinition
      use eN_eventDefinition
      use leptonicID
      use constants, only: g_A
      use FF_QE_nucleonScattering, only: MA_in,MV2,mup,mun

      type(electronNucleon_event), intent(in) :: eN

      integer :: nuswitch= 0 ! switch for neutrino/antineutrino in structure function
      ! nuswitch = 0 for em, = +1 for neutrino, -1 for antineutrino

      real :: sinsqthetahalf,cossqthetahalf
      real :: omega, Q2   ! energy transfer, four-momentum transfer
      integer :: IP
      real :: GMV0,GA0,MA,GMV,GA,GMV2,GA2


      Q2=eN%Q2                     ! Q^2
      omega=eN%boson%mom(0)              ! omega = energy transfer

      !    vector and axial coupling constants and cutoff masses
      GMV0 = mup - mun           ! vector component of magnetic coupling GM
      GA0 = - g_A
      MA = MA_in

      !   vector and axial coupling formfactors

      GMV = GMV0 * 1./(1 + Q2/MV2)**2
      GA = GA0 * 1./(1 + Q2/MA**2)**2
      GMV2 = GMV**2
      GA2 = GA**2

      IP = abs(eN%IdProcess)
      if (IP==1) then
         nuswitch = 0
         GA2=0
      else
         nuswitch = sign(1,eN%IdProcess)
      end if

      sinsqthetahalf = 0.5*(1. - en%lepton_out%mom(3)/k1)
      cossqthetahalf = 1 - sinsqthetahalf

      ME_W1W2W3 =  sinsqthetahalf * 2*W1(Q2,omega,GMV2,GA2) &
           &  + cossqthetahalf * W2(Q2,omega,GMV2,GA2) &
           &  - nuswitch*(en%lepton_in%mom(0)+en%lepton_out%mom(0))/mN &
           &  * W3(Q2,omega,GMV,GA) * sinsqthetahalf

      ME_W1W2W3 = ME_W1W2W3 * ME_Norm_QE(IP)

      if (nuswitch /= 0 ) ME_W1W2W3 = ME_W1W2W3 * 2. * (T + 1)
      ! multiplication with isospin factor for neutrinos only, T from
      ! module neutrinoParms
      ! Enforce positivity constraint on structure functions by keeping
      ! matrixelement positive

      if (ME_W1W2W3 < 0) then
!!$         write(*,*) 'enforce positivity constraint in 2p2h, set to 0:',ME_W1W2W3
!!$         write(*,*) 'Q2 =',Q2,'omega =', omega
         ME_W1W2W3 = 0.
      end if

    end function ME_W1W2W3

    !**************************************************************************
    !****if* ME_W1W2W3/W1
    ! NAME
    ! real function W1(Q2,omega,GM2,GA2)
    ! PURPOSE
    ! Structure function W1 (for electrons: W1E, for neutrinos: W1NU)
    !**************************************************************************
    real function W1(Q2,omega,GM2,GA2)

      real, intent(in) :: Q2,omega,GM2,GA2

      IP = abs(eN%idProcess)
      select case (IP)

      case (1)                                    ! electron
         if (DoQE) then !=== N N final state ===
            W1= W1E(Q2,omega)
         else           !=== N Delta final state ===
            W1= W1E(Q2,omega)
         end if

      case (2:)                                    ! CC and NC

         if (DoQE) then !=== N N final state ===
            W1= W1NU(Q2,omega,GM2,GA2)
         else           !=== N Delta final state ===
            W1= W1NU(Q2,omega,GM2,GA2)
         end if

      case default
         write(*,*) 'ProcessID Error in 2p2h'
      end select

    end function W1


    !**************************************************************************
    !****if* W1/W1E
    ! NAME
    ! real function W1E(Q2,omega)
    ! PURPOSE
    ! Structure function for electrons, parametrizations for MEC term only
    !**************************************************************************
    real function W1E(Q2,omega)

      use constants, only: mN
      use nucleus, only: getTarget
      use nucleusdefinition
      use minkowski, only: abs3,abs4sq

      real, intent(in) :: Q2,omega

      type(tnucleus), pointer :: targetNuc
      integer :: Atarget

      real :: a1,b1,c1,t1,dw2,Wrecsq
      real :: ENfree,pNplusQ
      real :: p18,f

      ! parameters from Bosted arXiV:1203.2262:
      real, parameter :: p0=0.005138
      real, parameter :: p1=0.980710
      real, parameter :: p2=0.046379
      real, parameter :: p3=1.643300
      real, parameter :: p4=6.982600
      real, parameter :: p5=-0.226550
      real, parameter :: p19=-0.045536

      real, save :: q02
      !real, save :: a1
      real, save :: a2
      !real, save :: b1
      !real, save :: c1
      real, save :: b2
      real, save :: c2
      real :: numin, Wrecsqmin, Y, F1A, F1B, F1MEC

      targetNuc => getTarget()
      Atarget = targetNuc%mass


      ! select Fermi-smearing
      select case (inmedW)
      case (1)
         Wrecsq = mN**2 + 2*mN*omega - Q2
      case (2)
         Wrecsq = abs4sq(eN%boson%mom+eN%nucleon%mom)
      case (3)
         ENfree = sqrt(abs3(eN%nucleon%mom)**2  +mN**2)
         pNplusQ = abs3(eN%nucleon%mom + eN%boson%mom)
         Wrecsq = (omega + ENfree)**2 - pNplusQ**2
      case default
         write(*,*) 'wrong case for Wrecsq in 2p2h'
      end select

      if (Wrecsq <= 0.0) then
         W1E = 0.0
         return
      end if


      select case (ME_Version)

      case (4)
         !   This case returns the value for the structure function W1(MEC)
         !   fitted by E.Christy to inclusive electron scattering data,
         !   E. Christy, priv. comm., August 2015, fitted for C12

         a1 = 6.049*Q2**2 * exp(-1.0*Q2/1.298)/(0.314+Q2)**5.708
         b1 = 0.791 + 0.154*Q2
         c1 = 0.290
         t1 = (Wrecsq - b1)**2/c1**2/2.
         dw2 = Wrecsq + Q2 * (1. - 1./2.2) - 1.0*mN*mN
         if (dw2 < 0.0) dw2 = 0.

         if (DoQE) then !=== N N final state ===
            W1E = a1 * (exp(-1.*t1)*sqrt(dw2))/mN
         else           !=== N Delta final state ===
            W1E = a1 * (exp(-1.*t1)*sqrt(dw2))/mN
         end if

         !    A-dependence of structure function
         W1E = W1E * FAdep_2p2h(Atarget)

      case (5)
         !   This case returns the value for the structure function W1(MEC)
         !   fitted by Mamyan and Bosted to inclusive electron scattering data,
         !   http://arxiv.org/abs/arXiv:1203.2262, 2012

         select case (Atarget)
         case (5:20)
            p18 = 215
         case (21:50)
            p18 = 235
         case (51:)
            p18 =230
         end select

         f=(1. + max(0.3,Q2)/p3 )**p4 /( omega**p5 &
              * (1.+p18*Atarget**(1.+p19*Q2/2./mN/omega)))

         !     In the Mamyan-Bosted paper eq.(10) is wrong, corrected here,
         !     following Bosted code

         if (DoQE) then !=== N N final state ===
            W1E = p0*exp( -(Wrecsq-p1)**2/p2 )/f/mN
            !        write (*,*) 'MEC=',MEC, '    Q2=',Q2,'   omega=',omega
         else           ! === N Delta final state ===
            W1E = p0*exp( -(Wrecsq-p1)**2/p2 )/f /mN
         end if

      case(6)
         ! this case contains the parametrization of Bodek and Christy in
         ! Phys.Rev.C 106 (2022) 6, L061305 • e-Print: 2208.14772 [hep-ph]
         ! Original paper contains typos, corrected here through privat comm.
         ! with Eric Christy
         numin = 0.0165
         if (omega < numin) then
            W1E = 0.0
         else
            q02 = 1.e-4
            a1 = 0.091648
            a2 = 0.01045
            b1 = 0.77023
            c1 = 0.077051 + 0.26795*Q2
            b2 = 1.275
            c2 = 0.265
            Wrecsqmin = mN**2 + 2*mN*numin - Q2
            Y = Atarget * exp(-Q2**2/12.715) * (Q2 + q02)**2/(0.13380 + Q2)**6.90679
            F1A = a1 * Y * (Wrecsq - Wrecsqmin)**1.5 * exp(-(Wrecsq - b1)**2/(2*c1**2))
            F1B = a2 * Y * (Q2 + q02)**1.5 * exp(-(Wrecsq - b2)**2/(2*c2**2))
            F1MEC = max((F1A + F1B),0.0)
            W1E = F1MEC/mN
         end if

            W1E = W1E * FAdep_2p2h(Atarget) * 12./Atarget

      case default
         write(*,*) 'ME_Version does not exist: ',ME_Version
         call TRACEBACK()

      end select

    end function W1E

    !**************************************************************************
    !****if* W1E/FAdep_2p2h
    ! NAME
    ! real function FAdep_2p2h(Atarget)
    ! PURPOSE
    !**************************************************************************
    real function FAdep_2p2h(Atarget)

      integer, intent(in) :: Atarget

      select case (Adep)

      case (1)
         FAdep_2p2h =  (0.145 - 0.147*Atarget**(-1./3.))/ &
              &    (0.145 - 0.147*12.**(-1./3.)) * Atarget/12.
         ! The A-dependence here is taken from Mosel&Gallmeister,Phys. Rev. C94 (2016) 3, 035502
         ! normalized to C12, since Christy-Bosted fit was for this nucleus

      case (2)
         FAdep_2p2h = Atarget/12.

      case default
         write (*,*) 'error in 2p2h mass dependence'

      end select

    end function FAdep_2p2h

    !**************************************************************************
    !****is* W1/Transverse_resp
    ! NAME
    ! subroutine Transverse_resp(Q2,omega,GM2,RT,kinfact)
    ! PURPOSE
    ! transverse response = reduced structure function
    ! NOTES
    ! cf. ME_ODW
    !**************************************************************************
    subroutine Transverse_resp(Q2,omega,GM2,RT,kinfact)

      use constants, only: mN
      Use neutrinoParms, only : VAfact

      real, intent(in) :: Q2,omega,GM2
      real, intent(out) :: RT, kinfact
      real :: qvec2,W

      qvec2 = Q2 + omega**2
      W = sqrt(mN**2 + 2*mN*omega - Q2)
     
      call VAfact(ME_ODW,qvec2,omega,kinfact)

      RT = 1./(kinfact*GM2) * W1E(Q2,omega)

    end subroutine Transverse_resp


    !**************************************************************************
    !****if* W1/W1NU
    ! NAME
    ! real function W1NU(Q2,omega,GMV2,GA2)
    ! PURPOSE
    ! structure function W1 for neutrino-induced CC and NC MEC process
    !**************************************************************************
    real function W1NU(Q2,omega,GMV2,GA2)

      real, intent(in) :: Q2,omega,GMV2,GA2
      real :: RT,kinfact

      call Transverse_resp(Q2,omega,GMV2,RT,kinfact)

      W1NU = (kinfact*GMV2  + GA2) * RT

    end function W1NU

    !**************************************************************************
    !****if* ME_W1W2W3/W2
    ! NAME
    ! real function W2(Q2,omega,GMV2,GA2)
    ! PURPOSE
    ! Structure function W2
    !**************************************************************************
    real function W2(Q2,omega,GMV2,GA2)

      use constants, only: mN

      real, intent(in) :: Q2,omega,GMV2,GA2
      real, parameter :: MDelta = 1.232
      real :: qvec2
      integer :: IP

      IP = abs(eN%idProcess)

      qvec2 = Q2 + omega**2
      !   vector and axial coupling constants and cutoff masses

      W2 = ME_Transversity(IP) * Q2/qvec2 * W1(Q2,omega,GMV2,GA2)
      ! W2: term necessary for purely transverse interaction,
      ! could be turned off by setting ME_Transversity = 0, default = 1

      if (ME_Long(IP) > 0) &
           W2 = W2  + GA2*(MDelta - mN)**2/(2.*(Q2 + omega**2))* 1./(1 + Q2/0.3**2)**2&
           & * ME_LONG(iP) * 1.e-5
      ! Structure of longitudinal term follows Martini et al (PRC 2009)
      ! ME_LONG: allows to turn off longitudinal contribution to 2p2h,
      ! default = 0
      ! W2: 2nd term for longitudinal response, strength function untested
    end function W2

    !**************************************************************************
    !****if* ME_W1W2W3/W3
    ! NAME
    ! real function W3(Q2,omega,GMV,GA)
    ! PURPOSE
    ! Structure function W3, relevant only for neutrinos
    ! W3 is directly related to W1, either according to Martini and Ericsson,
    ! or to O'Connell, Donnelly, Walecka
    !**************************************************************************
    real function W3(Q2,omega,GMV,GA)

      real, intent(in) :: Q2,omega,GMV,GA
      real :: RT,kinfact,GMV2

      integer :: IP

      IP = abs(eN%idProcess)
      if (IP==1) then
         W3 = 0.0
         return
      end if

      GMV2 = GMV**2

      call Transverse_resp(Q2,omega,GMV2,RT,kinfact)

      W3 = 2*GA*GMV*RT

    end function W3

  end function lepton2p2h_XS


end module lepton2p2h
