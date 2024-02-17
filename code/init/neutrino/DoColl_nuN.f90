!******************************************************************************
!****m* /Coll_nuN
! NAME
! module Coll_nuN
!
! PURPOSE
! neutrino Nucleon scatterings with PYTHIA
!******************************************************************************
module Coll_nuN

  use particleDefinition
  use CollTools
  use hadronFormation
  use eN_eventDefinition

  implicit none
  private

  public :: DoColl_nuN_Py
  public :: AnaEstimate
  public :: AnaEstimatePythia
!!$   PUBLIC :: GetPythiaXY
  public :: CalcXY

  integer, save :: MSEL0,MSTP21,MSTP23,MSTP32,MSTP61,MSTP71,MSTP11
  real, save :: CKIN01,CKIN05,CKIN06,CKIN39,PARP91
  real, save :: PMAS1,PMAS2,PMAS150,PMAS152,PMAS153,PMAS156
  integer, save :: maxMSTI52


contains
  !****************************************************************************
  !****s* Coll_nuN/DoColl_nuN_Py
  ! NAME
  ! subroutine DoColl_nuN_Py(eNev,outPart,flagOK,cross,respectHad,doMSTP23)
  !
  ! PURPOSE
  ! generate a high energy neutrino event with PYTHIA
  !
  ! returned cross section is dsigma/dE'dcost in mb/GeV
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev     -- kinematics
  ! * logical :: respectHad -- flag to respect hadronization effects
  ! * logical :: doMSTP23 -- flag to set MSTP(23)=1
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: outPart  -- outgoing particles
  ! * logical                     :: flagOK   -- .TRUE., if everything okay
  ! * real                        :: cross    -- cross section
  !
  ! NOTES
  ! With the PYTHIA-Option MSTP(23)=1, the routine PYREMN tries to
  ! keep the chosen x and Q2 values of the partonic process also for the
  ! outer hadronic process.
  !****************************************************************************
  subroutine DoColl_nuN_Py(eNev,outPart,flagOK,cross,respectHad,doMSTP23)

    use GetLeading
    use EventInfo_HiLep
    use LorentzTrafo
    use rotation
    use PythiaSpecFunc
    use eN_event
    use CallStack
    use output
    use constants, only: mN

    type(electronNucleon_event), intent(in)   :: eNev
    type(particle),dimension(:), intent(inout):: outPart
    logical,                     intent(out)  :: flagOK
    real,                        intent(out)  :: cross
    logical,                     intent(in)   :: respectHad
    logical,                     intent(in)   :: doMSTP23

    ! common blocks:

    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    integer MSTU,MSTJ
    double precision PARU,PARJ
    SAVE /PYDAT1/

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    COMMON/PYINT1/MINT(400),VINT(400)
    integer MINT
    double precision VINT
    SAVE /PYINT1/

    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
    integer NGENPD,NGEN
    double precision XSEC
    SAVE /PYINT5/

    ! internal variables:

    real :: Wfree, Whad
    character*20 :: cTarget,cBeam,cNu(3) = (/'nu_e  ','nu_mu ','nu_tau'/)
    character*20 :: charged_lepton(3) = (/'e-  ','mu- ','tau-'/)
    character*20 :: charged_antilepton(3) = (/'e+  ','mu+ ','tau+'/)
    real, dimension(0:3) :: pB,pT
    real :: phi,theta,phiLepton,phi2
    real, dimension(1:3) :: betacm
    integer :: iEv, nEV=10, iTry,i
    integer :: iL_lep
    integer :: MSTP23def

    real :: x_wish,y_wish
    logical :: flagXY

    ! switch on/off ISR,FSR and intrinsic kT:
    logical, parameter :: doCollinear = .TRUE.

    ! switch on/off masses of quarks and diquarks:
    logical, parameter :: doMassless = .TRUE.

    ! special treatment of events with cluster -> 1 particle fragmentation
    ! 0: ignore
    ! 1: redo
    ! 2: return with failure flag; this modifies the cross section!
    integer :: treat1particle = 2

    ! W_hadron is strange, but check whether it is too strange:
    logical :: doCheckW = .TRUE.


    !...reset flags according input:
    doCheckW = .not.doMSTP23
    if (doMSTP23) then
       MSTP23def = 1
    else
       MSTP23def = 0
    end if

    if (.not.respectHad) then
       treat1particle = 1
       doCheckW = .FALSE.
    end if

    !...set default output

    outPart%ID = 0
    flagOK = .FALSE.
    cross = 0.0

    !...save original PYTHIA values:
    call StorePyVars

    !...set up PYTHIA:
    call SetSomeDefaults_PY

!    MSTP(122) = 2          ! write out maximum finding

    Wfree = eNev%W_free

    MSEL = 0       ! selection of subprocesses
    MSUB     = 0   ! all subprocesses off
    MSUB(10) = 1   ! selection of subprocesses

    MSTP(127) = 1  ! continue PYINIT after error

    select case (eNev%idProcess)
    case (-1,1) !=== EM
       MSTP(21) = 4  ! gamma*/Z  neutral current only
       MSTP(11) = 0  ! initial lepton w/o initial radiation, also unresolved because MSTP(12)=0 default
    case (-2,2) !=== CC
       MSTP(21) = 5 ! charged current only
    case (-3,3) !=== NC
       MSTP(21) = 3 ! neutral current, Z0 only
    case default
       write(*,*) 'idProcess=',eNev%idProcess
       call TRACEBACK('wrong idProcess!')
    end select

    MSTP(32)= 5 ! Q2 definition: -t

    CKIN( 1) = 0.1
    CKIN(39) = 0.1  ! W2_min in DIS (loose cut)
    CKIN( 6) = 0.0  ! singular if masses below this cut
    ! adjust the pThat cut off according W:
    if (2*CKIN(5).ge.Wfree) CKIN(5) = Wfree/2-1e-3

    if (doCollinear) then
       PARP(91)=0.00 ! width intrinsic kT
       MSTP(61)=0    ! master: ISR (QCD/QED)  switch initial-state showers off
       MSTP(71)=0    ! master: FSR (QCD/QED)  switch final-state showers off
    end if

    if (doMassless) then
       PMAS(1,1)=0.001 ! Mass of d quark
       PMAS(2,1)=0.001 ! Mass of u quark
       PMAS(150,1)=0.002 ! Mass of diquark
       PMAS(152,1)=0.002 ! Mass of diquark
       PMAS(153,1)=0.002 ! Mass of diquark
       PMAS(156,1)=0.002 ! Mass of diquark
    end if

    call Init_VM_Mass(Wfree,eNev%nucleon%pos)

    if (eNev%nucleon_free%charge.gt.0) then
       call PYNAME(2212,cTarget)          ! Target is proton
    else
       call PYNAME(2112,cTarget)          ! Target is neutron
    end if

    select case (eNev%idProcess)
    case (1)
       cBeam = trim(charged_lepton(eNev%idFamily))
    case (-1)
       cBeam = trim(charged_antilepton(eNev%idFamily))
    case (2,3)
       cBeam = trim(cNu(eNev%idFamily))
    case (-2,-3)
       cBeam = trim(cNu(eNev%idFamily))//'bar'
    end select

    !...set incoming vectors

    call CalcXY(eNev%lepton_in%mom,eNev%lepton_out%mom,&
         & eNev%nucleon_free%mom, x_wish, y_wish, flagXY)
    if (.not.flagXY) call TRACEBACK('Error in CalcXY')

!    write(*,*) 'wish: ', x_wish, eNev%Q2

    CKIN(21) = 1.0
    CKIN(23) = x_wish
    CKIN(24) = x_wish*1.001

    CKIN(35) = eNev%Q2
    CKIN(36) = eNev%Q2*1.001

    pB = eNev%lepton_in%mom
    pT = eNev%nucleon_free%mom

    call eNeV_GetLeptonCM(eNev, betacm,phi,theta,phiLepton)
!    theta = 0.0

    call lorentz(betacm,pB)
    call lorentz(betacm,pT)

    pB(1:3) = rotateZY(theta, phi, pB(1:3))
    pT(1:3) = rotateZY(theta, phi, pT(1:3))

    pB(1:3) = rotateZY(0.0, eNev%phiLepton, pB(1:3))
    pT(1:3) = rotateZY(0.0, eNev%phiLepton, pT(1:3))

    P(1,1:3) = pB(1:3)
    P(1,4)   = pB(0)
    P(1,5)   = eNev%lepton_in%mass

    P(2,1:3) = pT(1:3)
    P(2,4)   = pT(0)
    P(2,5)   = eNev%nucleon_free%mass

    !...Initialize Pythia

    MSTP(111) = 0          ! master switch fragmentation/decay
    MSTP(23)  = 0          ! fix x and Q2 for DIS

    call PYINIT('5MOM', cBeam, cTarget, 99.9)

    if (MSTI(53).ne.0) then
       write(*,*) 'Problems in PYINIT (x=',CKIN(23),&
            ' Q2=',CKIN(35),')'
       return ! -> FAILURE
    end if

    do iEv=1,nEv
       if (useJetSetVec) call GetJetsetVecINIT
       call PYEVNT

!       write(124,*) MSTI(52),MINT(51)
!       write(*,*) MSTI(52),MINT(51)

!!$          call PYGIVE('VINT( 1)=') ! sqrt(s)
!!$          write(*,*) '== tau:'
!!$          call PYGIVE('VINT(11)=') ! taumin
!!$          call PYGIVE('VINT(21)=') ! tau
!!$          call PYGIVE('VINT(31)=') ! taumax
!!$          write(*,*) '== y:'
!!$          call PYGIVE('VINT(12)=') ! ystmin
!!$          call PYGIVE('VINT(22)=') ! yst
!!$          call PYGIVE('VINT(32)=') ! ystmax
!!$          write(*,*) '== cos(theta):'
!!$          call PYGIVE('VINT(13)=') ! theta: ctnmin
!!$          call PYGIVE('VINT(14)=') ! theta: ctpmin
!!$          call PYGIVE('VINT(23)=') ! theta: cth
!!$          call PYGIVE('VINT(33)=') ! theta: ctnman
!!$          call PYGIVE('VINT(34)=') ! theta: ctpmax
!!$          write(*,*) '== x_1,2, Q2:'
!!$          call PYGIVE('VINT(41)=') ! x_1
!!$          call PYGIVE('VINT(42)=') ! x_2
!!$          call PYGIVE('VINT(52)=') ! inner Q2
!!$          call PYGIVE('VINT(54)=') ! outer Q2
!!$          write(*,*) '== shat,that,uhat:'
!!$          call PYGIVE('VINT(44)=')
!!$          call PYGIVE('VINT(45)=')
!!$          call PYGIVE('VINT(46)=')
!!$
!!$          call PYGIVE('MINT(43)=')
!!$          stop

    end do

    MSTP(111) = 1          ! master switch fragmentation/decay
    MSTP(23)  = MSTP23def  ! fix x and Q2 for DIS

    if (nEV.gt.1) then
       cross=XSEC(0,3) ! cross section of all subprocesses in mb
       if (cross.eq.0.0) then
          N = 0
          call ResetPyVars
          return ! -> FAILURE
       end if
    end if

    iTry = 0
    TryLoop: do
       iTry = iTry+1

       if (iTry.ge.20) then
          if (DoPr(-1)) then
             write(*,'(A,i4)')       'DoColl_nuN: itry=',iTry
             write(*,'(A,1P,3e13.5)')'          : ',eNev%lepton_in%mom(0),&
                  & eNev%lepton_out%mom(0),eNev%Q2
          end if

          N = 0
          call ResetPyVars
          return ! -> FAILURE
       end if
       call PYINIT('5MOM', cBeam, cTarget, 99.9)
       if (useJetSetVec) call GetJetsetVecINIT
       call PYEVNT
!       write(124,*) MSTI(52),MINT(51)
!       write(*,*) MSTI(52),MINT(51)
       !    call PYLIST(2)

       !...check success:

       if (MSTU(24).ne.0) cycle ! retry
       if (MINT(51).eq.2) cycle ! retry <---- Main Error !

       call MarkLepton(eNev,iL_lep)
       !write(*,*) 'iL_lep=',iL_lep

!       call PYLIST(2)

       if (useJetSetVec) then
          call GetJetsetVec(.TRUE.)
          call GetJetsetVecCheckT(-1d-5)
          call GetJetsetVecPYEDIT
       end if

       call GetLeading_PY         ! find leading particles

       ! check, whether it was a "cluster->1" deacy. Then the outgoing
       ! lepton was modified by the fragmentation

       if (treat1particle>0) then
          do i=N,1,-1
             if (K(i,1).ge.10) cycle ! not final particle
             if ((K(i,4).eq.3).and.(K(i,5).eq.1)) then
                if (treat1particle.eq.1) then
                   cycle TryLoop
                else
                   call ResetPyVars
                   return ! --> failure
                end if
             end if
          end do
       end if

       exit ! success, leave the Try-Loop

    end do TryLoop


    if (nEV.le.1) cross=XSEC(0,3) ! cross section of all subprocesses in mb

    cross = cross / ((CKIN(36)-CKIN(35))*(CKIN(24)-CKIN(23)))
    ! Transform from dsigma/dxdQ2 into dsigma/dxdy:
    cross = cross * eNev%Q2*eNev%lepton_in%mom(0)/eNev%boson%mom(0)
    ! Transform from dsigma/dxdy into dsigma/dcostdEprime:
    cross = cross * eNev%lepton_out%mom(0)/(mN*eNev%boson%mom(0))

    call ResetPyVars

    phi2 = atan2(P(iL_lep,2),P(iL_lep,1))

    call PYROBO(1,N, 0.0,phiLepton-phi2, 0.0,0.0,0.0)
    call PYROBO(1,N, theta,phi, betacm(1),betacm(2),betacm(3))

    if (useJetSetVec) then
       call GetJetsetVecPYROBO(0.0,  eNev%phiLepton, 0.0,0.0,0.0)
       call GetJetsetVecPYROBO(theta,phi, eNev%betacm(1),eNev%betacm(2),eNev%betacm(3))
    end if

!    call PYLIST(2)
    call PYEDIT(1)            ! clean up event list

    !...Copy Particles to ouput-vector:
    call SetVectorFromPYJETS(outPart, eNev%Q2)

    ! A final check:
    if (doCheckW) then
       if (doMassless .and. doCollinear) then

          ! in this case, the invariant mass of the hadronic system is always
          ! ~200 MeV below the desired W. Therefore we throw away all events,
          ! where the hadronic W is larger than the W at the photon vertex.

          Whad = sqrtS(outPart)
          if (Whad>Wfree) then
             outPart%ID = 0
             cross = 0.0
             return ! --> failure
          end if

       else

          ! do not yet know how to do this here !

       end if
    end if


    flagOK = .TRUE.

  end subroutine DoColl_nuN_Py


!!$  !*************************************************************************
!!$  !****s* Coll_nuN/DoColl_nuN_Py_fixMC
!!$  ! NAME
!!$  ! subroutine DoColl_nuN_Py_fixMC(eNev,outPart,flagOK, doMassless,cross)
!!$  !
!!$  ! PURPOSE
!!$  ! generate a high energy neutrino event with PYTHIA
!!$  !
!!$  ! Here the outgoing kinematics is fixed from the very beginning.
!!$  !
!!$  ! INPUTS
!!$  ! * type(electronNucleon_event) :: eNev     -- kinematics
!!$  ! * logical                     :: doMassless -- neglect masses of quarks and
!!$  !   diquarks
!!$  !
!!$  ! OUTPUT
!!$  ! * type(electronNucleon_event) :: eNev     -- kinematics
!!$  ! * type(particle),dimension(:) :: outPart  -- outgoing particles
!!$  ! * logical                     :: flagOK   -- .TRUE., if everything okay
!!$  ! * real                        :: cross    -- cross section
!!$  !
!!$  ! NOTES
!!$  ! We give as input the neutrino and the nucleon in medium and free nucleon
!!$  ! kinematics. We are trying to get as close as possible to the kinematics
!!$  ! of the outgoing lepton.
!!$  !
!!$  ! TODO: find a mechanism to detect unphysical kinematics/regions of x and y
!!$  !
!!$  ! TODO: We have to calculate the corresponding cross section
!!$  !
!!$  !*************************************************************************
!!$  subroutine DoColl_nuN_Py_fixMC(eNev,outPart,flagOK, doMassless, cross)
!!$
!!$    use constants, only: pi
!!$    use GetLeading
!!$    use EventInfo_HiLep
!!$    use LorentzTrafo
!!$    use rotation
!!$    use VMMassPythia
!!$    use eN_event
!!$    use minkowski, only: abs4,abs4sq,SP
!!$    use CallStack
!!$
!!$    implicit none
!!$
!!$    type(electronNucleon_event), intent(inout):: eNev
!!$    logical,                     intent(in)   :: doMassless
!!$    type(particle),dimension(:), intent(inout):: outPart
!!$    logical,                     intent(out)  :: flagOK
!!$    real,                        intent(out)  :: cross
!!$
!!$    ! common blocks
!!$
!!$    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
!!$    integer MSTU,MSTJ
!!$    double precision PARU,PARJ
!!$    SAVE /PYDAT1/
!!$
!!$    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
!!$    integer KCHG
!!$    double precision PMAS,PARF,VCKM
!!$    SAVE /PYDAT2/
!!$
!!$    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
!!$    integer N,NPAD,K
!!$    double precision P,V
!!$    SAVE /PYJETS/
!!$
!!$    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!!$    integer MSTP,MSTI
!!$    double precision PARP,PARI
!!$    SAVE /PYPARS/
!!$
!!$    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
!!$    integer MSEL,MSELPD,MSUB,KFIN
!!$    double precision CKIN
!!$    SAVE /PYSUBS/
!!$
!!$    COMMON/PYINT1/MINT(400),VINT(400)
!!$    integer MINT
!!$    double precision VINT
!!$    SAVE /PYINT1/
!!$
!!$    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
!!$    integer NGENPD,NGEN
!!$    double precision XSEC
!!$    SAVE /PYINT5/
!!$
!!$
!!$    ! internal variables:
!!$
!!$    real :: Wfree
!!$    character*20 :: cTarget,cNu(3) = (/'nu_e  ','nu_mu ','nu_tau'/)
!!$    real, dimension(0:3) :: pB,pT
!!$    real :: phi,phi2,theta,phiLepton
!!$    real, dimension(1:3) :: betacm
!!$    integer :: iEv, nEV=10000, iTry,i
!!$    integer :: iL_lep
!!$
!!$    real, dimension(0:3) :: pP,pQ,pK,pKprime
!!$    real :: pPpQ, pPpK, pQpQ
!!$
!!$    real :: x_wish,y_wish, x_try, y_try
!!$    real :: DistMin, Dist
!!$    logical :: flagXY
!!$
!!$  integer, save :: iEv_max = 100000
!!$  real, save :: DistMinCut = 1e-4
!!$
!!$    !...set default output
!!$
!!$    outPart%ID = 0
!!$    flagOK = .FALSE.
!!$    cross = 0.0
!!$
!!$    !...save original PYTHIA values
!!$
!!$    MSEL0 = MSEL
!!$    MSUB10 = MSUB(10)
!!$    MSTP21 = MSTP(21)
!!$    MSTP32 = MSTP(32)
!!$    MSTP61 = MSTP(61)
!!$    MSTP71 = MSTP(71)
!!$    CKIN01 = CKIN( 1)
!!$    CKIN05 = CKIN( 5)
!!$    CKIN06 = CKIN( 6)
!!$    CKIN39 = CKIN(39)
!!$    PARP91 = PARP(91)
!!$
!!$    PMAS1   = PMAS(1,1)
!!$    PMAS2   = PMAS(2,1)
!!$    PMAS150 = PMAS(150,1)
!!$    PMAS152 = PMAS(152,1)
!!$    PMAS153 = PMAS(153,1)
!!$    PMAS156 = PMAS(156,1)
!!$
!!$
!!$    !...set up PYTHIA
!!$
!!$    call SetSomeDefaults_PY
!!$
!!$    Wfree = eNev%W_free
!!$
!!$    MSEL = 0       ! selction of subprocesses
!!$    MSUB(10) = 1   ! selction of subprocesses
!!$
!!$    select case(eNev%idProcess)
!!$    case (-1,1) !=== EM
!!$       write(*,*) 'not prepared for EM!'
!!$       call TRACEBACK('wrong idProcess!')
!!$    case (-2,2) !=== CC
!!$       MSTP(21) = 5 ! charged current only
!!$    case(-3,3) !=== NC
!!$       MSTP(21) = 3 ! neutral current, Z0 only
!!$    case default
!!$       write(*,*) 'idProcess=',eNev%idProcess
!!$       call TRACEBACK('wrong idProcess!')
!!$    end select
!!$
!!$    MSTP(32)= 5 ! Q2 definition: -t
!!$
!!$    CKIN( 1) = 0.1
!!$    CKIN(39) = 0.1  ! W2_min in DIS (loose cut)
!!$    CKIN( 6) = 0.0  ! singular if masses below this cut
!!$    ! adjust the pThat cut off according W:
!!$    if (2*CKIN(5).ge.Wfree) CKIN(5) = Wfree/2-1e-3
!!$
!!$!    PARP(91)=0.00 ! width intrinsic kT
!!$!    MSTP(61)=0    ! master: ISR (QCD/QED)
!!$!    MSTP(71)=0    ! master: FSR (QCD/QED)
!!$
!!$    if (doMassless) then
!!$       PMAS(1,1)=0.010 ! Mass of d quark
!!$       PMAS(2,1)=0.010 ! Mass of u quark
!!$       PMAS(150,1)=0.020 ! Mass of diquark
!!$       PMAS(152,1)=0.020 ! Mass of diquark
!!$       PMAS(153,1)=0.020 ! Mass of diquark
!!$       PMAS(156,1)=0.020 ! Mass of diquark
!!$    endif
!!$
!!$    call Init_VM_Mass(Wfree,eNev%nucleon%pos)
!!$
!!$    if (eNev%nucleon_free%charge.gt.0) then
!!$       call PYNAME(2212,cTarget)          ! Target is proton
!!$    else
!!$       call PYNAME(2112,cTarget)          ! Target is neutron
!!$    endif
!!$
!!$    !...set incoming vectors (free kinematics)
!!$
!!$    pB = eNev%lepton_in%mom
!!$    pT = eNev%nucleon_free%mom
!!$
!!$    call eNeV_GetLeptonCM(eNev, betacm,phi,theta,phiLepton)
!!$!    theta = 0.0
!!$
!!$    call lorentz(betacm,pB)
!!$    call lorentz(betacm,pT)
!!$
!!$    call rotateZY(theta,phi,pB(1:3),pB(1:3))
!!$    call rotateZY(theta,phi,pT(1:3),pT(1:3))
!!$
!!$    call rotateZY(0.0,phiLepton,pB(1:3),pB(1:3))
!!$    call rotateZY(0.0,phiLepton,pT(1:3),pT(1:3))
!!$
!!$    pKprime = eNev%lepton_out%mom
!!$    pQ = eNev%boson%mom
!!$
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pK     ',pB
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pT     ',pT
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pKprime',pKprime
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pQ     ',pQ
!!$
!!$    call lorentz(betacm,pKprime)
!!$    call lorentz(betacm,pQ)
!!$
!!$    call rotateZY(theta,phi,pKprime(1:3),pKprime(1:3))
!!$    call rotateZY(theta,phi,pQ(1:3),pQ(1:3))
!!$
!!$    call rotateZY(0.0,phiLepton,pKprime(1:3),pKprime(1:3))
!!$    call rotateZY(0.0,phiLepton,pQ(1:3),pQ(1:3))
!!$
!!$!    write(*,*)
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pK     ',pB
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pT     ',pT
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pKprime',pKprime
!!$!    write(*,'(A," = ",1P,4e13.5)') 'pQ     ',pQ
!!$
!!$
!!$    P(1,1:3) = pB(1:3)
!!$    P(1,4)   = pB(0)
!!$    P(1,5)   = eNev%lepton_in%mass
!!$
!!$    P(2,1:3) = pT(1:3)
!!$    P(2,4)   = pT(0)
!!$    P(2,5)   = eNev%nucleon_free%mass
!!$
!!$
!!$    call CalcXY(eNev%lepton_in%mom,eNev%lepton_out%mom,&
!!$         & eNev%nucleon_free%mom, x_wish, y_wish, flagXY)
!!$    if (.not.flagXY) call TRACEBACK('Error in CalcXY')
!!$
!!$!    write(*,*) 'wish: x,y:',x_wish,y_wish
!!$
!!$
!!$    !...Initialize Pythia
!!$
!!$    MSTP(111) = 0          ! master switch fragmentation/decay
!!$    call PYINIT('5MOM', cNu(eNev%idFamily), cTarget, eNev%lepton_in%mom(0))
!!$
!!$    DistMin = 99.9 ! DUMMY
!!$
!!$    iEv = 0
!!$    do
!!$       iEv = iEv+1
!!$
!!$       if (iEv.ge.iEv_max) then
!!$          write(*,*) 'WARNING: iEv > ',iEv_max,' !'
!!$          return ! ---> FAILURE
!!$       end if
!!$
!!$
!!$       if (useJetSetVec) call GetJetsetVecINIT
!!$       call PYEVNT
!!$       if (MINT(51).eq.2) then
!!$          cycle ! retry
!!$       endif
!!$
!!$!       call PYLIST(2)
!!$
!!$       call MarkLepton(eNev,iL_lep)
!!$
!!$!       write(*,*) 'iL_lep=',iL_lep
!!$
!!$       pKprime = (/P(iL_lep,4),P(iL_lep,1:3)/)
!!$
!!$       call CalcXY((/P(1,4),P(1,1:3)/), pKprime, &
!!$            & (/P(2,4),P(2,1:3)/),x_try, y_try, flagXY)
!!$       if (.not.flagXY) call TRACEBACK('Error in CalcXY (try)')
!!$
!!$       Dist = (x_try-x_wish)**2+(y_try-y_wish)**2
!!$
!!$       if (Dist.lt.DistMin) then
!!$
!!$          ! now do the fragmentation:
!!$
!!$!          call PYLIST(2)
!!$
!!$          MSTP(111) = 1          ! master switch fragmentation/decay
!!$          call PYEXEC
!!$          MSTP(111) = 0          ! master switch fragmentation/decay
!!$
!!$          if (MSTU(24).ne.0) cycle ! retry
!!$          if (MINT(51).eq.2) cycle ! retry
!!$
!!$          if (useJetSetVec) then
!!$             call GetJetsetVec(.TRUE.)
!!$             call GetJetsetVecCheckT(-1d-5)
!!$             call GetJetsetVecPYEDIT
!!$          end if
!!$
!!$          call GetLeading_PY         ! find leading particles
!!$
!!$          ! check, whether it was a "cluster->1" deacy. Then the outgoing
!!$          ! lepton was modified by the fragmentation
!!$
!!$          do i=N,1,-1
!!$             if (K(i,1).ge.10) cycle ! no final particle
!!$             if ((K(i,4).eq.3).and.(K(i,5).eq.1)) stop
!!$             if ((K(i,4).eq.3).and.(K(i,5).eq.1)) cycle
!!$          end do
!!$
!!$!          call PYLIST(2)
!!$
!!$          DistMin = Dist
!!$!          write(*,*) 'x,y:',x_try,y_try
!!$!          write(*,*) 'DistMin = ',sqrt(DistMin),iEv
!!$
!!$
!!$          if (DistMin.lt.DistMinCut) exit
!!$!          stop
!!$
!!$       endif
!!$
!!$
!!$    enddo
!!$
!!$!    write(*,*) 'DistMin = ',sqrt(DistMin)
!!$
!!$    ! now we store the last event in the particle vector:
!!$    cross=XSEC(0,3) ! cross section of all subprocesses in mb
!!$    ! we have to scale the cross section by the number of events
!!$    ! that were necessary to get to this event !!!!
!!$
!!$    cross = cross / (iEv*DistMinCut*4./pi)
!!$
!!$    !...Copy Particles to ouput-vector
!!$
!!$!    write(*,*) '####'
!!$!    call PYLIST(2)
!!$
!!$    phi2 = atan2(pKprime(2),pKprime(1))
!!$
!!$    call PYROBO(1,N, 0.0,phiLepton-phi2, 0.0,0.0,0.0)
!!$    call PYROBO(1,N, theta,phi, betacm(1),betacm(2),betacm(3))
!!$
!!$    if (useJetSetVec) then
!!$       call GetJetsetVecPYROBO(0.0,  phiLepton, 0.0,0.0,0.0)
!!$       call GetJetsetVecPYROBO(theta,phi, betacm(1),betacm(2),betacm(3))
!!$    end if
!!$
!!$!    call PYLIST(2)
!!$
!!$    call PYEDIT(1)            ! clean up event list
!!$
!!$    call SetVectorFromPYJETS(outPart, eNev%Q2)
!!$    call ResetPyVars
!!$    flagOK = .TRUE.
!!$
!!$    write(*,*) iEv
!!$
!!$    return
!!$    stop
!!$
!!$
!!$  end subroutine DoColl_nuN_Py_fixMC


!!$  !****************************************************************************
!!$  !****s* Coll_nuN/DoColl_nuN_Py_Free
!!$  ! NAME
!!$  ! subroutine DoColl_nuN_Py_Free(eNev,outPart,flagOK, doMassless,cross)
!!$  !
!!$  ! PURPOSE
!!$  ! generate a high energy neutrino event with PYTHIA
!!$  !
!!$  ! Here the outgoing kinematics is free at input and is set in this routine
!!$  !
!!$  ! INPUTS
!!$  ! * type(electronNucleon_event) :: eNev     -- kinematics
!!$  ! * logical                     :: doMassless -- neglect masses of quarks and
!!$  !   diquarks
!!$  !
!!$  ! OUTPUT
!!$  ! * type(electronNucleon_event) :: eNev     -- kinematics
!!$  ! * type(particle),dimension(:) :: outPart  -- outgoing particles
!!$  ! * logical                     :: flagOK   -- .TRUE., if everything okay
!!$  ! * real                        :: cross    -- cross section
!!$  !
!!$  ! NOTES
!!$  ! We give as input the neutrino and the nucleon in medium and free nucleon
!!$  ! kinematics. Since we can not ask for any cuts on the outgoing
!!$  ! lepton/neutrino, we set the momentum of this particle as output
!!$  !****************************************************************************
!!$  subroutine DoColl_nuN_Py_Free(eNev,outPart,flagOK, doMassless, cross)
!!$
!!$    use GetLeading
!!$    use EventInfo_HiLep
!!$    use LorentzTrafo
!!$    use rotation
!!$    use PythiaSpecFunc, only: Init_VM_Mass
!!$    use eN_event
!!$    use minkowski, only: abs4,abs4sq
!!$    use CallStack
!!$
!!$    type(electronNucleon_event), intent(inout):: eNev
!!$    logical,                     intent(in)   :: doMassless
!!$    type(particle),dimension(:), intent(inout):: outPart
!!$    logical,                     intent(out)  :: flagOK
!!$    real,                        intent(out)  :: cross
!!$
!!$    ! common blocks
!!$
!!$    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
!!$    integer KCHG
!!$    double precision PMAS,PARF,VCKM
!!$    SAVE /PYDAT2/
!!$
!!$    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
!!$    integer N,NPAD,K
!!$    double precision P,V
!!$    SAVE /PYJETS/
!!$
!!$    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
!!$    integer MSTP,MSTI
!!$    double precision PARP,PARI
!!$    SAVE /PYPARS/
!!$
!!$    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
!!$    integer MSEL,MSELPD,MSUB,KFIN
!!$    double precision CKIN
!!$    SAVE /PYSUBS/
!!$
!!$    COMMON/PYINT1/MINT(400),VINT(400)
!!$    integer MINT
!!$    double precision VINT
!!$    SAVE /PYINT1/
!!$
!!$    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
!!$    integer NGENPD,NGEN
!!$    double precision XSEC
!!$    SAVE /PYINT5/
!!$
!!$
!!$    ! internal variables:
!!$
!!$    real :: Wfree
!!$    character*20 :: cTarget,cNu(3) = (/'nu_e  ','nu_mu ','nu_tau'/)
!!$    real, dimension(0:3) :: pB,pT
!!$    real :: phi,theta
!!$    real, dimension(1:3) :: beta
!!$    integer :: iEv, nEV=10, iTry,i
!!$    integer :: iL_lep
!!$
!!$    !...set default output
!!$
!!$    outPart%ID = 0
!!$    flagOK = .FALSE.
!!$    cross = 0.0
!!$
!!$    !...save original PYTHIA values
!!$    call StorePyVars
!!$
!!$    !...set up PYTHIA
!!$    call SetSomeDefaults_PY
!!$
!!$    Wfree = eNev%W_free
!!$
!!$    MSEL = 0       ! selction of subprocesses
!!$    MSUB     = 0   ! all subprocesses off
!!$    MSUB(10) = 1   ! selction of subprocesses
!!$
!!$    select case (eNev%idProcess)
!!$    case (-1,1) !=== EM
!!$       write(*,*) 'not prepared for EM!'
!!$       call TRACEBACK('wrong idProcess!')
!!$    case (-2,2) !=== CC
!!$       MSTP(21) = 5 ! charged current only
!!$    case (-3,3) !=== NC
!!$       MSTP(21) = 3 ! neutral current, Z0 only
!!$    case default
!!$       write(*,*) 'idProcess=',eNev%idProcess
!!$       call TRACEBACK('wrong idProcess!')
!!$    end select
!!$
!!$    MSTP(32)= 5 ! Q2 definition: -t
!!$
!!$    CKIN( 1) = 0.1
!!$    CKIN(39) = 0.1  ! W2_min in DIS (loose cut)
!!$    CKIN( 6) = 0.0  ! singular if masses below this cut
!!$    ! adjust the pThat cut off according W:
!!$    if (2*CKIN(5).ge.Wfree) CKIN(5) = Wfree/2-1e-3
!!$
!!$!    PARP(91)=0.00 ! width intrinsic kT
!!$!    MSTP(61)=0    ! master: ISR (QCD/QED)
!!$!    MSTP(71)=0    ! master: FSR (QCD/QED)
!!$
!!$    if (doMassless) then
!!$       PMAS(1,1)=0.010 ! Mass of d quark
!!$       PMAS(2,1)=0.010 ! Mass of u quark
!!$       PMAS(150,1)=0.020 ! Mass of diquark
!!$       PMAS(152,1)=0.020 ! Mass of diquark
!!$       PMAS(153,1)=0.020 ! Mass of diquark
!!$       PMAS(156,1)=0.020 ! Mass of diquark
!!$    end if
!!$
!!$    call Init_VM_Mass(Wfree,eNev%nucleon%pos)
!!$
!!$    if (eNev%nucleon_free%charge.gt.0) then
!!$       call PYNAME(2212,cTarget)          ! Target is proton
!!$    else
!!$       call PYNAME(2112,cTarget)          ! Target is neutron
!!$    end if
!!$
!!$    !...set incoming vectors
!!$
!!$    ! please note: we need the electron vectors in the x-z plane,
!!$    ! therefore we have to append a additional rotation with phi
!!$    ! around the z-axis in order to set the y-component to 0 !!!
!!$
!!$    ! BlaBlaBla!!! electron_out=electron_in at the moment, therefore
!!$    ! all angles arbitrary!!!
!!$
!!$    pB = eNev%lepton_in%mom
!!$    pT = eNev%nucleon_free%mom
!!$
!!$    phi = atan2(eNev%pcm(2),eNeV%pcm(1))
!!$    theta = atan2(sqrt(eNev%pcm(1)**2+eNev%pcm(2)**2),eNev%pcm(3))
!!$    beta = eNev%betacm
!!$
!!$!    theta = 0.0
!!$
!!$    call lorentz(beta,pB)
!!$    call lorentz(beta,pT)
!!$
!!$    pB(1:3) = rotateZY (theta, phi, pB(1:3))
!!$    pT(1:3) = rotateZY (theta, phi, pT(1:3))
!!$
!!$    pB(1:3) = rotateZY (0.0, eNev%phiLepton, pB(1:3))
!!$    pT(1:3) = rotateZY (0.0, eNev%phiLepton, pT(1:3))
!!$
!!$    P(1,1:3) = pB(1:3)
!!$    P(2,1:3) = pT(1:3)
!!$
!!$    !...Initialize Pythia
!!$
!!$    MSTP(111) = 0          ! master switch fragmentation/decay
!!$    call PYINIT('3MOM', cNu(eNev%idFamily), cTarget, eNev%lepton_in%mom(0))
!!$
!!$    do iEv=1,nEv
!!$       if (useJetSetVec) call GetJetsetVecINIT
!!$       call PYEVNT
!!$    end do
!!$
!!$    cross=XSEC(0,3) ! cross section of all subprocesses in mb
!!$
!!$    MSTP(111) = 1          ! master switch fragmentation/decay
!!$    iTry = 0
!!$    TryLoop: do
!!$       iTry = iTry+1
!!$       if (iTry.ge.100) then
!!$          write(*,*) 'DoColl_nuN: itry=',iTry
!!$          N = 0
!!$          call ResetPyVars
!!$          return ! -> FAILURE
!!$       end if
!!$
!!$       call PYINIT('3MOM', cNu(eNev%idFamily), cTarget, eNev%lepton_in%mom(0))
!!$       if (useJetSetVec) call GetJetsetVecINIT
!!$       call PYEVNT
!!$
!!$       !    call PYLIST(2)
!!$
!!$       !...check success:
!!$
!!$       if (MINT(51).eq.2) then
!!$          cycle TryLoop ! retry
!!$       end if
!!$
!!$       call MarkLepton(eNev,iL_lep)
!!$!       write(*,*) 'iL_lep=',iL_lep
!!$
!!$!       call PYLIST(2)
!!$
!!$       if (useJetSetVec) then
!!$          call GetJetsetVec(.TRUE.)
!!$          call GetJetsetVecCheckT(-1d-5)
!!$          call GetJetsetVecPYEDIT
!!$       end if
!!$
!!$       call GetLeading_PY         ! find leading particles
!!$
!!$       ! check, whether it was a "cluster->1" deacy. Then the outgoing
!!$       ! lepton was modified by the fragmentation
!!$
!!$       do i=N,1,-1
!!$          if (K(i,1).ge.10) cycle ! no final particle
!!$!          if ((K(i,4).eq.3).and.(K(i,5).eq.1)) stop
!!$          if ((K(i,4).eq.3).and.(K(i,5).eq.1)) cycle TryLoop
!!$       end do
!!$
!!$       exit ! success, leave the Try-Loop
!!$
!!$    end do TryLoop
!!$
!!$
!!$    if (nEV.le.1) cross=XSEC(0,3) ! cross section of all subprocesses in mb
!!$
!!$    !...restore original PYTHIA values
!!$
!!$    call ResetPyVars
!!$
!!$    call PYROBO(1,N, 0.0,eNev%phiLepton, 0.0,0.0,0.0)
!!$    call PYROBO(1,N, theta,phi, eNev%betacm(1),eNev%betacm(2),eNev%betacm(3))
!!$
!!$    if (useJetSetVec) then
!!$       call GetJetsetVecPYROBO(0.0,  eNev%phiLepton, 0.0,0.0,0.0)
!!$       call GetJetsetVecPYROBO(theta,phi, eNev%betacm(1),eNev%betacm(2),eNev%betacm(3))
!!$    end if
!!$
!!$!    call PYLIST(2)
!!$
!!$    !...Set initial event kinematics
!!$
!!$    eNev%lepton_out%mom(1:3) = P(iL_lep,1:3)
!!$    eNev%lepton_out%mom(0)   = P(iL_lep,4)
!!$    eNev%boson%mom = eNev%lepton_in%mom-eNev%lepton_out%mom
!!$    eNev%Q2     = -abs4Sq(eNev%boson%mom)
!!$    eNev%W            = abs4(eNev%boson%mom+eNev%nucleon%mom)
!!$    eNev%W_free       = abs4(eNev%boson%mom+eNev%nucleon_free%mom)
!!$
!!$!    call write_electronNucleon_event(eNev,.FALSE.,.FALSE.)
!!$
!!$    call PYEDIT(1)            ! clean up event list
!!$
!!$    !...Copy Particles to ouput-vector
!!$
!!$    call SetVectorFromPYJETS(outPart, eNev%Q2)
!!$
!!$    flagOK = .TRUE.
!!$
!!$  end subroutine DoColl_nuN_Py_Free

  !****************************************************************************
  !****is* Coll_nuN/MarkLepton
  ! NAME
  ! subroutine MarkLepton(eNev, iL_lep)
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- event information
  !
  ! OUTPUT
  ! * integer :: iL_lep -- line number of lepton
  ! * Pythia array entry K(iL_lep,1) changed
  !
  ! PURPOSE
  ! find the outgoing lepton and mark it a documentation line
  !****************************************************************************
  subroutine MarkLepton(eNev, iL_lep)

    type(electronNucleon_event), intent(in):: eNev
    integer, intent(out) :: iL_lep

    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    integer N,NPAD,K
    double precision P,V
    SAVE /PYJETS/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    integer KF0,i

    KF0 = 10+eNev%idFamily*2    ! KF code of the outgoing lepton
    if (abs(eNev%idProcess).eq.1) KF0=KF0-1 !=== el-m
    if (abs(eNev%idProcess).eq.2) KF0=KF0-1 !=== CC
    if (eNev%idProcess.lt.0) KF0=-KF0

    !      call PYGIVE('MSTI(4)=')

    do i=MSTI(4)+1,N          ! skip documentation
       if (K(i,2).eq.KF0) then
          if (K(i,3).gt.MSTI(4)) then
             write(*,*) 'Ooops, what an event. Stop!'
             call PYLIST(2)
             call PYGIVE('MSTI(4)=')
             stop
          end if
          K(i,1) = 21
          iL_lep = i
          return ! --> success
       end if
    end do
  end subroutine MarkLepton


  !****************************************************************************
  !****is* Coll_nuN/StorePyVars
  ! NAME
  ! subroutine StorePyVars
  !
  ! PURPOSE
  ! store Pytha variables, which are changed during the routine
  !****************************************************************************
  subroutine StorePyVars

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    integer :: SetMaxMSTI52 ! prototype

    MSEL0 = MSEL
    MSTP11 = MSTP(11)
    MSTP21 = MSTP(21)
    MSTP23 = MSTP(23)
    MSTP32 = MSTP(32)
    MSTP61 = MSTP(61)
    MSTP71 = MSTP(71)
    CKIN01 = CKIN( 1)
    CKIN05 = CKIN( 5)
    CKIN06 = CKIN( 6)
    CKIN39 = CKIN(39)
    PARP91 = PARP(91)

    PMAS1   = PMAS(1,1)
    PMAS2   = PMAS(2,1)
    PMAS150 = PMAS(150,1)
    PMAS152 = PMAS(152,1)
    PMAS153 = PMAS(153,1)
    PMAS156 = PMAS(156,1)

    maxMSTI52 = SetMaxMSTI52(5)

  end subroutine StorePyVars

  !****************************************************************************
  !****is* Coll_nuN/ResetPyVars
  ! NAME
  ! subroutine ResetPyVars
  !
  ! PURPOSE
  ! reset Pytha variables to its stored values
  !****************************************************************************
  subroutine ResetPyVars

    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    integer KCHG
    double precision PMAS,PARF,VCKM
    SAVE /PYDAT2/

    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    integer MSTP,MSTI
    double precision PARP,PARI
    SAVE /PYPARS/

    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    integer MSEL,MSELPD,MSUB,KFIN
    double precision CKIN
    SAVE /PYSUBS/

    integer :: dummy, SetMaxMSTI52

    MSEL = MSEL0
    MSTP(11) = MSTP11
    MSTP(21) = MSTP21
    MSTP(23) = MSTP23
    MSTP(32) = MSTP32
    MSTP(61) = MSTP61
    MSTP(71) = MSTP71
    CKIN( 1) = CKIN01
    CKIN( 5) = CKIN05
    CKIN( 6) = CKIN06
    CKIN(39) = CKIN39
    PARP(91) = PARP91

    PMAS(1,1) = PMAS1
    PMAS(2,1) = PMAS2
    PMAS(150,1) = PMAS150
    PMAS(152,1) = PMAS152
    PMAS(153,1) = PMAS153
    PMAS(156,1) = PMAS156

    dummy = SetMaxMSTI52(maxMSTI52)

  end subroutine ResetPyVars

  !****************************************************************************
  !****is* Coll_nuN/CalcXY
  ! NAME
  ! subroutine CalcXY(pK,pKprime,pP, x,y, flagOK)
  !
  ! PURPOSE
  ! calculate x and y from given momenta
  !****************************************************************************
  subroutine CalcXY(pK,pKprime,pP, x,y, flagOK)
    use minkowski, only: SP

    real, dimension(0:3),intent(in) :: pP,pK,pKprime
    real, intent(out)               :: x,y
    logical, intent(out)            :: flagOK

    real, dimension(0:3)            :: pQ
    real :: pPpQ, pPpK, pQpQ

    flagOK = .false.

    pQ = pK-pKprime

!!$    write(*,'(A," = ",1P,4e13.5)') 'pK     ',pK
!!$    write(*,'(A," = ",1P,4e13.5)') 'pP     ',pP
!!$    write(*,'(A," = ",1P,4e13.5)') 'pKprime',pKprime
!!$    write(*,'(A," = ",1P,4e13.5)') 'pQ     ',pQ

    pPpQ = SP(pP,pQ)
    pPpK = SP(pP,pK)
    pQpQ = SP(pQ,pQ)

!!$    write(*,'(A," = ",1P,e13.5)') 'pPpQ',pPpQ
!!$    write(*,'(A," = ",1P,e13.5)') 'pPpK',pPpK
!!$    write(*,'(A," = ",1P,e13.5)') 'pQpQ',pQpQ

    if (pPpK.eq.0.0) return
    if (pPpQ.eq.0.0) return

    y = pPpQ/pPpK
    x = -pQpQ/(2*pPpQ)

    flagOK = .true.

  end subroutine CalcXY

  !****************************************************************************
  !****f* Coll_nuN/AnaEstimate
  ! NAME
  ! real function CalcAnaEstimate(eNev)
  ! PURPOSE
  ! Calculate the cross sections according eqs. (C.37)-(C.39) of Phys.Rep.
  ! OUTPUT
  ! * function value = dsigma/dE'dcost in mb/GeV
  !****************************************************************************
  real function AnaEstimate(eNev) result(sigma)

    use CallStack, only: TRACEBACK
    use eN_event
!    use minkowski, only : abs4Sq
    use constants, only: pi, GF, alphaQED, hbarc, mN

    type(electronNucleon_event), intent(inout):: eNev

    real:: x,y
    logical :: flagXY
    real, dimension(-25:25) :: xq
    real :: srts2,sigma0
    real :: nu,Q2,W,Wfree,eps,fT
    real, parameter :: hc2 = 10*hbarc**2 ! = 0.389 mb GeV^2

    call CalcXY(eNev%lepton_in%mom,eNev%lepton_out%mom,&
         & eNev%nucleon_free%mom, x, y, flagXY)
    if (.not.flagXY) call TRACEBACK('Error in CalcXY')

    call eNeV_GetKinV(eNev, nu,Q2,W,Wfree,eps,fT)

    srts2 = 2*mN*eNev%lepton_in%mom(0)
!    srts2 = abs4sq(eNev%lepton_in%mom+eNev%nucleon_free%mom)

    call pypdfl(2212, x,Q2,xq)

    select case ( eNev%idProcess )
    case (-1,1) !=== EM
       sigma0 = 2*pi*alphaQED**2*srts2*hc2/Q2**2*(1+(1-y)**2) ! in mb !!!
!       sigma0 = 2*pi*alphaQED**2*srts2*hc2/Q2**2*y**2/(1-eps) ! in mb !!!
    case (-2,2) !=== CC
       sigma0 = GF**2*srts2*hc2/(2*pi)  ! in mb !!!
    case (-3,3) !=== NC
       sigma0 = 0. ! === DIS not yet implemented here === !
    case default
       write(*,*) 'idProcess=',eNev%idProcess
       call TRACEBACK('wrong idProcess!')
    end select

    select case ( eNev%idProcess )
    case (-1,1) !=== EM
       if (eNev%nucleon%charge .gt. 0) then ! *** proton ***
          sigma = sigma0 * ( 1 * (xq(1) + xq(-1)) &
               &            +4 * (xq(2) + xq(-2)) &
               &            +1 * (xq(3) + xq(-3)) )/9.
       else
          sigma = sigma0 * ( 4 * (xq(1) + xq(-1)) &
               &            +1 * (xq(2) + xq(-2)) &
               &            +4 * (xq(3) + xq(-3)) )/9.
       end if

    case (-2) !=== CC, antiparticle

       if (eNev%nucleon%charge .gt. 0) then ! *** proton ***
          sigma = 2*sigma0 * (xq(-1)+ (1-y)**2*xq(2)) ! x(dbar + (1-y)**2u)
       else
          sigma = 2*sigma0 * (xq(-2)+ (1-y)**2*xq(1)) ! x(ubar + (1-y)**2d)
       end if

    case (2) !=== CC, particle

       if (eNev%nucleon%charge .gt. 0) then ! *** proton ***
          sigma = 2*sigma0 * (xq(1)+ (1-y)**2*xq(-2)) ! x(d + (1-y)**2ubar)
       else
          sigma = 2*sigma0 * (xq(2)+ (1-y)**2*xq(-1)) ! x(u + (1-y)**2dbar)
       end if

    end select

    ! Transform from dsigma/dxdy into dsigma/dcostdEprime:
    sigma = sigma * eNev%lepton_out%mom(0)/(mN*eNev%boson%mom(0))


  end function AnaEstimate

  !****************************************************************************
  !****f* Coll_nuN/AnaEstimatePythia
  ! NAME
  ! real function CalcAnaEstimatePythia(eNev)
  ! PURPOSE
  ! Calculate the cross sections according Pythia for the em case
  ! OUTPUT
  ! * function value = dsigma/dE'dcost in mb/GeV
  !****************************************************************************
  real function AnaEstimatePythia(eNev) result(sigma)

    use CallStack, only: TRACEBACK
    use eN_event
    use constants, only: pi, alphaQED, hbarc

    type(electronNucleon_event), intent(inout):: eNev

    real:: x,y
    logical :: flagXY
    real, dimension(-25:25) :: xq
    real :: nu,Q2,W,Wfree,eps,fT
    real, parameter :: hc2 = 10*hbarc**2 ! = 0.389 mb GeV^2
    real :: r2GA

    call CalcXY(eNev%lepton_in%mom,eNev%lepton_out%mom,&
         & eNev%nucleon_free%mom, x, y, flagXY)
    if (.not.flagXY) call TRACEBACK('Error in CalcXY')

    call eNeV_GetKinV(eNev, nu,Q2,W,Wfree,eps,fT)

!    eps = (1-y)/(1-y+0.5*y**2)

    sigma = eNev%lepton_out%mom(0)*alphaQED*hc2/pi &
         & * y*(1-x)/Q2 * 1/(1-eps) & !         & * (1-x)/Q2 * (1+(1-y)**2)/y &
         & * 4*pi*alphaQED/(Q2*(1-x))

    ! respect longitudinal contribution:

    if (Q2.gt.0.35) then
       r2GA = 0.0635*(1+12*Q2*0.125**2/((Q2+1)* &
            & (0.125**2+X**2)))/log(Q2/0.04) &
            & + 0.5747/Q2 - 0.3534/(Q2**2+0.09)
    else
       r2GA = 0.0635*(1+12*0.35*0.125**2/((0.35+1)* &
            & (0.125**2+X**2)))/log(0.35/0.04) &
            & + 0.5747/0.35 - 0.3534/(0.35**2+0.09)
    end if

!    r2GA = Q2/nu**2
    sigma = sigma*(1+eps*r2GA)


    call pypdfl(2212, x,Q2,xq)
    if (eNev%nucleon%charge .gt. 0) then ! *** proton ***
       sigma = sigma * ( 1 * (xq(1) + xq(-1)) &
            &           +4 * (xq(2) + xq(-2)) &
            &           +1 * (xq(3) + xq(-3)) )/9.
    else
       sigma = sigma * ( 4 * (xq(1) + xq(-1)) &
            &           +1 * (xq(2) + xq(-2)) &
            &           +4 * (xq(3) + xq(-3)) )/9.
    end if

    sigma = sigma * Q2**2/(Q2+0.776**2)**2 ! additional form factor

  end function AnaEstimatePythia


!!$   subroutine GetPythiaXY(x,y)
!!$
!!$     real, intent(out) :: x,y
!!$
!!$     COMMON/PYINT1/MINT(400),VINT(400)
!!$     integer MINT
!!$     double precision VINT
!!$     SAVE /PYINT1/
!!$
!!$     x = VINT(21)
!!$     y = VINT(22)
!!$
!!$   end subroutine GetPythiaXY

end module Coll_nuN
