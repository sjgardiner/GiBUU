!******************************************************************************
!****m* /initNeutrino
! NAME
! module initNeutrino
!
! PURPOSE
! This module is the main module for neutrino-induced reactions. It provides
! the calculation of various inclusive neutrino cross sections, depending
! only on the exchanged gauge boson momentum.
! It also sets the produced outgoing particles into the pertubativeVector.
!******************************************************************************
module initNeutrino

  use neutrino_IDTable
  use Electron_origin, only: origin_singlePi, origin_doublePi, &
       origin_DIS, origin_2p2hQE, origin_2p2hDelta
  use Minkowski, only : abs3
  use CallStack

  implicit none

  private

  public :: init_neutrino
  public :: getFirstEventRange
  public :: getNeutrinoInfo
  public :: nuXsectionMode, process_ID, flavor_ID
  public :: nuExp
  public :: includeQE, includeDELTA, includeRES, include1pi, includeDIS, &
       include2p2hQE, include2p2hDelta, include2pi
  public :: neutrinoInit_getRealRun
  public :: cleanup

  public :: get_init_namelist
  public :: get_runtime_vars
  public :: max_finalstate_ID, max_Hist,includeHist,K2Hist,   &
       numberOfExperiments, OscLength,isOsc
  public :: cost_min,cost_max,delta_cost, Elept_min,Elept_max,delta_Elept
  public :: pL_min,pL_max,pT_min,pT_max, delta_pL,delta_pT



  !****************************************************************************
  !****g* initNeutrino/process_ID
  ! SOURCE
  integer, save :: process_ID = 2
  ! PURPOSE
  ! Determine the process (cf. module leptonicID):
  ! * 1 = EM
  ! * 2 = CC
  ! * 3 = NC
  ! * -1 = antiEM
  ! * -2 = antiCC
  ! * -3 = antiNC
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/flavor_ID
  ! SOURCE
  integer, save :: flavor_ID = 2
  ! PURPOSE
  ! Determine the lepton flavor:
  ! * 1 = electron
  ! * 2 = muon
  ! * 3 = tau
  !****************************************************************************

  character*(*), dimension(3), parameter :: sProcess = (/"EM","CC","NC"/)
  character*(*), dimension(3), parameter :: sFamily  = (/"e  ","mu ","tau"/)

  !****************************************************************************
  !****g* initNeutrino/nuXsectionMode
  ! SOURCE
  integer, save :: nuXsectionMode = 0
  ! PURPOSE
  ! To choose which kind of Xsection is calculated. All values set in
  ! module neutrino_IDTable.f90
  !
  ! possible values:
  ! * 0 = integratedSigma: required input: enu
  ! * 1 = dSigmadCosThetadElepton: required input: enu, costheta, elepton
  ! * 2 = dSigmadQ2dElepton: required input: enu, Q2, elepton
  ! * 4 = dSigmadCosTheta: required input: enu, costheta
  ! * 5 = dSigmadElepton: required input: enu, elepton
  ! * 6 = dSigmaMC: required input: enu
  ! * 7 = dSigmaMC_dW: required input: enu, W
  ! * 3 = dSigmaMC_dQ2: required input: enu, Q2
  !
  ! calculation for specific experiments taking into account the flux
  ! (choose your favorite experiment with flag nuExp):
  ! * 10 = EXP_dSigmadEnu
  ! * 11 = EXP_dSigmadCosThetadElepton
  ! * 12 = EXP_dSigmadQ2dElepton
  ! * 14 = EXP_dSigmadCosTheta
  ! * 15 = EXP_dSigmadElepton
  ! * 16 = EXP_dSigmaMC
  ! * 17 = EXP_dSigmaMC_dW
  ! * 13 = EXP_dSigmaMC_dQ2
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/nuExp
  ! SOURCE
  integer, save :: nuExp = 0
  ! PURPOSE
  ! * 0 = no specific experiment
  ! * 1 = MiniBooNE neutrino flux (in neutrino mode = positive polarity)
  ! * 2 = ANL
  ! * 3 = K2K
  ! * 4 = BNL
  ! * 5 = MiniBooNE antienutrino flux (in antineutrino mode = negative polarity)
  ! * 6 = MINOS muon-neutrino  in neutrino mode
  ! * 7 = MINOS muon-antineutrino  in neutrino mode
  ! * 8 = NOVA neutrino (medium energy NuMI, 14 mrad off-axis), FD
  ! * 9 = T2K neutrino off-axix 2.5 degrees ( at ND280 detector )
  ! * 10 = uniform distribution from Eflux_min to Eflux_max
  !       (see namelist nl_neutrino_energyFlux in the module expNeutrinoFluxes)
  ! * 11 = MINOS muon-neutrino  in antineutrino mode
  ! * 12 = MINOS muon-antineutrino  in antineutrino mode
  ! * 13 = MINERvA muon neutrino, old flux
  ! * 14 = MINERvA muon antineutrino, old flux
  ! * 15 = LBNF/DUNE in neutrino mode
  ! * 16 = LBNF/DUNE in antineutrino mode
  ! * 17 = LBNO neutrino in neutrino mode
  ! * 18 = NOMAD
  ! * 19 = BNB nue          BNB= Booster Neutrino Beam
  ! * 20 = BNB nuebar
  ! * 21 = BNB numu
  ! * 22 = BNB numubar
  ! * 23 = NOvA ND
  ! * 24 = T2K on axis
  ! * 25 = MINERvA, 2016 flux
  ! * 99 = user provided input file
  !****************************************************************************

  integer, parameter :: numberOfExperiments=25



  character*(*), dimension(0:numberOfExperiments), parameter ::  sExp  = (/ &
        "no specific experiment ", &
        "MiniBooNE nu           ", "ANL                    ", &
        "K2K                    ", &
        "BNL                    ", "MiniBooNE barnu        ", &
        "MINOS nu numode        ", "MINOS barnu numode     ", &
        "NOvA FD                ", &
        "T2K OffAxis 2.5deg     ", "uniform distribution   ", &
        "MINOS nu barnumode     ", "MINOS barnu barnumode  ", &
        "MINERvA nu numode      ", "MINERvA barnu barnumode", &
        "LBNF-DUNE nu           ", "LBNF-DUNE barnu        ", &
        "LBNO nu numode         ", "NOMAD                  ", &
        "BNB nue                ", "BNB nuebar             ", &
        "BNB numu               ", "BNB numubar            ", &
        "NOvA ND                ", "T2K on axis            ", &
        "MINERvA, 2016 flux     " /)


  real, dimension(0:numberOfExperiments), parameter :: OscLength = &
  (/ 0., 0.541, 0., 250., 0., 0.541, 735., 735., 810., 295., 0., 735., 735., &
     0.5, 0.5, 1300., 1300., 2300., 0.6262,0.,0.,0.,0.,0.,0.,0. /)
  ! oscillation length for various experiments in kilometers

  logical, dimension(0:numberOfExperiments), parameter:: Osc = &
  (/ .FALSE.,.FALSE.,.FALSE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,&
     .TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,&
     .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE. /)
  ! OSC is true for oscillation experiments, false otherwise
  !


  !****************************************************************************
  !****g* initNeutrino/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag = .false.
  ! PURPOSE
  ! To switch on debugging information
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeQE
  ! SOURCE
  logical, save :: includeQE = .true.
  ! PURPOSE
  ! include QE scattering
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeDELTA
  ! SOURCE
  logical, save :: includeDELTA = .true.
  ! PURPOSE
  ! include Delta excitation
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeRES
  ! SOURCE
  logical, save :: includeRES = .true.
  ! PURPOSE
  ! include excitation of higher resonances
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include1pi
  ! SOURCE
  logical, save :: include1pi = .false.
  ! PURPOSE
  ! include one-pion cross section
  ! see neutrinoXsection.f90 for details: there
  ! one might choose between different models and
  ! also whether it is taken as background or as
  ! total cross section
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include2pi
  ! SOURCE
  logical, save :: include2pi = .false.
  ! PURPOSE
  ! include 2 pion background channel
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/includeDIS
  ! SOURCE
  logical, save :: includeDIS = .false.
  ! PURPOSE
  ! include DIS contribution
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include2p2hQE
  ! SOURCE
  logical, save :: include2p2hQE = .false.
  ! PURPOSE
  ! include 2p2h QE contribution
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/include2p2hDelta
  ! SOURCE
  logical, save :: include2p2hDelta = .false.
  ! PURPOSE
  ! include 2p2h Delta contribution
  !****************************************************************************


  !****************************************************************************
  !****g* initNeutrino/realRun
  ! SOURCE
  logical, save :: realRun = .false.
  ! PURPOSE
  ! Do not initialize the final state particles as perturbative particles but
  ! as real ones.
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/printAbsorptionXS
  ! SOURCE
  logical, save :: printAbsorptionXS = .false.
  ! PURPOSE
  ! flag to produce output about inclusive (absorption) cross sections
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/printInclHist
  ! SOURCE
  logical, save :: printInclHist = .true.
  ! PURPOSE
  ! flag to produce additional output about inclusive cross sections
  !
  ! only checked, if printAbsorptionXS = T
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/storeNucleon
  ! SOURCE
  integer, save :: storeNucleon = 2
  ! PURPOSE
  ! indicate which kind of struck nucleon to save:
  ! * 1: free Nucleon (i.e. potential removed)
  ! * 2: bound nucleon
  !
  ! NOTES
  ! real check of energy and momentum conservation only possible with '2'
  !****************************************************************************

  real, save :: raiseVal


  !****************************************************************************
  !****g* initNeutrino/max_finalstate_ID
  ! SOURCE
  integer,parameter :: max_finalstate_ID=37
  ! * Parameter determines the reaction mechanism and kind of final states
  ! * Final states are numbered (often IP) by
  ! * 1: nucleon (QE)
  ! * 2-31: non-strange baryon resonance (as in IdTable)
  ! * 32: pi neutron-background  (e.g. nu + n -> mu + pi+ + n)
  ! * 33: pi proton-background   (e.g. nu + n -> mu + pi0 + p)
  ! * 34: DIS
  ! * 35: 2p2h QE
  ! * 36: 2p2h Delta
  ! * 37: two pion background
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/sigmacut
  ! SOURCE
  real, save :: sigmacut=10e-4
  ! in usual units (10^(-33) cm^2 for electrons, 10^(-38) cm^2 for neutrinos)
  ! PURPOSE
  ! events with a cross section smaller than this value are skipped.
  !****************************************************************************



  !****************************************************************************
  !****g* neutrinoAnalysis/cost_min
  ! SOURCE
  real, save :: cost_min= -1.0
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Minimal cos(theta) of outgoing leptons, used in 2D dsigma/dEdcos(theta)
  ! This cut affects *only* the outgoing lepton
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/cost_max
  ! SOURCE
  real, save :: cost_max= +1.0
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Maximal cos(theta) of outgoing leptons, used in 2D dsigma/dEdcos(theta)
  ! This cut affects *only* the outgoing lepton
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/delta_cost
  ! SOURCE
  real, save :: delta_cost = 0.1      ! for MB,T2K, for higher energies smaller
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! stepsize of cos(theta) of outgoing leptons, used in 2D dsigma/dEdcos(theta)
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/Elept_min
  ! SOURCE
  real, save :: Elept_min = 0.0
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! minimal energy of outgoing leptons, used in 2D dsigma/dEdcos(theta)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Elept_max
  ! SOURCE
  real, save :: Elept_max = 2.0     ! for MB,T2K, for higher energies larger
  ! PURPOSE
  ! if detailed_diff_output or printAbsorption are TRUE:
  ! maximal energy of outgoing leptons, used in 2D dsigma/dEdcos(theta)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/delta_Elept
  ! SOURCE
  real, save :: delta_Elept = 0.01   ! for MB,T2K
  ! PURPOSE
  ! if detailed_diff_output or printAbsorption are TRUE:
  ! stepsize of energy of outgoing leptons, used in 2D dsigma/dEdcos(theta)
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/pL_min
  ! SOURCE
  real, save :: pL_min = 0.0
  ! PURPOSE
  ! if detailed_diff or printAbsorption are TRUE:
  ! minimal long. momentum  of outgoing leptons, used in 2D dsigma/dpLdpT
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/pL_max
  ! SOURCE
  real, save :: pL_max = 20.0
  ! PURPOSE
  ! if detailed_diff_output or printAbsorption are TRUE:
  ! maximal long. momentum  of outgoing leptons, used in 2D dsigma/dpLdpT
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/delta_pL
  ! SOURCE
  real, save :: delta_pL = 0.25
  ! PURPOSE
  ! if detailed_diff_output or printAbsorption are TRUE:
  ! stepsize of long. momentum of outgoing leptons, used in 2D dsigma/dpLdpT
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/pT_min
  ! SOURCE
  real, save :: pT_min = 0.0
  ! PURPOSE
  ! if detailed_diff or printAbsorption are TRUE:
  ! minimal transv. momentum  of outgoing leptons, used in 2D dsigma/dpLdpT
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pT_max
  ! SOURCE
  real, save :: pT_max = 2.5
  ! PURPOSE
  ! if detailed_diff_output or printAbsorption are TRUE:
  ! maximal transv. momentum  of outgoing leptons, used in 2D dsigma/dpLdpT
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/delta_pT
  ! SOURCE
  real, save :: delta_pT = 0.1
  ! PURPOSE
  ! if detailed_diff_output or printAbsorption are TRUE:
  ! binwidth of transv. momentum of outgoing leptons, used in 2D dsigma/dpLdpT
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Q2_Max
  ! SOURCE
  real, save :: Q2_Max = 100.
  ! PURPOSE
  ! maximal value of Q2 in Q2-distribution
  !****************************************************************************

  ! !****************************************************************************
  ! !****g* neutrinoAnalysis/delta_Q2
  ! ! SOURCE
  ! real, save :: delta_Q2 = 0.1
  ! ! PURPOSE
  ! ! binwidth in Q2-distribution
  ! !****************************************************************************

  ! !****************************************************************************
  ! !****g* neutrinoAnalysis/numax
  ! ! SOURCE
  ! real, save :: numax = 1000
  ! ! PURPOSE
  ! ! max value in nu-distribution (energy-transfer distribution)
  ! !****************************************************************************

  ! !****************************************************************************
  ! !****g* neutrinoAnalysis/delta_nu
  ! ! SOURCE
  ! real, save :: delta_nu = 10
  ! ! PURPOSE
  ! ! bin width in nu-distribution (energy-transfer distribution)
  ! !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/FileNameFlux
  ! PURPOSE
  ! The absolute filename of the file containing flux info, if user supplied
  !
  ! possible values:
  ! * if given, but does not contain '/':
  !   default is '[path_To_Input]/[FileNameFlux]'
  ! * otherwise: filename is absolute, including path ('~' is okay)
  !
  ! NOTE
  ! if you want to use the file 'XXX.dat' in the actual directory,
  ! give it as './XXX.dat'
  !
  ! SOURCE
  !
  character(1000), save :: FileNameFlux = ''
  !****************************************************************************

  !****************************************************************************
  !****g* initNeutrino/equalWeights_Mode
  ! SOURCE
  integer, save :: equalWeights_Mode = 0
  !
  ! PURPOSE
  ! possible values are:
  ! * 0: default perweight mode is used (default)
  ! * 1: default perweight mode is used, but max is printed
  ! * 2: MC rejection method is used.
  !
  ! In the default mode, the perweights of the final particles are given by
  ! cross section/(A * numEnsembles)
  !
  ! If equalWeightsMode==2, then the perweights are given by
  ! equalWeights_Max/(A * numEnsembles)
  !
  ! Please check in the output the line "numberOfSuccess = ..." for the
  ! number of events actually generated.
  !****************************************************************************


  !****************************************************************************
  !****g* initNeutrino/equalWeights_Max
  ! SOURCE
  real, save :: equalWeights_Max = -1e99
  !
  ! PURPOSE
  ! The maximum value the MC-rejection method is done against.
  !****************************************************************************



  logical, save :: initFlag = .true.
  logical, save :: readinputflag = .true.
  integer, save :: first=0

  !resulting cross sections that can be accessed by analysis routines
  real,dimension(0:max_finalstate_ID)      :: sigabsArr
  real,dimension(0:max_finalstate_ID),save :: sigabsArrFinal

  logical, dimension(1:max_finalstate_ID) :: includeK
  integer, dimension(1:max_finalstate_ID),parameter :: K2Hist = (/&
  ! 1=QE, 2=Delta, 3=highRES, 4=1piBG, 5=DIS, 6=2p2hQE, 7=2p2hDelta, 8=2pi
       & 1,2,3,3,3,3,3,3,3,3,&
       & 3,3,3,3,3,3,3,3,3,3,&
       & 3,3,3,3,3,3,3,3,3,3,&
       & 3,4,4,5,6,7,8/)
  ! K2Hist controls output of various reaction components.
  ! The first 31 elements stand for the number of nucleon resonances taken
  ! into account. The many '3' have the effect that all higher-lying
  ! resonances beyond the Delta are for the analysis being summed into
  ! one effective higher-lying resonance.

  integer, parameter :: max_Hist = 8 ! number of different histograms

  logical, dimension(0:max_Hist), save :: includeHist

  ! needed only for checkEvent:
  integer,dimension(1:max_finalstate_ID),parameter :: eOrigin = (/&
       1,2,3,4,5,6,7,8,9,10, &
       11,12,13,14,15,16,17,18,19,20, &
       21,22,23,24,25,26,27,28,29,30, &
       31,origin_singlePi,origin_singlePi,origin_DIS, &
       origin_2p2hQE, origin_2p2hDelta, origin_doublePi /)
  ! origin_singlePi appears twice because of 2 possible final pion charges

  real, save :: ww = -99.9 ! for equalWeights mode

contains

  !****************************************************************************
  !****s* initNeutrino/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! This subroutine reads input out of jobcard
  ! from namelist 'neutrino_induced'.
  !****************************************************************************
  subroutine readInput
    use output
    use esample
    use inputGeneral, only: path_To_Input, ExpandPath
    use neutrinoParms, only: setDefault_Res

    integer :: ios

    !**************************************************************************
    !****n* initNeutrino/neutrino_induced
    ! NAME
    ! NAMELIST /neutrino_induced/
    ! PURPOSE
    ! This Namelist includes:
    ! * process_ID
    ! * flavor_ID
    ! * nuXsectionMode
    ! * nuExp
    ! * includeQE
    ! * includeDELTA
    ! * includeRES
    ! * include1pi
    ! * include2pi
    ! * includeDIS
    ! * include2p2hQE
    ! * include2p2hDelta
    ! * sigmacut
    ! * realRun
    ! * printAbsorptionXS
    ! * printInclHist
    ! * FileNameFlux
    ! * storeNucleon
    ! * equalWeights_Mode
    ! * equalWeights_Max
    !**************************************************************************
    NAMELIST /neutrino_induced/  process_ID,flavor_ID,nuXsectionMode,nuExp, &
         includeQE,includeDELTA,includeRES,include1pi,includeDIS,&
         include2p2hQE, include2p2hDelta, include2pi, storeNucleon, &
         sigmacut, realRun, printAbsorptionXS, printInclHist, FileNameFlux, &
         equalWeights_Mode, equalWeights_Max

    !**************************************************************************
    !****n* initNeutrino/lepton_bin
    ! NAME
    ! NAMELIST /lepton_bin/
    ! PURPOSE
    ! This Namelist includes:
    ! * cost_min
    ! * cost_max
    ! * delta_cost
    ! * Elept_min
    ! * Elept_max
    ! * delta_Elept
    ! * pL_min
    ! * pL_max
    ! * delta_pL
    ! * pT_min
    ! * pT_max
    ! * delta_pT
    ! * Q2_Max
    !**************************************************************************
    NAMELIST /lepton_bin/ cost_min,cost_max,delta_cost, &
         Elept_min,Elept_max,delta_Elept,pL_min,pL_max,delta_pL, &
         pT_min,pT_max,delta_pT,Q2_Max !,delta_Q2,numax,delta_nu

    if (.not.readinputflag) return

    call Write_ReadingInput('neutrino_induced',0)
    rewind(5)
    read(5,nml=neutrino_induced,IOSTAT=ios)
    call Write_ReadingInput("neutrino_induced",0,ios)

    select case (flavor_ID)
    case (electron,muon,taulepton)
       write(*,'(A,A,i5)') ' ...flavor = ',sFamily(flavor_ID)
    case default
       write(*,'(A,A,i5)') ' ...flavor = ','***unknown*** ', flavor_ID
       call TRACEBACK()
    end select

    select case (process_ID)
    case (1:3)
       write(*,'(A,A,i5)') ' ...Process: ',sProcess(process_ID)
    case (-3:-1)
       write(*,'(A,A,i5)') ' ...Process: anti-',sProcess(-process_ID)
    case default
       write(*,'(A,A,i5)') ' ...Process: ','***unknown*** ',process_ID
       call TRACEBACK()
    end select
    call setDefault_Res(process_ID)

    write(*,'(a,2I3)') ' nuXsectionMode,nuExp:',nuXsectionMode,nuExp

    write(*,*)
    write(*,'(a,L2)')' ...include QE        : ',includeQE
    write(*,'(a,L2)')' ...include Delta     : ',includeDELTA
    write(*,'(a,L2)')' ...include higher RES: ',includeRES
    write(*,'(a,L2)')' ...include 1pi       : ',include1pi
    write(*,'(a,L2)')' ...include 2pi       : ',include2pi
    write(*,'(a,L2)')' ...include DIS       : ',includeDIS
    write(*,'(a,L2)')' ...include 2p2h QE   : ',include2p2hQE
    write(*,'(a,L2)')' ...include 2p2h Delta: ',include2p2hDelta
    write(*,*)

    if (include2p2hDelta) call notInRelease("2p2p Delta")

    if (nuXsectionMode.ge.10) then
       select case (nuExp)
       case (1:numberOfExperiments)
          write(*,*) '##### calculation is done for the ',trim(sExp(nuExp)),&
               & ' experiment #####'
       case (99)
          write(*,*) '##### calculation is done for a user given flux #####'
          if (len_trim(FileNameFlux)>0) then
             if (index(FileNameFlux,"/")>0) then
                call ExpandPath(FileNameFlux)
                FileNameFlux = trim(FileNameFlux)
             else
                FileNameFlux = trim(path_to_Input)//'/'//trim(FileNameFlux)
             end if
             write(*,*) 'Flux.dat : ', trim(FileNameFlux)
          else
             write(*,*) 'filename "FileNameFlux" not provided. Stop!'
             call TRACEBACK()
          end if
       case default
          write(*,*) 'combination nuXsectionMode.ge.10.and.nuExp makes no sense -> STOP', &
             & nuexp,nuXsectionmode
          call TRACEBACK()
       end select
    end if

    if (nuExp.gt.0.and.nuXsectionMode.lt.10) then
       write(*,*) 'combination nuExp.gt.0.and.nuXsectionMode.lt.10 makes no sense -> STOP', &
          & nuexp,nuXsectionmode
       call TRACEBACK()
    end if

    select case (nuExp)
    case (1)
       if (process_ID.le.0) then
          write(*,*) 'MiniBooNE neutrino flux and antineutrino run makes no sense &
               &  -> STOP', nuexp, process_ID
          call TRACEBACK()
       end if
    case (5)
       if (process_ID.ge.0) then
          write(*,*) 'MiniBooNE antineutrino flux and neutrino run makes no sense &
               & -> STOP', nuexp, process_ID
          call TRACEBACK()
       end if
    case (6,11)
       if (process_ID.le.0) then
          write(*,*) 'MINOS neutrino flux and antineutrino run makes no sense &
               & -> STOP', nuexp, process_ID
          call TRACEBACK()
       end if
    case (7,12)
       if (process_ID.ge.0) then
          write(*,*) 'MINOS antineutrino flux and neutrino run makes no sense &
               & -> STOP', nuexp, process_ID
          call TRACEBACK()
       end if
    end select

    select case (storeNucleon)
    case (1)
       write(*,*) 'store nucleon for output: 1 = free nucleon'
    case (2)
       write(*,*) 'store nucleon for output: 2 = bound nucleon'
    case default
       write(*,*) 'store nucleon for output:', storeNucleon
       call TRACEBACK('wrong number for storeNucleon')
    end select

    select case (equalWeights_Mode)
    case (0)
       write(*,*) 'equalWeights:',equalWeights_Mode
    case (1:2)
       write(*,*) 'equalWeights:',equalWeights_Mode, equalWeights_Max
       if (equalWeights_Max <= 0.) then
          call TRACEBACK("You have to give 'equalWeights_Max'.")
       end if
    case default
       call TRACEBACK('wrong number for equalWeights_Mode')
    end select

    if (realRun) write(*,*) '#### REAL RUN ####'

    call Write_ReadingInput('neutrino_induced',1)

    ! ===================

    call Write_ReadingInput('lepton_bin',0)
    rewind(5)
    read(5,nml=lepton_bin,IOSTAT=ios)
    call Write_ReadingInput('lepton_bin',0,ios)
    write (*,*) 'cost_min =', cost_min,'cost_max =',cost_max, &
         'delta_cost =', delta_cost
    write (*,*) 'Elept_min=',Elept_min,'Elept_max=',Elept_max, &
         'delta_Elept=', delta_Elept

    call Write_ReadingInput('lepton_bin',1)

    readinputflag = .false.


  end subroutine readInput

  !****************************************************************************
  !****************************************************************************


  subroutine cleanUp
    use formfactors_A_main, only: cleanupMAID => cleanup

    call cleanupMAID
  end subroutine cleanUp


  !****************************************************************************
  !****s* initNeutrino/init_neutrino
  ! NAME
  ! subroutine init_neutrino(realParticles,pertParticles,raiseFlagIn,
  ! num_runs_sameEnergy,targetNuc)
  !
  ! PURPOSE
  ! This subroutine initializes a neutrino event on each nucleon in the
  ! realparticles vector (given to the routine). The resulting particles
  ! are set into the pertParticles vector. The reaction process is
  ! determined by the values read in by 'readInput':
  ! * process_ID
  ! * flavor_ID
  ! * nuXsectionMode
  !
  ! The user might also choose which contributions should be included
  ! (QE, DELTA and/or higher resonances) and whether the calculation should
  ! be done for a specific experiment.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParticles
  ! * integer                        :: num_runs_sameEnergy
  ! * logical                        :: raiseFlagIn -- if .true. then the
  !   energy etc is raised by raiseValue
  ! * type(tnucleus), pointer        :: targetNuc
  !
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: pertParticles
  !
  !****************************************************************************
  subroutine init_neutrino(realParticles,pertParticles,raiseFlagIn,&
       & num_runs_sameEnergy,targetNuc)
    use particleDefinition
    use random, only: rn
    use idtable
    use pauliBlockingModule, only: checkPauli
    use neutrinoSigma
    use propagation, only: gradients,updateVelocity
    use collisionNumbering, only: pert_numbering,real_numbering
    use insertion, only: setIntoVector
    use inputGeneral, only: fullEnsemble
    use output, only: WriteParticleVector,WriteParticle_debug,IntToChar
    use hist
    use hist2D
    use neutrinoInfoStorage
    use neutrinoProdInfo, only: neutrinoProdInfo_Init,neutrinoProdInfo_Store
    use offShellPotential, only: setOffShellParameter
    use expNeutrinoFluxes, only: getFluxEnu
    use esample
    use eN_eventDefinition
    use eN_event
    use eventGenerator_eN_lowEnergy, only: checkEvent
    use Coll_nuN, only: CalcXY
    use nucleusDefinition
    use NuclearPDF, only: SetNuclearPDFA
    use monteCarlo, only: MonteCarloChoose
    use constants, only: mN
    use residue, only: InitResidue, ResidueAddPH, ResidueSetWeight

    type(particle), dimension(:,:),intent(inOut) :: realParticles
    type(particle), dimension(:,:),intent(inOut) :: pertParticles
    integer, intent(in) :: num_runs_sameEnergy
    logical, intent(in) :: raiseFlagIn
    type(tnucleus), pointer :: targetNuc

    logical :: raiseFlag

    integer :: numtry=0
    integer,dimension(1:max_finalstate_ID)       :: numberofsuccess=0
    integer,dimension(1:max_finalstate_ID), save :: numberofsuccessfinal=0
    real :: sigtot!,sig
    real, dimension(1:max_finalstate_ID) ::  sigma=0.

    type(particle),dimension(1:max_finalstate_ID,1:2),target :: OutPart
    type(particle),dimension(1:20),target :: OutPartDIS
    type(particle),dimension(1:3), target :: OutPart2pi
    type(particle),dimension(:), pointer :: pOutPart
    type(particle),dimension(:), allocatable :: finalstate
    integer :: i, j, k

    real :: totalWeight
    logical :: setflag
    integer, save :: numberofcalls=0
    integer :: countfluxcutoff
    integer :: countsigmacutoff
    type(histogram), save :: energyInit

    real :: flux_enu
    integer :: firstEvent
    integer :: whichReal_nucleon,numNucleons
    integer :: failuresV=0,failuresO=0
    logical :: success
    real,dimension(1:3)  :: grad_P
    real :: fak1


    type(electronNucleon_event) :: eNeV0,eNev1
    type(electronNucleon_event),dimension(1:max_finalstate_ID) :: eNev
    integer :: number
    logical :: NumbersAlreadySet
    logical,save :: MCmode

    type(histogram), save :: hSigmaMC_qz(0:max_Hist)
    type(histogram), save :: hSigmaMC_nu(0:max_Hist)
    type(histogram), save :: hSigmaMC_Q2(0:max_Hist)
    type(histogram), save :: hSigmaMC_X(0:max_Hist)
    type(histogram), save :: hSigmaMC_Xrec(0:max_Hist)
    type(histogram2D), save :: hSigmaMC_nuQ2(0:max_Hist)
    type(histogram2D), save :: hSigmaMC_EprimeCost(0:max_Hist)
    type(histogram2D), save :: hsigmaMC_pLpT(0:max_Hist)
    type(histogram2D), save :: hsigmaMC_pLW(0:max_Hist)
    type(histogram2D), save :: hsigmaMC_pTW(0:max_Hist)
    type(histogram2D), save :: hSigmaMC_XY(0:max_Hist)
!
!  There are 3 different definitions for the invariant mass W:
!  1. a boson with four-momentum q hits Fermi-moving nucleon in potential well: W
!  2. a boson with four-momentum q hits Fermi-moving nucleon without potential: W_free
!  3. a boson with four-momentum q hits free nucleon: Wrec
!
    type(histogram), save :: hSigmaMC_W(0:max_Hist)
    type(histogram), save :: hSigmaMC_Wrec(0:max_Hist)
    type(histogram), save :: hSigmaMC_Wfree(0:max_Hist)

    real :: Q2max,numax,W2max,Eprime,cost,pL,pT,plept,Wrec

    logical :: flagDUMMY

    write(*,*)
    write(*,*) '################### NEUTRINO INIT STARTS #######################'

    if (initFlag) then
       call readInput
       call DoInit
       initFlag=.false.
    end if

    if (fullEnsemble.and.realRun) then
       write(*,*) 'FullEnsemble+Real particles in final state not yet implemented!!! STOP'
       call TRACEBACK()
    end if


    numberofcalls=numberofcalls+1
    if (numberofcalls.gt.num_runs_sameEnergy) then
       numberofcalls=1
       sigabsArrFinal = 0.
       numberofsuccessfinal=0
    end if

    sigabsArr = 0.
    numberofsuccess=0
    countsigmacutoff=0
    countfluxcutoff=0
    failuresV=0
    failuresO=0
    first=0

    !loop to determine numtry (number of testteilchen)
    numtry=0
    do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
          if (realParticles(i,j)%ID.ne.nucleon) cycle
          numtry=numtry+1
       end do
    end do

    call InitResidue(numtry,1,targetNuc%mass,targetNuc%charge)

    raiseflag = raiseflagin
    call neutrinoProdInfo_Init(numtry)

    call SetNuclearPDFA(targetNuc%mass)

    select case (nuExp)
    case (1)
       if (nuXsectionMode.eq.EXP_dSigmaMC_dQ2 &
            & .or.nuXsectionMode.eq.EXP_dSigmadEnu &
            & .or.nuXsectionMode.eq.EXP_dSigmaMC) then
          call neutrinoInfoStorage_Init(numtry)
       end if
    case (3)
       if (nuXsectionMode.eq.EXP_dSigmadEnu &
            & .or. nuXsectionMode.eq.EXP_dSigmaMC) then
          call neutrinoInfoStorage_Init(numtry)
       end if
    end select

    ! set the overall kinematics (most of it as dummy):
    call eNev_SetProcess(eNev0, process_ID,flavor_ID)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  ENSEMBLE LOOP  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Now a loop over all ensembles starts
! realParticles: in runs with frozen nuclear configuration these are the
! target nucleons, 1st index gives ensemble, 2nd index gives particle
! identity, p or n
!

    loopEnsemble: do i = lbound(realParticles,dim=1),ubound(realParticles,dim=1)
       if (realRun) then
          numNucleons=0
          do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)
             if (realParticles(i,j)%ID.ne.nucleon) cycle
             numNucleons=numNucleons+1
          end do
          if (numNucleons.gt.0) then
             ! Choose randomly one nucleon in the ensemble to make
             ! the collsion with:
             whichReal_nucleon=1+int(rn()*numNucleons)
          else
             cycle loopEnsemble
          end if
          numNucleons=0
       end if

       ! print out status info;
       if (mod(i,50)==0) write(*,*) 'now starting ensemble ',i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   BOUND NUCLEON LOOP   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       loopNucleon: do j = lbound(realParticles,dim=2),ubound(realParticles,dim=2)

          if (realParticles(i,j)%ID.ne.nucleon) cycle
          if (realRun) then
             ! Only allow for a collision with the chosen nucleon.
             numNucleons=numNucleons+1
             if (.not.(numNucleons.eq.whichReal_nucleon)) cycle
          end if

          flux_enu = -99.9
          if (nuExp.gt.0) then
             ! return sampled neutrino energy for given flux distribution:
             flux_enu = getFluxEnu(nuExp,&
                  Flavor_ID,Process_ID,FileNameFlux)

             ! take possible flux-cuts into account:
             if (flux_enu.lt.Enu_lower_cut .or. flux_enu.gt.Enu_upper_cut) then
                countfluxcutoff=countfluxcutoff+1
                cycle
             end if

             ! for statistics:
             ! The flux distribution is reflected in the density of sampled
             ! energies
             call AddHist(energyInit,flux_enu,1.)

          end if

          ! Step 1 of "neutrino init": Set the target nucleon:
          eNev1 = eNev0
          call eNev_init_nuStep1(eNev1,realParticles(i,j))

          sigma=0.

          ! prepare calculation of cross section if in MC mode:
          if (MCmode) call SetXsecMC(eNev1,flux_enu,nuXsectionMode)

          !          call write_electronNucleon_event(eNev1)

          call resetNumberGuess()

!!$          ! workaround for Christine (Hang Qi):
!!$          if (eNev1%Q2 < 2.0) then
!!$             write(*,*) "Q2 too small:", eNev1%Q2
!!$             cycle loopNucleon
!!$          else
!!$             write(*,*) "Q2 okay:", eNev1%Q2
!!$          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!   LOOP OVER REACTION MECHANISMS and FINAL STATES   !!!!!!!!!!!!!!
          !
          ! The next loop runs over the IDs of all possible final states
          ! max_finalstate_ID set earlier in this file
          !
          ! This loop calculates the inclusive cross section for all processes
          ! (aka final state) and generates an array.
          !
          ! This array will be used below for the MC decision.
          !
          ! Also the possible final state particles for every process are stored

          loopFinalState: do k=1, max_finalstate_ID ! = ..., 37

             if (.not.includeK(k)) cycle

             ! set event info for this channel to the default
             eNev(k) = eNev1

             ! Set pointer to array, where outgoing particles to store:
             ! (2pi needs larger array, DIS even larger)
             select case (k)
             case (chDIS)
                pOutPart => OutPartDIS(:)
             case (chTwoPion)
                pOutPart => OutPart2pi(:)
             case default
                pOutPart => OutPart(k,:)
             end select


             select case (nuXsectionMode)

             case (integratedSigma)
                call Xsec_integratedSigma(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadCosThetadElepton)
                call Xsec_dSigmadCosThetadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadQ2dElepton)
                call Xsec_dSigmadQ2dElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadCosTheta)
                call Xsec_dSigmadCosTheta(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmadElepton)
                call Xsec_dSigmadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmaMC)
                call Xsec_SigmaMC(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmaMC_dW)
                call Xsec_SigmaMC_W(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

             case (dSigmaMC_dQ2)
                call Xsec_SigmaMC_Q2(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k))

            !cross sections for specific experiments with given flux

             case (EXP_dSigmadEnu)
                call Xsec_integratedSigma(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadCosThetadElepton)
                call Xsec_dSigmadCosThetadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadQ2dElepton)
                call Xsec_dSigmadQ2dElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadCosTheta)
                call Xsec_dSigmadCosTheta(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmadElepton)
                call Xsec_dSigmadElepton(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmaMC)
                call Xsec_SigmaMC(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmaMC_dW)
                call Xsec_SigmaMC_W(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)

             case (EXP_dSigmaMC_dQ2)
                call Xsec_SigmaMC_Q2(eNev(k), k, &
                     &  raiseFlag,raiseVal,pOutPart,sigma(k), flux_enu)


             case default
                write(*,*) 'nuXsectionMode=',nuXsectionMode
                call TRACEBACK('error in case nuXsectionMode')

             end select

             raiseflag=.false.
!!!!!!!!!!!!!!!!!!!!!!
!  Pauli blocking

             if (sigma(k).ne.0.) then
                if (.not.checkPauli(pOutPart,realParticles)) sigma(k)=0.
             end if

!  End Pauli blocking
!!!!!!!!!!!!!!!!!!!!!!

          end do loopFinalState
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          if (realRun) sigma=sigma*real(size(realParticles,2))/11.
          ! number 11 comes from main.f90 where lengthReal is
          ! set to lengthReal=targetNuc%mass+targetNuc%mass*10

          sigtot=sum(sigma)
          if (debugflag) write(*,'(A,g13.5)') ' sigtot= ', sigtot

!!$   write(*,'(A,g13.5)') ' sigtot= ', sigtot
!!$   write(*,'(2g13.5)') sigma(1:2)

          ! check if total cross section smaller than fixed cutoff,
          ! if so, then throw cross section away

          if (sigtot.lt.sigmacut) then
             countsigmacutoff=countsigmacutoff+1
             if (debugflag) write(*,'(A,i7.3)') ' countsigmacutoff=',&
                  & countsigmacutoff
             if (debugflag) write(*,'(2(A,g12.5))') &
                  & 'In initNeutrino.f90: sigtot=',sigtot,&
                  &' is less than sigmacut=',&
                  & sigmacut
             cycle loopNucleon
          end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!   DO THE MC-DECISION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          k = MonteCarloChoose(sigma,totalWeight)
          ! write(*,*) 'MonteCarlo k =', k,'sigma(k)=',sigma(k)
          ! The field 'sigma' contains the cross sections of the various
          ! subprocesses, as calculated above, 'total weight' gives the total
          ! cross section, summed over all subprocesses.
          ! 'k' is then the statistically chosen final state

          ! In the 'equal-weights-mode', we have to do an additional rejection:
          if(equalWeights_Mode>0) then

             if (abs(totalWeight)> ww) then
                ww = abs(totalWeight)
                write(87,*) ww
                flush(87)
             end if

             if(equalWeights_Mode>1) then
                if (totalWeight > equalWeights_Max) then
                   write(*,*) 'totalWeight > equalWeights_Max: ',&
                        totalWeight,equalWeights_Max
                   call TRACEBACK("You have to increase 'equalWeights_Max'.")
                end if

                if (totalWeight < rn()*equalWeights_Max) cycle ! event rejected

                sigma(1:) = sigma(1:)*equalWeights_Max/abs(totalWeight) ! --> statistics
                totalWeight = sign(equalWeights_Max, totalWeight)   ! --> perweight

             end if

          end if

          firstEvent=getFirstEvent()


          select case (k)
          case (0)
             write(*,*) 'Problem initNeutrino: no event generated:', sigma
             call TRACEBACK()
          case (chDIS)
             pOutPart => OutPartDIS(:)
          case (chTwoPion)
             pOutPart => OutPart2pi(:)
          case default
             pOutPart => OutPart(k,:)
          end select

          allocate(finalstate(size(pOutPart)))
          finalstate = pOutPart

          finalState%pert=(.not.realRun)
          finalState%firstEvent=FirstEvent         ! Number of first event,
                                                   ! stays constant during run
          finalState%perweight=totalWeight/float(numtry)
          finalState%history=0

          numberofsuccess(k)=numberofsuccess(k)+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          fak1 = totalWeight/float(numtry)
          call doStoreNeutrinoInfo(firstEvent, k, fak1, eNev(k))

          call fillHist(k, fak1, eNev(k))

          sigabsArr(0)  = sigabsArr(0)  + fak1
          sigabsArr(1:) = sigabsArr(1:) + sigma(1:)/float(numtry)


          ! calculate the finalstate offshellness

          select case (k)
          case default ! ==== single-pi, DIS, 2p2h, two-pi-backg

             finalstate%offshellPar=0.
             call updateVelocity(finalstate)

          case (:31)    ! ==== RES or QE

             call setOffShellParameter(finalstate(1),success)
             if (.not.success) then
                if (debugflag) write(*,*) &
                     'offshell parameter > max_offshellparameter &
                     & in initNeutrino => sig=0'
                failuresO=failuresO+1
                deallocate(finalstate)
                cycle loopNucleon
             end if

             call gradients(finalstate(1),grad_P) ! Evaluate dH/dp
             finalstate(1)%vel=grad_P
             if (1. - Dot_Product(grad_P(1:3),grad_P(1:3)) .le. 0.) then
                write(*,'(A,5G13.5)')'problems in initNeutrino: &
                     & velocity**2 greater or equal 1.', &
                     & k, Dot_Product(grad_P(1:3) &
                     & ,grad_P(1:3)), finalstate(1)%mass,&
                     & finalstate(1)%offshellPar,sigma(k)/float(numtry)
                if (debugflag) then
                   call WriteParticle_debug(finalstate(1))
                   write(*,*)
                   write(*,*) '...this funny particle is now deleted'
                end if
                failuresV=failuresV+1
                deallocate(finalstate)
                cycle loopNucleon
             end if

          end select



          ! now check events for charge and momentum conservation,
          ! if violated throw away

          if (.not.checkEvent(eNev(k),finalstate,eOrigin(k))) then
             deallocate(finalstate)
             write(*,*) 'conservations violated'
             cycle loopNucleon
 !            call TRACEBACK('conservations violated.')
          end if

          ! now we can give them the first-event-number:
          if (realRun) then
             number = real_numbering()
          else
             number = pert_numbering(realParticles(i,j))
          end if
          finalState%event(1) = number
          finalState%event(2) = number

          ! Add a hole in the target nucleus (analysis only):
          call ResidueAddPH(firstEvent,realParticles(i,j))
          call ResidueSetWeight(firstEvent,finalstate(1)%perweight)

          if (k==chQE2p2h) &
               call ResidueAddPH(firstEvent,eNev(k)%nucleon2)

          ! set the particles in the particle vector:

          if (k.ne.chDIS) call resetNumberGuess()
          NumbersAlreadySet = AcceptGuessedNumbers()

          if (fullEnsemble) then
             if (realRun) &
                  call TRACEBACK('RealRun+fullEnsemble not yet implemented: initNeutrino')
             call setIntoVector(finalState,pertParticles,setFlag,NumbersAlreadySet)
          else
             if (realRun) then
                call TRACEBACK('RealRun seems to be corrupted: initNeutrino')

!!! ATTENTION !!! The following seems to be untested !!!!!!!!

                ! Delete the scattering partner
                realParticles(i,j)%ID=0
                ! Set new particles into real particle vector
                realParticles(i:i,:)%perweight=finalState(1)%perweight      ! ????????
                realParticles(i:i,:)%firstEvent=finalState(1)%firstEvent    ! ????????
                call setIntoVector(finalState,realParticles(i:i,:),setFlag,NumbersAlreadySet)
                ! ????????
             else
                call setIntoVector(finalState,pertParticles(i:i,:),setFlag,NumbersAlreadySet)
             end if
          end if

          deallocate(finalstate)

          if (.not.setFlag) then
             write(*,*) 'error setIntoVector in initNeutrino'
             if (.not.realRun) then
                call WriteParticleVector('pert',pertParticles)
             else
                call WriteParticleVector('real',realParticles)
             end if
             call TRACEBACK()
          end if

       end do loopNucleon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!   END OF BOUND NUCLEON LOOP   !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do loopEnsemble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!   END OF LOOPS OVER ENSEMBLE   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    call SetNuclearPDFA(1)

    write(*,'(40("-"))')
    write(*,'(a,I10)') ' failures v>c            :', failuresV
    write(*,'(a,I10)') ' failures due to offshell:', failuresO
    write(*,'(a,I10)') ' numtry                  =', numtry
    write(*,'(a,I10)') ' numberOfSuccess         =', sum(numberofsuccess)
    write(*,'(a,I10)') ' countSigmaCutoff        =', countsigmacutoff
    write(*,'(a,I10)') ' countFluxCutoff         =', countfluxcutoff
    write(*,'(a,F12.4)')         ' sigabs       =', sigabsArr(0)
    write(*,'(a,F12.4,a,I10,a)') ' sigabsQE     =', sigabsArr(1), '(#',numberofsuccess(1),')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabsDelta  =', sigabsArr(2), '(#',numberofsuccess(2),')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabshighRes=', sum(sigabsArr(3:31)), &
                                                     & '(#',sum(numberofsuccess(3:31)), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs1pi (n)=', sigabsArr(32), '(#',numberofsuccess(32), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs1pi (p)=', sigabsArr(33), '(#',numberofsuccess(33), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs2pi    =', sigabsArr(37), '(#',numberofsuccess(37), ')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabsDIS    =', sigabsArr(34), '(#',numberofsuccess(34),')'
    write(*,'(a,F12.4,a,I10,a)') ' sigabs2p2h   =', sum(sigabsArr(35:36)), &
                                                     & '(#',sum(numberofsuccess(35:36)),')'
    write(*,'(40("-"))')

    sigabsArrFinal = sigabsArrFinal+sigabsArr
    numberofsuccessfinal=numberofsuccess+numberofsuccessfinal

    call DoWrite

    !**************************************************************************
    !****o* initNeutrino/neutrino_initialized_energyFlux.dat
    ! NAME
    ! file neutrino_initialized_energyFlux.dat
    ! PURPOSE
    ! This file provides the incoming neutrino flux
    ! Only energies at which actual sampling has taken place are in this file
    ! The density of sampled energies reflects the flux distribution which in
    ! turn is reflected in this histogram of sampled energies.
    !
    ! Columns:
    ! * #1: Enu (in GeV)
    ! * #2: flux (in 1/GeV), for the calculations of X-sections relevant is only
    !   the normalized flux
    ! * #3: number of events
    ! * #4: should be ignored
    !**************************************************************************
    if (nuExp.gt.0) then
       call writeHist(energyInit,10, &
            & file='neutrino_initialized_energyFlux.dat')
    end if


    if (printAbsorptionXS) then      ! inclusive X-sections are printed
       if (printInclHist) call writeInclHist()
    end if

    write(*,*) '################### NEUTRINO INIT FINISHED #######################'

  contains
    !**************************************************************************
    !****s* init_neutrino/DoInit
    ! NAME
    ! subroutine DoInit
    ! PURPOSE
    ! Doing the actual initialzation steps
    !**************************************************************************
    subroutine DoInit()

      integer :: iHist

      !---------------------------------------------------------------------
      ! setting up some arrays for switching on/off the channels
      !---------------------------------------------------------------------
      includeHist = (/.true.,includeQE, includeDELTA, includeRES, &
           include1pi, includeDIS, include2p2hQE, include2p2hDelta, &
           include2pi /)

      includeK = .false.
      do k=1,max_finalstate_ID

         select case (k) ! === check for inclusion of process
         case (nucleon)
            if (.not.includeQE) cycle
         case (delta)
            if (.not.includeDELTA) cycle
         case (P11_1440:F37_1950)
            if (.not.includeRES) cycle
         case (chOnePionN,chOnePionP)
            if (.not.include1pi) cycle
         case (chDIS)
            if (.not.includeDIS) cycle
         case (chQE2p2h)
            if (.not.include2p2hQE) cycle
         case (chDelta2p2h)
            if (.not.include2p2hDelta) cycle
         case (chTwoPion)
            if (.not.include2pi) cycle
         end select

         select case (k) ! === check for inclusion of resonance
         case (S11_2090,D13_2080,G17_2190,P11_2100,P13_1900,F15_2000,S31_1900,&
              D33_1940,D35_1930,D35_2350,P31_1750,F35_1750)
            cycle
         end select

         includeK(k) = .true.
      end do

      !---------------------------------------------------------------------
      ! initialize some output files
      !---------------------------------------------------------------------
      open(10,File='neutrino_info.dat')
      if (process_ID.lt.0) then
         write(10,*)'# process_ID : anti-',sProcess(-process_ID)
      else
         write(10,*)'# process_ID : ',sProcess(process_ID)
      end if
      write(10,*)'# flavor_ID  : ',sFamily(flavor_ID)
      write(10,*)'# '
      write(10,'(a,L2)')' #...include QE        : ',includeQE
      write(10,'(a,L2)')' #...include Delta     : ',includeDELTA
      write(10,'(a,L2)')' #...include higher RES: ',includeRES
      write(10,'(a,L2)')' #...include 1pi       : ',include1pi
      write(10,'(a,L2)')' #...include 2pi       : ',include2pi
      write(10,'(a,L2)')' #...include DIS       : ',includeDIS
      write(10,'(a,L2)')' #...include 2p2h QE   : ',include2p2hQE
      write(10,'(a,L2)')' #...include 2p2h Delta: ',include2p2hDelta
      write(10,*)'# '
      write(10,*)'# nuXsectionMode: ',sXsectionMode(MOD(nuXsectionMode,10))
      write(10,*)'# '
      if (nuExp.gt.0) write(10,*)'##### calculation is done for the ',&
           & trim(sExp(nuExp)),' experiment #####'
      close(10)

      if (printAbsorptionXS) then

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section.dat
      ! NAME
      ! file neutrino_absorption_cross_section.dat
      ! PURPOSE
      ! The file is produced in the runs with eventtype=5=neutrino .
      !
      ! The file shows the absorption cross section for lepton
      ! ( neutrino or charged lepton) scattering for the sum of all channels
      ! which were set to TRUE in the namelist "neutrino_induced"
      ! ( QE+Delta+highRES+1pi+DIS if includeQE, includeDELTA, includeRES,
      ! include1pi, includeDIS were TRUE)
      !
      ! For process_ID=CC and NC  the units 10^{-38} cm^2 for integrated xsec
      ! (10^{-38)cm^2/GeV  for  dsigma/dElepton,  10^{-38)cm^2/GeV^2  for
      ! dsigma/dQ^2, and so on)
      ! For process_ID=EM the units are nanobars=10^{-33}cm^2
      ! All cross sections are given per nucleon (1/A)
      !
      ! Columns:
      ! * #1: variable which was raised
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQ2 mode, Elepton for
      !   nuXsectionMode=2=dSigmadQ2dElepton  and so on)
      ! * #2: cross section
      !************************************************************************
      open(10,File='neutrino_absorption_cross_section.dat')
      write(10,*) '#  raiseVal, xsection'
      close(10)

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_ALL.dat
      ! NAME
      ! file neutrino_absorption_cross_section_ALL.dat
      ! PURPOSE
      ! More detailed information than in
      ! neutrino_absorption_cross_section.dat:
      !
      ! All cross sections are given per nucleon (1/A)
      ! Columns:
      ! * #1: variable which was raised (the same as
      !   neutrino_absorption_cross_section.dat)
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQ2 mode, Elepton for
      !   nuXsectionMode=2=dSigmadQ2dElepton )
      ! * #2: cross section, sum over all channels  (the same as
      !   neutrino_absorption_cross_section.dat)
      ! * #3: cross section for QE events (the same as column 2 in
      !   neutrino_absorption_cross_section_QE.dat)
      ! * #4: cross section for DELTA events (the same as column 2 in
      !   neutrino_absorption_cross_section_Delta.dat)
      ! * #5: cross section for hihgRES events (the same as column 2 in
      !   neutrino_absorption_cross_section_highRES.dat)
      ! * #6: cross section for 1pi events (the same as column2+column3 in
      !   neutrino_absorption_cross_section_1pi.dat)
      ! * #7: cross section for DIS events (the same as column 2 in
      !   neutrino_absorption_cross_section_DIS.dat)
      ! * #8: cross section for 2p2h QE events
      ! * #9: cross section for 2p2h Delta events
      ! * #10: cross section for 2 pion background events
      !
      !************************************************************************
      open(10,File='neutrino_absorption_cross_section_ALL.dat')
      write(10,*) '# 1:var 2:sum 3:QE 4:Delta 5:highRES 6:1pi 7:DIS 8:2p2h-QE &
                  & 9:2p2h-Delta 10:2pi'
      close(10)

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_QE.dat
      ! NAME
      ! file neutrino_absorption_cross_section_QE.dat
      ! PURPOSE
      ! The same structure as neutrino_absorption_cross_section.dat,
      ! but only for QE events (=the first interaction act was quasielastic
      ! or elastis scattering)
      !************************************************************************
      if (includeQE) then
         open(10,File='neutrino_absorption_cross_section_QE.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_2p2h.dat
      ! NAME
      ! file neutrino_absorption_cross_section_2p2h.dat
      ! PURPOSE
      ! The same structure as neutrino_absorption_cross_section.dat,
      ! but only for 2p2h events (=the first interaction was 2p2h)
      !************************************************************************
      if (include2p2hQE) then
         open(10,File='neutrino_absorption_cross_section_2p2h.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_Delta.dat
      ! NAME
      ! file neutrino_absorption_cross_section_Delta.dat
      ! PURPOSE
      ! The Delta production events (=the first interaction was
      ! production of the Delta resonance)
      !************************************************************************
      if (includeDELTA) then
         open(10,File='neutrino_absorption_cross_section_Delta.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_highRES.dat
      ! NAME
      ! file neutrino_absorption_cross_section_highRES.dat
      ! PURPOSE
      ! For highRES production events (=the first interaction was
      ! production of any resonance beyond Delta)
      !
      ! Columns:
      ! * #1: variable which was raised (the same as
      !   neutrino_absorption_cross_section.dat)
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQ2 mode, Elepton for
      !   nuXsectionMode=2=dSigmadQ2dElepton )
      ! * #2: cross section, sum over all higher resonances beyond the Delta
      ! * #3 - 31: contribution of individual nucleon resonances beyond Delta
      !   Individual resonance numbers in Module IdTable
      !   3: P11(1440), 4: S11(1535), ..., 7: D13(1520), ... etc
      !
      !************************************************************************
      if (includeRES) then
         open(10,File='neutrino_absorption_cross_section_highRES.dat')
         write(10,*) '#  raiseVal, sum, single contribution xsection'
         write(10,*) '#  columns 3 - 31 give results for N* > Delta &
              &resonances from module IdTable'
         write(10,'(A9,30I13)') '#       ',(k,k=2,31)
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_1pi.dat
      ! NAME
      ! file neutrino_absorption_cross_section_1pi.dat
      ! PURPOSE
      ! Nearly the same structure as neutrino_absorption_cross_section.dat,
      ! but only for nonresonant 1-pion production events (=the first
      ! interaction was production of 1-pion final state)
      !
      ! Columns:
      ! * #2: cross section  for the channel with neutron in the final state
      !   ("final" here is "after the first interaction act", that is before
      !   final state interactions,
      !   e.g. nu n \to mu- n pi^0)
      ! * #3: cross section  for the channel with proton in the final state,
      !   e.g. nu n \to mu- p pi^-)
      !************************************************************************
      if (include1pi) then
         open(10,File='neutrino_absorption_cross_section_1pi.dat')
         write(10,*) '#  raiseVal, pi neutron-backgr, pi proton-backgr'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_2pi.dat
      ! NAME
      ! file neutrino_absorption_cross_section_2pi.dat
      ! PURPOSE
      ! The same structure as neutrino_absorption_cross_section.dat,
      ! but only for 2piBG production events (=the first interaction was 2piBG)
      !************************************************************************
      if (include2pi) then
         open(10,File='neutrino_absorption_cross_section_2pi.dat')
         write(10,*) '#  raiseVal, 2pi BG'
         close(10)
      end if


      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_DIS.dat
      ! NAME
      ! file neutrino_absorption_cross_section_DIS.dat
      ! PURPOSE
      ! The same structure as neutrino_absorption_cross_section.dat,
      ! but only for DIS production events (=the first interaction was DIS)
      !************************************************************************
      if (includeDIS) then
         open(10,File='neutrino_absorption_cross_section_DIS.dat')
         write(10,*) '#  raiseVal, xsection'
         close(10)
      end if

      !************************************************************************
      !****o* initNeutrino/neutrino_absorption_cross_section_numbers.dat
      ! NAME
      ! file neutrino_absorption_cross_section_numbers.dat
      ! PURPOSE
      ! Shows the number of test particles which interacted (here we mean
      ! the first interaction) via various channels
      !
      ! Columns:
      ! * #1: variable which was raised (the same as
      !   neutrino_absorption_cross_section.dat)
      !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQ2 mode, Elepton for
      !   nuXsectionMode=2=dSigmadQ2dElepton )
      ! * #2: QE channel
      ! * #3: Delta production
      ! * #4: P_{11}(1440) production
      ! * #5: S_{11}(1535) production
      ! * #6: D_{13}(1520) production
      ! * #7: neutron and pi+  in the final state
      ! * #8: proton and pi0 in the final state
      ! * #9: number of test particles that underwent any kind of interaction
      ! * #10: number of test particles (= number of tries)
      !************************************************************************
      open(10,File='neutrino_absorption_cross_section_numbers.dat')
      write(10,*) '#  raiseVal, numberofQE, numberofDELTA, numberP11,', &
           &' numberS11, numberD13,', &
           & '# number_npi+, number_ppi0, numberTOT, numberTestParticles'
      close(10)

      end if

      !---------------------------------------------------------------------
      ! Allocate some histograms
      !---------------------------------------------------------------------

      MCmode = ((nuXsectionMode.eq.dSigmaMC)&
           & .or.(nuXsectionMode.eq.EXP_dSigmaMC) &
           & .or.(nuXsectionMode.eq.dSigmaMC_dQ2) &
           & .or.(nuXsectionMode.eq.EXP_dSigmaMC_dQ2) &
           & .or.(nuXsectionMode.eq.dSigmaMC_dW) &
           & .or.(nuXsectionMode.eq.EXP_dSigmaMC_dW) )



      if (nuExp.gt.0) then
         ! create a histogram for the sampled energies
         ! the histogram reflects the distribution of sampled energies;
         ! this distribution is determined by the original flux distribution
         call CreateHist(energyInit, 'initialized energy',0.,150.,0.02)
      end if

      ! this is to estimate the "reasonable" energy transfer in the neutrino
      ! experiments in order to use it for initialisation of histograms
      ! for differential absorption cross sections
      if (nuExp.gt.0) then
         numax = getNuMaxExp()
      else
      ! if nuExp = 0 then numax is set to enu
         call get_sigma_namelist(XsectionMode=nuXsectionMode,Genu=numax)

      end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! In these 2D histograms the last 3 2d arrays determine the binning.
! For example, in the EprimeCost histogram the input
! (/0.,-1./),(/numax,1.0/),(/0.05,0.01/) stands for 0 < Eprime < numax with a
! binwidth of 0.05 and -1. cost < +1. with a bindwidth of 0.01
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! now set parameters for binning
      numax = max(numax,1e-10)
      W2max = mN**2 + 2*mN*numax     ! rough estimate
      Q2max = 2*mN*numax             ! rough estimate
      Q2max = min(Q2_Max,Q2max)

      write(*,*) 'In initialization of the histograms for absorption xsec:'
      write(*,*) '    numax=', numax, '   Q2max=', Q2max, '   W2max=',W2max

      do iHist=0,8
         if (.not.includeHist(iHist)) cycle
         call CreateHist2D(hSigmaMC_nuQ2(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs nu and Q2 ', &
              & (/0.,0./),(/numax,Q2max/),(/numax/100.,Q2max/500.0/))
         call CreateHist2D(hSigmaMC_EprimeCost(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Eprime and cost ', &
              & (/Elept_min,cost_min/),(/Elept_max,cost_max/),(/delta_Elept,delta_cost/))
         call CreateHist2D(hSigmaMC_pLpT(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs pL and pT ', &
              & (/pL_min,pT_min/),(/pL_max,pT_max/),(/delta_pL,delta_pT/))
         call CreateHist2D(hSigmaMC_pLW(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs pL and Wrec ', &
              & (/pL_min,0.8/),(/pL_max,6.0/),(/delta_pL,0.1/))
         call CreateHist2D(hSigmaMC_pTW(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs pT and Wrec ', &
              & (/pT_min,0.8/),(/pT_max,6.0/),(/delta_pT,0.1/))
         call CreateHist2D(hSigmaMC_XY(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs x and y ', &
              & (/0.,0./),(/2.1,2.1/),(/1.0/100,1.0/100/))
         call CreateHist(hSigmaMC_nu(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs nu ', &
              & 0.,numax,numax/100.)
         call CreateHist(hSigmaMC_qz(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs qz=sqrt(nu^2+Q2) ', &
              & 0.,numax,numax/100.)
         call CreateHist(hSigmaMC_Q2(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Q2 ', &
              & 0.,Q2max,Q2max/500.)
         call CreateHist(hSigmaMC_X(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs X_Bjorken ', &
              & 0.,2.0,0.01)
         call CreateHist(hSigmaMC_Xrec(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs X_Bjorken reconstructed = &
              & Q2/2./mN/nu', 0.,2.0,0.01)
         call CreateHist(hSigmaMC_Wrec(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Wrec = &
              & sqrt(mN^2 +2*mN*nu - Q2) ', &
              & 0.,sqrt(W2max),sqrt(W2max)/100.)
         call CreateHist(hSigmaMC_Wfree(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs Wfree', &
              & 0.,sqrt(W2max),sqrt(W2max)/100.)
         call CreateHist(hSigmaMC_W(iHist), &
              & 'sigma ('//trim(intToChar(iHist))//') vs W', &
              & 0.,sqrt(W2max),0.02)
      end do

    end subroutine DoInit


    !**************************************************************************
    !****s* init_neutrino/DoWrite
    ! NAME
    ! subroutine DoWrite
    ! PURPOSE
    ! Doing the write out of most data
    !**************************************************************************
    subroutine DoWrite()

      if (.not.printAbsorptionXS) return

      open(10,File='neutrino_absorption_cross_section.dat',&
           & position='append')
      if (numberofcalls.ne.1) backspace(10)
      write(10,'(10g13.5)') raiseVal, sigabsArrFinal(0)/real(numberofcalls)
      close(10)

      open(10,File='neutrino_absorption_cross_section_ALL.dat',&
           & position='append')
      if (numberofcalls.ne.1) backspace(10)
      write(10,'(30g13.5)') raiseVal, sigabsArrFinal(0)/real(numberofcalls),&
           & sigabsArrFinal(1)/real(numberofcalls), &
           & sigabsArrFinal(2)/real(numberofcalls), &
           & sum(sigabsArrFinal(3:31))/real(numberofcalls), &
           & sum(sigabsArrFinal(32:33))/real(numberofcalls), &
           & sigabsArrFinal(34)/real(numberofcalls), &
           & sigabsArrFinal(35)/real(numberofcalls), &
           & sigabsArrFinal(36)/real(numberofcalls), &
           & sigabsArrFinal(37)/real(numberofcalls)
      close(10)

!!$      write(113,'(10g13.5)') raiseVal, sigabsArrFinal(0)/real(numberofcalls),&
!!$           & sigabsArrFinal(1)/real(numberofcalls), &
!!$           & sigabsArrFinal(2)/real(numberofcalls), &
!!$           & sum(sigabsArrFinal(3:31))/real(numberofcalls), &
!!$           & sum(sigabsArrFinal(32:33))/real(numberofcalls), &
!!$           & sigabsArrFinal(34)/real(numberofcalls)

      if (includeQE) then
         open(10,File='neutrino_absorption_cross_section_QE.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(1)/real(numberofcalls)
         close(10)
      end if

      if (include2p2hQE) then
         open(10,File='neutrino_absorption_cross_section_2p2h.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(35)/real(numberofcalls)
         close(10)
      end if

      if (includeDELTA) then
         open(10,File='neutrino_absorption_cross_section_Delta.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(2)/real(numberofcalls)
         close(10)
      end if

      if (includeRES) then
         open(10,File='neutrino_absorption_cross_section_highRES.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(31g13.5)') raiseVal, &
              & sum(sigabsArrFinal(3:31))/real(numberofcalls),&
              & sigabsArrFinal(3:31)/real(numberofcalls)
         close(10)
      end if

      if (include1pi) then
         open(10,File='neutrino_absorption_cross_section_1pi.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, &
              & sigabsArrFinal(32:33)/real(numberofcalls)
         close(10)
      end if

      if (includeDIS) then
         open(10,File='neutrino_absorption_cross_section_DIS.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(34)/real(numberofcalls)
         close(10)
      end if

      if (include2pi) then
         open(10,File='neutrino_absorption_cross_section_2pi.dat',&
              & position='append')
         if (numberofcalls.ne.1) backspace(10)
         write(10,'(10g13.5)') raiseVal, sigabsArrFinal(37)/real(numberofcalls)
         close(10)
      end if

      open(10,File='neutrino_absorption_cross_section_numbers.dat',&
           & position='append')
      if (numberofcalls.ne.1) backspace(10)
      write(10,'(10g13.5)') raiseVal, numberofsuccessfinal(1)/numberofcalls,&
           & numberofsuccessfinal(2)/numberofcalls, &
           & numberofsuccessfinal(3)/numberofcalls,&
           & numberofsuccessfinal(4)/numberofcalls, &
           & numberofsuccessfinal(7)/numberofcalls,&
           & numberofsuccessfinal(32)/numberofcalls, &
           & numberofsuccessfinal(33)/numberofcalls,&
           & sum(numberofsuccessfinal)/numberofcalls,numtry
      close(10)

    end subroutine DoWrite


    !**************************************************************************
    !****is* init_neutrino/doStoreNeutrinoInfo
    ! PURPOSE
    ! store the neutrino and the nucleon info
    !**************************************************************************
    subroutine doStoreNeutrinoInfo(firstEvent, k, fak, eNev)
      use eN_eventDefinition

      integer, intent(in) :: firstEvent, k
      real, intent(in) :: fak
      type(electronNucleon_event),intent(in) :: eNeV

      select case (storeNucleon)
      case (1)
         call neutrinoProdInfo_Store(firstEvent, k, fak,&
              eNev%lepton_in%mom,&
              eNev%lepton_out%mom,&
              eNev%boson%mom,&
              eNev%nucleon_free%mom,&
              eNev%nucleon_free%charge)
      case (2)
         call neutrinoProdInfo_Store(firstEvent, k, fak,&
              eNev%lepton_in%mom,&
              eNev%lepton_out%mom,&
              eNev%boson%mom,&
              eNev%nucleon%mom,&
              eNev%nucleon%charge)
      end select

    end subroutine doStoreNeutrinoInfo

    !**************************************************************************
    !****is* init_neutrino/fillHist
    !**************************************************************************
    subroutine fillHist(k,fak, eNev)
      use eN_eventDefinition

      integer, intent(in) :: k
      real, intent(in) :: fak
      type(electronNucleon_event),intent(in) :: eNeV

      integer :: iHist

      iHist = K2Hist(k)

      call AddHist2D(hSigmaMC_nuQ2(0),hSigmaMC_nuQ2(iHist),&
           (/eNev%boson%mom(0),eNev%Q2/), fak)
      call AddHist(hSigmaMC_nu(0),hSigmaMC_nu(iHist),&
           eNev%boson%mom(0), fak)
      call AddHist(hSigmaMC_qz(0),hSigmaMC_qz(iHist),&
           sqrt(eNev%boson%mom(0)**2+eNev%Q2), fak)
      call AddHist(hSigmaMC_Q2(0),hSigmaMC_Q2(iHist),&
           eNev%Q2, fak)
      call AddHist(hSigmaMC_Xrec(0),hSigmaMC_Xrec(iHist),&
           eNev%Q2/(2*mN*eNev%boson%mom(0)), fak)

      if (eNev%W > 0) then
         call AddHist(hSigmaMC_W(0),hSigmaMC_W(iHist),&
              eNev%W, fak)
      end if

      if (eNev%W_free > 0) then
         !        Write(*,*) 'Wfree', eNev%W_free
         call AddHist(hSigmaMC_Wfree(0),hSigmaMC_Wfree(iHist),&
              eNev%W_free, fak)
      end if

      if (eNev%W_rec > 0) then
         call AddHist(hSigmaMC_Wrec(0),hSigmaMC_Wrec(iHist),&
              eNev%W_rec, fak)
      end if

      Eprime = eNev%lepton_out%mom(0)
      cost = eNeV_Get_CostLepton(eNev)
      plept = abs3(eNev%lepton_out%mom)
      pL = plept * cost
      pT = plept * sqrt(1. - cost**2)
      Wrec = eNev%W_rec

      call AddHist2D(hSigmaMC_EprimeCost(0),hSigmaMC_EprimeCost(iHist),&
           (/Eprime,cost/), fak)

      call AddHist2D(hSigmaMC_pLpT(0),hSigmaMC_pLpT(iHist),&
           (/pL,pT/), fak)

      call AddHist2D(hSigmaMC_pLW(0),hSigmaMC_pLW(iHist),&
           (/pL,Wrec/), fak)

      call AddHist2D(hSigmaMC_pTW(0),hSigmaMC_pTW(iHist),&
           (/pT,Wrec/), fak)

      if (.not.MCmode) then
         call CalcXY(eNev%lepton_in%mom, &
              eNev%lepton_out%mom,&
              eNev%nucleon%mom, MC_x,MC_y, flagDUMMY)
      end if

      call AddHist2D(hSigmaMC_XY(0),hSigmaMC_XY(iHist),&
           (/MC_x,MC_y/), fak)
      call AddHist(hSigmaMC_X(0),hSigmaMC_X(iHist),&
           MC_x, fak)



    end subroutine fillHist

    !**************************************************************************
    !****o* initNeutrino/neutrino.NuQ2planeXS.ZZZ.dat
    ! NAME
    ! file neutrino.NuQ2planeXS.ZZZ.dat,   ZZZ=000 - 008
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! nu(transferred energy)-Q^2 plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2/GeV^3 for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV^3  for EM
    ! All cross sections are given per nucleon (1/A)
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! ZZZ = ...:
    ! * 000: sum over all channels
    ! * 001: QE cross section
    ! * 002: Delta
    ! * 003: highRES
    ! * 004: 1pi
    ! * 005: DIS
    ! * 006: 2p2h QE
    ! * 007: 2p2h Delta
    ! * 008: 2pi background
    !
    ! Columns:
    ! * #1: nu (transfered energy),  GeV
    ! * #2: Q^2, GeV^2
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.Q2NuplaneXS.ZZZ.dat
    ! NAME
    ! file neutrino.Q2NuplaneXS.ZZZ.dat
    ! PURPOSE
    ! Nearly the same as  neutrino.NuQ2planeXS.ZZZ.dat, but:
    ! * #1: Q^2, GeV^2
    ! * #2: nu,  GeV
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.EprimeCostplaneXS.ZZZ.dat
    ! NAME
    ! file neutrino.EprimeCostplaneXS.ZZZ.dat,   ZZZ=000 - 008
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! Eprime(energy of the outgoing lepton)-costheta(cos of the angle between
    ! the incoming and outgoing leptons) plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2/GeV for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV  for EM
    ! All cross sections are given per nucleon (1/A)
    !
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! ZZZ is the origin (the first interaction vertex) of the event:
    ! * 000: sum over all origins
    ! * 001: QE cross section
    ! * 002: Delta
    ! * 003: highRES
    ! * 004: 1pi
    ! * 005: DIS
    ! * 006: 2p2h QE
    ! * 007: 2p2h Delta
    ! * 008: 2pi background
    !
    ! Columns:
    ! * #1: Eprime (energy of the outgoing lepton),  GeV
    ! * #2: costheta(cos of the angle between the incoming and outgoing leptons)
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.pLpTXS.ZZZ.dat
    ! NAME
    ! file neutrino.EprimepLpTXS.ZZZ.dat,   ZZZ=000 - 008
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! longitudinal momentum  of the outgoing lepton
    ! - transverse momentum of the outgoing lepton plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2/GeV for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV  for EM
    ! All cross sections are given per nucleon (1/A)
    !
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! ZZZ is the origin (the first interaction vertex) of the event:
    ! * 000: sum over all origins
    ! * 001: QE cross section
    ! * 002: Delta
    ! * 003: highRES
    ! * 004: 1pi
    ! * 005: DIS
    ! * 006: 2p2h QE
    ! * 007: 2p2h Delta
    ! * 008: 2pi background
    !
    ! Columns:
    ! * #1: longitudinal momentum of the outgoing lepton,  in GeV
    ! * #2: transverse momentum of the outgoing lepton,  in GeV
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************


    !**************************************************************************
    !****o* initNeutrino/neutrino.XYplaneXS.ZZZ.dat
    ! NAME
    ! file neutrino.XYplaneXS.ZZZ.dat
    ! ZZZ=000 - 008  is the origin (the first interaction vertex) of the event:
    !                (see description in  neutrino.EprimeCostplaneXS.ZZZ.dat)
    ! PURPOSE
    ! This files provides absorption cross section in the
    ! x_{Bjorken}-y(=nu/Enu) plane.
    ! Format is suitable as input to gnuplots "splot" command
    !
    ! The units are 10^{-38} cm^2 for process_ID=CC, NC and
    ! 10^{-33} cm^2  for EM
    ! All cross sections are given per nucleon (1/A)
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! Columns:
    ! * #1: x_Bjorken
    ! * #2: y (=nu/Enu)
    ! * #3: cross section
    ! * #4: number of entries that lead to cross section
    ! * #5: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.NuXS.ZZZ.dat
    ! NAME
    ! file neutrino.NuXS.ZZZ.dat
    ! ZZZ=000 - 008  is the origin (the first interaction vertex) of the event:
    !                (see description in  neutrino.EprimeCostplaneXS.ZZZ.dat)
    ! PURPOSE
    ! This files provides absorption cross section versus nu(=transfered energy)
    ! Format is suitable as input to gnuplots "splot" command
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! The units are 10^{-38} cm^2/GeV for process_ID=CC, NC and
    ! 10^{-33} cm^2/GeV  for EM
    ! All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: nu (transfered energy), GeV
    ! * #2: cross section
    ! * #3: number of entries that lead to cross section
    ! * #4: should be ignored
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.Q2XS.ZZZ.dat
    ! NAME
    ! file neutrino.NuQ2.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus Q2
    !
    ! Coulmns:
    ! * #1: Q2 , GeV^2
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.XXS.ZZZ.dat
    ! NAME
    ! file neutrino.XXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus x_Bjorken
    !
    ! Columns:
    ! * #1: x_Bjorken
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.XrecXS.ZZZ.dat
    ! NAME
    ! file neutrino.XrecXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus reconstructed x_Bjorken,
    ! defined for nucleon at rest: Xrec = Q^2/(2*mN*nu)
    !
    ! Columns:
    ! * #1: free x_Bjorken
    !**************************************************************************

    !**************************************************************************
    !****o* initNeutrino/neutrino.WrecXS.ZZZ.dat
    ! NAME
    ! file neutrino.WrecXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus Wrec=qsrt(mN^2+2*mN*nu-Q2) ,
    ! where mN=0.938 GeV = nucleon mass (see constants.f90), i.e. W
    ! reconstructed for nucleon at rest
    !
    ! Columns:
    ! * #1: Wrec , GeV
    !**************************************************************************


    !**************************************************************************
    !****o* initNeutrino/neutrino.WfreeXS.ZZZ.dat
    ! NAME
    ! file neutrino.WfreeXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus Wfree, i.e. W value for
    ! boson with four-momentum q and Fermi-moving nucleon, without potential
    !
    ! Columns:
    ! * #1: Wfree , GeV
    !**************************************************************************


    !**************************************************************************
    !****o* initNeutrino/neutrino.WXS.ZZZ.dat
    ! NAME
    ! file neutrino.WXS.ZZZ.dat
    ! PURPOSE
    ! Similar to neutrino.NuXS.ZZZ.dat, but versus W, i.e. W value for
    ! boson with four-momentum q and Fermi-moving nucleon in the potential
    !
    ! Columns:
    ! * #1: Wfree , GeV
    !**************************************************************************

    !**************************************************************************
    !****is* init_neutrino/writeInclHist
    !**************************************************************************
    subroutine writeInclHist

      integer :: iHist

      do iHist=0,max_Hist
         if (.not.includeHist(iHist)) cycle
         call WriteHist2D_Gnuplot(hSigmaMC_nuQ2(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.NuQ2planeXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist2D_Gnuplot(hSigmaMC_nuQ2(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.Q2NuplaneXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false., SwapXY=.true.)
         call WriteHist2D_Gnuplot(hSigmaMC_EprimeCost(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.EprimeCostplaneXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist2D_Gnuplot(hSigmaMC_pLpT(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.pLpTXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist2D_Gnuplot(hSigmaMC_pLW(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.pLWrecXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist2D_Gnuplot(hSigmaMC_pTW(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.pTWrecXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist2D_Gnuplot(hSigmaMC_XY(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.XYplaneXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_nu(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.NuXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_Q2(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.Q2XS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_qz(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.qzXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_X(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.XXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_Xrec(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.XrecXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_Wrec(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.WrecXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_Wfree(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.WfreeXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
         call WriteHist(hSigmaMC_W(iHist),10, &
              & mul=1.0/numberofcalls, add=1e-20,&
              & file='neutrino.WXS.'//trim(intToChar(iHist))//'.dat',&
              & dump=.false.)
      end do


    end subroutine writeInclHist

  end subroutine init_neutrino

  !****************************************************************************
  !****************************************************************************
  integer function getFirstEvent()

    first=first+1
    getFirstEvent=first
  end function getFirstEvent

  !****************************************************************************
  !****************************************************************************
  function getFirstEventRange()

    integer, dimension (1:2) :: getFirstEventRange
    getFirstEventRange=(/1,first/)
  end function getFirstEventRange

  !****************************************************************************
  !****************************************************************************
  subroutine getNeutrinoInfo(raiseVal_)

    real, intent(out) :: raiseVal_
    raiseVal_=raiseVal
  end subroutine getNeutrinoInfo

  !****************************************************************************
  !****************************************************************************
  logical function neutrinoInit_getRealRun()

    if (readinputflag) then
       call readInput
    end if
    neutrinoInit_getRealRun=realRun
  end function neutrinoInit_getRealRun

  !****************************************************************************
  !****s* initNeutrino/get_init_namelist
  ! NAME
  ! subroutine get_init_namelist(process_ID, flavor_ID,
  ! nuXsectionMode, nuExp, debugflag, includeQE, includeDELTA, includeRES,
  ! include1pi, realRun)
  !
  ! PURPOSE
  ! This subroutine returns any entry of the neutrino init namelist.
  !
  ! OUTPUT
  ! * logical, optional :: debugflag,includeQE,includeDELTA,
  !   includeRES,include1pi,realRun
  ! * integer, optional :: process_ID,flavor_ID,nuXsectionMode,nuExp
  !
  !****************************************************************************
  subroutine get_init_namelist(Gprocess_ID,Gflavor_ID, &
       & GnuXsectionMode, &
       & GnuExp,Gdebugflag,GincludeQE,GincludeDELTA,GincludeRES,Ginclude1pi,&
       & GrealRun, outLepton_ID, outLepton_charge)

    logical, optional, intent(out) :: Gdebugflag, &
         & GincludeQE,GincludeDELTA,GincludeRES,Ginclude1pi,GrealRun
    integer, optional, intent(out) :: Gprocess_ID,Gflavor_ID, &
         & GnuXsectionMode, GnuExp, outLepton_ID, outLepton_charge

    integer :: leptonout_ID, leptonout_charge

    if (present(Gflavor_ID)) Gflavor_ID=flavor_ID
    if (present(Gprocess_ID)) Gprocess_ID=process_ID
    if (present(GnuXsectionMode)) GnuXsectionMode=nuXsectionMode
    if (present(GnuExp)) GnuExp=nuExp
    if (present(Gdebugflag)) Gdebugflag=debugflag
    if (present(GincludeQE)) GincludeQE=includeQE
    if (present(GincludeDELTA)) GincludeDELTA=includeDELTA
    if (present(GincludeRES)) GincludeRES=includeRES
    if (present(Ginclude1pi)) Ginclude1pi=include1pi
    if (present(GrealRun)) GrealRun=realRun


     if (abs(process_ID).eq.2 .or. abs(process_ID).eq.1) then
          leptonout_ID= 900+flavor_ID
          leptonout_charge=-sign(1, process_ID)
     else if (abs(process_ID).eq.3) then
          leptonout_ID= sign(910+flavor_ID, process_ID)
          leptonout_charge=0
     end if

    if (present(outLepton_ID)) outLepton_ID=leptonout_ID
    if (present(outLepton_charge)) outLepton_charge=leptonout_charge


  end subroutine get_init_namelist

  !****************************************************************************
  !****s* initNeutrino/get_runtime_vars
  ! NAME
  ! subroutine get_runtime_vars(sigabsArrFinal,sigabsArr)
  !
  ! PURPOSE
  ! This subroutine returns variables that are changed with every init.
  !
  ! OUTPUT
  ! * real,dimension(0:max_finalstate_ID),optional :: sigabsArrFinal, sigabsArr
  !
  !****************************************************************************
  subroutine get_runtime_vars(GsigabsArrFinal,GsigabsArr)

    real,dimension(0:max_finalstate_ID),optional,intent(out) :: GsigabsArrFinal,GsigabsArr
    if (present(GsigabsArrFinal)) GsigabsArrFinal=sigabsArrFinal
    if (present(GsigabsArr)) GsigabsArr=sigabsArr
  end subroutine get_runtime_vars

  !****************************************************************************
  !****f* initNeutrino/isOsc
  ! NAME
  ! logical function isOsc
  !
  ! PURPOSE
  ! return true, if nuExp points to a neutrino oscllation experiment
  !****************************************************************************
  logical function isOsc()
    if (nuExp.le.numberOfExperiments) then
       isOsc = Osc(nuExp)
    else
       isOsc = .FALSE.
    end if
  end function isOsc

  !****************************************************************************
  !****f* initNeutrino/getNuMax
  ! NAME
  ! real function getNuMax
  !
  ! PURPOSE
  ! return some estimate of a "reasonable" numax
  !****************************************************************************
  real function getNuMaxExp()

    use expNeutrinofluxes, only: userFlux, userFluxEMax

    ! numbers given in array nuMaxArr are upper boundaries for energy-transfer
    real, dimension(1:numberOfExperiments), parameter :: nuMaxArr = (/&
         & 2.5, 3.0, 4.0, 5.0, 2.5, & ! MiniBooNE-nu,ANL,K2K,BNL,MiniBooNE-barnu
         & 30.0, 30.0, 15.0, 5.0, &   ! MINOS-nu-numode, MINOS-barnu-numode, NOvA, T2K OA2.5,
         & 3.0, 30.0,  &              ! uniform distr, MINOS-nu-barnumode,
         & 30.0,      &               ! MINOS-barnu-barnumode
         & 30.0, 30.0, 30.0, 30.0, &  ! MINERvA-numu, MINERvA-antinumu,DUNE-nu,DUNE-barnu
         & 30.0, 300.0,      &        ! LBNO-nu, NOMAD
         & 7.5, 7.5, 7.5, 7.5, &      ! BNB-nue,BNB-nuebar,BNBnumu,BNBnumubar
         & 15., 20., 20. /)           ! NOvA, T2Koa, MINERvA 2016



    getNuMaxExp = 0.
    if (nuExp.gt.0) then
       if (nuExp.le.numberOfExperiments) then
          getNuMaxExp = nuMaxArr(nuExp)
       else
          getNuMaxExp = userFlux(FileNameFlux) ! dummy
          getNuMaxExp = userFluxEmax()
       end if
    end if

  end function getNuMaxExp


end module initNeutrino
