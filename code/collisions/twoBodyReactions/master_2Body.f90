!******************************************************************************
!****m* /master_2Body
! NAME
! module master_2Body
! PURPOSE
! Implements all 2-body processes (a b -> X).
!******************************************************************************
module master_2Body

  use CallStack, only: Traceback

  implicit none
  private

  !****************************************************************************
  !****g* master_2Body/debug
  ! SOURCE
  logical, parameter :: debug=.false.
  ! PURPOSE
  ! Switch for debug information
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/useHiEnergy
  ! SOURCE
  !
  logical, save :: useHiEnergy = .true.
  ! PURPOSE
  ! Switch to turn HiEnergy on/off. Formerly known as "useFritiof".
  !
  ! NOTES
  ! Please be very sure what you are doing when setting this parameter
  ! to .false.!
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/HiEnergyThresholdBarMes
  ! SOURCE
  !
  real, save ::  HiEnergyThresholdBarMes = 2.2
  ! PURPOSE
  ! Sqrt(s) threshold for HiEnergy in Baryon-Meson Reactions
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/HiEnergyThresholdBarMesDelta
  ! SOURCE
  !
  real, save ::  HiEnergyThresholdBarMesDelta = 0.2
  ! PURPOSE
  ! width for the Sqrt(s) threshold for HiEnergy in Baryon-Meson Reactions
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/HiEnergyThresholdBarBar
  ! SOURCE
  !
  real, save ::  HiEnergyThresholdBarBar = 3.4
  ! PURPOSE
  ! Sqrt(s) threshold for HiEnergy in Baryon-Baryon Reactions
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/HiEnergyThresholdBarBarDelta
  ! SOURCE
  !
  real, save :: HiEnergyThresholdBarBarDelta = 0.1
  ! PURPOSE
  ! width for the Sqrt(s) threshold for HiEnergy in Baryon-Baryon Reactions
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/HiEnergyThresholdBarAntibar
  ! SOURCE
  !
  real, save ::  HiEnergyThresholdBarAntibar = 2.38
  ! PURPOSE
  ! Sqrt(s) threshold for HiEnergy in Baryon-Antibaryon Reactions
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/HiEnergyThresholdBarAntibarDelta
  ! SOURCE
  !
  real, save :: HiEnergyThresholdBarAntibarDelta = 0.0
  ! PURPOSE
  ! width for the Sqrt(s) threshold for HiEnergy in Baryon-Antibaryon Reactions
  !****************************************************************************


  !****************************************************************************
  !****g* master_2Body/xSectionCutOff
  ! SOURCE
  !
  real, parameter, public :: xSectionCutOff=1.E-10
  ! PURPOSE
  ! No reaction takes place if Xsection is below this value! Units: mB
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/correctEnergy
  ! SOURCE
  !
  logical,save :: correctEnergy=.true.
  ! PURPOSE
  ! Scale final state momenta to fulfill energy and momentum conservation.
  ! If .false. energy conservation is violated
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/correctEnergy_message
  ! SOURCE
  !
  logical,save :: correctEnergy_message=.true.
  ! PURPOSE
  ! Switch off error message for energy correction failures.
  !****************************************************************************


  !****************************************************************************
  !****g* master_2Body/usePythia
  ! SOURCE
  !
  integer,save :: usePythia = 1
  ! PURPOSE
  ! This flag decides whether to use Fritiof or Pythia for high-energy
  ! collisions:
  ! * 0: use Fritiof
  ! * 1: use Pythia
  ! NOTES
  ! * This flag is not used in the baryon-antibaryon channel
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/usePythia_BaB
  ! SOURCE
  !
  integer,save :: usePythia_BaB = 0
  ! PURPOSE
  ! This flag decides whether to use Fritiof or Pythia for high-energy
  ! baryon-antibaryon collisions:
  ! * 0: use Fritiof
  ! * 1: use Pythia
  !****************************************************************************


  !****************************************************************************
  !****g* master_2Body/useManni
  ! SOURCE
  !
  logical,save :: useManni = .true.
  ! PURPOSE
  ! Flag, whether to use meson-baryon annhilation as proposed by
  ! Markus Wagner (Diploma, Giessen 2004), but with some enhanced treatment
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/ElastAngDist
  ! SOURCE
  !
  integer,save :: ElastAngDist = 3
  ! PURPOSE
  ! Choice of angular distribution in (high-energy) elastic collisions
  ! (cf. DoColl_Elast):
  ! * 1 = isotropic
  ! * 2 = J. Cugnon et al., NPA 352, 505 (1981)
  ! * 3 = Pythia (default)
  !****************************************************************************


  !----------------------------------------------------------------------------
  ! Switches to turn on/off channels:
  !----------------------------------------------------------------------------

  !****************************************************************************
  !****g* master_2Body/baryonBaryonScattering
  ! SOURCE
  !
  logical, save :: baryonBaryonScattering=.true.
  ! PURPOSE
  ! Switch to turn off baryon-baryon-Scattering
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/baryonMesonScattering
  ! SOURCE
  !
  logical, save :: baryonMesonScattering=.true.
  ! PURPOSE
  ! Switch to turn off baryon-Meson-Scattering
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/mesonMesonScattering
  ! SOURCE
  !
  logical, save :: mesonMesonScattering=.true.
  ! PURPOSE
  ! Switch to turn off meson-Meson-Scattering
  !****************************************************************************

  !----------------------------------------------------------------------------
  ! Cut-Off parameter
  !----------------------------------------------------------------------------

  !****************************************************************************
  !****g* master_2Body/coarse
  ! SOURCE
  !
  real, dimension(1:3), save :: coarse=(/3.,4.,4./) !baryonBaryon, baryonMeson, mesonMeson
  ! PURPOSE
  ! coarse maximal impact parameter (in fm)
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/bmax_nucleonNucleon
  ! SOURCE
  !
  real, save :: bmax_nucleonNucleon=2.52
  ! PURPOSE
  ! Real maximal impact parameter for nucleon-nucleon-scattering.
  ! Maximal crossection is
  !   bMax**2 * pi * 10 mb/fm**2 = (2.52**2*pi*10) mb  = 199.5 mb
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/bmax_nucleonResonance
  ! SOURCE
  !
  real, save :: bmax_nucleonResonance=1.60
  ! PURPOSE
  ! Real maximal impact parameter for nucleon-resonance scattering.
  ! Maximal crossection is
  !   bMax**2 * pi * 10 mb/fm**2 = (1.60**2*pi*10) mb  = 80.4 mb
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/bmax_hyperonNucleon
  ! SOURCE
  !
  real, save :: bmax_hyperonNucleon=2.52
  ! PURPOSE
  ! Real maximal impact parameter for hyperon-nucleon-scattering.
  ! Maximal crossection is
  !   bMax**2 * pi * 10 mb/fm**2 = (2.52**2*pi*10) mb  = 199.5 mb
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/bmax_baryonPion
  ! SOURCE
  !
  real, save :: bmax_baryonPion=2.52
  ! PURPOSE
  ! real maximal impact parameter for baryon pion scattering
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/bmax_baryonMeson
  ! SOURCE
  !
  real, save :: bmax_baryonMeson=2.52
  ! PURPOSE
  ! real maximal impact parameter for baryon-Meson-scattering
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/bmax_mesonMeson
  ! SOURCE
  !
  real, save :: bmax_mesonMeson=2.
  ! PURPOSE
  ! real maximal impact parameter for meson-meson-scattering
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/flagElastBB
  ! SOURCE
  !
  logical, save :: flagElastBB = .false.
  ! PURPOSE
  ! If .true., use a constant elastic baryon-baryon cross section of 40 mb
  ! and no inelastic baryon-baryon scattering.
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/OverideSigma_PiN
  ! SOURCE
  !
  real, save :: OverideSigma_PiN=-99.9
  ! PURPOSE
  ! Parameter to replace the calculated cross section for pi+N collision
  ! by a fixed value (in mb). Only in use if >= 0.
  !
  ! The elastic cross section is assumed to be 1/10 of the given value.
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/OverideSigma_RhoN
  ! SOURCE
  !
  real, save :: OverideSigma_RhoN=-99.9
  ! PURPOSE
  ! Parameter to replace the calculated cross section for rho+N collision
  ! by a fixed value (in mb). Only in use if >= 0.
  !
  ! The elastic cross section is assumed to be 1/10 of the given value.
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/OverideSigma_PiPi
  ! SOURCE
  !
  real, save :: OverideSigma_PiPi=-99.9
  ! PURPOSE
  ! Parameter to replace the calculated cross section for pi+pi collision
  ! by a fixed value (in mb). Only in use if >= 0.
  !
  ! We set sigma_elast = sigma_tot
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/Overide_PiPi_ResIsElast
  ! SOURCE
  !
  logical, save :: Overide_PiPi_ResIsElast=.false.
  ! PURPOSE
  ! Flag to replace the calculated cross section for pi+pi collision;
  ! The calculated resonant cross section will be transformed into the
  ! elastic cross section. Thus no resonances will be propagated explicitely,
  ! but they show up in the cross section
  !
  ! We set sigma_elast = sigma_Res, sigma_Res = 0, sigma_tot = sigma_elast
  !
  ! please note: background processes as pi pi <-> K K~ are *not* affected
  ! by this switch. You have to disable those additionally by hand,
  ! see mesMes_do2to2
  !****************************************************************************


  !****************************************************************************
  !****g* master_2Body/omega_K_factor
  ! SOURCE
  !
  real, save :: omega_K_factor = 2.
  ! PURPOSE
  ! Modification factor for the inelastic omega-nucleon cross section.
  ! Necessary to describe transpacrency ratio data measured by CBELSA/TAPS,
  ! see: http://arxiv.org/abs/1210.3074
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/mesMes_do2to2
  ! SOURCE
  !
  logical, save :: mesMes_do2to2 = .true.
  ! PURPOSE
  ! flag whether to do m m' <-> K K~, K K*~ etc.
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/mesMes_useWidth
  ! SOURCE
  !
  logical, save :: mesMes_useWidth = .false.
  ! PURPOSE
  ! flag whether to use the width in m m' <-> K K~, K K*~ etc.
  ! This is needed to enforce detailed balance. Otherwise only pole
  ! masses are used.
  !****************************************************************************

  !****************************************************************************
  !****g* master_2Body/doScaleResidue
  ! SOURCE
  !
  logical, save :: doScaleResidue = .true.
  ! PURPOSE
  ! scale the cross section of real-pert collisions by a factor N'/N or Z'/Z
  ! for a scattering on a neutron or proton, where N' and Z' are calculated
  ! via the residuum.
  !****************************************************************************


  logical, save :: initFlag=.true.

  public :: collide_2Body
  public :: XsectionMaster
  public :: generateFinalState
  public :: check_veryRough, checkKodama_Rough
  public :: setKinematics, setKinematicsHiEnergy
  public :: HiEnergyContrib


  ! Just to make the code more readable:

  integer,parameter :: scenarioBarBar=1
  integer,parameter :: scenarioBarMes=2
  integer,parameter :: scenarioMesMes=3

  integer,parameter :: scenarioBarAntiBar=2
  integer,parameter :: scenarioAntiBarAntiBar=3

  integer,parameter :: fehler_max = 100000

  integer,parameter :: isigTot   = 0  ! total XS
  integer,parameter :: isigElast = 1  ! elastic XS
  integer,parameter :: isigCEX   = 2  ! charge exchange XS
  integer,parameter :: isigAnni  = 3  ! annihilation XS
  integer,parameter :: isigLbar  = 4  ! Lambda LambdaBar excl. XS
  integer,parameter :: isigSbar  = 5  ! Lambda SigmaBar + LambdaBar Sigma
  integer,parameter :: isigXibar = 6  ! Xi XiBar excl. XS
  integer,parameter :: isigJPsi  = 7  ! JPsi excl. XS

  integer,parameter :: isigRes   = 2  ! resonant XS
  integer,parameter :: isigBG    = 3  ! background XS




contains

  !****************************************************************************
  !****s* master_2Body/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "master_2Body"
  !****************************************************************************
  subroutine readInput

    use output, only: Write_ReadingInput
    use hadronFormation, only: forceInitFormation
    use twoBodyPhaseSpace, only: setMaxSqrts_2bps => setMaxSqrts
    use dimi, only: setMaxSqrts_dimi => setMaxSqrts
    use mesonMeson_Xsections, only: mesMes_Tabulate

    !**************************************************************************
    !****n* master_2Body/master_2body
    ! NAME
    ! NAMELIST master_2Body
    ! PURPOSE
    ! Includes the switches:
    ! * correctEnergy
    ! * baryonBaryonScattering
    ! * baryonMesonScattering
    ! * mesonMesonScattering
    ! * usePythia
    ! * usePythia_BaB
    ! * useHiEnergy
    ! * HiEnergyThresholdBarMes
    ! * HiEnergyThresholdBarMesDelta
    ! * HiEnergyThresholdBarBar
    ! * HiEnergyThresholdBarBarDelta
    ! * HiEnergyThresholdBarAntibar
    ! * HiEnergyThresholdBarAntibarDelta
    ! * useManni
    ! * ElastAngDist
    ! * flagElastBB
    ! * coarse
    ! * bmax_nucleonNucleon
    ! * bmax_nucleonResonance
    ! * bmax_hyperonNucleon
    ! * bmax_baryonPion
    ! * bmax_baryonMeson
    ! * bmax_mesonMeson
    ! * correctEnergy_message
    ! * OverideSigma_PiN
    ! * OverideSigma_RhoN
    ! * OverideSigma_PiPi
    ! * Overide_PiPi_ResIsElast
    ! * omega_K_factor
    ! * mesMes_do2to2
    ! * mesMes_useWidth
    ! * doScaleResidue
    !**************************************************************************
    NAMELIST /master_2Body/ &
         correctEnergy, baryonBaryonScattering, baryonMesonScattering, &
         mesonMesonScattering, usePythia, usePythia_BaB, &
         useHiEnergy, HiEnergyThresholdBarMes, HiEnergyThresholdBarMesDelta, &
         HiEnergyThresholdBarBar, &
         HiEnergyThresholdBarBarDelta, HiEnergyThresholdBarAntibar, &
         HiEnergyThresholdBarAntibarDelta, &
         useManni, ElastAngDist, flagElastBB, coarse, &
         bmax_nucleonNucleon, bmax_nucleonResonance, bmax_hyperonNucleon, &
         bmax_baryonPion, &
         bmax_baryonMeson, bmax_mesonMeson, correctEnergy_message, &
         OverideSigma_PiN, OverideSigma_RhoN, OverideSigma_PiPi,&
         Overide_PiPi_ResIsElast, &
         omega_K_factor, &
         mesMes_do2to2, mesMes_useWidth, doScaleResidue

    character(60), dimension(1:3), parameter :: NNe = (/ &
         'isotropic                            ', &
         'J. Cugnon et al., NPA 352, 505 (1981)', &
         'Pythia (default)                     ' /)

    character(12), dimension(0:1), parameter :: NNp = (/ &
         'use Fritiof', &
         'use Pythia ' /)
    integer :: ios

    call Write_ReadingInput('master_2Body',0)
    rewind(5)
    read(5,nml=master_2Body,IOSTAT=ios)
    call Write_ReadingInput('master_2Body',0,ios)

    write(*,*) 'Correct energy by momentum scaling after collisions= ', &
         correctEnergy
    write(*,*) 'Baryon-Baryon-Scattering = ',baryonBaryonScattering
    write(*,*) 'Baryon-Meson-Scattering  = ',baryonMesonScattering
    write(*,*) 'Meson-Meson-Scattering   = ',mesonMesonScattering

    if (.not.correctEnergy_message) write(*,'(/,A,/)') &
         ' WARNING: You have switched off energy correction failure message in case of threshold problems!'
    write(*,*) 'use HiEnergy                     : ',useHiEnergy
    write(*,'(A,f5.2," +-",f5.2," GeV")') ' srts threshold bar-mes    :', &
         HiEnergyThresholdBarMes, HiEnergyThresholdBarMesDelta
    write(*,'(A,f5.2," +-",f5.2," GeV")') ' srts threshold bar-bar    :', &
         HiEnergyThresholdBarBar, HiEnergyThresholdBarBarDelta
    write(*,'(A,f5.2," +-",f5.2," GeV")') ' srts threshold bar-antibar:', &
         HiEnergyThresholdBarAntibar, HiEnergyThresholdBarAntibarDelta

    write(*,*) 'use Manni                        : ',useManni

    if (usePythia<0 .or. usePythia>1) then
       write(*,'(A,i3," = ",A)') ' use PYTHIA directly              : ',&
            usePythia,' !!!! STOP'
       call Traceback()
    end if
    write(*,'(A,i3," = ",A)') ' use PYTHIA directly              : ',&
         usePythia,trim(NNp(usePythia))

    if (ElastAngDist<1.or.ElastAngDist>3) then
       write(*,'(A,i3," = ",A)') ' Elastic Angular distribution     : ',&
            ElastAngDist,' !!!! STOP'
       call Traceback()
    end if
    write(*,'(A,i3," = ",A)') ' Elastic Angular distribution     : ',&
         ElastAngDist,trim(NNe(ElastAngDist))

    if (usePythia_BaB<0 .or. usePythia_BaB>1) then
       write(*,'(A,i3," = ",A)') ' BaB: use PYTHIA directly         : ',&
            usePythia_BaB,' !!!! STOP'
       call Traceback()
    end if
    write(*,'(A,i3," = ",A)') ' BaB: use PYTHIA directly         : ',&
         usePythia_BaB,trim(NNp(usePythia_Bab))

    write(*,*)
    write(*,'(A,3(1x,f5.2))') &
         ' coarse maximum distances for bar-bar, bar-mes, mes-mes, fm : ',&
         coarse
    write(*,'(A,f5.2)') &
         ' max. impact parameter for nucleon-nucleon scatt.,  fm :', &
         bmax_nucleonNucleon
    write(*,'(A,f5.2)') &
         ' max. impact parameter for nucleon-resonace scatt., fm :', &
         bmax_nucleonResonance
    write(*,'(A,f5.2)') &
         ' max. impact parameter for hyperon-nucleon scatt.,  fm :', &
         bmax_hyperonNucleon
    write(*,'(A,f5.2)') &
         ' max. impact parameter for baryon-pion scatt.,      fm :', &
         bmax_baryonPion
    write(*,'(A,f5.2)') &
         ' max. impact parameter for baryon-meson scatt.,     fm :', &
         bmax_baryonMeson
    write(*,'(A,f5.2)') &
         ' max. impact parameter for meson-meson scatt.,      fm :', &
         bmax_mesonMeson
    write(*,*) 'flagElastBB: ', flagElastBB

    if (OverideSigma_PiN.ge.0.0) then
       write(*,'(/,A,g12.3,A,/)') ' !!!!! YOU HAVE SET PION NUCLEON CROSS SECTION BY HAND: sigma = ',&
          & OverideSigma_PiN,' mb !!!!!'
    end if
    if (OverideSigma_RhoN.ge.0.0) then
       write(*,'(/,A,g12.3,A,/)') ' !!!!! YOU HAVE SET RHO NUCLEON CROSS SECTION BY HAND:  sigma = ',&
          &  OverideSigma_RhoN,' mb !!!!!'
    end if

    if (OverideSigma_PiPi.ge.0.0) then
       write(*,'(/,A,g12.3,A,/)') ' !!!!! YOU HAVE SET PION PION CROSS SECTION BY HAND: sigma = ',&
          &  OverideSigma_PiPi,' mb !!!!!'
    end if

    if (Overide_PiPi_ResIsElast) then
       write(*,'(/,A,g12.3,A,/)') ' !!!!! YOU HAVE SET RESONANT PION PION CROSS SECTION TO BE ELASTIC !!!!!'
    end if

    write(*,*) 'omega_K_factor = ', omega_K_factor
    write(*,*)
    write(*,*) "meson-meson: do m m' <-> K K~, K K*~ etc.: ",mesMes_do2to2
    write(*,*) "meson-meson: use width                   : ",mesMes_useWidth
    write(*,*) "rescale of real-pert collisions via residuum        :",doScaleResidue

    call Write_ReadingInput('master_2Body',1)


    initFlag=.false.

    call forceInitFormation
    if (mesMes_useWidth) then
       call mesMes_Tabulate
    end if

    ! set tabulation bound of two-body phase space to highest energy reachable
    ! in low-energy mode (resonance model)
    call setMaxSqrts_2bps (HiEnergyThresholdBarBar + HiEnergyThresholdBarBarDelta)
    call setMaxSqrts_dimi (HiEnergyThresholdBarBar + HiEnergyThresholdBarBarDelta)

  end subroutine readInput


  !****************************************************************************
  !****s* master_2Body/collide_2body
  ! NAME
  ! subroutine collide_2body(partIn, partOut, time, collFlag,
  ! HiEnergyFlag, HiEnergyType,
  ! collTime, weightLocal, sigTot_out, pauliIncluded_out)
  !
  ! PURPOSE
  ! Evaluates for two given particles the final states. First it checks the
  ! collision criteria, then it evaluates the final state if the collision
  ! criteria is fulfilled.
  !
  ! INPUTS
  ! * type(particle),dimension(:) :: partIn  -- incoming pair of particles
  ! * type(particle),dimension(:) :: partOut -- produced final state particles
  ! * real                        :: time    -- time in fm
  !
  ! OUTPUT
  ! * type(particle),dimension(:) :: partOut -- produced final state particles
  ! * logical :: collFlag     -- true if collision takes place
  ! * logical :: HiEnergyFlag -- true if HiEnergy was used
  ! * integer :: HiEnergyType -- 0:LowEnergy, 1:Fritiof, 2:Pythia
  ! * real, optional  :: collTime -- Time instant in the
  !   computational frame when the two particles collide
  !   (with respect to the current "BUU time" used for time stepping)
  ! * integer,  optional  :: weightLocal -- Weight used for local ensemble
  !   method for reweighting the probability in the collision criteria.
  !   Not necessary in case of localEnsemble=.false.
  ! * real, optional    :: sigTot_out  -- total cross section
  ! * logical, optional :: pauliIncluded_out  -- true if Pauli blocking was
  !   already included in the cross section evaluation
  !
  ! NOTES
  ! "partOut" has to be given as input with some values set properly!
  !****************************************************************************
  subroutine collide_2body(partIn, partOut, time, collFlag, &
       HiEnergyFlag, HiEnergyType, &
       collTime, weightLocal, sigTot_out, pauliIncluded_out)

    use IdTable, only: isMeson
    use particleDefinition
    use collisionCriteria, only: kodama_time
    use inputGeneral, only: delta_T,localEnsemble,fullensemble,numEnsembles
    use twoBodyTools, only: sqrtS_free, LogicMatrix
    use history, only: setHistory
    use residue, only: ResidueGetScale

    type(particle),dimension(:),intent(in)    :: partIn
    type(particle),dimension(:),intent(inout) :: partOut
    real, intent(in)                          :: time

    logical,intent(out)                       :: collFlag, HiEnergyFlag
    integer,intent(out)                       :: HiEnergyType
    real, optional,intent(out)                :: collTime
    integer, optional,intent(in)              :: weightLocal
    real, optional,intent(out)                :: sigTot_out
    logical, optional, intent(out)            :: pauliIncluded_out

    ! scaling factor of cross section due to formation time effects:
    real :: stringFactor
    real :: sqrtS_vacuum, XShelp
    integer :: scenario
    logical :: pauliIncluded

    collFlag=.false.
    pauliIncluded=.false.

    if (present(pauliIncluded_out)) then
       pauliIncluded_out=.false.
    end if

    ! Define Sqrt(s) in the vacuum
    sqrtS_vacuum = sqrtS_free(partIn)
    if (sqrtS_vacuum < Sum(partIn%mass)) then
       write(*,*) '***************************'
       write(*,*) 'Error in master_2body: sqrtS_vacuum.lt.Sum(partIn%mass): ',&
            'sqrtS_vacuum',  sqrtS_vacuum,Sum(partIn%mass)
       write(*,*) partIn%mass, partIn%ID, partIn%charge
       write(*,*) partIn%pert, partIn%anti
       write(*,*) partIn(1)%mom, partIn(1)%vel
       write(*,*) partIn(2)%mom, partIn(2)%vel
       write(*,*) 'stopping'
       call Traceback()
    else if (sqrtS_vacuum-Sum(partIn%mass) < 1.e-03) then
       return
    end if

    ! (1) Initialize switches at first call
    if (initFlag) call ReadInput

    ! (2) Divide in different scattering scenarios:
    scenario = LogicMatrix(isMeson(partIn(1)%ID), isMeson(partIn(2)%ID))

    select case (scenario)
    case (scenarioBarBar)
      if (.not.baryonBaryonScattering) return
    case (scenarioBarMes)
      if (.not.baryonMesonScattering) return
    case (scenarioMesMes)
      if (.not.mesonMesonScattering) return
    end select

    ! (3) Evaluate scaling factor for leading hadrons
    stringFactor = 1.
    if (partIn(1)%inF) stringFactor=stringFactor*partIn(1)%scaleCS
    if (partIn(2)%inF) stringFactor=stringFactor*partIn(2)%scaleCS

    if (doScaleResidue) then ! we abuse stringfactor for the scaling
       if (partIn(2)%pert) then ! particle 1 is always real
          if (partIn(1)%ID==1) then
             stringFactor = stringFactor &
                  * ResidueGetScale(partIn(2)%firstEvent,partIn(1)%charge)
          end if
       end if
    end if

    if (stringFactor<1E-8) return

    ! (4) In case of the "naive" geometrical collision criteria, we here
    !     exclude collisions based on rough assumptions about the locality
    !     in space and time of these collisions:
    if (.not.(localEnsemble.and.fullEnsemble)) then
       ! (4a) Use very rough collision Criteria with cut-off's which are not
       !      relativistic
       if (.not.check_veryRough(partIn,stringFactor,numEnsembles,scenario)) return

       ! (4b) Apply time criteria to ensure that collision happens at minimal
       !      distance (similar to Kodama)
       if (.not.kodama_time(partIn,delta_T,collTime))  return
       if (debug) write(*,*) ' After time criteria in master_2body/collide_2Body'

       ! (4c) Use rough collision Criteria with cut-off's to ensure that
       !      collisions are local in space
       if (.not.checkKodama_rough(partIn,stringFactor,numEnsembles,scenario)) return
    end if

    if (debug) write(*,*) 'Velos in master_2_body=', &
         partIn(1)%vel, '##'  ,partIn(2)%vel

    call generateFinalState(partIn, partOut, stringFactor, numEnsembles, &
         time, collFlag, HiEnergyFlag, HiEnergyType, weightLocal, scenario, &
         XShelp, PauliIncluded_out=PauliIncluded, srtS_vacuum_in=sqrtS_vacuum)

    if (present(sigTot_out)) then
       sigTot_out = XShelp
    end if

    if (present(pauliIncluded_out)) then
       pauliIncluded_out=pauliIncluded
    end if

    call setHistory(partIn, partOut)

  end subroutine collide_2body


  !****************************************************************************
  !****f* master_2Body/check_veryRough
  ! NAME
  ! logical function check_veryRough(partIn,stringFactor,numEnsembles,scenario)
  !
  ! PURPOSE
  ! Very Rough collision criterium
  !
  ! INPUTS
  ! * type(particle), dimension(1:2) :: partIn
  ! * real                           :: stringFactor
  ! * integer                        :: numEnsembles -- number of ensembles
  ! * integer                        :: scenario -- type of scattering scenario
  !
  ! OUTPUT
  ! * logical                        :: flag
  !****************************************************************************
  function check_veryRough(partIn,stringFactor,numEnsembles,scenario) &
       Result(flag)

    use particleDefinition
    use inputGeneral, only: delta_T,fullEnsemble
    use barbar_barbar, only:delta2Body_inMedium_treatment
    use idTable, only: delta, nucleon

    type(particle), dimension(1:2), intent(in) :: partIn
    real, intent(in) :: stringFactor
    integer, intent(in) :: numEnsembles
    integer, intent(in) :: scenario

    logical :: flag
    real :: bMax_coarse    ! maximal impact parameters
    real, dimension(1:3) :: dx   ! DeltaX of both Particles

    flag=.true.

    ! First coarse approximate collision criteria for spacial distance:
    ! (To quote the old code:) The following prescription is not covariant,
    ! but useful to save cpu time :

    if (((partIn(1)%ID.eq.delta).or.(partIn(2)%ID.eq.delta)).and. &
         & delta2Body_inMedium_treatment()) then
       ! NO maximal impact parameter for the case of Oset model
       ! in ND->NN and ND->ND collisions:
       if ((partIn(1)%ID.eq.nucleon).or.(partIn(2)%ID.eq.nucleon)) then
          return
       end if
    end if

    if (fullEnsemble) then
       bmax_coarse=sqrt( delta_T**2 &
            +coarse(scenario)**2*stringFactor/real(numEnsembles) )
    else
       bmax_coarse=sqrt( delta_T**2 &
            +coarse(scenario)**2*stringFactor )
    end if

    dx=partIn(1)%pos-partIn(2)%pos
    if (abs(dx(1)) > bmax_coarse) then
       flag=.false.
       return
    else if (abs(dx(2)) > bmax_coarse) then
       flag=.false.
       return
    else if (abs(dx(3)) > bmax_coarse) then
       flag=.false.
       return
    else if (Dot_Product(dx,dx) > bmax_coarse**2) then
       flag=.false.
       return
    end if
  end function check_veryRough


  !****************************************************************************
  !****f* master_2Body/checkKodama_Rough
  ! NAME
  ! logical function checkKodama_Rough(partIn,stringFactor,numEnsembles,scenario)
  !
  ! PURPOSE
  ! Checks in a rough way wether the "partIn" of particles is so close in
  ! space that a collision should take place. Therefore we check with the
  ! kodama criteria whether a collision would take place, if we assume that
  ! the Xsection is maximal, which means that
  !                     xSection = 2 pi cutOff**2.
  ! If "kodama_position" returns .true. than the "flag" is set .true. .
  !
  ! NOTES
  ! This check is meant to save computation time. If a collision is already
  ! excluded assuming the maximal cross section, then we do not even have to
  ! check for the real cross section and the correct collision criteria!
  !
  ! INPUTS
  ! * type(particle), dimension(1:2):: partIn         -- incoming particles
  ! * real        :: stringFactor -- XS reduction factor
  ! * integer     :: numEnsembles -- number of ensembles used
  ! * integer     :: scenario     -- kind of scattering scenario
  !
  ! OUTPUT
  ! * logical :: flag
  !****************************************************************************
  function checkKodama_Rough(partIn,stringFactor,numEnsembles,scenario) &
       Result(flag)

    use particleDefinition
    use particleProperties, only: hadron
    use collisionCriteria, only: kodama_position
    use inputGeneral, only: fullEnsemble
    use idTable, only: delta, nucleon,pion
    use barbar_barbar, only:delta2Body_inMedium_treatment

    type(particle), dimension(1:2),intent(in) :: partIn
    real,   intent(in) :: stringFactor
    integer,intent(in) :: numEnsembles
    integer,intent(in) :: scenario

    logical :: flag
    real :: bMax    ! maximal impact parameter

    flag=.true.

    ! Evaluate maximum impact parameter
    select case (scenario)
    case (scenarioBarBar) ! baryon baryon
       if ( (hadron(partIn(1)%ID)%strangeness.ne.0).or.(hadron(partIn(2)%ID)%strangeness.ne.0))&
          & then
          bmax=bmax_hyperonNucleon*sqrt(stringfactor)     ! hyperon-nucleon
       else if ((partIn(1)%ID.eq.nucleon).and.(partIn(2)%ID.eq.nucleon)) then
          bmax=bmax_nucleonNucleon*sqrt(stringfactor)     ! nucleon-nucleon
       else
          bmax=bmax_nucleonResonance*sqrt(stringfactor)
          if (((partIn(1)%ID.eq.delta).or.(partIn(2)%ID.eq.delta)).and. &
               & delta2Body_inMedium_treatment()) then
             ! NO maximal impact parameter for the case of Oset model
             ! in ND->NN and ND->ND collisions:
             ! We choose b_max=5fm we should be numerically the same as
             ! "no maximal b_max".
             if ((partIn(1)%ID.eq.nucleon).or.(partIn(2)%ID.eq.nucleon)) then
                bmax=5.
             end if
          end if
       end if
    case (scenarioBarMes) ! meson baryon
       if ((partIn(1)%ID.eq.pion).or.(partIn(2)%ID.eq.pion)) then
          bmax=bmax_baryonPion*sqrt(stringfactor)         ! pion baryon
       else
          bmax=bmax_baryonMeson*sqrt(stringfactor)        ! meson baryon
       end if
    case (scenarioMesMes) ! meson meson
       bmax=bmax_mesonMeson*sqrt(stringfactor)
    end select

    if (fullEnsemble) bMax=bMax/sqrt(real(numEnsembles))

    ! use original space-criteria for scattering (Kodama) with bMax:
    if (.not.kodama_position(partIn,bMax)) then
       flag=.false.
       return
    end if
  end function checkKodama_Rough

  !****************************************************************************
  !****s* master_2Body/generateFinalState
  ! NAME
  ! subroutine generateFinalState(partIn, partOut, stringFactor,
  ! numEnsembles, time, collFlag, HiEnergyFlag, HiEnergyType,
  ! weightLocal, scenario0, sigTot_out, pauliIncluded_out, srts_vacuum_in)
  !
  ! PURPOSE
  ! This is the routine, which really does the collision
  !
  ! INPUTS
  ! * type(particle), dimension(1:2) :: partIn         -- incoming particles
  ! * real                :: stringFactor -- rescaling of XS
  ! * integer             :: numEnsembles -- number of ensembles
  ! * real                :: time         -- actual time step
  ! * integer, optional   :: weightLocal  -- weight which is used for the
  !   local ensemble method.
  !   Used for reweighting the probability in the collision criteria.
  !   Not necessary in case of localEnsemble=.false.
  ! * integer, optional   :: scenario0    -- if given : BarBar,BarMes,MesMes
  ! * real, optional      :: srtS_vacuum_in -- free sqrt(s)
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: partOut     -- outgoing particles
  ! * logical            :: collFlag  -- true if collisions was okay
  ! * logical            :: HiEnergyFlag   -- true if HiEnergy was used
  ! * integer            :: HiEnergyType   -- 0:LowEnergy, 1:Fritiof, 2:Pythia
  ! * real, optional     :: sigTot_out     -- total cross section
  ! * logical, optional  :: pauliIncluded_out  -- true if Pauli blocking was
  !   already included in the cross section evaluation
  !
  !****************************************************************************
  subroutine generateFinalState(partIn, partOut, stringFactor, &
       numEnsembles, time, &
       collFlag, HiEnergyFlag, HiEnergyType, weightLocal, scenario0, &
       sigTot_out, PauliIncluded_out, srts_vacuum_in)

    use densitymodule, only: densityAt
    use collisionCriteria, only: kodama_position,localCollisionCriteria
    use twoBodyTools, only: sqrtS_free, LogicMatrix
    use particleDefinition, only: particle, setToDefault, resetNumberGuess,&
         setnumberguess, sqrtS, useProductionPos, setProductionPos
    use mediumDefinition
    use inputGeneral, only: fullEnsemble,localEnsemble,delta_t
    use lorentzTrafo, only: lorentzCalcBeta, lorentz
    use minkowski, only: abs4Sq
    use propagation, only: updateVelocity, checkVelo
    use dichtedefinition
    use preEventDefinition
    use idTable, only: pion, delta, nucleon, isMeson, isBaryon
    use RMF, only: getRMF_flag
    use energyCalc, only: energyCorrection, energyDetermination
    use XsectionRatios, only: getSigmaScreened
    use PIL_mesonMom, only: PIL_mesonMom_PUT
    use offShellPotential, only: setOffShellParameter
    use mediumModule, only: mediumAt,getMediumCutOff
    use output, only: writeParticle_debug, writeParticle
    use deuterium_PL, only: deuteriumPL_inUse
    use Annihilation, only: annihilate
    use Dilepton_Analysis, only: Dilep_Brems
    use constants, only: pi

    type(particle), dimension(1:2), intent(in) :: partIn
    type(particle), dimension(:), intent(inout) :: partOut
    real, intent(in) :: stringFactor
    integer, intent(in) :: numEnsembles
    real, intent(in) :: time

    logical, intent(out) :: collFlag
    logical, intent(out) :: HiEnergyFlag
    integer, intent(out) :: HiEnergyType
    integer, optional, intent(in)  :: weightLocal
    integer, optional, intent(in)  :: scenario0
    real, optional, intent(out)    :: sigTot_out
    logical, optional, intent(out) :: pauliIncluded_out
    real, optional, intent(in)     :: srtS_vacuum_in

    logical :: PauliIncluded, success, successFlag, potFailure
    real, dimension(0:3) :: momentum_calc, momentum_LRF
    real :: srtS,srtS_vacuum,srtS_corr,srtS_XS,srtS_final,ratio
    type(dichte) :: density
    type(medium) :: mediumAtColl
    real, dimension(1:3)  :: position, betaToLRF, betaToCM
    integer :: maxID, k, i, nAttempts, scenario
    real, dimension(1:2) :: mstar
    real, dimension(0:7) :: sigs ! holds sigma, cf. isigXXX
    real :: bmax, sigmaTot_screened, sigmaTot_rescal
    type(preEvent), dimension(1:4) :: chosenEvent
    real, dimension(0:3) :: mesonMomentum_CM

    if (present(pauliIncluded_out)) then
       pauliIncluded_out=.false.
    end if

    pauliIncluded=.false.

    ! (0) assume failure at return:
    collFlag=.false.

    ! (1) Evaluate Sqrt(s) using the calculation frame
    momentum_calc(0:3)=partIn(1)%mom(0:3)+partIn(2)%mom(0:3)
    srtS = sqrtS(partIn,"generateFinalState, srtS")

    ! (2) Define Sqrt(s) in the vacuum
    if (present(srtS_vacuum_in)) then
       srtS_vacuum = srtS_vacuum_in
    else
       srtS_vacuum = sqrtS_free(partIn)
    end if
    if (srtS_vacuum.lt.Sum(partIn%mass)) then
       write(*,*) 'Error in master_2body: srtS_vacuum.lt.Sum(partIn%mass): ',&
            'srtS_vacuum',Sum(partIn%mass), srtS_vacuum
       write(*,*) partIn%mass, partIn%ID, partIn%charge
       write(*,*) partIn%pert, partIn%anti
       write(*,*) partIn(1)%mom, partIn(1)%vel
       write(*,*) partIn(2)%mom, partIn(2)%vel
       write(*,*) 'stopping'
       call Traceback()
    else if (srtS_vacuum-Sum(partIn%mass).lt.1.e-03) then
       return
    end if

    ! (3) Read out medium at collision point in LRF
    position=(partIn(1)%pos+partIn(2)%pos)/2.
    density=densityAt(position)
    mediumAtColl=mediumAt(density,position)
    if (deuteriumPL_inUse()) mediumAtColl%useMedium = .true.

    ! (4) Define total momentum in LRF.
    ! If density not negligible then boost is needed.
    momentum_LRF=momentum_Calc
    if (density%baryon(0).gt.getMediumCutOff()/100.) then
       if(.not.getRMF_flag()) then
          betaToLRF = lorentzCalcBeta(density%baryon)
       else
          if(abs4Sq(density%baryon).gt.0.) then
             betaToLRF(1:3) = density%baryon(1:3)/density%baryon(0)
          else
             betaToLRF = 0.
          end if
       end if
       call lorentz(betaToLRF, momentum_LRF)
    else
       betaToLRF=0.
    end if

    if (debug) write(*,*) srtS_vacuum,mediumAtColl%useMedium, momentum_LRF

    ! (5) Evaluate cross sections and set some special chosen channel for
    ! the final state.

    if (present(scenario0)) then
       scenario = scenario0
    else
       scenario = LogicMatrix(isMeson(partIn(1)%ID), isMeson(partIn(2)%ID))
    end if

    if ( .not.getRMF_flag() ) then
       srtS_XS = srtS_vacuum ! This is the srtS value the XS should be calculated with
    else
       mstar(1) = sqrtS(partIn(1),'generateFinalState, mstar(1)')
       mstar(2) = sqrtS(partIn(2),'generateFinalState, mstar(2)')

       srtS_corr = srtS - mstar(1) - mstar(2) + partIn(1)%mass + partIn(2)%mass
       srtS_XS = srtS_corr
    end if

    call XsectionMaster(srtS_XS, partIn, mediumAtColl, momentum_LRF, &
         chosenEvent, sigs, HiEnergyFlag, PauliIncluded_out=pauliIncluded)

    if (present(pauliIncluded_out)) pauliIncluded_out = pauliIncluded
    if (present(sigTot_out)) then
       sigTot_out = sigs(isigTot)
       if (sigs(isigTot)<0.0001) return
    end if
    if (sigs(isigTot)<XsectionCutOff) return


    ! In-medium screening:
    sigmaTot_screened = getSigmaScreened(partIn,sigs(isigTot),srtS_XS,mediumAtColl)
    ratio = sigmaTot_screened/sigs(isigTot)
    sigs = sigs * ratio

    if (debug) then
       write(*,*) ' Xsections:'
       write(*,*) ' sigmaTot, sigmaElast, sigmaCEX: ', &
            sigs(isigTot), sigs(isigElast), sigs(isigCEX)
       write(*,*) ' sigmaAnni, sigmaLbar, sigmaSbar: ', &
            sigs(isigAnni), sigs(isigLbar), sigs(isigSbar)
       write(*,*) ' sigmaXiBar, sigmaJPsi, HiEnergyFlag: ', &
            sigs(isigXiBar), sigs(isigJPsi), HiEnergyFlag
    end if

    ! (6) Intialize output
    call setToDefault(partOut)
    partOut%pert = (partIn(1)%pert.or.partIn(2)%pert)

    partOut(1:4)%ID=chosenEvent(1:4)%ID
    partOut(1:4)%charge=chosenEvent(1:4)%charge
    partOut(1:4)%anti=chosenEvent(1:4)%anti
    partOut(1:4)%mass=chosenEvent(1:4)%mass

    partOut%prodTime = time
    partOut%formTime  = time

    ! (7) Use exact space-criteria for scattering
    !
    ! * There are two different schemes for full ensemble runs
    ! a) Kodama (e.g. Effenberger Phd thesis)
    ! b) Local Ensemble (e.g. Lang Phd thesis)
    !
    ! * In parallel ensemble mode, only Kodama can be used.
    !
    if (.not.(localEnsemble.and.fullensemble)) then
       ! (7a) Use Kodama criteria with the real sigma :
       bmax=SQRT(sigs(isigTot)/10./pi*stringfactor)
       if (fullEnsemble) then
          bMax=bMax/sqrt(real(numEnsembles))
          if (debug) &
               write(*,*) 'fullensemble: bmax=bmax/ ', sqrt(real(numEnsembles))
       end if
       if (debug) write(*,*) 'bmax=',bmax

       if (.not.kodama_position(partIn,bMax)) return

    else
       ! (7b) Use local ensemble Collision Criteria
       ! This is only possible in the full ensemble mode!
       sigmaTot_rescal=sigs(isigTot)*stringfactor
       if(.not.localCollisionCriteria(partIn,sigmaTot_rescal,weightLocal,numEnsembles,delta_T))&
         & return

    end if

    collFlag=.true. ! now assume success at return...
    if (debug) write(*,*) 'set kinematics'

    ! (7.5) Calculate dilepton production via Bremsstrahlung for NN and piN
    ! collisions
    ! (only if dilepton analysis is enabled).
    call Dilep_Brems(partIn, srtS_XS, sigs(isigElast), sigs(isigTot))


    ! (9) Set positions and kinematics of outgoing particles

    betaToCM = lorentzCalcBeta(momentum_calc)

    if (.not.HiEnergyFlag) then  !**** NOT HIGH ENERGY ****

       if ( isBaryon(partIn(1)%Id) .and. isBaryon(partIn(2)%Id) .and. &
            & (partIn(1)%anti.neqv.partIn(2)%anti) .and. &
            & partOut(1)%Id.eq.pion .and. partOut(2)%Id.eq.pion .and.&
            & partOut(3)%Id.eq.pion ) then

          ! Simulate annihilation:

          if (partIn(1)%anti) then
             call annihilate(partIn(1),partIn(2),time,partOut,&
                  collFlag,HiEnergyType)
          else
             call annihilate(partIn(2),partIn(1),time,partOut,&
                  collFlag,HiEnergyType)
          end if

          ! Bad trick to avoid charge check after annihilation:
          ! (unfortunately charge is not always conserved by annihilate)
          HiEnergyFlag=.true.
          return

       else
          ! Usual low-energy mode:

          call ResetPosition ! this also sets maxID
          HiEnergyType=0
          call setKinematics(srtS, srtS_XS, betaToLRF, betaToCM, &
               mediumAtColl, partIn, partOut(1:maxID), collFlag)
          if (.not. collFlag) return

       end if

    else                       !**** USE HIGH ENERGY ****

       ! write(*,*) ' Before setKinematicsHiEnergy:',srtS,sigmaTot,sigmaElast

       nAttempts = 0
       do
          nAttempts = nAttempts + 1

          if (nAttempts>10) then
             write(*,*) ' In generateFinalState: energy correction FAILED:'
             write(*,*) ' Colliding particles:', &
                  partIn(1:2)%Id, partIn(1:2)%anti
             if (getRMF_flag()) write(*,*) ' Corrected srtS:', srtS_corr
             write(*,*) ' srtS for XS:   ', srtS_XS
             write(*,*) ' wished srtS:   ', srtS
             write(*,*) ' bare masses:', partIn(1:2)%mass
             if (getRMF_flag()) write(*,*) ' Effective masses:', mstar(1:2)
             write(*,*) ' final particles:', &
                  partOut(1:maxId)%Id, partOut(1:maxId)%anti
             write(*,*) ' final bare masses:',  partOut(1:maxId)%mass
             write(*,*) ' sum of final bare masses:', &
                  sum(partOut(1:maxId)%mass)

             srtS_final = sqrtS(partOut(1:maxId), 'srtS_final')

             write(*,*) ' final srtS:', srtS_final
             collFlag = .false.
             return
          end if

          call setKinematicsHiEnergy(srtS, srtS_XS, sigs, &
               betaToCM, partIn, time, partOut, collFlag, &
               HiEnergyType, scenario)

          if (.not.collFlag) return

          call ResetPosition

          ! Perform energy correction for all events
          ! excepting baryon-antibaryon annihilation,
          ! where energy correction has been already done
          if (HiEnergyType.eq.-3) exit

          ! This is needed because energyCorrection takes finalState in the
          ! CM frame and returns it in the calculational frame:
          do i=1,maxId
             call lorentz(betaToCM,partOut(i)%mom)
          end do

          call energyCorrection(srtS, betaToCM, &
               partOut(1:maxId), successFlag, potFailure)

          if (potFailure) then
             if (correctEnergy_message) then
                write(*,'(2A)') 'generateFinalState: ',&
                     'Energy correction failed. (potentialFailure)'
             end if
             collFlag = .false.
             return
          end if
          if (successFlag) exit

       end do

    end if


    if (maxID==1) then
       call resetNumberGuess()
       call setNumberGuess(partOut(1))

       ! (10) store additional information for Delta
       if (partOut(1)%ID==delta .and. &
            & ((partIn(1)%ID==pion .and. partIn(2)%ID==nucleon) .or. &
            &  (partIn(2)%ID==pion .and. partIn(1)%ID==nucleon)) ) then

          if (partIn(1)%ID==pion) then
             mesonMomentum_CM=partIn(1)%mom
          else
             mesonMomentum_CM=partIn(2)%mom
          end if
          call lorentz(betaToCM,mesonMomentum_CM)
          call PIL_mesonMom_PUT(partOut(1)%number, mesonMomentum_CM(1:3))
       end if
    end if

    ! (11) Update velocities
    do i = 1,maxId
        if (partOut(i)%id.eq.0) cycle
        if (.not.getRMF_flag()) call energyDetermination(partOut(i))
        call setOffShellParameter(partOut(i), success)
        if (.not.success) then
            collFlag = .false.
            return
        end if
        if (.not.getRMF_flag()) then
            call updateVelocity(partOut(i),success)
        else
           partOut(i)%vel = partOut(i)%mom(1:3)/partOut(i)%mom(0)
           success=checkVelo(partOut(i))
        end if
        if (.not.success) then
            write(*,*) 'Master_2Body(1): velocity ge 1: collFlag -> FALSE'
            write(*,*) 'initial particles:',partIn(1)%ID,partIn(2)%ID
            if (debug) call WriteParticle_debug(partIn(1))
            if (debug) call WriteParticle_debug(partIn(2))
            collFlag = .false.
            return
        end if
    end do

    ! (12) set statistical factors
    partOut(:)%lastCollTime=time

    if (debug) then
       write(*,*) 'Id1,Id2,Charge1,Charge2:', partIn%ID,partIn%charge
       write(*,*) 'Antiparticles?:', partIn%anti
       write(*,*) 'srtS_vacuum:', srtS_vacuum

       do i=1,maxId
          if (partOut(i)%id <= 0) exit
          write(*,*) 'Id, antiparticle?, charge:',&
               & partOut(i)%Id,partOut(i)%anti,&
               & partOut(i)%charge
       end do
    end if

    if (useProductionPos) call setProductionPos(partOut)

  contains
    !**************************************************************************
    !****s* generateFinalState/ResetPosition
    ! NAME
    ! subroutine ResetPosition
    ! PURPOSE
    ! * calculate maxID, i.e. the number of finalState-particles
    ! * set position of finalState-particles:
    !   It loops over the FS particles vector and sets for the first
    !   FS-particle, which is baryon/meson as the first IS-particle,
    !   the spatial position as the first IS particle.
    !   Looping continues and the same is done for the next FS particle,
    !   which is of the same type as the secon IS particle.
    !   All other particles get as position the position of the interaction.
    !
    ! NOTES
    ! It is even a little bit more involved; get the power, read the source ;)
    !**************************************************************************
    subroutine ResetPosition

      use deuterium_pl, only: deuteriumPL_inUse

      logical :: flag1,flag2
      maxId= ubound(partOut,dim=1)
      do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
         if (partOut(k)%Id.eq.0) then
            maxID=k-1
            exit
         end if
      end do

      if (maxID==1) then
         ! In processes where 2->1 (i.e. resonance production), the
         ! resonance is created in the middle
         ! For Deuterium it is not the middle but the position of the
         ! old baryon if it's a baryon.
         if (.not. deuteriumPL_inUse()) then
            partOut(lbound(partOut,dim=1))%pos = position
         else if (isBaryon(partIn(1)%ID) .and. isBaryon(partOut(1)%ID)) then
            partOut(lbound(partOut,dim=1))%pos = partIn(1)%pos
         else if (isBaryon(partIn(2)%ID) .and. isBaryon(partOut(1)%ID)) then
            partOut(lbound(partOut,dim=1))%pos = partIn(2)%pos
         else
            partOut(lbound(partOut,dim=1))%pos = position
         end if
         return
      end if

      flag1= .false.
      flag2= .false.
      do k=lbound(partOut,dim=1),maxID
         if ( .not.flag1 .and. &
              ( (isBaryon(partIn(1)%ID) .and. isBaryon(partOut(k)%ID) .and. &
              (partIn(1)%anti.eqv.partOut(k)%anti))&
              .or. &
              (isMeson(partIn(1)%ID) .and. isMeson(partOut(k)%ID)) ) ) then
            partOut(k)%pos = partIn(1)%pos
            flag1= .true.
         else if ( .not.flag2 .and. &
              ( (isBaryon(partIn(2)%ID) .and. isBaryon(partOut(k)%ID) .and. &
              (partIn(2)%anti.eqv.partOut(k)%anti))&
              .or. &
              (isMeson(partIn(2)%ID) .and. isMeson(partOut(k)%ID)) ) ) then
            partOut(k)%pos = partIn(2)%pos
            flag2 = .true.
         else
            partOut(k)%pos = position
         end if
      end do

    end subroutine ResetPosition

  end subroutine generateFinalState



  !****************************************************************************
  !****s* master_2Body/setKinematics
  ! NAME
  ! subroutine setKinematics(srts, srts_vacuum_in, betaToLRF, betaToCM,
  ! mediumAtColl, partIn, partOut, collFlag, verbose)
  !
  ! PURPOSE
  ! Evaluates the kinematics for the "finalState" particles assuming vacuum
  ! kinematics (-> no potentials).
  !
  ! INPUTS
  ! * type(particle),dimension(1:2):: partIn      -- incoming initial particles
  ! * type(particle),dimension(:)  :: partOut     -- outgoing particles
  ! * real                         :: srts        -- sqrt(s)
  ! * real                         :: srts_vacuum_in -- sqrt(s) in the vacuum
  ! * real, dimension(1:3)         :: betaToLRF
  !   -- beta of calc frame to LRF frame
  ! * real, dimension(1:3)         :: betaToCM
  !   -- beta of calc frame to CM frame
  ! * type(medium)                 :: mediumAtColl -- medium information
  ! * logical, OPTIONAL            :: verbose
  !   -- flag print messages on/off (default: on)
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: partOut      -- outgoing particles
  ! * logical                      :: collFlag
  !   -- "true" if kinematics could be set
  !
  ! NOTES
  ! * It's important that the Id's, positions, perturbative flags and charges
  !   of the "partOut" particles are already set when calling this
  !   subroutine.
  ! * Only the kinematics (including masses of the finalState) will be set.
  !****************************************************************************
  subroutine setKinematics(srts, srts_vacuum_in, betaToLRF, betaToCM, &
       mediumAtColl, partIn, partOut, collFlag, verbose)

    use mediumDefinition
    use particleDefinition
    use finalStateModule, only: assMass, massAss
    use energyCalc, only: energyCorrection
    use lorentzTrafo, only: lorentz
    use output, only: WriteParticle
    use idTable, only: pion,nucleon,photon,isBaryon
    use particleProperties, only: hadron
    use RMF, only: getRMF_flag, flagCorThr
    use densitymodule, only: getGridIndex, SelfEnergy_scalar
    use nBodyPhaseSpace, only: momenta_in_4BodyPS
    use baryonWidthMedium, only: get_MediumSwitch_coll

    real,                         intent(in)    :: srts, srts_vacuum_in
    real, dimension(1:3),         intent(in)    :: betaToLRF, betaToCM
    type(medium),                 intent(in)    :: mediumAtColl
    type(particle),dimension(1:2),intent(in)    :: partIn
    type(particle), dimension(:), intent(inOut) :: partOut
    logical,                      intent(out)   :: collFlag
    logical, OPTIONAL,            intent(in)    :: verbose

    integer :: i,j
    real, parameter, dimension (1:3) ::spotOut=0. ! Neglecting possible
                                                  ! scalar potentials !!!
    logical :: flag, successFlag, verb
    logical :: potFailure ! set by energycorrection if there is a threshold
                          ! violation due to potentials which are neglected
    integer, parameter :: maxCorrectLoop=20  ! maximal number of iterations
                                             ! for energy correction
    real :: srts_vacuum

    ! (1) Initialize switches at first call
    if (initFlag) call ReadInput

    collFlag=.true.
    verb = .true.
    if (present(verbose)) verb = verbose

    ! CHECK INPUT
    do i=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (partOut(i)%Id.eq.0) then
          write(*,*) 'Critical Error in setKinematics'
          write(*,*) 'Particle ID equals 0!!!'
          write(*,*) 'stop'
       end if
    end do

    flag=.true.

    srts_vacuum = srts_vacuum_in

    ! Set Energy and momenta

    select case (size(partOut,dim=1))
    case (0)
       write(*,'(A)') 'Master_2Body: Error in setKinematics: ',&
            'No finalstate particles! '

    case (1)  ! ===== one body final state =====
       partOut(1)%mom = partIn(1)%mom(0:3) + partIn(2)%mom(0:3)

    case (2)  ! ===== two body final state =====
       call Do2Body

    case (3)  ! ===== three body final state =====
       call Do3Body

    case (4)  ! ===== four body final state =====
       call Do4Body

    case default
       write(*,'(A)') 'Master_2Body: Error in setKinematics: ',&
            'No treatment of final state with 4+ particles implemented yet!'
       write(*,*) size(partOut), partOut(:)%ID
       call Traceback()
    end select

    ! set velocities in vacuum approximation

    !    Do i=lBound(partOut,dim=1), uBound(partOut, dim=1)
    !       partOut(i)%vel=partOut(i)%mom(1:3)/partOut(i)%mom(0)
    !    End do

    partOut%scaleCS=1.
    partOut%inF=.false.

  contains
    !**************************************************************************
    subroutine Do2Body
      use constants, only: mN, mPi
      integer, save :: fehlerZaehler=0
      real,    dimension(1:2) :: spotOut
      integer,    dimension(1:3) :: ind

      if(isBaryon(partOut(1)%id) .and. isBaryon(partOut(2)%id) .and. getRMF_flag() .and. flagCorThr)  then
          spotOut(1:2)=0.
          if(getGridIndex(partOut(1)%pos,ind,0)) &
              & spotOut(1) = SelfEnergy_scalar(ind(1),ind(2),ind(3),partOut(1)%id,.false.)
          if(getGridIndex(partOut(2)%pos,ind,0)) &
              & spotOut(2) = SelfEnergy_scalar(ind(1),ind(2),ind(3),partOut(2)%id,.false.)
          srts_vacuum = srtS - spotOut(1) - spotOut(2)
      end if

      energyCorrectLoop : do i=1, maxCorrectLoop
         if (debug) write(*,*) 'Vor massAss'
         call massass(srts_vacuum,mediumAtColl,partIn,partOut(1:2),&
              betaToLRF,betaToCM,0,flag)

         if (.not.flag) then

            do j=1,2
               select case (partIn(j)%ID)
               case (nucleon)
                  if (partIn(j)%mass<mN .and. .not.get_MediumSwitch_coll()) then
                     write(*,*) '... it was a nucleon with strange mass! (2) ',&
                          partIn(j)%mass
                     collFlag=.false. ! Kill Event
                     partOut(1:2)%ID=0
                     return ! --> Exit
                  end if
               case (pion)
                  if ((partIn(j)%mass < mPi) .and. (srts_vacuum-partIn(3-j)%mass < mPi)) then
                     write(*,*) '... it was a Jetset-Pion. threshold violated. ok!'
                     collFlag=.false. ! Kill Event
                     partOut(1:2)%ID=0
                     return ! --> Exit
                  end if
               case (photon)
                  collFlag=.false. ! Kill Event
                  partOut(1:2)%ID=0
                  return ! --> Exit

               end select

               if ((partIn(j)%mass .lt. hadron(partIn(j)%ID)%minmass+0.005 )) then
                  write(*,*) '... it was a resonance with strange mass! (1)', &
                             & partIn(j)%ID,partIn(j)%mass
                  write(*,*) partIn%ID,'->',partOut%ID
                  collFlag=.false. ! Kill Event
                  partOut(1:2)%ID=0
                  return ! --> Exit
               end if
            end do

            write(*,*) 'master_2Body: Impossible to find final state (1)'
            write(*,*) partIn%Id, partOut%ID
            write(*,*) srts, srts_vacuum

            call WriteParticle(6,99,1,partIn(1))
            call WriteParticle(6,99,2,partIn(2))

            call WriteParticle(6,-99,1,partOut(1))
            call WriteParticle(6,-99,2,partOut(2))

            fehlerZaehler=fehlerzaehler+1
            if (fehlerZaehler>fehler_max) then
               write(*,*) 'too many errors in master_2Body: massass'
               call Traceback()
            end if

            collFlag=.false. ! Kill Event
            partOut(1:2)%ID=0
            return ! --> Exit
         end if

         if ( correctEnergy .or. getRMF_flag() ) then
            if (debug) write(*,*) '*************************'
            if (debug) write(*,*) 'wished srts=', srts
            if (debug) write(*,*) 'initial srts=', sqrts(partIn(1),partIn(2))
            call energyCorrection(srts,betaToCM, &
                 & partOut, successFlag,potFailure, verb)
            if (debug) write(*,*) 'final srts=', sqrts(partOut(1),partOut(2))
            if (debug) write(*,*) 'initial srts=', sqrts(partIn(1),partIn(2))
            if (successFlag) exit energyCorrectLoop
            if (potFailure) exit energyCorrectLoop ! Failed due to threshold violation
         else
            ! just boost to Calculation frame
            successFlag=.true.
            do j=1,2
               call lorentz(-betaToCM, partOut(j)%mom)
            end do
            exit energyCorrectLoop
         end if
      end do energyCorrectLoop

      if (.not.successFlag) then

         if (potFailure) then
            if (correctEnergy_message) then
               write(*,'(A)') 'SetKinematics[2]: Energy correction failed.' &
                    //' Threshold is violated since potentials are neglected.' &
                    //' ok!'
            end if
            collFlag=.false. ! Kill Event
            partOut(1:2)%ID=0
            return ! --> Exit
         end if

         if (verb) then
            write(*,'(A,2I4,A,2I4,A,2ES12.4)') &
                 'SetKinematics[2]: Energy correction failed.', &
                 partIn(1:2)%Id,' ->',partOut(1:2)%Id,' @',srts,srts_vacuum
         end if

!!$         call WriteParticle(6,99,1,partIn(1))
!!$         call WriteParticle(6,99,2,partIn(2))
!!$
!!$         call WriteParticle(6,-99,1,partOut(1))
!!$         call WriteParticle(6,-99,2,partOut(2))


         collFlag=.false. ! Kill Event
         partOut(1:2)%ID=0
         return ! --> Exit
      end if

    end subroutine Do2Body

    !**************************************************************************
    subroutine Do3Body

      integer, save :: fehlerZaehler=0

      do i=1, maxCorrectLoop
         ! Neglecting possible scalar potentials, since spotOut=0
         call assMass(srts_vacuum,mediumAtColl,partIn,partOut,spotOut(1:3), &
              &betaToLRF,betaToCM,flag)
         if (.not.flag) then
            write(*,*) 'master_2Body, after assmass: ',&
                 'Impossible to find final state (2)'
          write(*,*) partIn%Id, partOut%ID, srts, srts_vacuum
          write(*,*) "partIn(1)"
          write(*,*) partIn(1)%mass, partIn(1)%mom,&
               partIn(1)%charge, partIn(1)%anti
          write(*,*) "partIn(2)"
          write(*,*) partIn(2)%mass, partIn(2)%mom,&
               partIn(2)%charge, partIn(2)%anti
          successFlag=.false.
          cycle
         end if

         if (correctEnergy) then
            call energyCorrection(srts,betaToCM, &
                 & partOut, successFlag,potFailure)
            if (successFlag) exit
            if (potFailure)  exit  ! Failed due to threshold violation
         else
            ! just boost to Calculation frame
            successFlag=.true.
            do j=1,3
               call lorentz(-betaToCM,partOut(j)%mom)
            end do
         end if
      end do

      if (.not.successFlag) then
         write(*,'(A,2I4,A,3I4,A,2ES12.4)') &
              'SetKinematics[3]: Energy correction failed.', &
              partIn(1:2)%Id,' ->',partOut(1:3)%Id,' @',srts,srts_vacuum

         fehlerZaehler=fehlerzaehler+1
         if (fehlerZaehler>fehler_max) then
            write(*,*) 'too many errors in master_2Body: Do3Body'
            call Traceback()
         end if

         collFlag=.false. ! Kill Event
         partOut(1:3)%ID=0
         return ! --> Exit
      end if

    end subroutine Do3Body

    !**************************************************************************
    subroutine Do4Body
      use IdTable, only: isHadron
      use particleProperties, only: hadron
      use propagation, only: checkVelo
      !
      ! This routine is only prepared for 'stable' outgoing particles,
      ! as e.g. in pi N -> pi pi pi N.
      !
      real, dimension(4)   :: mass
      real, dimension(3,4) :: p4
      integer, save :: fehlerZaehler=0

      do j=1,4
         if (isHadron(partOut(j)%ID)) then
            mass(j) = hadron(partOut(j)%ID)%mass
         else
            write(*,*) 'setKinematics,Do4Body: unknown ID:',partOut(j)%ID
            call Traceback()
         end if
         partOut(j)%mass = mass(j)
      end do
      do i=1, maxCorrectLoop

         p4 = momenta_in_4BodyPS (srts_vacuum, mass)

         do j=1,4
            partOut(j)%mom(1:3) = p4(1:3,j)
            partOut(j)%mom(0)   = sqrt(mass(j)**2 + &
                 DOT_PRODUCT(partOut(j)%mom(1:3),partOut(j)%mom(1:3)))

            partOut(j)%vel = partOut(j)%mom(1:3)/partOut(j)%mom(0)
            if (.not. checkVelo(partOut(j))) then
               write(*,*) 'SetKinematics,do4Body: checkVelo failed!'
               collFlag=.false.
               partOut(1:3)%ID=0
               return
            end if
         end do


!!$         call WriteParticle(6)
!!$         call WriteParticle(6,99,1,partIn(1))
!!$         call WriteParticle(6,99,2,partIn(2))
!!$         write(*,*) !!$         do j=1,4
!!$            call WriteParticle(6,99,j,partOut(j))
!!$         end do
!!$
!!$         write(*,*) 'Boost:',betaToCM



         if (correctEnergy) then
            call energyCorrection(srts,betaToCM, &
                 & partOut, successFlag,potFailure)
            if (successFlag) exit
            if (potFailure)  exit  ! Failed due to threshold violation
         else
            ! just boost to Calculation frame
            successFlag=.true.
            do j=1,4
               call lorentz(-betaToCM,partOut(j)%mom)
            end do
         end if

      end do


!!$      write(*,*) '****** final ******'
!!$      do j=1,4
!!$         call WriteParticle(6,99,j,partOut(j))
!!$      end do
!!$      call Traceback()

      if (.not.successFlag) then
         write(*,'(A,2I4,A,4I4,A,2ES12.4)') &
              'SetKinematics[4]: Energy correction failed.', &
              partIn(1:2)%Id,' ->',partOut(1:4)%Id,' @',srts,srts_vacuum

         fehlerZaehler=fehlerzaehler+1
         if (fehlerZaehler>fehler_max) then
            write(*,*) 'too many errors in master_2Body: Do4Body'
            call Traceback()
         end if

         collFlag=.false. ! Kill Event
         partOut(1:4)%ID=0
         return ! --> Exit
      end if

    end subroutine Do4Body

  end subroutine setKinematics


  !****************************************************************************
  !****s* master_2Body/setKinematicsHiEnergy
  ! NAME
  ! subroutine setKinematicsHiEnergy(srts, srts_vacuum, sigs,
  ! betaToCM, partIn, time, partOut, collFlag,
  ! HiEnergyType, scenario0)
  !
  ! PURPOSE
  ! Evaluates the kinematics for the "finalState" particles assuming
  ! vacuum kinematics (-> no potentials) using FRITIOF or PYTHIA.
  !
  ! INPUTS
  ! * type(particle),dimension(1:2):: partIn      -- incoming initial particles
  ! * type(particle), dimension(:) :: partOut     -- outgoing particles
  ! * real                         :: srts           -- sqrt(s)
  ! * real                         :: srts_vacuum    -- sqrt(s) in the vacuum
  ! * real, dimension(1:3)         :: betaToCM       -- beta of calc frame to CM frame
  ! * real, dimenssion(0:7)        :: sigs           -- XS, cf. isigXXX
  ! * real                         :: time           -- time of event
  ! * integer, OPTIONAL            :: scenario0      -- if given : BarBar,BarMes,MesMes
  !
  ! OUTPUT
  ! * type(particle), dimension(:) :: partOut
  ! * logical                      :: collFlag  -- "true" if kinematics could be set
  ! * integer                      :: HiEnergyType    -- (see below)
  !
  ! possible values for "HiEnergyType":
  ! *  -3: BaB    (Baryon-Antibaryon-Annihilation)
  ! *  -2: Manni  (Meson-Baryon-Annihilation)
  ! *  -1: Elastic or charge exchange
  ! *   0: == Low Energy ==
  ! *   1: FRITIOF
  ! *   2: PYTHIA
  !
  ! NOTES
  ! 0) Formerly known as 'setKinematicsFritiof'.
  ! 1) there are some "hardwired" switches:
  !    * DoFritiofCor (=.true.) -- Switch to vary the srts-description
  ! 2) The cross sections CEX, Anni, Lbar, Sbar, XiBar and JPsi are
  !    relevant only for baryon-antibaryon collisions.
  !    Otherwise they are set to zero.
  !****************************************************************************
  subroutine setKinematicsHiEnergy(srts, srts_vacuum, sigs, &
       betaToCM, partIn, time, partOut, collFlag, &
       HiEnergyType, scenario0)

    use random, only: rn
    use particleDefinition
    use PythiaSpecFunc, only: Init_VM_Mass
    use LorentzTrafo, only: lorentz, lorentzCalcBeta
    use IdTable, only: JPsi, isMeson, OmegaResonance
!    use output, only: WriteParticle
    use RMF, only: getRMF_flag
    use Annihilation, only: annihilate, annihilate_to_meson
    use monteCarlo, only: MonteCarloChoose
    use CallStack, only: Traceback
    use twoBodyTools, only: LogicMatrix

    use Coll_Elastic, only: DoColl_Elast
    use Coll_ChEx, only: DoColl_ChEx
    use Coll_Manni, only: DoColl_ManniProb, DoColl_Manni, DoColl_ManniCheck
    use Coll_Pythia, only: DoColl_Pythia
    use Coll_Fritiof, only: DoColl_Fritiof
    use Coll_HypAntiHyp, only: DoColl_YYbar

    real,                         intent(in)    :: srts, srts_vacuum
    real, dimension(0:7),         intent(in)    :: sigs
    real, dimension(1:3),         intent(in)    :: betaToCM
    type(particle),dimension(1:2),intent(in)    :: partIn
    real,                         intent(in)    :: time
    type(particle), dimension(:), intent(inOut) :: partOut
    logical,                      intent(out)   :: collFlag
    integer,                      intent(out)   :: HiEnergyType
    integer, OPTIONAL,            intent(in)    :: scenario0

    logical, parameter :: DoFritiofCor=.true.

    real :: srtsCor
    real, dimension(1:3) :: betaCor
    real, dimension(0:3) :: ptot,pcm
    logical, parameter :: debugger=.false.
    integer :: scenario
    real, dimension(1:8) :: sigsTemp ! isigXXX + 1


    ! Initialize switches at first call
    if (initFlag) call ReadInput

    HiEnergyType  = 0
    collFlag=.FALSE.

    if (present(scenario0)) then
       scenario = scenario0
    else
       scenario = LogicMatrix(isMeson(partIn(1)%ID), isMeson(partIn(2)%ID))
    end if


    ! Define CM - Momentum :
    pcm=partIn(1)%mom
    call lorentz(betaToCM, pcm)


    srtscor=srts_vacuum
    betacor=betaToCM

    if ( .not.getRMF_flag() ) then
       if (srts.lt.sum(partIn%mass)  .or. DoFritiofCor) then

          !  the following must be done because the potential cannot be
          !  accounted for in the FRITIOF/PYTHIA event generation
          !  (as a result energy and momentum are not exactly conserved)

          pTot(0) = SQRT(partIn(1)%mass**2+Sum(partIn(1)%mom(1:3)**2)) &
               &  + SQRT(partIn(2)%mass**2+Sum(partIn(2)%mom(1:3)**2))
          pTot(1:3) = partIn(1)%mom(1:3)+partIn(2)%mom(1:3)
          srtsCor=SQRT(pTot(0)**2-Dot_Product(pTot(1:3),pTot(1:3)))

          betaCor = lorentzCalcBeta(pTot)
       else
          srtscor=srts
       end if
    else

       srtscor=srtS_vacuum
       betacor=betaToCM

    end if

    if (debugger) &
         write(*,*) 'In setKinematicsHiEnergy, SRTS, srtS_vacuum, srtscor:',&
         & srtS, srtS_vacuum, srtscor

    ! save srtscor & position for calls of VM_Mass
    call Init_VM_Mass(srtscor,(partIn(1)%pos+partIn(2)%pos)/2.)

    ! Choose whether to make elastic, charge exchange or inelastic event :
    !
    ! At this place we only subtract elastic and charge exchange cross
    ! sections since this is common for all possible collisions.
    ! Further subtractions might be done later-on specifically for every
    ! channel (baryon-meson, baryon-baryon or baryon-antibaryon).
    ! Of course one should very carefully set the corresponding  cross
    ! sections at the input to this routine.

    sigsTemp(1:2) = sigs(1:2)
    sigsTemp(3) = sigs(0) - sum(sigs(1:2)) ! = sigmaInel

    select case (MonteCarloChoose(sigsTemp(1:3)))
    case (1) ! === Elast ===
       if (debugger) write(*,*) 'Elastic event',srtsCor,sigs(0:2)

       call DoColl_Elast(partIn,partOut,collFlag, &
            srtscor,pcm,betacor, ElastAngDist)
       HiEnergyType=-1
       return

    case (2) ! === CEX ===
       if (debugger) write(*,*) 'Charge exchange event',srtsCor,sigs(0:2)

       call DoColl_ChEx(partIn,partOut,collFlag, srtscor,pcm,betacor, 3)
       HiEnergyType=-1
       return

    case (3) ! === Inelastic ===

       ! do nothing, just follow below

    case default

       call TRACEBACK("wrong value in MC choice 1")

    end select

    ! === Inelastic ===
    ! Do inelastic events...

    if (debugger) write(*,*) 'InElastic event', srtsCor, sigs(0:1)

    ! Set the production time for all finalState:

    partOut%prodTime = time
    partOut%inF = .TRUE.

    select case (scenario)

    case (scenarioBarMes)
       if (debugger) write(*,*) 'Baryon--Meson'

       if (useManni) then
          if (DoColl_ManniProb(partIn,srtscor) > rn()) then
             if (debugger) write(*,*) 'Vor Manni'
             call DoColl_Manni(partIn,partOut,collFlag, &
                  srtscor,pcm,betacor)
             if (collFlag) &
                  call DoColl_ManniCheck(partIn,partOut,collFlag)
             if (debugger) write(*,*) 'Nach Manni'
             HiEnergyType = -2
             return
          end if
       end if

       ! Select whether Fritiof or PYTHIA should be called:
       if (usePythia==1) then
          HiEnergyType = 2
          if (debugger) write(*,*) 'Vor Coll_Py'
          call DoColl_Pythia(partIn, partOut, collFlag, &
               srtscor, pcm, betacor)
          if (debugger) write(*,*) 'Nach Coll_Py'
       else
          HiEnergyType = 1
          if (debugger) write(*,*) 'Vor Coll_Fr'
          call DoColl_Fritiof(partIn,partOut,collFlag, &
               srtscor,pcm,betacor)
          if (debugger) write(*,*) 'Nach Coll_Fr'
       end if
       return

    case (scenarioBarBar)

       select case (LogicMatrix(partIn(1)%anti,partIn(2)%anti))

       case (scenarioBarBar)

          if (debugger) write(*,*) 'BaryonBaryon'

          ! Select whether Fritiof or PYTHIA should be called:
          if (usePythia==1) then
             HiEnergyType = 2
             if (debugger) write(*,*) 'Vor Coll_Py'
             call DoColl_Pythia(partIn, partOut, collFlag, srtscor, pcm, betacor)
             if (debugger) write(*,*) 'Nach Coll_Py'
          else
             HiEnergyType = 1
             if (debugger) then
                write(*,*) ' Vor Coll_Fr:'
                write(*,*) ' Ids, charges:', partIn%Id, partIn%charge
                write(*,*) ' masses:', partIn%mass
                write(*,*) ' srtscor:',  srtscor
             end if
             call DoColl_Fritiof(partIn,partOut,collFlag, srtscor,pcm,betacor)
             if (debugger) write(*,*) 'Nach Coll_Fr'
          end if
          return

       case (scenarioBarAntiBar)

          if (debugger) write(*,*) 'BaryonAntiBaryon'

          ! Do Pythia/Fritiof or BaB event:

          sigsTemp(1:7) = sigs(1:7)
          sigsTemp(8) = sigs(0) - sum(sigs(1:7)) ! = sigmaProd
          sigsTemp(1:2) = 0.0 ! elast and CEX already handled, here only inelast

          select case (MonteCarloChoose(sigsTemp(3:8))+2)

          case (3) ! === Annihilation ===
             if (partIn(1)%anti) then
                call annihilate(partIn(1),partIn(2),time,partOut,collFlag,HiEnergyType)
             else
                call annihilate(partIn(2),partIn(1),time,partOut,collFlag,HiEnergyType)
             end if

          case (4) ! === LambdaBar+Lambda ===
             if (debugger) write(*,*) 'Vor DoColl_YYbar,1'
             call DoColl_YYbar(partIn,partOut,collFlag,srtscor,pcm,betacor,1)
             if (debugger) write(*,*) 'Nach DoColl_YYbar,1'

          case (5) ! === Lambda+SigmaBar, LambdaBar+Sigma ===
             if (debugger) write(*,*) 'Vor DoColl_YYbar,2'
             call DoColl_YYbar(partIn,partOut,collFlag,srtscor,pcm,betacor,2)
             if (debugger) write(*,*) 'Nach DoColl_YYbar,2'

          case (6) ! === XiBar+Xi ===
             if (debugger) write(*,*) 'Vor DoColl_YYbar,3'
             call DoColl_YYbar(partIn,partOut,collFlag,srtscor,pcm,betacor,3)
             if (debugger) write(*,*) 'Nach DoColl_YYbar,3'

          case (7) ! === Jpsi ===
             if (debugger) write(*,*) 'Vor J/Psi prod'
             partOut(1)%id=JPsi
             partOut(1)%charge=0
             call annihilate_to_meson(partIn,time,partOut,collFlag,HiEnergyType)
             if (debugger) write(*,*) 'Nach J/Psi prod'

          case (8) ! === Prod ===
             ! Select whether Fritiof or PYTHIA should be called:
             ! (for Omegas we always call Pythia, because Fritiof produces errors)
             if (usePythia_BaB==1 .or. partIn(1)%ID==OmegaResonance .or. &
                & partIn(2)%ID==OmegaResonance) then
                HiEnergyType = 2
                if (debugger) write(*,*) 'Vor Coll_Py'
                call DoColl_Pythia(partIn, partOut, collFlag, srtscor, pcm, betacor)
                if (debugger) write(*,*) 'Nach Coll_Py'
             else
                HiEnergyType = 1
                if (debugger) write(*,*) 'Vor Coll_Fr'
                call DoColl_Fritiof(partIn,partOut,collFlag, srtscor,pcm,betacor)
                if (debugger) write(*,*) 'Nach Coll_Fr'
             end if

          case default
             call TRACEBACK("wrong value in MC choice 2")
          end select

          return

       case (scenarioAntiBarAntiBar)

          if (debugger) write(*,*) 'AntiBaryonAntiBaryon'
          ! --- not yet implemented
          return

       end select

    end select


  end subroutine setKinematicsHiEnergy


  !****************************************************************************
  !****s* master_2Body/XsectionMaster
  ! NAME
  ! subroutine XsectionMaster(srts, partIn, mediumAtColl, momLRF,&
  ! partOut, sigs, HiEnergyFlag, plotFlag, scenario0, ForceHiEnergy,&
  ! pauliIncluded_out)
  !
  ! PURPOSE
  ! Gets as input two colliding particles and returns the total, elastic and
  ! annihilation cross sections.
  ! In case of resonance model processes it also already
  ! returns the final state particles.
  !
  ! INPUTS
  ! * real                          :: srts         -- sqrt(s) of the process
  ! * type(particle),dimension(1:2) :: partIn       -- colliding particles
  ! * type(medium)                  :: mediumAtColl -- Medium informations
  ! * real,dimension(0:3)           :: momLRF       -- Total Momentum in LRF
  ! * logical, optional             :: plotFlag     -- Switch plotting the XS
  ! * integer, optional             :: scenario0    -- if given: BarBar,BarMes,MesMes
  ! * logical, optional             :: ForceHiEnergy -- if given: force scenario selection
  !
  ! OUTPUT
  ! * type(preEvent),dimension(:)   :: partOut       -- produced particles
  ! * real, dimension(0:7)          :: sigs          -- XS, cf. isigXXX
  ! * logical                       :: HiEnergyFlag  -- true if HiEnergy is
  !   used to generate the final state
  ! * logical,optional              :: pauliIncluded_out -- true if Pauli
  !   blocking was already included in the cross section evaluation
  ! NOTES
  ! The cross sections for CEX, Anni, Lbar, Sbar and XiBar are relevant only for
  ! baryon-antibaryon collisions. Otherwise they are set to zero.
  !****************************************************************************
  subroutine XsectionMaster(srts, partIn, mediumAtColl, momLRF, &
       partOut, sigs, HiEnergyFlag, plotFlag, scenario0, ForceHiEnergy,&
       pauliIncluded_out)

    use mediumDefinition
    use particleDefinition
    use preEventDefinition
    use IdTable, only: nucleon, JPsi, isMeson
    use twoBodyTools, only: LogicMatrix

    real,                           intent(in)  :: srts
    type(particle), dimension(1:2), intent(in)  :: partIn
    type(medium),                   intent(in)  :: mediumAtColl
    real, dimension(0:3),           intent(in)  :: momLRF

    type(preEvent), dimension(:),   intent(out) :: partOut
    real, dimension(0:7),           intent(out) :: sigs
    logical,                        intent(out) :: HiEnergyFlag

    logical, optional,              intent(in)  :: plotFlag
    integer, optional,              intent(in)  :: scenario0
    logical, optional,              intent(in)  :: ForceHiEnergy
    logical, optional,              intent(out) :: pauliIncluded_out

    logical :: plotIt, pauliIncluded
    integer :: scenario

    ! Initialize switches at first call
    if (initFlag) call ReadInput

    if (present(plotFlag)) then
       plotIT = plotFlag
    else
       plotIT = .false.
    end if

    pauliIncluded = .false.
    if (present(pauliIncluded_out)) pauliIncluded_out = .false.

    ! Initialize output (preEvent)
    partOut(:)%ID=0
    partOut(:)%charge=0
    partOut(:)%anti=.false.
    partOut(:)%mass=0.

    sigs = 0.

    ! Check kind of collision

    ! This is the original code, but leads to segmentation errors:
!!$    if ((present(scenario0)).and.(scenario0.gt.-1000)) then
!!$       scenario = scenario0
!!$    else
!!$       scenario = LogicMatrix(isMeson(partIn(1)%ID), isMeson(partIn(2)%ID))
!!$    end if
    ! This is the replacement:
    scenario = -1000
    if (present(scenario0)) scenario = scenario0
    if (scenario.le.-1000) then
       scenario = LogicMatrix(isMeson(partIn(1)%ID), isMeson(partIn(2)%ID))
    end if



    select case (scenario)
    case (scenarioBarBar)
       if (present(ForceHiEnergy)) then
          HiEnergyFlag = ForceHiEnergy
       else if (partIn(1)%anti.neqv.partIn(2)%anti) then
          HiEnergyFlag = DecideHiEnergy(srts, HiEnergyThresholdBarAntibar, &
               &                              HiEnergyThresholdBarAntibarDelta)
       else
          HiEnergyFlag = DecideHiEnergy(srts, HiEnergyThresholdBarBar, &
               &                              HiEnergyThresholdBarBarDelta)
       end if
       call BarBar(srts, partIn(1:2), mediumAtColl, &
            partOut(1:4), sigs, PauliIncluded)

    case (scenarioBarMes)
       if (present(ForceHiEnergy)) then
          HiEnergyFlag = ForceHiEnergy
       else if (    max(partIn(1)%id,partIn(2)%id).eq.JPsi &
            & .and. min(partIn(1)%id,partIn(2)%id).eq.nucleon ) then
          HiEnergyFlag = .false.
       else
          HiEnergyFlag = DecideHiEnergy(srts, HiEnergyThresholdBarMes, &
               &                              HiEnergyThresholdBarMesDelta)
       end if
       call BarMes(srts, partIn(1:2), mediumAtColl, momLRF, &
            partOut(1:4), sigs(0:1) )

    case (scenarioMesMes)
       HiEnergyFlag = .false.
       call MesMes(srts, partIn(1:2), mediumAtColl, momLRF, &
            partOut(1:3), sigs(0:1) )

    end select


    if (present(pauliIncluded_out)) pauliIncluded_out = pauliIncluded


  contains

    !**************************************************************************
    !****if* XsectionMaster/DecideHiEnergy
    ! NAME
    ! logical function DecideHiEnergy(srts, srts0, dsrts)
    ! PURPOSE
    ! Decide whether to use high-energy prescription or not.
    ! For the smooth transition, the logic is the following:
    ! * srts < srts0-dsrts : .false.
    ! * srts > srts0+dsrts : .true.
    ! * in between: linear increase of probability
    ! INPUTS
    ! * real :: srts -- sqrt(s) value
    ! * real :: srts0 -- threshold
    ! * real :: dsrts -- width of transition region
    ! OUTPUT
    ! function value
    !**************************************************************************
    logical function DecideHiEnergy(srts, srts0, dsrts)
      use random, only: rn

      real, intent(in) :: srts, srts0, dsrts

      if (.not. useHiEnergy) then
        DecideHiEnergy = .false.
        return
      end if

      if (srts < srts0-dsrts) then
        DecideHiEnergy = .false.
        return
      end if

      if (srts > srts0+dsrts) then
        DecideHiEnergy = .true.
        return
      end if

      DecideHiEnergy = ((srts-(srts0-dsrts)) > rn()*2.*dsrts)

    end function DecideHiEnergy

    !**************************************************************************
    !****if* XsectionMaster/mesMesResMass
    ! NAME
    ! real function mesMesResMass(partIn, ID)
    !
    ! PURPOSE
    ! Calculate the vacuum mass of the resonance
    ! INPUTS
    ! * type(particle), dimension(1:2) :: partIn -- Incoming particles
    ! * integer :: ID -- ID of outgoing resonance
    ! NOTES
    ! * we have to subtract potentials before calculating the mass; This is
    !   also true for mesons, if one considers Coulomb
    ! * this is analogous to 'resonanceCrossSections/resonanceMass' for barMes
    ! * sqrts_free gives different results
    ! * The Coulomb potential of the final resonance is subtracted, not that
    !   of the incoming particles. This is the potential, which is added again
    !   in 'energyDetermination'. (Coulomb for pi and rho may be different due
    !   to different masses!)
    ! * without Coulomb, it just returns 'sqrtS(partIn)'
    !**************************************************************************
    real function mesMesResMass(partIn, ID)

      use coulomb, only: emfoca

      type(particle),dimension(1:2), intent(in)  :: partIn
      integer, intent(in) :: ID

      integer :: IQ
      real, dimension(0:3) :: mom
      real, dimension(1:3) :: pos

      IQ = Sum(partIn%charge)
      mom =  partIn(1)%mom + partIn(2)%mom
      pos = (partIn(1)%pos + partIn(2)%pos)/2

      if (IQ.ne.0) then
         mom(0) = mom(0) - emfoca(pos, mom(1:3), IQ, ID)
      end if

      mesMesResMass = mom(0)**2 - sum(mom(1:3)**2)
      if (mesMesResMass > 0.) then
         mesMesResMass = sqrt(mesMesResMass)
      else
         mesMesResMass = 0.
      end if

    end function mesMesResMass


    !**************************************************************************
    !****is* XsectionMaster/mesMes
    ! NAME
    ! subroutine mesMes(srts, partIn, mediumAtColl, momLRF, partOut, sigs)
    !
    ! PURPOSE
    ! Calculate the total and elastic cross section of meson-meson collisions.
    !
    ! INPUTS
    ! * real                          :: srts  -- sqrt(s) of the process
    !   (in vacuum!)
    ! * type(particle),dimension(1:2) :: partIn       -- colliding particles
    ! * type(medium)                  :: mediumAtColl -- Medium informations
    ! * real,dimension(0:3)           :: momLRF       -- Total Momentum in LRF
    !
    ! OUTPUT
    ! * type(preEvent),dimension(:)   :: partOut      -- produced particles
    ! * real,dimension(0:1)           :: sigs         -- XS, cf. isigXXX
    !
    ! NOTES
    ! original pi pi <-> K Kbar from Martin Effenberger,
    ! other channels added by Markus Wagner
    !**************************************************************************
    subroutine mesMes(srts, partIn, mediumAtColl, momLRF, partOut, sigs)

      use DecayChannels, only: Decay2BodyMeson, nDecay2bodyMeson
      use clebschGordan, only: clebschSquared
      use IdTable, only: pion, kaon, kaonBar, kaonStar, kaonStarBar, rho, &
           omegaMeson, phi, sigmaMeson, etaPrime, nmes
      use mesonWidthMedium, only: decayWidthMesonMedium, WidthMesonMedium
      use random, only: rn
      use mesonMeson_Xsections, only: sig_pipi2kkbar, kkbar_cross, &
           kkbar_out, kstarkbar_cross
      use particleProperties, only: hadron, nDecays
      use constants, only: pi, mK, GeVSquared_times_mb
      use mesonPotentialMain, only: vecMes_massShift
      use twoBodyTools, only: pCM_sqr
      use monteCarlo, only: MonteCarloChoose
      use CallStack, only: Traceback

      real,                          intent(in)  :: srts
      type(particle),dimension(1:2), intent(in)  :: partIn
      type(medium),                  intent(in)  :: mediumAtColl
      real,dimension(0:3),           intent(in)  :: momLRF

      type(preEvent),dimension(1:3), intent(out) :: partOut
      real, dimension(0:1),          intent(out) :: sigs


      real, dimension(1:3) :: sigsTemp ! 1:3 = (Elast, Res, BG)

      integer :: izt,ichannel,k,idRes
      logical :: hasResonant
      real :: is1,is2,is3,iz1,iz2, multipl, isofac
      real :: pinitial2,m0,mR,gamtot, srts0
      real, dimension(pion:pion+nmes-1) :: sigRes, mRes

      real, dimension(1:nDecays) :: width_meson
      real, dimension(1:21) :: msigbg

      real, parameter :: const = 2.  ! assumed constant cross-section (mbarn)
                                     ! for KKbar production

      sigs = 0.
      sigsTemp = 0.
      sigRes = 0.

      if (debug) then
         write(*,*) 'In mesmes:'
         write(*,*) partIn%ID,partIn%mass,srts
      end if

      if ((OverideSigma_PiPi > 0).and.(partIn(1)%ID==pion) &
           .and. (partIn(2)%ID==pion)) goto 111

      !    Determine input channel, i.e. decay channel for a resonance:
      hasResonant = .false.
      do ichannel=1,nDecay2bodyMeson
         if ((partIn(1)%ID==Decay2BodyMeson(ichannel)%ID(1) .and. &
            & partIn(2)%ID==Decay2BodyMeson(ichannel)%ID(2)) .or. &
             (partIn(1)%ID==Decay2BodyMeson(ichannel)%ID(2) .and. &
            & partIn(2)%ID==Decay2BodyMeson(ichannel)%ID(1))) then
            hasResonant = .true.
            exit
         end if
      end do

      !    c.m. momentum of colliding particles:
      pinitial2 =  pCM_sqr( srts**2, partIn(1)%mass**2, partIn(2)%mass**2 )

      !    Isospins of incoming mesons:
      is1 = 0.5*float(hadron(partIn(1)%ID)%isospinTimes2)
      is2 = 0.5*float(hadron(partIn(2)%ID)%isospinTimes2)

      !    z-component of the isospins of incoming mesons
      !    by a Gell-Mann - Nishijima formula:
      iz1 = 0.5*float(2*partIn(1)%charge &
           -hadron(partIn(1)%ID)%strangeness &
           -hadron(partIn(1)%ID)%charm)
      iz2 = 0.5*float(2*partIn(2)%charge &
           -hadron(partIn(2)%ID)%strangeness &
           -hadron(partIn(2)%ID)%charm)

      !    Total charge:
      izt = Sum(partIn%charge)

      !    Resonance contribution to the cross section:
      if (hasResonant) then

         multipl = (2.*hadron(partIn(1)%ID)%spin+1.) &
              &  * (2.*hadron(partIn(2)%ID)%spin+1.)

         do idRes = pion,pion+nmes-1

            if (.not.hadron(idRes)%propagated) cycle

            mR = mesMesResMass(partIn, idRes)
            mRes(idRes) = mR

            if (idRes==rho .or. idRes==omegaMeson .or. idRes==phi) then
               ! possible in-medium mass shift of the vector mesons:
               m0 = hadron(idRes)%mass &
                    + vecMes_massShift(idRes, mediumAtColl%density)
            else
               m0 = hadron(idRes)%mass
            end if
            do k=1,nDecays
               if (hadron(idRes)%decaysID(k)/=ichannel) cycle
               if (hadron(idRes)%decays(k)<=0.001) cycle

               is3 = 0.5*float(hadron(idRes)%isospinTimes2) ! isospin of meson
               isofac = clebschSquared(is1,is2,is3,iz1,iz2) ! isospin factor

               ! In the case of two identical incoming particles
               ! the isospin factor must be corrected:
               if (partIn(1)%ID==partIn(2)%ID) isofac = 2.*isofac

               width_meson = decayWidthMesonMedium(idRes, mR, izt)
               ! get the in-width
               gamtot = WidthMesonMedium(idRes, mR, momLRF, mediumAtColl)
               ! total in-medium width (incl. collisional width)

               sigRes(idRes) = 4.*pi/pinitial2*isofac &
                    * (2.*hadron(idRes)%spin+1.)/(multipl*GeVSquared_times_mb) &
                    * width_meson(k) * gamtot * mR**2 &
                    / ((mR**2-m0**2)**2+(gamtot*mR)**2)

               sigsTemp(isigRes) = sigsTemp(isigRes) + sigRes(idRes)

               exit ! --> next idRes
            end do

         end do ! loop over resonances

      end if

      ! Background
      ! (sum over final isospin states is assumed in computing sigbgt):
      srts0 = 2.*mK

      if ((partIn(1)%ID==pion .or. partIn(1)%ID==rho) .and. &
          (partIn(2)%ID==pion .or. partIn(2)%ID==rho)) then

         ! === pi pi, rho rho --> K Kbar  or pi rho --> K^* Kbar, K Kbar^*

         if (srts>srts0) then
            !        determine isospin factor:
            if (abs(izt)==1) then
               isofac=0.5
            else if (izt==0) then
               if (partIn(1)%charge==0) then
                  isofac=1./3.
               else
                  isofac=5./6.
               end if
            else
               isofac=0.
            end if
            if (isofac>0.) then
               if (partIn(1)%ID/=partIn(2)%ID) then
                  ! for pi rho channel there are 2 times more outgoing
                  ! states because of K^* or Kbar^*:
                  srts0 = mK+hadron(kaonStar)%minmass
                  sigsTemp(isigBG) = 2. * isofac*9./4.*sig_pipi2kkbar(srts,srts0)
               else
                  sigsTemp(isigBG) = isofac*9./4.*sig_pipi2kkbar(srts,srts0)
               end if
            end if
         end if

      else if (partIn(1)%ID>=pion .and. partIn(1)%ID<=etaPrime .and. &
               partIn(2)%ID>=pion .and. partIn(2)%ID<=etaPrime) then

         ! === nonstrange + nonstrange (other than pi pi, pi rho, rho rho)
         ! --> K + Kbar, K^* Kbar, K Kbar^*:

         if (srts>srts0 .and. abs(izt)<2) sigsTemp(isigBG) = const

      else if ((partIn(1)%ID==kaon .and. partIn(2)%ID==kaonBar) .or. &
               (partIn(2)%ID==kaon .and. partIn(1)%ID==kaonBar)) then

         ! === K Kbar --> pi pi, pi rho ...

         call kkbar_cross(izt,srts,pinitial2,const,msigbg,mesMes_useWidth)

         sigsTemp(isigBG) = Sum(msigbg)

      else if ((partIn(1)%ID==kaon     .and. partIn(2)%ID==kaonStarBar) .or. &
               (partIn(1)%ID==kaonStar .and. partIn(2)%ID==kaonBar)) &
           then

         ! === K Kbar^*, K^* Kbar --> nonstrange mesons

         call kstarkbar_cross(izt,srts,pinitial2,const,msigbg,mesMes_useWidth)

         sigsTemp(isigBG) = Sum(msigbg)

      end if


      ! Correction because of positive parity of sigma:
      if ((partIn(1)%ID==sigmaMeson .or. partIn(2)%ID==sigmaMeson) &
          .and. partIn(1)%ID/=partIn(2)%ID) then
         sigsTemp(isigBG) = 0.
      end if

      ! here you can switch off the 2 to 2 interactions:
      ! (could be optimized: do not calculate the cross sections,
      ! if you know you will throw them anyhow)
      if (.not.mesMes_do2to2) then
         sigsTemp(isigBG) = 0.
      end if

      ! Correction for one incoming meson with spin=1:
      if (srts<mK+hadron(kaonStar)%mass &
           .and. hadron(partIn(1)%ID)%spin+hadron(partIn(2)%ID)%spin==1. &
           .and. hadron(partIn(1)%ID)%strangeness==0 &
           .and. hadron(partIn(2)%ID)%strangeness==0) then
         sigsTemp(isigBG) = 0.
      end if

111   continue

      if ((OverideSigma_PiPi > 0) &
           .and. (partIn(1)%ID==pion) .and. (partIn(2)%ID==pion)) then
         sigsTemp(isigRes)   = 0.0
         sigsTemp(isigBG)    = 0.0
         sigsTemp(isigElast) = OverideSigma_PiPi
      end if

      if ((Overide_PiPi_ResIsElast) &
           .and. (partIn(1)%ID==pion) .and. (partIn(2)%ID==pion)) then
         sigsTemp(isigElast) = sigsTemp(isigRes)
         sigsTemp(isigRes)   = 0.0
      end if

      sigs(isigElast) = sigsTemp(isigElast)

      select case (MonteCarloChoose( sigsTemp, sigs(isigTot) ))
      case (0) ! === failure ===
         if (sigs(isigTot)==0.) then
            partOut(1:3)%ID = 0
            return
         end if
         call TRACEBACK("strange failure")

      case (1) ! === elastic ===

         partOut(1:2)%ID = partIn(1:2)%ID
         partOut(1:2)%charge = partIn(1:2)%charge
         partOut(1:2)%anti = partIn(1:2)%anti
         partOut(1:2)%mass = partIn(1:2)%mass

      case (2) ! === res ===

         idRes = MonteCarloChoose(sigRes) ! select resonance
         if (idRes==0) call TRACEBACK("iRes = 0")
         idRes = idRes + pion-1  ! set it to a meson resonace number!

         !      Set final state:
         partOut(1:3)%ID = (/ idRes, 0, 0 /)
         if (idRes==rho .or. idRes==omegaMeson .or. idRes==phi) then
            partOut(1)%mass = mRes(idRes)&
                 - vecMes_massShift(idRes, mediumAtColl%density)
         else
            partOut(1)%mass = mRes(idRes)
         end if
         partOut(1)%charge = izt
         partOut(1)%anti = .false.

      case (3) ! === background ===

         partOut(3)%ID = 0

         if (partIn(1)%ID>=pion .and. partIn(1)%ID<=etaPrime .and. &
             partIn(2)%ID>=pion .and. partIn(2)%ID<=etaPrime) then

            !        nonstrange + nonstrange --> K + Kbar, K^* Kbar, K Kbar^*:

            partOut(1:2)%ID = (/ kaon, kaonBar /)

            !        Correction for one incoming vector meson:
            if (hadron(partIn(1)%ID)%spin+hadron(partIn(2)%ID)%spin==1.) then
               if (rn()<=0.5) then
                  partOut(1)%ID = kaonStar
               else
                  partOut(2)%ID = kaonStarBar
               end if
            end if

            !        Charge assignment:
            if (izt==0) then
               if (rn()<0.5) then
                  partOut(1:2)%charge = (/0,0/)
               else
                  partOut(1:2)%charge = (/1,-1/)
               end if
            else if (izt==1) then
               partOut(1:2)%charge = (/1,0/)
            else if (izt==-1) then
               partOut(1:2)%charge = (/0,-1/)
            else
               write(*,*) 'error in mesmes: izt=', izt
               write(*,*) partIn%ID
               write(*,*) partIn%charge
               call Traceback()
            end if

         else if ( (partIn(1)%ID==kaon .and. partIn(2)%ID==kaonBar) &
              .or. (partIn(2)%ID==kaon .and. partIn(1)%ID==kaonBar)) then
            !        K Kbar --> pi pi, pi rho ...
            call kkbar_out(izt,msigbg,sigsTemp(isigBG),partOut)

         else if ( (partIn(1)%ID==kaon     .and. partIn(2)%ID==kaonStarBar) &
              .or. (partIn(1)%ID==kaonStar .and. partIn(2)%ID==kaonBar)) then
            !        K Kbar^*, K^* Kbar --> nonstrange mesons
            call kkbar_out(izt,msigbg,sigsTemp(isigBG),partOut)

         end if

         partOut(1:2)%mass = hadron(partOut(1:2)%ID)%mass
         partOut(1:2)%anti = .false.

      end select


    end subroutine mesMes


    !**************************************************************************
    !****is* XsectionMaster/barMes
    ! NAME
    ! subroutine barMes(srts, partIn, mediumAtColl, momLRF, partOut, sigs)
    !
    ! PURPOSE
    ! Calculate the total and elastic cross section of baryon-meson collisions.
    !
    ! INPUTS
    ! * real                          :: srts     -- sqrt(s) of the process
    ! * type(particle),dimension(1:2) :: partIn       -- colliding particles
    ! * type(medium)                  :: mediumAtColl -- Medium informations
    ! * real,dimension(0:3)           :: momLRF       -- Total Momentum in LRF
    !
    ! OUTPUT
    ! * type(preEvent),dimension(:)   :: partOut       -- produced particles
    ! * real,dimension(0:1)           :: sigs          -- XS, cf. isigXXX
    !
    ! NOTES
    ! the internal variable "HiEnergyFlag" is input.
    !**************************************************************************
    subroutine barMes(srts, partIn_Orig, mediumAtColl, momLRF, partOut, sigs)

      use mediumDefinition
      use idTable
      use particleDefinition
      use parametrizationBarMes_HighEnergy, only: paramBarMesHE
      use preEventDefinition
      use output, only: writeParticle
      ! The cross sections for the individual channels:
      use pionNucleon, only: pionNuc
      use omegaNucleon, only: omegaNuc
      use phiNucleon, only: phiNuc
      use rhoNucleon, only: rhoNuc
      use etaNucleon, only: etaNuc
      use sigmaNucleon, only: sigmaNuc
      use kaonNucleon, only: kaonNuc
      use antiKaonNucleon, only: kaonBarNuc
      use JPsiNucleon, only: JPsiNuc
      use mesonHyperon, only: mesonY
      use pionP11_1440_resonance, only: pionP11_1440
      use etaDelta_resonance, only: etaDelta
      use rhoDelta_resonance, only: rhoDelta
      use kaonSigma_resonance, only: kaonSigma
      use pionDelta_resonance, only: pionDelta
      use kaonLambda_resonance, only: kaonLambda

      real,                          intent(in)  :: srts
      type(particle),dimension(1:2), intent(in)  :: partIn_Orig
      type(medium),                  intent(in)  :: mediumAtColl
      real, dimension(0:3),          intent(in)  :: momLRF

      type(preEvent),dimension(1:4), intent(out) :: partOut
      real, dimension(0:1),          intent(out) :: sigs


      integer :: mesonID, mesonCharge, BaryonID, baryonCharge
      logical :: antiBaryon = .false.
      type(particle),dimension(1:2)   :: partIn        ! colliding particles

      sigs = 0.

      ! We copy the input, since for the ease of calculation we want to
      ! transform the meson to a particle if it's an antiparticle in the input.
      ! This could not be done otherwise since partIn_Original has
      ! intent(in).
      partIn =partIn_Orig


      if ( isMeson(partIn(1)%ID) ) then
         if (partIn(1)%anti) then
            call WriteParticle(6,99,1,partIn(1))
            call TRACEBACK('There is a meson declared as anti-meson!')
         end if
         mesonID=partIn(1)%ID
         mesonCharge=partIn(1)%Charge

         if (.not. isBaryon(partIn(2)%ID) ) then
            write(*,*) 'id 2 not a baryon!!!', partIn(1:2)%ID
            call TRACEBACK()
         end if

         baryonID=partIn(2)%ID
         antiBaryon=partIn(2)%anti
         baryonCharge=partIn(2)%charge


      else if ( isMeson(partIn(2)%ID) ) then
         if (partIn(2)%anti) then
            call WriteParticle(6,99,2,partIn(2))
            call TRACEBACK('There is a meson declared as anti-meson!')
         end if
         mesonID=partIn(2)%ID
         mesonCharge=partIn(2)%Charge

         if (.not. isBaryon(partIn(1)%ID) ) then
            write(*,*) 'id 1 not a baryon!!!', partIn(1:2)%ID
            call TRACEBACK()
         end if

         baryonID=partIn(1)%ID
         antiBaryon=partIn(1)%anti
         baryonCharge=partIn(1)%charge
      else
         write(*,*) 'Error in barMes. No meson!!!', partIn(1:2)%ID
      end if

      if ((OverideSigma_PiN>=0.) &
           .and. (baryonID==nucleon) .and. (mesonID==pion)) then
         if (.not.HiEnergyFlag) then
            call pionNuc(srts, partIn(1:2), mediumAtColl, &
                 momLRF, partOut(1:4), sigs(isigTot), sigs(isigElast),&
                 plotIT)
         else
            sigs(isigTot) = 99.9
         end if
         if (sigs(isigTot)>0.) then
            sigs(isigTot) = OverideSigma_PiN
            sigs(isigElast) = sigs(isigTot)/10.
         end if
         return
      end if

      if ((OverideSigma_RhoN>=0.) &
           .and. (baryonID==nucleon) .and. (mesonID==rho)) then
         if (.not.HiEnergyFlag) then
            call rhoNuc(srts, partIn(1:2), mediumAtColl, momLRF, &
                 partOut(1:3), sigs(isigTot), sigs(isigElast), &
                 useHiEnergy, HiEnergyThresholdBarMes, plotIT)
         else
            sigs(isigTot) = 99.9
         end if
         if (sigs(isigTot)>0.) then
            sigs(isigTot) = OverideSigma_RhoN
            sigs(isigElast) = sigs(isigTot)/10.
         end if
         return
      end if


      if (HiEnergyFlag) then
         if (antiBaryon) baryonID=-baryonID     ! Setting for "paramBarMesHE"
         ! (old BUU scheme = antiparticles have minus sign in ID)
         call paramBarMesHE(srts,mesonID,baryonID,mesonCharge,baryonCharge, &
              mediumAtColl, sigs(isigTot),sigs(isigElast))
         if (mesonID == omegaMeson .and. omega_K_factor>1.)   &
              sigs(isigTot) = (sigs(isigTot)-sigs(isigElast))*omega_K_factor &
              &                + sigs(isigElast)
         return
      end if


      select case (baryonID)
      case (nucleon)
         select case (mesonID)
         case (pion)
            call pionNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:4), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         case (omegaMeson)
            call omegaNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), &
                 useHiEnergy, HiEnergyThresholdBarMes, plotIT, omega_K_factor)
         case (rho)
            call rhoNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), &
                 useHiEnergy, HiEnergyThresholdBarMes, plotIT)
         case (sigmaMeson)
            call sigmaNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), &
                 useHiEnergy, HiEnergyThresholdBarMes, plotIT)
         case (eta)
            call etaNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), &
                 useHiEnergy, HiEnergyThresholdBarMes, plotIT)
         case (phi)
            call phiNuc(srts, partIn, mediumAtColl, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), &
                 useHiEnergy, HiEnergyThresholdBarMes, plotIT)
         case (Kaon)
            if (.not.antiBaryon) then
               call kaonNuc(srts, partIn, partOut(1:3), &
                    sigs(isigTot), sigs(isigElast), plotIT)
            else
               call kaonBarNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                               & sigs(isigTot), sigs(isigElast), plotIT)
            end if
         case (KaonBar)
            if (.not.antiBaryon) then
               call kaonBarNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                               & sigs(isigTot), sigs(isigElast), plotIT)
            else
               call kaonNuc(srts, partIn, partOut(1:3), sigs(isigTot), sigs(isigElast), plotIT)
            end if
         case (KaonStar)
            if (.not.antiBaryon) then
               ! Not implemented.
!!$          call KaonStarNuc(srts,partIn,mediumAtColl,momLRF,partOut(1:3), &
!!$                          &  sigs(isigTot),sigs(isigElast),plotIT)
            else
               call kaonBarNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                    & sigs(isigTot), sigs(isigElast), plotIT)
            end if
         case (KaonStarBar)
            if (.not.antiBaryon) then
               call kaonBarNuc(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                    & sigs(isigTot), sigs(isigElast), plotIT)
            else
               ! Not implemented.
            end if
         case (JPsi)
            call JPsiNuc(srts,partIn,partOut(1:3),sigs(isigTot),sigs(isigElast),plotIT)
         end select
      case (delta)
         select case (mesonID)
         case (pion)
            call pionDelta(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         case (rho)
            call rhoDelta(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         case (eta)
            call etaDelta(srts, partIn, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         end select
      case (P11_1440)
         select case (mesonID)
         case (pion)
            call pionP11_1440(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                              &sigs(isigTot), sigs(isigElast), plotIT)
         end select
      case (SigmaResonance)
         select case (mesonID)
         case (pion)
            call mesonY(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         case (kaon)
            call kaonSigma(srts, partIn, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         end select
      case (Lambda)
         select case (mesonID)
         case (pion,eta)
            call mesonY(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         case (kaon)
            call kaonLambda(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         end select
      case (Sigma_1385:Sigma_1915)
         select case (mesonID)
         case (pion)
            call mesonY(srts, partIn, mediumAtColl, momLRF, partOut(1:3), &
                 sigs(isigTot), sigs(isigElast), plotIT)
         end select
!!$    case (Lambda_CPlus)
!!$       select case (mesonID)
!!$       case (pion)
!!$          call pionLambda_CPlus(srts,partIn,mediumAtColl,momLRF,  &
!!$          &  partOut(1:3),sigs(isigTot),sigs(isigElast),plotIT)
!!$       End Select
!!$
!!$    case (Xi_C)
!!$       select case (mesonID)
!!$       case (pion)
!!$          call pionXi_c(srts,partIn,mediumAtColl,momLRF,partOut(1:3), &
!!$                       & sigs(isigTot),sigs(isigElast),plotIT)
!!$       End Select

      end select


    end subroutine barMes

    !**************************************************************************
    !****is* XsectionMaster/barBar
    ! NAME
    ! subroutine barBar(srts, partIn, mediumAtColl, partOut, sigs,
    ! PauliIncluded)
    !
    ! PURPOSE
    ! Calculate the total and elastic cross section of baryon-baryon collisions.
    !
    ! INPUTS
    ! * real                          :: srts      -- sqrt(s) of the process
    ! * type(particle),dimension(1:2) :: partIn       -- colliding particles
    ! * type(medium)                  :: mediumAtColl -- Medium informations
    !
    ! OUTPUT
    ! * type(preEvent),dimension(1:4) :: partOut       -- produced particles
    ! * real,dimension(0:7)           :: sigs          -- Xsection, cf. isigXXX
    ! * logical  :: pauliIncluded     -- true if Pauli blocking is included
    !   when calculating the Xsection
    !
    ! NOTES
    ! the internal variable "HiEnergyFlag" is input.
    !**************************************************************************
    subroutine barBar(srts, partIn, mediumAtColl, partOut, &
         sigs, PauliIncluded)

      use mediumDefinition
      use particleDefinition
      use particleProperties, only: hadron
      use barAntiBar, only: sigmaBarAntiBar
      use barBar_Main, only: XsectionBarBar
      use barbar_barbar, only: nukNuk_nukNuk
      use antiBarBar_main, only: XsectionAntiBarBar
      use preEventDefinition
      use constants, only: mN

      real,                          intent(in)  :: srts
      type(particle),dimension(1:2), intent(in)  :: partIn
      type(medium),                  intent(in)  :: mediumAtColl

      type(preEvent),dimension(1:4), intent(out) :: partOut
      real, dimension(0:7),          intent(out) :: sigs
      logical,                       intent(out) :: pauliIncluded

      !real :: facts
      real :: dummy

      partOut%ID=0
      sigs = 0.
      pauliIncluded=.false.

      !*** test -- remove collisions of strange baryons: ********************************************
      !if ( hadron(partIn(1)%ID)%strangeness.ne.0 .or. hadron(partIn(2)%ID)%strangeness.ne.0 ) return
      !**********************************************************************************************

      if (flagElastBB) then
         sigs(isigElast) = 40. ! mb
         sigs(isigTot) = sigs(isigElast)
         partOut(1:2)%Id=partIn(1:2)%Id
         partOut(1:2)%anti=partIn(1:2)%anti
         partOut(1:2)%charge=partIn(1:2)%charge
         return
      end if

HiEn: if (HiEnergyFlag) then
         if (debug) write(*,*) ' barBar: high energy', srts

         ! no collisions of charmed baryons
         if (max(hadron(partIn(1)%ID)%charm,hadron(partIn(2)%ID)%charm)>0) return

 NNbar:  if (partIn(1)%anti.neqv.partIn(2)%anti) then

            ! Baryon-Antibaryon annihilation:
            !  facts=0.
            !  ! We count the number of possible q qbar pairs that could
            !  ! annihilate, assuming that u and dBar are also allowed to
            !  ! annihilate. Maximal number of pairs is 9
            !  if ((max(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.0)   &
            !       &     .and. &
            !       & (min(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.0)) then
            !     ! S=0 and S=0
            !     facts=1.       ! =3*3 / 9
            !  else if ((max(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.0)   &
            !       &     .and.(min(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-1)) then
            !     ! S=0 and S=-1
            !     facts=2./3.    ! =(3*2)  / 9
            !  else if ((max(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.0)   &
            !       &     .and.(min(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-2)) then
            !     ! S=0 and S=-2
            !     facts=1./3.    ! =(3*1) / 9
            !  else if ( (min(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-1) &
            !       .and. (max(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-1)) then
            !     ! S=-1 and S=-1
            !     facts=5./9.    !=(2*2+1) / 9
            !  else if ( (min(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-2) &
            !       .and. (max(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-1)) then
            !     ! S=-1+S=-2
            !     facts=4./9.    !=(2+2) / 9
            !  else if ( (min(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-2) &
            !       .and. (max(hadron(partIn(1)%ID)%strangeness,hadron(partIn(2)%ID)%strangeness).eq.-2)) then
            !     ! S=-2+S=-2  !=(2*2+1) / 9
            !     facts=5./9.
            !  end if

            call sigmaBarAntiBar(srts, partIn, mediumAtColl, &
                 sigs(isigTot), &
                 sigs(isigElast), &
                 sigCEX=sigs(isigCEX), &
                 sigAnni=sigs(isigAnni), &
                 sigLambdaBar=sigs(isigLbar), &
                 sigSigmaBar=sigs(isigSbar), &
                 sigXiBar=sigs(isigXiBar), &
                 sigOmegaBar=dummy, &
                 sigJPsi=sigs(isigJPsi) )

            if (debug) then
               write(*,*) 'sigmaTot:   ', sigs(isigTot)
               write(*,*) 'sigmaElast: ', sigs(isigElast)
               write(*,*) 'sigmaCEX:   ', sigs(isigCEX)
               write(*,*) 'sigmaAnni:  ', sigs(isigAnni)
               write(*,*) 'sigmaLbar:  ', sigs(isigLbar)
               write(*,*) 'sigmaSbar:  ', sigs(isigSbar)
               write(*,*) 'sigmaXiBar: ', sigs(isigXiBar)
               write(*,*) 'sigmaJPsi:  ', sigs(isigJPsi)
            end if

            !            sigmaTot=facts*sigmaTot
            !            sigmaElast=facts*sigmaElast
            !            sigmaAnni=facts*sigmaAnni
            !            sigmaLbar=facts*sigmaLbar
            !            sigmaSbar=facts*sigmaSbar
            !            sigmaXiBar=facts*sigmaXiBar

         else NNbar
            ! baryon Baryon collisions
            sigs(isigTot)   = sigtot_BB(srts, partIn%ID, sum(partIn(1:2)%charge))
            sigs(isigElast) = nukNuk_nukNuk(srts, (/mN,mN/), sum(partIn(1:2)%charge))
         end if  NNbar
      else  HiEn
         if (debug) write(*,*) ' barBar: low energy', srts

         if ( partIn(1)%anti .and. partIn(2)%anti ) then
            ! At low energies there is no antibaryon-antibaryon
            ! collision-Xsection yet
            partOut%ID=0
            !            sigs(0:1) = 0.
         else if ( partIn(1)%anti .or. partIn(2)%anti ) then
            ! Low energy cross sections for antibaryon-baryon-collisions
            call XsectionAntiBarBar(srts, partIn, mediumAtColl, &
                 partOut, &
                 sigs(isigTot), sigs(isigElast), sigs(isigCEX), sigs(isigAnni),&
                 pauliIncluded, plotIT)

         else
            ! Low energy cross sections for baryon-baryon-collisions
            if (plotIT) then
               call XsectionBarBar(srts, partIn, mediumAtColl, &
                    partOut, &
                    sigs(isigTot), sigs(isigElast), pauliIncluded, &
                    "BaryonBaryon_CrossSection")
            else
               call XsectionBarBar(srts, partIn, mediumAtColl, &
                    partOut, &
                    sigs(isigTot), sigs(isigElast), pauliIncluded)
            end if
         end if
      end if  HiEn

    end subroutine barBar

  end subroutine XsectionMaster


  !****************************************************************************
  !****f* master_2Body/sigtot_BB
  ! NAME
  ! real function sigtot_BB(srts, ID, totalcharge)
  !
  ! PURPOSE
  ! This routine provides parametrizations of the total baryon-baryon
  ! cross sections (possibly consisting of different pieces for different
  ! energy regimes).
  !
  ! INPUTS
  ! * real,    intent(in) :: srts         -- sqrt(s) of the process
  ! * integer, intent(in) :: ID(1:2)      -- IDs of colliding particles
  ! * integer, intent(in) :: totalcharge  -- total charge of colliding particles
  !
  ! References:
  ! * PDG, J. Beringer et al., Phys. Rev. D 86 (2012) 010001
  ! * J. Cugnon, J. Vandermeulen, D. L'Hote,
  !   Nuclear Instrum. Methods B 111 (1996) 215–220.
  !****************************************************************************
  real function sigtot_BB(srts, ID, totalcharge)

    use idTable, only: isHyperon
    use constants, only: mN, pi, GeVSquared_times_mb
    use particleProperties, only: hadron
    use barbar_barbar, only: nukNuk_nukNuk

    real, intent(in) :: srts
    integer, intent(in) :: ID(1:2)
    integer, intent(in) :: totalcharge

    real :: p,x
    real, parameter :: B = pi/GeVSquared_times_mb/2.15**2

    if (isHyperon(ID(1)) .or. isHyperon(ID(2))) then

      ! high-energy parametrization for Sigma-p according to PDG 2012
      x = (hadron(ID(1))%mass + hadron(ID(2))%mass + 2.15)**2 / srts**2
      sigtot_BB = 34.9 + B*log(1./x)**2 - 55.*x**0.462 + 57.*x**0.550

    else

      p = sqrt(srts**4/(4.*mN**2)-srts**2)  ! lab. momentum

      if (totalcharge == 1) then
        ! p+n
        if (p>2.245) then
          ! high-energy parametrization for p+n according to PDG 2012
          x = (hadron(ID(1))%mass + hadron(ID(2))%mass + 2.15)**2 / srts**2
          sigtot_BB = 35.0 + B*log(1./x)**2 + 12.19*x**0.462 - 6.62*x**0.550
        else if (p>0.9) then
          ! Cugnon (this part is rather crude, needs improvement!)
          sigtot_BB = 24.2 + 8.9*p
        else
          ! elastic
          sigtot_BB = nukNuk_nukNuk(srts, (/mN,mN/), totalcharge)
        end if
      else
        ! p+p, n+n
        if (p>4.78) then
          ! high-energy parametrization for p+p according to PDG 2012
          x = (hadron(ID(1))%mass + hadron(ID(2))%mass + 2.15)**2 / srts**2
          sigtot_BB = 34.71 + B*log(1./x)**2 + 12.72*x**0.462 - 7.35*x**0.550
        else if (p>1.5) then
          ! Cugnon
          sigtot_BB = 41. + 60.*(p-0.9)*exp(-1.2*p)
        else if (p>0.7) then
          ! Cugnon
          sigtot_BB = 23.5 + 24.6/(1.+exp((1.2-p)/0.1))
        else
          ! elastic
          sigtot_BB = nukNuk_nukNuk(srts, (/mN,mN/), totalcharge)
        end if
      end if

    end if

  end function sigtot_BB


  !****************************************************************************
  !****f* master_2Body/HiEnergyContrib
  ! NAME
  ! real function HiEnergyContrib(srts, ID)
  !
  ! PURPOSE
  ! This routine returns the relative importance of the HiEnergy part,
  ! i.e. the return value varies between 0 (full low energy resonance model)
  ! and 1 (full HiEnergy model)
  !
  ! This routine is not used by XsectionMaster and related routines.
  ! But it may be helpful for plotting the XS.
  !
  ! INPUTS
  ! * real,    intent(in)  :: srts     -- sqrt(s) of the process
  ! * integer, intent(in)  :: ID(1:2)  -- IDs of colliding particles
  !
  ! OUTPUT
  ! function value as described above
  !
  ! NOTES
  ! This routine assumes a smooth transition.
  ! If you want to have the (old) sharp transition, just compare the function
  ! value against 0.5: <0.5=resonance, >0.5=HiEnergy
  !****************************************************************************
  real function HiEnergyContrib(srts, ID, Anti)
    use IdTable, only: isMeson
    use twoBodyTools, only: LogicMatrix

    real,    intent(in) :: srts
    integer, intent(in) :: ID(1:2)
    logical, intent(in) :: Anti(1:2)

    integer :: scenario
    real :: schwelle = 0.0, delta = 0.0

    ! Initialize switches at first call
    if (initFlag) call ReadInput

    HiEnergyContrib = 0.0

    scenario = LogicMatrix(isMeson(ID(1)), isMeson(ID(2)))

    select case (scenario)
    case (scenarioBarBar)
       if (Anti(1).neqv.Anti(2)) then
         schwelle = HiEnergyThresholdBarAntibar
         delta    = HiEnergyThresholdBarAntibarDelta
       else
         schwelle = HiEnergyThresholdBarBar
         delta    = HiEnergyThresholdBarBarDelta
       end if
    case (scenarioBarMes)
       schwelle = HiEnergyThresholdBarMes
       delta    = HiEnergyThresholdBarMesDelta
    case (scenarioMesMes)
       return
    end select

    if (srts<schwelle-delta) return

    HiEnergyContrib = 1.0
    if (srts>schwelle+delta) return

    if (delta>0) HiEnergyContrib = (srts - (schwelle-delta)) / (2.*delta)

  end function HiEnergyContrib


end module master_2Body
