!******************************************************************************
!****m* /LowPhotonAnalysis
! NAME
! module LowPhotonAnalysis
! PURPOSE
! Includes analysis routines for photon induced runs
!******************************************************************************
module LowPhotonAnalysis
  implicit none
  private

  public :: analyze_Photon, twopi_output, kruscheAnalyse


  !****************************************************************************
  !****g* lowPhotonAnalysis/outputEvents
  ! SOURCE
  !
  logical,save  :: outputEvents=.false.
  ! PURPOSE
  ! If .true. then all events are printed to file.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/outputEvents_onlyFree
  ! SOURCE
  !
  logical,save  :: outputEvents_onlyFree=.false.
  ! PURPOSE
  ! If outputEvents=.true. then only particles which may leave the nucleus,
  ! i.e. may become "free", are printed to file.
  !****************************************************************************


  !****************************************************************************
  !****g* lowPhotonAnalysis/FissumOutput
  ! SOURCE
  !
  logical,save  ::  FissumOutput=.false.
  ! PURPOSE
  ! If .true. then we perform an analysis as in PRC 53,#3 pages 1278 ff. (1996)
  ! Produces dsigma/dOmega/dT_pi for pi^+
  !****************************************************************************


  !****************************************************************************
  !****g* lowPhotonAnalysis/KruscheOutput
  ! SOURCE
  !
  logical,save  ::  KruscheOutput=.false.
  ! PURPOSE
  ! If .true. then we perform an analysis as in EPJA22 347-351 (2004)
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/KruscheAnalyse_cut
  ! SOURCE
  !
  real, save :: KruscheAnalyse_cut=0. ! GeV
  ! PURPOSE
  ! Value of the cut for the deltaE cut in EPJA22 347-351 (2004).
  !****************************************************************************


  !****************************************************************************
  !****g* lowPhotonAnalysis/TwoPiOutput
  ! SOURCE
  !
  logical,save :: TwoPiOutput = .false.
  ! PURPOSE
  ! If .true. then we perform an analysis for 2pi production,
  ! including statistics for the mass of the pi-pi pair.
  !****************************************************************************


  !****************************************************************************
  !****g* lowPhotonAnalysis/photonAnalyse
  ! SOURCE
  !
  logical,save :: photonAnalyse=.false.
  ! PURPOSE
  ! Special analysis for final state photons
  !****************************************************************************


  !****************************************************************************
  !****g* lowPhotonAnalysis/pi0gamma_analysis
  ! SOURCE
  !
  logical, save :: pi0gamma_analysis = .false.
  ! PURPOSE
  ! Do analysis of pi0 gamma pairs (dsigma/dm), to reconstruct invariant mass
  ! spectrum of omega mesons.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/pi0gamma_momcut
  ! SOURCE
  !
  real, save :: pi0gamma_momcut = 0.5
  ! PURPOSE
  ! Cut on the absolute omega three momentum in GeV,
  ! being applied to the pi0 gamma spectrum.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/pi0gamma_masscut
  ! SOURCE
  !
  real, dimension(1:2), save :: pi0gamma_masscut = (/0.,2./)
  ! PURPOSE
  ! Cuts on the pi0-gamma invariant mass in GeV, being applied to all
  ! pi0-gamma spectra (except the mass spectrum).
  ! First component is lower limit, second component is upper limit.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/pi0gamma_mombin
  ! SOURCE
  !
  real, save :: pi0gamma_mombin = 0.050
  ! PURPOSE
  ! Bin size for the pi0 gamma momenentum spectrum in GeV.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/pi0gamma_massres_sigma
  ! SOURCE
  !
  real, save :: pi0gamma_massres_sigma = 0.025
  ! PURPOSE
  ! Sigma parameter for the exp. resolution smearing (width of the Gauss or
  ! Novosibirsk function in GeV).
  ! See also pi0gamma_massres_tau.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/pi0gamma_massres_tau
  ! SOURCE
  !
  real, save :: pi0gamma_massres_tau = -0.090
  ! PURPOSE
  ! Skewness parameter tau of the Novosibirsk function (for exp. resolution
  ! smearing).
  ! See also pi0gamma_massres_sigma.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/Ekin_pi0_cut
  ! SOURCE
  !
  real, save :: Ekin_pi0_cut = 0.
  ! PURPOSE
  ! Cut on the kinetic energy of neutral pions in the pi0gamma_analysis.
  ! Only pions with kinetic energies larger than this cutoff are used for the
  ! analysis.
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/eta_analysis
  ! SOURCE
  !
  logical, save :: eta_analysis = .false.
  ! PURPOSE
  ! Do analysis of eta mesons
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/MissingMass_analysis
  ! SOURCE
  !
  logical, save :: MissingMass_analysis = .false.
  ! PURPOSE
  ! Do analysis of missing mass distributions for 1 nucleon and 2 nucleon
  ! final states
  !****************************************************************************

  !****************************************************************************
  !****g* lowPhotonAnalysis/DoOutChannels
  ! SOURCE
  !
  logical, save :: DoOutChannels = .false.
  ! PURPOSE
  ! switch on/off: reporting of all final state channels
  !****************************************************************************


  logical, save :: initFlag=.true.

contains

  !****************************************************************************
  !****s* lowPhotonAnalysis/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads out the jobcard "photonAnalysis".
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput
    use output, only: Write_ReadingInput

    integer :: IOS
    !**************************************************************************
    !****n* lowPhotonAnalysis/lowPhotonAnalysis
    ! NAME
    ! NAMELIST /lowPhotonAnalysis/
    ! PURPOSE
    ! Includes the switches:
    ! * outputEvents
    ! * outputEvents_onlyFree
    ! * KruscheOutput
    ! * KruscheAnalyse_cut
    ! * FissumOutput
    ! * photonAnalyse
    ! * TwoPiOutput
    ! * pi0gamma_analysis
    ! * pi0gamma_momcut
    ! * pi0gamma_masscut
    ! * pi0gamma_mombin
    ! * pi0gamma_massres_sigma
    ! * pi0gamma_massres_tau
    ! * Ekin_pi0_cut
    ! * eta_analysis
    ! * MissingMass_analysis
    ! * DoOutChannels
    !**************************************************************************
    NAMELIST /lowPhotonAnalysis/ outputEvents, outputEvents_onlyFree, &
         KruscheOutput, kruscheAnalyse_cut, &
         fissumOutput, photonAnalyse, twoPiOutput, &
         pi0gamma_analysis, pi0gamma_momcut, pi0gamma_masscut, &
         pi0gamma_mombin, &
         pi0gamma_massres_sigma, pi0gamma_massres_tau, Ekin_pi0_cut, &
         eta_analysis, MissingMass_analysis, DoOutChannels

    call Write_ReadingInput('lowPhotonAnalysis',0)
    rewind(5)
    read(5,nml=lowPhotonAnalysis,IOSTAT=IOS)
    if (outputEvents) then
       if (outputEvents_onlyFree) then
          write(*,*) 'All events are printed to files!                   '
       else
          write(*,*) 'All events are printed to files (but only free particles)'
       end if
    else
       write(*,*) 'Events are NOT printed to files!                   '
    end if
    write(*,*) 'Perform output for Fissum  pi^+ production ?', fissumOutput
    write(*,*) 'Perform output for Krusche pi^0 production ?', kruscheOutput
    if (kruscheOutput) &
         write(*,*) '  =>Cut for Krusche analysis=', kruscheAnalyse_cut
    write(*,*) 'Perform output for 2pi production          ?', TwoPiOutput
    write(*,*) 'Perform analysis for final state photons   ?', photonAnalyse
    write(*,*) 'Perform pi0 gamma analysis                 ?', pi0gamma_analysis
    if (pi0gamma_analysis) then
       write(*,'(A,1F7.4)') ' pi0 gamma momentum cut [GeV]               ?', &
            pi0gamma_momcut
       write(*,'(A,2F7.4)') ' pi0 gamma mass cut [GeV]                   ?', &
            pi0gamma_masscut
       write(*,'(A,1F7.4)') ' pi0 gamma momentum bin size [GeV]          ?', &
            pi0gamma_mombin
       write(*,'(A,1F7.4)') ' pi0 gamma mass resolution: sigma [GeV]     ?', &
            pi0gamma_massres_sigma
       write(*,'(A,1F7.4)') ' pi0 gamma mass resolution: tau             ?', &
            pi0gamma_massres_tau
       write(*,'(A,1F7.4)') ' pi0 kinetic energy cut                     ?', &
            Ekin_pi0_cut
    end if
    write(*,*) 'Perform eta analysis                       ?', eta_analysis
    write(*,*) 'Perform missing mass  analysis             ?', MissingMass_analysis

    write(*,*) 'DoOutChannels    : ',DoOutChannels

    call Write_ReadingInput('lowPhotonAnalysis',0,IOS)

  end subroutine readinput



  !****************************************************************************
  !****s* LowPhotonAnalysis/analyze_Photon
  ! NAME
  !  subroutine analyze_Photon(Parts,finalFlag)
  ! INPUTS
  !  * type(particle), dimension(:,:) :: Parts -- particle vector
  !  * logical, :: finalFlag -- .true. if it is the last call for one
  !    specific energy, therefore final output must be made.
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! Main routine for the analysis. Calls all subroutines.
  !****************************************************************************
  subroutine analyze_Photon(Parts,finalFlag)

    use lowElectronAnalysis, only: lowElectron_Analyze
    use particleDefinition
    use initLowPhoton, only: getTwoPi,getResonances,getSinglePi, getPi0Eta

    type(particle), dimension(:,:), intent(in)  :: Parts
    logical, intent (in)                        :: finalFlag

    if (initFlag) then
       initFlag=.false.
       call readInput()
    end if

    if (getTwoPi().and.TwoPiOutput) then
       call twoPi_outPut(Parts,finalFlag)
       call twoPi_outPut_hist(Parts,finalFlag)
    end if

    if (getResonances().or.getSinglePi().or.getTwoPi().or.getpi0Eta()) &
         call lowElectron_Analyze(Parts,finalFlag)
    if ((getResonances().or.getSinglePi().or.getTwoPi()).and.KruscheOutput) &
         call kruscheAnalyse(Parts,finalFlag)
    if ((getResonances().or.getSinglePi().or.getTwoPi()).and.FissumOutput) &
         call fissumAnalyse(Parts,finalFlag)

    ! Metag Omega Experiment
    if (photonAnalyse .or. outputEvents .or. DoOutChannels .or. MissingMass_analysis) &
         call photonFS_Analyse(Parts,finalFlag)

    if (pi0gamma_analysis) &
         call pi0gamma_analyze(Parts)

    if (eta_analysis) &
         call etaAnalyse(Parts,finalFlag)


  end subroutine analyze_Photon



  !****************************************************************************
  !****s* LowPhotonAnalysis/photonFS_Analyse
  ! NAME
  ! subroutine photonFS_Analyse(Parts)
  ! INPUTS
  ! * type(particle), target, dimension(:,:) :: Parts --- particle vector
  ! * logical, finalFlag -- .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  ! OUTPUT
  ! if "outputEvents" is true:
  !  * Prints all events to a file called 'FinalEvents.dat'
  !  * Events stemming from omega, rho or phi production are printed to
  !    seperate files 'FinalEvents.rho.dat', ...
  ! if "photonAnalyse" is true:
  !  * Produces files called Photons_m2_vs_m3*.dat and Photons_dsigma_dm*.dat.
  !  * The first one contains the distribution of the invariant mass of all
  !    3 gamma combinations vs the inv. mass of all 2gamma combinations.
  !  * The second one contains the cross section dSigma/dm vs the mass of
  !    3 gammas.
  ! PURPOSE
  ! Groups all final state particles into events and writes them to file
  !****************************************************************************
  subroutine photonFS_Analyse(Parts,finalFlag)
    use inputGeneral, only: current_run_number
    use particleDefinition
    use initLowElectron, only: le_get_FirstEventRange
    use Electron_origin, only: le_isOmegaEvent,le_ispi0etaEvent, &
         & le_isRhoEvent, &
         & le_isPiPiEvent,le_isResEvent,le_isPhiEvent
    use AnaEventDefinition
    use AnaEvent, only: event_init, event_dump, event_add, event_clear, &
         event_pairPhotons
    use hist2D
    use hist
    use PreEvListDefinition
    use PreEvList, only: CreateSortedPreEvent, PreEvList_CLEAR, PreEvList_INIT,&
         PreEvList_INSERT, PreEvList_Print
    use callStack, only: traceback
    use output, only: inttochar4
    use initLowPhoton, only: getEnergyGamma

    type(particle), target, dimension(:,:) :: Parts
    logical, intent(in) :: finalFlag


    integer, dimension (1:2) :: firstEvents
    type(particle), POINTER :: pPart
    type(histogram2D), dimension(0:6),save :: Hist_massPhotons
    type(histogram), dimension(0:6,0:8),save :: dsigma_dm_threeGammas

    type(tAnaEvent), allocatable, dimension(:) :: events ! A list of all events
    type(tAnaEvent), allocatable, dimension(:) :: events_rho,events_phi,events_omega ! A list of all events coming from rho, phi, omega
    type(tAnaEvent), allocatable, dimension(:) :: events_pipi,events_pi0eta,events_res! A list of all events coming from pi pi, pi0 eta,resonance production

    type(tPreEvListEntry), save :: preEv
    type(tPreEvList), save :: ListPreEv

    integer        :: i,j,first
    logical, save  :: initFlag=.true.
    integer, save  :: nCall = 0, nCall0 = 0
    real :: mul0, Egamma


    ! (0) Global set up
    !
    if (initFlag) then
       call PreEvList_INIT(ListPreEv)
       nCall = 0
       nCall0 = nCall0+1
       initFlag = .false.
    end if

    nCall = nCall+1

    ! (1) Setting up the particles into the events
    ! This is done with the help of %firstIndex:
    ! Particles stemming from the same event get in the init the same
    ! %firstIndex entry. During the run %firstinit stays constant and is
    ! inherited during collisions.
    !
    firstEvents=le_get_FirstEventRange()
    write(*,*) 'firstEvents',firstEvents

    allocate(Events(firstEvents(1):firstEvents(2)))
    allocate(Events_rho(firstEvents(1):firstEvents(2)))
    allocate(Events_phi(firstEvents(1):firstEvents(2)))
    allocate(Events_omega(firstEvents(1):firstEvents(2)))
    allocate(Events_pipi(firstEvents(1):firstEvents(2)))
    allocate(Events_res(firstEvents(1):firstEvents(2)))
    allocate(Events_pi0eta(firstEvents(1):firstEvents(2)))

    do i=firstEvents(1),firstEvents(2)
       call event_init(events(i))
       call event_init(events_rho(i))
       call event_init(events_phi(i))
       call event_init(events_omega(i))
       call event_init(events_pipi(i))
       call event_init(events_pi0eta(i))
       call event_init(events_res(i))
    end do

    do i=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do j=lbound(Parts,dim=2),ubound(Parts,dim=2)
          if (Parts(i,j)%ID.le.0) cycle
          pPart => Parts(i,j)
          first = pPart%firstEvent

          ! Add particle to the event with its firstEvent index.
          call event_add(events(first),pPart)
          if (le_isOmegaEvent(first)) &
               call event_add(events_omega(first), pPart)
          if (le_isRhoEvent(first)) &
               call event_add(events_rho(first), pPart)
          if (le_isPhiEvent(first)) &
               call event_add(events_phi(first), pPart)
          if (le_isPiPiEvent(first)) &
               call event_add(events_pipi(first), pPart)
          if (le_isResEvent(first)) &
               call event_add(events_res(first), pPart)
          if (le_ispi0etaEvent(first)) &
               call event_add(events_pi0eta(first), pPart)
       end do
    end do


    ! (2) Do Photon Analysis:
    !     Plot m(3gamma) vs. m(2gamma) and dSigma/dm(3gamma)

    if (photonAnalyse) then
       call event_pairPhotons(events,        '', &
            dsigma_dm_threeGammas(0,:),Hist_massPhotons(0))
       call event_pairPhotons(events_rho,    '.rho', &
            dsigma_dm_threeGammas(1,:),Hist_massPhotons(1))
       call event_pairPhotons(events_omega,  '.omega', &
            dsigma_dm_threeGammas(2,:),Hist_massPhotons(2))
       call event_pairPhotons(events_phi,    '.phi', &
            dsigma_dm_threeGammas(3,:),Hist_massPhotons(3))
       call event_pairPhotons(events_pipi,   '.pipi', &
            dsigma_dm_threeGammas(4,:),Hist_massPhotons(4))
       call event_pairPhotons(events_res,    '.res', &
            dsigma_dm_threeGammas(5,:),Hist_massPhotons(5))
       call event_pairPhotons(events_pi0eta, '.pi0eta', &
            dsigma_dm_threeGammas(6,:),Hist_massPhotons(6))
    end if


    ! (3) Dump all events to file:
    !     print ID, charge, perweight, position and momentum of all produced
    !     particles

    if (outputEvents) then
       write(*,*) 'Writing events to file'

       call event_dump(current_run_number,events,&
            'FinalEvents.dat',outputEvents_onlyFree)
       call event_dump(current_run_number,events_rho,&
            'FinalEvents.rho.dat',outputEvents_onlyFree)
       call event_dump(current_run_number,events_omega,&
            'FinalEvents.omega.dat',outputEvents_onlyFree)
       call event_dump(current_run_number,events_phi,&
            'FinalEvents.phi.dat',outputEvents_onlyFree)
       call event_dump(current_run_number,events_pipi,&
            'FinalEvents.pipi.dat',outputEvents_onlyFree)
       call event_dump(current_run_number,events_res,&
            'FinalEvents.res.dat',outputEvents_onlyFree)
       call event_dump(current_run_number,events_pi0eta,&
            'FinalEvents.pi0eta.dat',outputEvents_onlyFree)

    end if

    ! (4) Group events by their final state

    if (DoOutChannels) then

       do i=firstEvents(1),firstEvents(2)

          if (CreateSortedPreEvent(events(i),PreEv%preE)) then
             PreEv%weight = events(i)%particleList%first%V%perweight
             call PreEvList_INSERT(ListPreEv,PreEv)
          end if

       end do

    end if

    ! (5) do missing mass analysis for 1- and 2-proton final states

    if (MissingMass_analysis) then

       call doMissingMass()

    end if

    ! (*) Clean Up:
    !     Clear event lists and deallocate memory

    do i=firstEvents(1),firstEvents(2)
       call event_clear(events(i))
       call event_clear(events_rho(i))
       call event_clear(events_phi(i))
       call event_clear(events_omega(i))
       call event_clear(events_pipi(i))
       call event_clear(events_pi0eta(i))
       call event_clear(events_res(i))
    end do

    deallocate(Events)
    deallocate(Events_rho)
    deallocate(Events_phi)
    deallocate(Events_omega)
    deallocate(Events_pipi)
    deallocate(Events_res)
    deallocate(Events_pi0eta)


    ! (*) Global clean up (prepare for next energy)
    if (finalFlag) then

       ! now we are doing finally some of the output...

       if (DoOutChannels) then
          write(*,*) 'nCall=',nCall,nCall0
          mul0 = 1.0/nCall

          Egamma = getEnergyGamma()

          open(141,file='OutChannels.FINAL.'//inttochar4(nCall0)//'.dat', status='unknown')
          rewind(141)
          Egamma = getEnergyGamma()
          write(141,"(A,f8.4)") "# E_gamma =",Egamma
          call PreEvList_Print(141,ListPreEv,mul0)
          close(141)

          call PreEvList_CLEAR(ListPreEv)

!          call traceback("mutliple energies not yet implemented")
       end if

       initFlag = .true.
    end if


  contains

    subroutine doMissingMass

      use AnaEvent, only: event_getParticle
      use hist_multipleRuns_mc
      use output, only: realToChar4

      integer :: i
      type(particle) :: p1,p2
      logical :: flagOK
      real, dimension(0:3) :: mom
      real :: missMass2
      integer :: nPi
      logical, save :: initFlag = .true.

      type(histogram_mr_mc),save :: dN1dMM2

      if (initFlag) then
         call CreateHist_mr_mc(dN1dMM2, 'dN/dM^2', -2.0,1.0,0.005,4)

         dN1dMM2%total%ydesc(1:4) = (/ ' 0pi', ' 1pi', ' 2pi', '>2pi' /)
         dN1dMM2%total%xdesc = 'miss mass^2 [GeV^2]'
         initFlag = .false.
      end if

      call startRunHist_mr_mc(dN1dMM2)

      Egamma = getEnergyGamma()

      do i=firstEvents(1),firstEvents(2)

         nPi = sum(events(i)%numberParticles(1,:))
         nPi = min(nPi,3)


         select case (sum(events(i)%numberParticles(9,:)))
         case (1)

            flagOK = event_getParticle(events(i), 1,1, 1, p1, .false.)
            if (.not.flagOK) &
                 flagOK = event_getParticle(events(i), 1,0, 1, p1, .false.)
            if (.not.flagOK) &
                 call traceback("could not find the nucleon")

            mom = (/ Egamma + 0.938, 0., 0., Egamma /) ! mass of free nucleon
            mom = mom - p1%mom

            missMass2 = mom(0)**2-sum(mom(1:3)**2)

            write(*,*) 'missMass [1] = ',missMass2,nPi
            call AddHist_mr_mc(dN1dMM2, missMass2, nPi+1, p1%perweight)

         case (2)

            flagOK = event_getParticle(events(i), 1,1, 1, p1, .false.)
            if (.not.flagOK) then
               flagOK = event_getParticle(events(i), 1,0, 1, p1, .false.)
               if (.not.flagOK) &
                    call traceback("could not find the first nucleon")
               flagOK = event_getParticle(events(i), 1,0, 2, p2, .false.)
               if (.not.flagOK) &
                    call traceback("could not find the second nucleon")
            else
               flagOK = event_getParticle(events(i), 1,1, 2, p2, .false.)
               if (.not.flagOK) then
                  flagOK = event_getParticle(events(i), 1,0, 1, p2, .false.)
                  if (.not.flagOK) &
                       call traceback("could not find the second nucleon")
               end if
            end if

            mom = (/ Egamma + 11.17, 0., 0., Egamma /) ! mass of C12
            mom = mom - p1%mom - p2%mom

            missMass2 = mom(0)**2-sum(mom(1:3)**2)

            write(*,*) 'missMass [2] = ',missMass2,nPi

         case default
         end select

      end do

      call endRunHist_mr_mc(dN1dMM2)

      call WriteHist_mr_mc(dN1dMM2,978,&
           & file='missmass1_'//realToChar4(ANINT(Egamma*1000.))//'_MeV.dat')

      if (finalFlag) then
         call removeHist_mr_mc(dN1dMM2)
         initFlag = .true.
      end if

    end subroutine doMissingMass

  end subroutine photonFS_Analyse


  !****************************************************************************
  !****s* LowPhotonAnalysis/pi0gamma_analyze
  ! NAME
  ! subroutine pi0gamma_analyze(Parts)
  ! INPUTS
  ! * type(particle), target, dimension(:,:) :: Parts -- particle vector
  ! OUTPUT
  ! * Produces an invariant mass spectrum of all pi0 gamma pairs in the final
  !   state and writes it to a file called 'pi0gamma_dsigma_dm.dat'
  !****************************************************************************
  subroutine pi0gamma_analyze(Parts)
    use particleDefinition
    use particlePointerListDefinition
    use AnaEventDefinition
    use IDTable, only: pion, photon
    use minkowski, only: abs4,abs3
    use anaEvent, only: event_init, event_add, event_clear
    use initLowElectron, only: le_get_FirstEventRange
    use hist, only: histogram, CreateHist, AddHist, WriteHist, &
         WriteHist_Gauss, WriteHist_Novo
    use inputGeneral, only: num_Runs_sameEnergy, num_Energies
    use initLowPhoton, only: getEnergyGamma, getDeltaEnergy, energyWeight
    use constants, only: rhoNull
    use PIL_omegaDec, only: PIL_omegaDec_GET

    type(particle), dimension(:,:), target :: Parts
    type(particle), pointer :: pPart
    type(tParticleListNode), pointer :: p_pi, p_gam
    type(tAnaEvent), Allocatable, dimension(:) :: pions, photons
    type(histogram),save :: dsigma_dm, dsigma_dm_momcut, dsigma_dp, dsigma_dE, &
         dsigma_dEkin, omega_excit_func, decDensity
    integer, dimension (1:2) :: firstEvents
    integer :: i,j,first
    real :: m,E_gamma,dE,pabs,E_tot,E_kin,dens,w
    integer, save :: num_pi=0,num_gam=0
    character(len=100) :: title

    firstEvents=le_get_FirstEventRange()
    write(*,*) 'pi0gamma_analyze: firstEvents',firstEvents

    if (.not. dsigma_dm%initialized) then
       call CreateHist(dsigma_dm, &
            'dsigma/dm (pi0 gamma) in microbarn/GeV/A', 0., 1., 0.001)
       write(title,'(A,F6.3,A)') &
            'dsigma/dm (pi0 gamma) in microbarn/GeV/A, |p| < ',&
            pi0gamma_momcut,' GeV'
       call CreateHist(dsigma_dm_momcut, title, 0., 1., 0.001)
       call CreateHist(dsigma_dp, &
            'dsigma/dp (pi0 gamma) in microbarn/GeV/A', 0., 2., pi0gamma_mombin)
       call CreateHist(dsigma_dE, &
            'dsigma/dE (pi0 gamma) in microbarn/GeV/A', 0., 2., pi0gamma_mombin)
       call CreateHist(dsigma_dEkin, &
            'dsigma/dE_kin (pi0 gamma) in microbarn/GeV/A', &
            -0.2, 2., pi0gamma_mombin)
       call CreateHist(decDensity, &
            'density at decay point [rho/rho0]', 0., 1., 0.01)
    end if

    E_gamma = getEnergyGamma()

    if (.not. omega_excit_func%initialized) then
       dE = getDeltaEnergy()
       call CreateHist(omega_excit_func, &
            'omega excitation function (including final state interaction) in microbarn/A', &
            E_gamma-dE/2, E_gamma+(num_Energies-0.5)*dE, dE)
    end if

    allocate(pions(firstEvents(1):firstEvents(2)))
    allocate(photons(firstEvents(1):firstEvents(2)))

    do i=firstEvents(1),firstEvents(2)
       call event_init(pions(i))
       call event_init(photons(i))
    end do

    ! loop over all particles and construct events
    do i=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do j=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(i,j)

!!! write out mass distribution to reproduce Muehlich fig. 9.14.:
!!!if (pPart%ID==omegaMeson) write(667,'(i12,3G12.7)') pPart%number, abs4(pPart%mom), abs3(pPart%mom), pPart%mom(0)

          ! keep only pions and photons:
          if (pPart%ID/=pion .and. pPart%ID/=photon) cycle

          first=pPart%firstEvent
          ! Add particle to the event with its firstEvent index.
          if (pPart%ID == pion .and. pPart%charge == 0 .and. &
               kineticEnergy(pPart)>Ekin_pi0_cut) then
             call event_add(pions(first), pPart)
             num_pi=num_pi+1
          else if (pPart%ID == photon) then
             call event_add(photons(first), pPart)
             num_gam=num_gam+1
          end if
       end do
    end do

    write(*,*) 'pi0gamma_analyze: #pions = ',num_pi
    write(*,*) 'pi0gamma_analyze: #gammas = ',num_gam

    ! loop over all pi0 gamma pairs, make histogram of inv. mass
    do i=firstEvents(1),firstEvents(2)
       ! loop over all pions
       p_pi => pions(i)%particleList%first
       do
          if (.not. associated(p_pi)) exit
          ! loop over all gammas
          p_gam => photons(i)%particleList%first
          do
             if (.not. associated(p_gam)) exit
             ! invariant mass of pi0-gamma pair:
             m     = abs4(p_pi%V%mom + p_gam%V%mom)
             ! absolute 3-momentum of pi0-gamma pair:
             pabs  = abs3(p_pi%V%mom + p_gam%V%mom)
             ! total energy of pi0-gamma pair:
             E_tot = p_pi%V%mom(0) + p_gam%V%mom(0)
             ! kinetic energy of pi0-gamma pair:
             E_kin = E_tot - m
             w = p_pi%V%perweight/float(num_Runs_sameEnergy*num_Energies)
             call addHist(dsigma_dm, m, w)
             if (m>pi0gamma_masscut(1) .and. m<pi0gamma_masscut(2)) then
                call addHist(dsigma_dp, pabs, w)
                call addHist(dsigma_dE, E_tot, w)
                call addHist(dsigma_dEkin, E_kin, w)
                ! divide out energy weight to get a proper excitation function
                call addHist (omega_excit_func, E_gamma, &
                     w*num_Energies*energyWeight(E_gamma))
                if (PIL_omegaDec_Get(p_pi%V%firstEvent, dens)) & ! decay density
                     call addHist(decDensity, dens/rhoNull, y=w)
                if (pabs<pi0gamma_momcut) then ! momentum cut
                   call addHist(dsigma_dm_momcut, m, w)
                end if
             end if
             p_gam => p_gam%next
          end do
          p_pi => p_pi%next
       end do
    end do

    !**************************************************************************
    !****o* LowPhotonAnalysis/pi0gamma_dsigma_dm.dat
    ! NAME
    ! file pi0gamma_dsigma_dm.dat
    ! PURPOSE
    ! Pi0-gamma invariant-mass spectrum in microbarn/GeV/A.
    !**************************************************************************
    call WriteHist(dsigma_dm, file="pi0gamma_dsigma_dm.dat")
    call WriteHist_Gauss(dsigma_dm, "pi0gamma_dsigma_dm_gauss.dat", &
         pi0gamma_massres_sigma)
    call WriteHist_Novo(dsigma_dm, "pi0gamma_dsigma_dm_novo.dat", &
         pi0gamma_massres_sigma, pi0gamma_massres_tau)

    call WriteHist(dsigma_dm_momcut, file="pi0gamma_dsigma_dm_momcut.dat")
    call WriteHist_Gauss(dsigma_dm_momcut, &
         "pi0gamma_dsigma_dm_momcut_gauss.dat", pi0gamma_massres_sigma)
    call WriteHist_Novo(dsigma_dm_momcut, &
         "pi0gamma_dsigma_dm_momcut_novo.dat",  pi0gamma_massres_sigma, &
         pi0gamma_massres_tau)
    !**************************************************************************
    !****o* LowPhotonAnalysis/pi0gamma_dsigma_dp.dat
    ! NAME
    ! file pi0gamma_dsigma_dp.dat
    ! PURPOSE
    ! Pi0-gamma momentum spectrum in microbarn/GeV/A.
    !**************************************************************************
    call WriteHist(dsigma_dp, file="pi0gamma_dsigma_dp.dat")
    !**************************************************************************
    !****o* LowPhotonAnalysis/pi0gamma_dsigma_dE.dat
    ! NAME
    ! file pi0gamma_dsigma_dE.dat
    ! PURPOSE
    ! Pi0-gamma energy spectrum in microbarn/GeV/A.
    !**************************************************************************
    call WriteHist(dsigma_dE, file="pi0gamma_dsigma_dE.dat")
    !**************************************************************************
    !****o* LowPhotonAnalysis/pi0gamma_dsigma_dEkin.dat
    ! NAME
    ! file pi0gamma_dsigma_dEkin.dat
    ! PURPOSE
    ! Pi0-gamma kinetic energy spectrum in microbarn/GeV/A.
    !**************************************************************************
    call WriteHist(dsigma_dEkin, file="pi0gamma_dsigma_dEkin.dat")
    !**************************************************************************
    !****o* LowPhotonAnalysis/OmegaExcitFunc_pi0gamma.dat
    ! NAME
    ! file OmegaExcitFunc_pi0gamma.dat
    ! PURPOSE
    ! Contains the omega excitation function, i.e. the energy-dependent
    ! inclusive omega photo-production cross section, as observed through
    ! the decay into pi0-gamma. In contrast to OmegaExcitFunc.dat, this does
    ! include FSI effects.
    !**************************************************************************
    ! multiplicative factor prevents division by bin size
    call WriteHist(omega_excit_func, mul=omega_excit_func%xBin, &
         file="OmegaExcitFunc_pi0gamma.dat")

    call WriteHist(decDensity, file='pi0gammaDecayDensity.dat')

    do i=firstEvents(1),firstEvents(2)
       call event_clear(pions(i))
       call event_clear(photons(i))
    end do

    deallocate(pions)
    deallocate(photons)

  end subroutine pi0gamma_analyze


  !****************************************************************************
  !****s* LowPhotonAnalysis/kruscheAnalyse
  ! NAME
  ! subroutine kruscheAnalyse(Parts,finalFlag)
  ! PURPOSE
  ! Evaluates the cross section for pi^0 production in the same way as Krusche
  ! in EPJA 22, 347 (2004)
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts -- particle vector
  ! * logical, finalFlag -- .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  ! OUTPUT
  ! * See 'KruscheAnalyse.dat' for Xsections and delta(E) histograms are given
  !   in 'deltaE_E_gamma_'//realToChar4(energy_gamma)//'_GeV.dat'
  !****************************************************************************
  subroutine kruscheAnalyse(Parts,finalFlag)
    use particleDefinition
    use output, only: realToChar4
    use initLowPhoton, only: energy_gamma
    use hist
    use IDTable, only: pion

    type(particle), dimension(:,:), intent(in), target :: Parts
    logical, intent(in) :: finalFlag


    type(particle), pointer :: pPart
    integer, save         :: numberOfCalls=0 ! number of calls per energy
    logical, save         :: firstCall=.true.
    logical, save         :: initFlag=.true.
    real,    save         :: sigma, sigma_cut, error, error_cut
    type(histogram),save  :: hist_deltaE, hist_dp

    real           :: new, new_cut
    integer        :: i,j
    real           :: deltaE
    integer, dimension (1:2) :: count

    count=0

    if (firstCall) then

       initFlag=.true.

       !***********************************************************************
       !****o* kruscheAnalyse/KruscheAnalyse.dat
       ! NAME
       ! file KruscheAnalyse.dat
       ! PURPOSE
       ! Evaluate the cross section for pi^0 production in the same way
       ! as Krusche in EPJA 22, 347 (2004).
       ! PURPOSE
       ! Columns:
       ! * 1 : E_gamma [Gev]
       ! * 2 : sigma pi0 [mb]
       ! * 3 : sigma pi0 after cut [mb]
       ! * 4 : error sigma pi0 [mb]
       ! * 5 : error sigma pi0 after cut [mb]
       !***********************************************************************
       open(11,file='KruscheAnalyse.dat')
       !***********************************************************************
       !****o* kruscheAnalyse/KruscheAnalyse_prelim.dat
       ! NAME
       ! file KruscheAnalyse_prelim.dat
       ! PURPOSE
       ! Same as KruscheAnalyse.dat, but output is written afer each run and
       ! not only the averaged result after each energy.
       !***********************************************************************
       open(12,file='KruscheAnalyse_prelim.dat')
       write(11,*) '# E_gamma, sigma pi0, sigma pi0 after cut, error sigma pi0, error sigma pi0 after cut'
       write(12,*) '# E_gamma, sigma pi0, sigma pi0 after cut, error sigma pi0, error sigma pi0 after cut'
       firstCall=.false.

    else

       open(11,file='KruscheAnalyse.dat',position='append')
       open(12,file='KruscheAnalyse_prelim.dat',position='append')

    end if

    if (initFlag) then
       sigma=0.
       sigma_cut=0.
       error=0.
       error_cut=0.
       numberOfCalls=0
       call CreateHist(hist_deltaE, 'delta E according to Krusche', &
            -0.4,0.4,0.005)
       call CreateHist(hist_dp, 'Pion momenta in lab',0.,0.4,0.005)
       initFlag=.false.
    end if

    numberOfCalls=numberOfCalls+1
    new=0.
    new_cut=0.

    do i=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do j=lbound(Parts,dim=2),ubound(Parts,dim=2)
          ! Select all pi^0
          if (Parts(i,j)%ID.ne.pion) cycle
          if (Parts(i,j)%Charge.ne.0) cycle

          pPart => Parts(i,j)

          new = new + pPart%perweight

          ! write momentum to histogram
          call AddHist(hist_dp,absMom(pPart),pPart%perweight,1.)

          deltaE=get_deltaE(pPart)
          if (deltaE.gt.0) then
             count(1)=count(1)+1
          else
             count(2)=count(2)+1
          end if
          ! Write deltaE to histogram
          call AddHist(hist_deltaE,deltaE,pPart%perweight,1.)
          ! Use deltaE to cut on the data
          if (deltaE.ge.kruscheAnalyse_cut) then
             new_cut=new_cut+pPart%perweight
          end if
       end do
    end do

    error    =error    +new**2
    error_cut=error_cut+new_cut**2

    sigma    =sigma    +new
    sigma_cut=sigma_cut+new_cut

    if (finalFlag) then
       ! Error analysis
       if (numberOfCalls.gt.1) then
          error_cut=sqrt(abs((error_cut-1./float(numberOfCalls)*sigma_cut**2) &
               &  /float(numberOfCalls-1))/float(numberOfCalls))
          error   =sqrt(abs((error     -1./float(numberOfCalls)*sigma**2    ) &
               &  /float(numberOfCalls-1))/float(numberOfCalls))
       else
          error    =sigma    /float(numberOfCalls)
          error_cut=sigma_cut/float(numberOfCalls)
       end if

       ! OUTPUT
       write(11,'(5E15.4)') energy_gamma, sigma/numberOfCalls, &
            & sigma_cut/numberOfCalls, error, error_cut
       write(12,'(5E15.4)') energy_gamma, sigma/numberOfCalls, &
            & sigma_cut/numberOfCalls, error, error_cut

       initFlag=.true.
       call WriteHist(hist_deltaE,978,mul=1./float(numberOfCalls),&
            & file='deltaE_E_gamma_'//realToChar4(energy_gamma*1000.)//'_MeV.dat')
       call WriteHist(hist_dp,979,mul=1./float(numberOfCalls),&
            & file='dp_piZero_'//realToChar4(energy_gamma*1000.)//'_MeV.dat')
       call removeHist(hist_deltaE)
       call removeHist(hist_dp)
    else
       call WriteHist(hist_deltaE,978,mul=1./float(numberOfCalls),&
            & file='deltaE_E_gamma_'//realToChar4(energy_gamma*1000.)//'_MeV.dat')
       call WriteHist(hist_dp,979,mul=1./float(numberOfCalls),&
            & file='dp_piZero_'//realToChar4(energy_gamma*1000.)//'_MeV.dat')
       write(12,'(5E15.4)') energy_gamma, sigma/numberOfCalls, &
            & sigma_cut/numberOfCalls
    end if

    close(11)
    close(12)
    write(*,*) 'positive deltaE:', count(1)
    write(*,*) 'negative deltaE:', count(2)

  contains

    !**************************************************************************
    !****s*  kruscheAnalyse/get_deltaE
    ! NAME
    ! real function get_deltaE(part)
    ! PURPOSE
    ! Evaluates for a given particle "delta E" according to eq. 4
    ! in arXiv:nucl-ex/0406002v1 (June 2004).
    ! INPUTS
    ! * type(particle), intent(in) :: part
    !**************************************************************************
    real function get_deltaE(part)
      use lorentzTrafo, only: lorentz
      use constants, only: mN

      type(particle), intent(in) :: part
      real, dimension (1:3) :: beta
      real, dimension (0:3) :: mom

      ! Define velocity of gamma N CM-System in a system where the nucleon is at rest
      beta(3)=energy_gamma/(energy_gamma+mN)
      beta(1:2)=0.

      ! Boost particle momentum to gamma N CM System
      mom=part%mom(0:3)
      call lorentz(beta,mom)

      ! Define Delta E according to eq. 4 in arXiv:nucl-ex/0406002v1 (June 2004)
      get_deltaE=mom(0) - get_EPion(energy_gamma)

    end function get_deltaE



    !**************************************************************************
    !****s*  kruscheAnalyse/get_EPion
    ! NAME
    ! real function get_EPion(E_gamma)
    ! PURPOSE
    ! Evaluates the energy of a pion in the CM-frame of a gamma Nucleon
    ! collision.
    ! INPUTS
    ! * real:: E_gamma -- Energy of gamma in GeV in frame where nucleon rests
    !**************************************************************************
    real function get_EPion(E_gamma)
      use twobodyTools, only: pcm
      use constants, only: mN, mPi

      real, intent(in) :: E_gamma
      real :: sqrtS,p
      sqrtS=sqrt((E_gamma+mN)**2-E_gamma**2)
      if (sqrtS.gt.mN+ mPi) then
         p=pCM(sqrtS, mN, mPi)
         get_EPion=sqrt(p**2+mPi**2)
      else
         get_EPion=-9999.
      end if
    end function get_EPion

  end subroutine kruscheAnalyse



  !****************************************************************************
  !****s* LowPhotonAnalysis/fissumAnalyse
  ! NAME
  ! subroutine fissumAnalyse(Parts,finalFlag)
  ! PURPOSE
  ! Evaluates the cross section for pi^+ production in the same binnings as
  ! Fissum PRC 53,3 pages 1278ff. (19996)
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts -- particle vector
  ! * logical :: finalFlag -- .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  ! OUTPUT
  !
  !****************************************************************************
  subroutine fissumAnalyse(Parts,finalFlag)
    use particleDefinition
    use output, only: realToChar4, realToChar
    use initLowPhoton, only: energy_gamma
    use degRad_conversion
    use hist_multipleRuns
    use IDTable, only: pion

    type(particle), dimension(:,:), intent(in), target :: Parts
    type(particle), pointer :: pPart

    logical, intent(in) :: finalFlag

    logical, save         :: initFlag=.true.

    ! Index 1: 51 deg, 2: 81 deg, 3: 109 deg, 4: 141 deg
    type(histogram_mr),save,dimension(1:4)  :: dsigmadThetadT_pi
    type(histogram_mr),save                 :: dsigmadTheta_greater17MeV

    integer        :: i,j,channel
    real           :: theta_lab
    real, dimension(1:4), parameter :: degs=(/51,81,109,141/)

    if (initFlag) then
       do i=1,4
          call CreateHist_mr(dsigmadThetadT_pi(i), &
               'dsigma/(dTheta dT_pion) [...b / degrees GeV] for piPlus&
               & according to Fissum PRC 53,3 (1996) at fixed &
               & angle'//realToChar(degs(i))//' degrees',0.,0.1,0.01)
       end do
       call CreateHist_mr(dsigmadTheta_greater17MeV, &
            'dsigma/dTheta  [...b / degrees] for piPlus according to &
            & Fissum PRC 53,3 (1996)', 0.,180.,3.)
       initFlag=.false.
    end if

    do i=1,4
       call startRunHist_mr(dsigmadThetadT_pi(i))
    end do
    call startRunHist_mr(dsigmadTheta_greater17MeV)

    do i=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do j=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(i,j)

          ! Select all pi^+
          if ((pPart%ID.ne.pion).or.(pPart%Charge.ne.1)) cycle

          ! Make histograms for pi^+
          theta_lab=degrees(acos(pPart%mom(3)/absMom(pPart)))
          ! write theta to histogram
          if (kineticEnergy(pPart).gt.0.017) &
               call AddHist_mr(dsigmadTheta_greater17MeV, &
               theta_lab,pPart%perweight,1.)

          ! dsigma/dOmega/dT_pi at fixed angles:
          channelLoop: do channel=1,4
             if (abs(theta_lab-degs(channel)).lt.5.) then
                call AddHist_mr(dsigmadThetadT_pi(channel),&
                     kineticEnergy(pPart),pPart%perweight/10.,1.)
                exit channelLoop
             end if
          end do channelLoop
       end do
    end do

    do i=1,4
       call endRunHist_mr(dsigmadThetadT_pi(i))
    end do
    call endRunHist_mr(dsigmadTheta_greater17MeV)

    ! Write histograms to file
    do i=1,4
       call WriteHist_mr(dsigmadThetadT_pi(i),978,&
            & file='fissum_dsigma_dThetadT_'//realToChar(degs(i))//'deg_Egamma_'//realToChar4(energy_gamma*1000.)//'_MeV.dat')
       ! Prepare error analysis
       if (finalFlag) call removeHist_mr(dsigmadThetadT_pi(i))
    end do

    call WriteHist_mr(dsigmadTheta_greater17MeV,978,&
         & file='fissum_dsigmadTheta_greater17MeV_Egamma_'//realToChar4(energy_gamma*1000.)//'_MeV.dat')
    if (finalFlag) call removeHist_mr(dsigmadTheta_greater17MeV)

    if (finalFlag) initFlag=.true.

  end subroutine fissumAnalyse

  !****************************************************************************
  !****s* LowPhotonAnalysis/etaAnalyse
  ! NAME
  ! subroutine etaAnalyse(Parts,finalFlag)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts -- particle vector
  ! * logical :: finalFlag -- .true. if it is the last call for one specific
  !   energy, therefore final output must be made.
  ! OUTPUT
  !
  ! NOTES
  ! the inclusive production cross sections are averaged over the different
  ! runs.
  !****************************************************************************
  subroutine etaAnalyse(Parts,finalFlag)
    use particleDefinition
    use constants, only: mN
    use output, only: realToChar4, realToChar, intToChar4
    use initLowPhoton, only: energy_gamma
    use IDTable, only: pion, eta
    use Averager
    use hist
    use hist2D
    use minkowski, only: abs4

    type(particle), dimension(:,:), intent(in), target :: Parts
    logical, intent(in) :: finalFlag

    logical, save :: first = .true.
    logical, save :: initFlag=.true.
    integer, save  :: nCall = 0, nCall0 = 0

    type(tAverager), save :: avePi, aveEta, aveEtaCut
    type(histogram), save :: h_T
    type(histogram2D), save :: h2D_Pos

    integer :: i,j
    real :: sumPi, sumEta, sumEtaCut, mul0, bT, missMass
    real, dimension(1:3) :: pos
    logical :: lMissMassCut
    real :: w1, w2

    real, dimension(0:3) :: pIn

    character*(200) :: Buf

    if (first) then

       open(146,file='EtaAna_XS.dat', status='unknown')
       write(146,*) '# E_gamma, sigma_pi, Err(sigma_pi), sigma_eta, Err(sigma_eta), sigma_eta^cut, Err(sigma_eta^cut)'
       close(146)

       if (useProductionPos) then
          call CreateHist2D(h2D_Pos,"Production pos",(/0.,-10./),(/10.,10./),(/.1,.1/),.true.)
       end if

       call CreateHist(h_T,"T_kin", 0., 2., 0.01)

       first = .false.
    end if

    if (initFlag) then
       call AveragerClear(avePi)
       call AveragerClear(aveEta)
       call AveragerClear(aveEtaCut)

       if (useProductionPos) then
          call ClearHist2D(h2D_Pos)
       end if

       call ClearHist(h_T)

       nCall = 0
       nCall0 = nCall0+1
       initFlag = .false.
    end if

    nCall = nCall+1

    sumPi=0
    sumEta = 0
    sumEtaCut = 0

    pIn = (/energy_gamma + mN, 0., 0., energy_gamma/)

    do i=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do j=lbound(Parts,dim=2),ubound(Parts,dim=2)

          select case(Parts(i,j)%ID)
          case (pion)
             sumPi = sumPi + Parts(i,j)%perWeight

          case (eta)

             missMass = abs4(pIn - Parts(i,j)%mom) - mN
             lMissMassCut = (missMass < 0.140)

             w1 = Parts(i,j)%perWeight
             w2 = 0.
             if (lMissMassCut) w2 = w1

             sumEta = sumEta + w1
             sumEtaCut = sumEtaCut + w2

             if (useProductionPos) then
                pos = getProductionPos(Parts(i,j))
                bT = sqrt(pos(1)**2+pos(2)**2)
                call AddHist2D(h2D_Pos, (/bT,pos(3)/),w1/bT)
             end if

             call AddHist(h_T, kineticEnergy(Parts(i,j)), w1,w2)

          end select

       end do
    end do

    call AveragerAdd(avePi,  sumPi)
    call AveragerAdd(aveEta, sumEta)
    call AveragerAdd(aveEtaCut, sumEtaCut)



    if (finalFlag) then

       ! write stuff to file:

       mul0 = 1.0/nCall

       open(146,file='EtaAna_XS.dat', position='append')
       write(146,'(f13.5,1P,6e13.5)') energy_gamma, &
            AveragerMean(avePi),  AveragerStdDev(avePi), &
            AveragerMean(aveEta), AveragerStdDev(aveEta), &
            AveragerMean(aveEtaCut), AveragerStdDev(aveEtaCut)
       close(146)

       if (useProductionPos) then
          call WriteHist2D_Gnuplot(h2D_Pos, &
               file='eta_Pos.'//inttochar4(nCall0)//'.dat', &
               mul=mul0, add=1e-20, dump=.true.)
       end if

       write(Buf,'(A," (nCall= ",i4,"; E_gamma=",f13.5,")")') &
            "dN/dT_kin", nCall0,energy_gamma
       call setNameHist(h_T,trim(Buf))
       call WriteHist(h_T, &
            file='eta_T.'//inttochar4(nCall0)//'.dat', &
            mul=mul0, add=1e-20, dump=.true.)


       initFlag=.true.
    end if

  end subroutine etaAnalyse



  !****************************************************************************
  !****s* LowPhotonAnalysis/twoPi_output
  ! NAME
  ! subroutine twoPi_output(Parts, finalFlag,elab)
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts -- particle vector
  ! * logical :: finalFlag -- .true. if it is the last call for one specific &
  !   energy, therefore final output must be made.
  ! * real, optional :: elab
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! Makes output for two-pion production processes. Writes to the following
  ! files:
  ! * Total cross sections for pi pi production : "twoPi_Total.dat"
  ! * dsigma/dm :   "sigmadm*.dat"
  ! * dsigma/dp     "sigmaMom*.dat"
  ! * dsigma/dm/dp  "sigmadm_mom*.dat"
  !
  ! where * is the given energy elab or the photon energy if elab is not
  ! present.
  !
  ! Only for photon induced events also dsigm/dr is available in file
  ! 'sigmaRadius*.dat'
  !
  ! Here r is the radius of the production point of the measured pi pi pair.
  ! NOTES
  ! If you want to use this routine for non photon-induced events, then you
  ! must provide "elab"!!!
  !****************************************************************************
  subroutine twoPi_output(Parts, finalFlag,elab)
    use particleDefinition
    use idTable, only: pion
    use output, only: paragraph, realTochar4, intTochar
    use initLowPhoton, only: energy_gamma,getCoordinate
    use densityModule, only: densityAt
    use dichteDefinition, only: dichte
    use inputGeneral, only: fullEnsemble

    type(particle), dimension(:,:), intent(in)  :: Parts
    logical, intent (in)                        :: finalFlag
    real, intent(in), optional :: elab

    !local
    integer, parameter :: dim=25            ! number of dm bins
    real, parameter    :: tresh=0.25        ! Two Pion Treshhold
    real, parameter    :: topvalue=0.52     ! Maximal mass


    real, parameter    :: cut1=0.144     ! Cuts for the momentum
    real, parameter    :: cut2=0.207     ! Cuts for the momentum

    integer, parameter :: dimMOM=25        ! number of dm bins
    integer, parameter :: dimRadius=50     ! number of radius bins
    real, parameter :: treshMOM=0.0        ! Two Pion Treshhold
    real, parameter :: topvalueMOM=0.5     ! Maximal momentum
    real, parameter :: topvalueRadius=10.0 ! Maximal radius

    real, save,dimension(1:6,0:dim) :: sigmadm           !d(Sigma)/d(mass)
    real, dimension(1:6,0:dim) :: sigmadm_new            !d(Sigma)/d(mass)

    real, save,dimension(1:6,0:dim) :: sigmadm_noColl    !d(Sigma)/d(mass) of pions which did not scatter
    real, dimension(1:6,0:dim) :: sigmadm_new_noColl     !d(Sigma)/d(mass) of pions which did not scatter


    real, save,dimension(1:6,0:dim) :: sigmadm_cut1      !d(Sigma)/d(mass) for total momentum less than cut1
    real, dimension(1:6,0:dim) :: sigmadm_new_cut1       !d(Sigma)/d(mass)              "

    real, save,dimension(1:6,0:dim) :: sigmadm_cut2      !d(Sigma)/d(mass) for total momentum greater than cut1 and less than cut2
    real, dimension(1:6,0:dim) :: sigmadm_new_cut2       !d(Sigma)/d(mass)              "

    real, save,dimension(1:6,0:dim) :: sigmadm_cut3      !d(Sigma)/d(mass) for total momentum greater than cut2
    real, dimension(1:6,0:dim) :: sigmadm_new_cut3       !d(Sigma)/d(mass)              "


    real, save,dimension(-1:1,0:dimMom) :: sigmaMomentum !d(Sigma)/dp
    real, save,dimension(1:6,0:dimRadius) :: sigmaRadius !d(Sigma)/dr sigma as function of original production point.

    real, save,dimension(1:6,0:dim) :: error_sigmadm     !d(Sigma)/d(mass)

    real, save,dimension(1:6,0:dim,0:dimMom) :: sigma_dm_Momentum !d(Sigma)/dm/dp(1) mit Momentum Information eines Teilchens

    real, save :: totalevents  !number of two-pi-events
    real :: deltaM,m,mom,deltaMom                  !delta_mass, invariant mass
    integer ::i,imass  ,k,j
    integer ::charge,kanal
    real, dimension(1:3) :: mom1,mom2

    logical ::  flag
    integer, save :: numCalls = 0 ! Counts how often the analysis routine is
    ! called before "finalFlag" is set. Therefore it counts the number of runs
    ! with the same energy.
    character(LEN=30) :: name           ! name of output files
    !character(LEN=4)   :: word

    logical,save  :: initFlag=.true.

    integer :: iEns,iPart, iPart2, iEns2
    integer :: starter, endEnsemble
    integer :: indexMomentum,indexRadius
    real :: radius, deltaRadius,energy,pairMomentum
    real, dimension(1:3) :: originalCoordinate
    type(dichte) :: density
    logical, save :: veryFirstTime =.true.

    if (present(elab)) then
       energy=elab
    else
       energy=energy_gamma
    end if

    deltaM=(topvalue-tresh)/float(dim)            ! bin size for delta(mass)
    deltaMOM=(topvalueMOM-treshMOM)/float(dimMOM) ! bin size for delta(momentum)
    deltaRadius=(topvalueRadius)/float(dimRadius) ! bin size for delta(radius)

    numCalls=numCalls+1

    ! Initialization

    write(*,paragraph) ' In TwoPi_OutPut'
    write(*,*) "DeltaM=",deltaM

    if (initFlag) then
       totalevents=0
       do i=lbound(sigmaDM,dim=2),ubound(sigmaDM,dim=2)
          sigmadm(:,i)=0.
          sigmadm_cut1(:,i)=0.
          sigmadm_noColl(:,i)=0.
          sigmadm_cut2(:,i)=0.
          sigmadm_cut3(:,i)=0.
          error_sigmadm(:,i)=0.
       end do
       do i=lbound(sigmaMomentum,dim=2),ubound(sigmaMomentum,dim=2)
          sigmaMomentum(:,i)=0.
       end do
       do i=lbound(sigmaRadius,dim=2),ubound(sigmaRadius,dim=2)
          sigmaRadius(:,i)=0.
       end do
       do i=lbound(sigma_dm_Momentum,dim=2),ubound(sigma_dm_Momentum,dim=2)
          do j=lbound(sigma_dm_Momentum,dim=3),ubound(sigma_dm_Momentum,dim=3)
             sigma_dm_Momentum(:,i,j)=0.
          end do
       end do
       initFlag=.false.
    end if


    ! ******************cross section******************************************

    do i=lbound(sigmadm_new,dim=1),ubound(sigmadm_new,dim=1)
       sigmadm_new(i,:)=0.
       sigmadm_new_noColl(i,:)=0.
       sigmadm_new_cut1(i,:)=0.
       sigmadm_new_cut2(i,:)=0.
       sigmadm_new_cut3(i,:)=0.
    end do

    ensembleLoop :do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       indexLoop : do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)
          flag=.false.
          if (Parts(iEns,iPart)%ID.eq.pion) then
             !Search second Pion
             if (fullEnsemble) then
                endEnsemble=ubound(Parts,dim=1)
             else
                endEnsemble=iEns
             end if

             ensembleLoop_2: do iEns2=iEns,endEnsemble
                if (iEns.eq.iEns2) then
                   if (iPart.eq.ubound(Parts,dim=2)) then
                      cycle ensembleLoop_2
                   else
                      starter=iPart+1
                   end if
                else
                   starter=1
                end if
                indexLoop_2: do iPart2=starter,ubound(Parts,dim=2)
                   if ((Parts(iEns,iPart)%firstEvent.eq.Parts(iEns2,iPart2)%firstEvent)&
                        &.and.(Parts(iEns2,iPart2)%ID.eq.pion)) then
                      if (Parts(iEns,iPart)%perweight-Parts(iEns2,iPart2)%perweight.gt.0.01) then
                         write(*,*) 'Problems twoPion_output'
                         write(*,*) Parts(iEns,iPart)%firstEvent,Parts(iEns2,iPart2)%firstEvent
                         write(*,*) Parts(iEns,iPart)%perweight, Parts(iEns2,iPart2)%perweight
                         stop
                      end if
                      totalevents=totalevents+1
                      flag=.true. !=found second pion
                      exit ensembleLoop_2
                   end if
                end do indexLoop_2
             end do ensembleLoop_2
          end if

          if (flag) then
             ! Kanäle belegen :
             ! 1 PiNull PiNull
             ! 2 PiNull PiPlus
             ! 3 piNull piMinus
             ! 4 piPlus piMinus
             ! 5 PiPlus PiPlus
             ! 6 piMinus piMinus
             charge= Parts(iEns,iPart)%charge+ Parts(iEns2,iPart2)%charge

             select case (charge)
             case (-2)   !piMinus piMinus
                kanal=6
             case (-1)
                kanal=3
             case (0)
                if (Parts(iEns,iPart)%charge.ne.0) then
                   kanal=4
                else
                   kanal=1
                end if
             case (1)
                kanal=2
             case (2)
                kanal=5
             end select
             mom1= Parts(iEns,iPart)%mom(1:3)
             mom2= Parts(iEns2,iPart2)%mom(1:3)


             ! ######### dsigma/dm
             m= massInv_twoPi(mom1,mom2 , Parts(iEns,iPart)%charge, &
                  Parts(iEns2,iPart2)%charge)

             imass=max(min(Nint((m-Tresh)/deltaM) ,dim),0)

             ! Sum all particles
             sigmadm_new(kanal,imass)=sigmadm_new(kanal,imass)+Parts(iEns,iPart)%perweight/deltaM


             ! Sum only these particles which didn't collide
             if (abs(Parts(iEns,iPart)%lastCollTime+Parts(iEns2,iPart2)%lastCollTime).lt.0.0001) then
                sigmadm_new_nocoll(kanal,imass)=sigmadm_new_nocoll(kanal,imass)+Parts(iEns,iPart)%perweight/deltaM
             end if


             ! Sum all particles with a specific total momentum
             pairMomentum=sqrt(Dot_product(mom1+mom2,mom1+mom2))
             if (pairMomentum.lt.cut1) then
                sigmadm_new_cut1(kanal,imass)=sigmadm_new_cut1(kanal,imass)+Parts(iEns,iPart)%perweight/deltaM

             else if (pairMomentum.ge.cut1.and.pairMomentum.le.cut2) then
                sigmadm_new_cut2(kanal,imass)=sigmadm_new_cut2(kanal,imass)+Parts(iEns,iPart)%perweight/deltaM

             else if (pairMomentum.gt.cut2) then
                sigmadm_new_cut3(kanal,imass)=sigmadm_new_cut3(kanal,imass)+Parts(iEns,iPart)%perweight/deltaM
             end if

             ! ######### dsigma/dp
             mom=AbsMom(Parts(iEns,iPart))
             indexMomentum=max(min(Nint((mom-TreshMom)/deltaMom) ,dimMom),0)

             sigmaMomentum(Parts(iEns,iPart)%charge,indexMomentum)=&
                  &sigmaMomentum(Parts(iEns,iPart)%charge,indexMomentum)&
                  & +Parts(iEns,iPart)%perweight/deltaMom

             ! ######### dsigma/dm/dp

             ! First Particle
             sigma_dm_momentum(kanal,imass,indexMomentum)=sigma_dm_momentum(kanal,imass,indexMomentum) &
                  & +Parts(iEns,iPart)%perweight/deltaMom/deltaM/2.

             mom=AbsMom(Parts(iEns2,iPart2))
             indexMomentum=max(min(int((mom-TreshMom)/deltaMom) ,dimMom),0)
             sigmaMomentum(Parts(iEns2,iPart2)%charge,indexMomentum)=&
                  &sigmaMomentum(Parts(iEns2,iPart2)%charge,indexMomentum)&
                  & +Parts(iEns2,iPart2)%perweight/deltaMom

             ! Second Particle
             sigma_dm_momentum(kanal,imass,indexMomentum)=sigma_dm_momentum(kanal,imass,indexMomentum)&
                  & +Parts(iEns2,iPart2)%perweight/deltaMom/deltaM/2.

             ! ######### dsigma/dr
             if (.not.present(elab)) then
                call getCoordinate(Parts(iEns,iPart)%firstEvent, OriginalCoordinate)

                radius=sqrt(dot_product(OriginalCoordinate,OriginalCoordinate))

                indexRadius=max(min(int(radius/deltaRadius) ,dimRadius),0)
                sigmaRadius(kanal,indexRadius)=sigmaRadius(kanal,indexRadius)&
                     & +Parts(iEns,iPart)%perweight/deltaRadius
             end if

          end if
       end do indexLoop
    end do ensembleLoop

    do i=lbound(sigmadm,dim=1),ubound(sigmadm,dim=1)
       sigmadm(i,:)=sigmadm(i,:)+sigmadm_new(i,:)
       sigmadm_noColl(i,:)=sigmadm_noColl(i,:)+sigmadm_new_noColl(i,:)
       sigmadm_cut1(i,:)=sigmadm_cut1(i,:)+sigmadm_new_cut1(i,:)
       sigmadm_cut2(i,:)=sigmadm_cut2(i,:)+sigmadm_new_cut2(i,:)
       sigmadm_cut3(i,:)=sigmadm_cut3(i,:)+sigmadm_new_cut3(i,:)
       error_sigmadm(i,:)=error_sigmadm(i,:)+sigmadm_new(i,:)**2
    end do


    ! *************************************Generate output*********************

    do i=lbound(sigmadm,dim=1),ubound(sigmadm,dim=1)

       ! Evaluate the error of the mean value:
       if (numcalls.gt.1) then
          error_sigmadm(i,:)=sqrt(abs((error_sigmadm(i,:)-1./float(numcalls)*sigmadm(i,:)**2)&
               & /float(numcalls-1))/float(numcalls))
       else
          error_sigmadm(i,:)=sigmadm(i,:)
       end if
    end do

    ! *write results to file***************************************************
    name='sigmadm'//realTochar4(ANINT(energy*1000.))//'.dat'
    open(10,File=name)
    write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
    write(10,'(A)') '### Integrated Cross sections [10^-6 barn/A]:'
    write(10,'(A)') '###  Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    write(10,'(A, 7F12.3)') ' ###', deltaM*Sum(sigmadm(1,:))/float(numCalls),deltaM*Sum(sigmadm(2,:))/float(numCalls),&
         & deltaM*Sum(sigmadm(3,:))/float(numCalls),deltaM*Sum(sigmadm(4,:))/float(numCalls),&
         & deltaM*Sum(sigmadm(5,:))/float(numCalls),deltaM*Sum(sigmadm(6,:))/float(numCalls),&
         & deltaM*(Sum(sigmadm(2,:))+Sum(sigmadm(3,:)))/float(numCalls)
    write(10,*) '### dSigma/dm: [10^-6 barn/GeV/A] '
    write(10,*) '###  mass, Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    do i=lbound(sigmadm,dim=2),ubound(sigmadm,dim=2)
       write(10,'(15F12.3)') i*deltam+Tresh, sigmadm(1:6,i)/float(numCalls),(sigmadm(2,i)+sigmadm(3,i))/float(numCalls),&
            & error_sigmadm(1:6,i)/float(numCalls),(error_sigmadm(2,i)+error_sigmadm(3,i))/float(numCalls)
    end do
    close(10)


    name='sigmadm'//realTochar4(ANINT(energy*1000.))//'_noColl.dat'
    open(10,File=name)
    write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
    write(10,'(A)') '### Integrated Cross sections [10^-6 barn/A]:'
    write(10,'(A)') '###  Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    write(10,'(A, 7F12.3)') ' ###', deltaM*Sum(sigmadm_noColl(1,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_noColl(2,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_noColl(3,:))/float(numCalls),deltaM*Sum(sigmadm_noColl(4,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_noColl(5,:))/float(numCalls),deltaM*Sum(sigmadm_noColl(6,:))/float(numCalls),&
         & deltaM*(Sum(sigmadm_noColl(2,:))+Sum(sigmadm_noColl(3,:)))/float(numCalls)
    write(10,*) '### dSigma/dm: [10^-6 barn/GeV/A] '
    write(10,*) '###  mass, Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    do i=lbound(sigmadm_noColl,dim=2),ubound(sigmadm_noColl,dim=2)
       write(10,'(15F12.3)') i*deltam+Tresh, sigmadm_noColl(1:6,i)/float(numCalls),&
            & (sigmadm_noColl(2,i)+sigmadm_noColl(3,i))/float(numCalls)
    end do
    close(10)


    name='sigmadm'//realTochar4(ANINT(energy*1000.))//'.cut1.dat'
    open(10,File=name)
    write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
    write(10,'(A)') '### Integrated Cross sections [10^-6 barn/A]:'
    write(10,'(A)') '###  Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    write(10,'(A, 7F12.3)') ' ###', deltaM*Sum(sigmadm_cut1(1,:))/float(numCalls),deltaM*Sum(sigmadm_cut1(2,:))&
         & /float(numCalls),&
         & deltaM*Sum(sigmadm_cut1(3,:))/float(numCalls),deltaM*Sum(sigmadm_cut1(4,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_cut1(5,:))/float(numCalls),deltaM*Sum(sigmadm_cut1(6,:))/float(numCalls),&
         & deltaM*(Sum(sigmadm(2,:))+Sum(sigmadm(3,:)))/float(numCalls)
    write(10,*) '### dSigma/dm: [10^-6 barn/GeV/A] '
    write(10,*) '###  mass, Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    do i=lbound(sigmadm,dim=2),ubound(sigmadm,dim=2)
       write(10,'(15F12.3)') i*deltam+Tresh, sigmadm_cut1(1:6,i)/float(numCalls),(sigmadm_cut1(2,i)+sigmadm_cut1(3,i))&
            & /float(numCalls)
    end do

    close(10)

    name='sigmadm'//realTochar4(ANINT(energy*1000.))//'.cut2.dat'
    open(10,File=name)
    write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
    write(10,'(A)') '### Integrated Cross sections [10^-6 barn/A]:'
    write(10,'(A)') '###  Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    write(10,'(A, 7F12.3)') ' ###', deltaM*Sum(sigmadm_cut2(1,:))/float(numCalls),deltaM*Sum(sigmadm_cut2(2,:))&
         & /float(numCalls),&
         & deltaM*Sum(sigmadm_cut2(3,:))/float(numCalls),deltaM*Sum(sigmadm_cut2(4,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_cut2(5,:))/float(numCalls),deltaM*Sum(sigmadm_cut2(6,:))/float(numCalls),&
         & deltaM*(Sum(sigmadm_cut2(2,:))+Sum(sigmadm_cut2(3,:)))/float(numCalls)
    write(10,*) '### dSigma/dm: [10^-6 barn/GeV/A] '
    write(10,*) '###  mass, Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    do i=lbound(sigmadm,dim=2),ubound(sigmadm,dim=2)
       write(10,'(15F12.3)') i*deltam+Tresh, sigmadm_cut2(1:6,i)/float(numCalls),(sigmadm_cut2(2,i)+sigmadm_cut2(3,i))&
            & /float(numCalls)
    end do

    close(10)

    name='sigmadm'//realTochar4(ANINT(energy*1000.))//'.cut3.dat'
    open(10,File=name)
    write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
    write(10,'(A)') '### Integrated Cross sections [10^-6 barn/A]:'
    write(10,'(A)') '###  Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    write(10,'(A, 7F12.3)') ' ###', deltaM*Sum(sigmadm_cut3(1,:))/float(numCalls),deltaM*Sum(sigmadm_cut3(2,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_cut3(3,:))/float(numCalls),deltaM*Sum(sigmadm_cut3(4,:))/float(numCalls),&
         & deltaM*Sum(sigmadm_cut3(5,:))/float(numCalls),deltaM*Sum(sigmadm_cut3(6,:))/float(numCalls),&
         & deltaM*(Sum(sigmadm_cut3(2,:))+Sum(sigmadm_cut3(3,:)))/float(numCalls)
    write(10,*) '### dSigma/dm: [10^-6 barn/GeV/A] '
    write(10,*) '###  mass, Pi0 Pi0, Pi0 Pi+,pi0 pi-,pi+ pi-, Pi+ Pi+, pi- pi-, pi0 +chargedPion'
    do i=lbound(sigmadm,dim=2),ubound(sigmadm,dim=2)
       write(10,'(15F12.3)') i*deltam+Tresh, sigmadm_cut3(1:6,i)/float(numCalls), &
            & (sigmadm_cut3(2,i)+sigmadm_cut3(3,i))/float(numCalls)
    end do

    close(10)



    name='twoPi_Total.dat'
    open(10,File=name,position='append')

    write(10,'(A)') '### 1:lab-energy  2:Pi0 Pi0   3:Pi0 Pi+  4:pi0 pi-  5:pi+ pi-  6:Pi+ Pi+  7:pi- pi-  8:pi0 +chargedPion'
    write(10,'( 8F15.6)') energy,deltaM*Sum(sigmadm(1,:))/float(numCalls),deltaM*Sum(sigmadm(2,:))/float(numCalls),&
         & deltaM*Sum(sigmadm(3,:))/float(numCalls),deltaM*Sum(sigmadm(4,:))/float(numCalls),&
         & deltaM*Sum(sigmadm(5,:))/float(numCalls),deltaM*Sum(sigmadm(6,:))/float(numCalls),&
         & deltaM*(Sum(sigmadm(2,:))+Sum(sigmadm(3,:)))/float(numCalls)
    close(10)

    if (finalFlag) then
       name='twoPi_Total_final.dat'
       if ( veryFirstTime ) then
          open(10,File=name)
          veryFirstTime=.false.
       else
          open(10,File=name,position='append')
       end if
       write(10,'(A)') '### 1:lab-energy  2:Pi0 Pi0   3:Pi0 Pi+  4:pi0 pi-  5:pi+ pi-  6:Pi+ Pi+  7:pi- pi-  8:pi0 +chargedPion'
       write(10,'( 8F15.6)') energy,deltaM*Sum(sigmadm(1,:))/float(numCalls),deltaM*Sum(sigmadm(2,:))/float(numCalls),&
            & deltaM*Sum(sigmadm(3,:))/float(numCalls),deltaM*Sum(sigmadm(4,:))/float(numCalls),&
            & deltaM*Sum(sigmadm(5,:))/float(numCalls),deltaM*Sum(sigmadm(6,:))/float(numCalls),&
            & deltaM*(Sum(sigmadm(2,:))+Sum(sigmadm(3,:)))/float(numCalls)
       close(10)
    end if


    !***************************************************************************
    !****o* LowPhotonAnalysis/sigmaMomEEEE.dat
    ! NAME
    ! file sigmaMomEEEE.dat,   EEEE= energy of the incoming photon  in MeV
    ! PURPOSE
    ! The files are produced in the runs with eventtype=3=LowPhoton
    ! They show  xsec (microbarn) for the  __two-pion events__   versus pion
    ! momentum (GeV)
    !
    ! Columns:
    ! * #1: momentum (absolute value of the pion 3-momentum)
    ! * #2: PiMinus (xsec fo pi- )
    ! * #3: PiNull  (xsec fo pi0 )
    ! * #4: PiPlus  (xsec fo pi+ )
    !***************************************************************************

    name='sigmaMom'//realTochar4(ANINT(energy*1000.))//'.dat'
    ! ANINT - round to nearest integer, return real;  introduced to preven output 199 for enrgy 200
    open(10,File=name)
    write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
    write(10,'(A)') '###  momentum, PiMinus, PiNull, PiPlus'
    do i=lbound(sigmaMomentum,dim=2),ubound(sigmaMomentum,dim=2)
       write(10,'(4F12.3)') i*deltaMom+TreshMom, sigmaMomentum(-1:1,i)/float(numCalls)
    end do
    close(10)


    !**************************************************************************
    !****o* LowPhotonAnalysis/sigmaRadiusEEEE.dat
    ! NAME
    ! file sigmaRadiusEEEE.dat,   EEEE= energy of the incoming photon  in MeV
    ! PURPOSE
    ! The files are produced in the runs with eventtype=3=LowPhoton
    ! They show  xsec (microbarn) for the  __two-pion events__   versus
    ! radius of the nucleus for various final states
    !
    ! Columns:
    ! * #1: radius[fm]
    ! * #2: baryon density[fm^-3]
    ! * #3: dsigma/dRadius (1:6) [10^-6 b/fm/A]
    ! * #4: xsec for   PiNull PiNull [microbarn]
    ! * #5  xsec for   PiNull PiPlus
    ! * #6  xsec for   piNull piMinus
    ! * #7  xsec for   piPlus piMinus
    ! * #8  xsec for   PiPlus PiPlus
    ! * #9  xsec for   piMinus piMinus
    !**************************************************************************

    if (.not.present(elab)) then
       name='sigmaRadius'//realTochar4(ANINT(energy*1000.))//'.dat'
       open(10,File=name)
       write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
       write(10,'(A)') '###  radius[fm], baryon density[fm^-3],dsigma/dRadius (1:6) [10^-6 b/fm/A]'
       write(10,'(A)') '# Kan�le :'
       write(10,'(A)') '# 1 PiNull PiNull'
       write(10,'(A)') '# 2 PiNull PiPlus'
       write(10,'(A)') '# 3 piNull piMinus'
       write(10,'(A)') '# 4 piPlus piMinus'
       write(10,'(A)') '# 5 PiPlus PiPlus'
       write(10,'(A)') '# 6 piMinus piMinus'

       do i=lbound(sigmaRadius,dim=2),ubound(sigmaRadius,dim=2)
          density=densityAt((/0.,0.,float(i)*deltaRadius/))
          write(10,'(8F12.3)') float(i)*deltaRadius, density%baryon(0),sigmaRadius(:,i)/float(numCalls)
       end do
       close(10)
    end if


    !**************************************************************************
    !****o* LowPhotonAnalysis/sigmadm_mom_EEEE_ZZZ.dat
    ! NAME
    ! file sigmadm_mom_EEEE_ZZZ.dat,
    ! EEEE= energy of the incoming photon  in MeV
    ! ZZZ = channel (predefined two-photon final state)
    ! PURPOSE
    ! The files are produced in the runs with eventtype=3=LowPhoton
    ! They show  xsec (microbarn) for the  __two-pion events__   versus pion
    ! momentum for various final states
    !
    ! Columns:
    ! * #1: mass [GeV]  = invariant mass of the two outgoing pions
    ! * #2: momentum [GeV]
    ! * #3: dsigma/dm/dp [10^-6 b/GeV^2/A]
    !**************************************************************************

    do k=1,6
       name='sigmadm_mom_'//realTochar4(ANINT(energy*1000.))//'_'//intTochar(k)//'.dat'
       open(10,File=name)
       write(10,'(A,F9.1,A,I12)') '### Counted two pion events:', totalEvents, 'number of runs:',numcalls
       write(10,'(A)') '###  mass, momentum, dsigma/dm/dp [10^-6 b/GeV^2/A]'
       write(10,'(A)') '# Filename denotes energy and channel: sigmadm_mom_"energy"_"channel".dat'
       write(10,'(A)') '# Kanaele :'
       write(10,'(A)') '# 1 PiNull PiNull'
       write(10,'(A)') '# 2 PiNull PiPlus'
       write(10,'(A)') '# 3 piNull piMinus'
       write(10,'(A)') '# 4 piPlus piMinus'
       write(10,'(A)') '# 5 PiPlus PiPlus'
       write(10,'(A)') '# 6 piMinus piMinus'

       write(10,*) '# Dimension der Matrix :', -lbound(sigma_dm_Momentum,dim=2)+ubound(sigma_dm_Momentum,dim=2),&
            &  'x',-lbound(sigma_dm_Momentum,dim=3)+ubound(sigma_dm_Momentum,dim=3)

       do i=lbound(sigma_dm_Momentum,dim=2),ubound(sigma_dm_Momentum,dim=2)
          do j=lbound(sigma_dm_Momentum,dim=3),ubound(sigma_dm_Momentum,dim=3)
             write(10,'(3F15.5)') i*deltam+Tresh, j*deltaMom+TreshMom, sigma_dm_Momentum(k,i,j)/float(numCalls)
          end do
       end do
       close(10)
    end do


    if (finalFlag) then
       initFlag=.true.
       numCalls=0
    end if

    return

  end subroutine twoPi_output


  !****************************************************************************
  !****s* LowPhotonAnalysis/twoPi_output_hist
  ! NAME
  ! subroutine twoPi_output_hist(Parts,finalFlag)
  ! PURPOSE
  ! Makes dsigma/dm output for two-pion production processes.
  !
  ! INPUTS
  !  * type(particle), dimension(:,:) :: Parts -- particle vector
  !  * logical :: finalFlag -- .true. if it is the last call for one specific
  ! energy, therefore final output must be made.
  ! OUTPUT
  ! Writes to the following files:
  ! * twopi_dsigma_dm_"energy_gamma"_MeV.dat'
  !
  ! "energy_gamma" is the photon energy
  !****************************************************************************
  subroutine twoPi_output_hist(Parts,finalFlag)
    use particleDefinition
    use output, only: realToChar4
    use initLowPhoton, only: energy_gamma
    use inputGeneral, only: fullEnsemble
    use hist_multipleRuns_mc
    use vector, only: absVec
    use IDTable, only: pion

    type(particle), dimension(:,:) :: Parts
    logical       , intent(in) :: finalFlag

    logical, save         :: initFlag=.true.

    type(histogram_mr_mc),save :: dsigmadm, dsigmadp_tot, dsigmadp_tot_noColl



    integer :: channel,charge!,i
    integer :: iEns, iPart,endEnsemble, iEns2, iPart2, starter
    logical :: flag
    real    :: m, ptot
    real, dimension(1:3) :: mom1,mom2

    if (initFlag) then
       call CreateHist_mr_mc(dsigmadm,  'dsigma/dm [mub / GeV] ',0.2,0.6,0.005,7)
       call CreateHist_mr_mc(dsigmadp_tot,  'dsigma/dp_tot [mub / GeV] ',0.,0.6,0.005,7)
       call CreateHist_mr_mc(dsigmadp_tot_noColl,  'dsigma/dp_tot [mub / GeV] of pions which did not rescatter ',0.,0.6,0.005,7)
       initFlag=.false.
       dsigmadm%total%ydesc(1:6)=(/ 'pi0 pi0', 'pi0 pi+','pi0 pi-','pi+ pi-','pi+ pi+' , 'pi- pi+' /)
       dsigmadm%total%ydesc(7)= 'pi-/+ pi0'
       dsigmadm%total%xdesc='mass(pi pi) in GeV'

       dsigmadp_tot_noColl%total%ydesc(1:6)=(/ 'pi0 pi0', 'pi0 pi+','pi0 pi-','pi+ pi-','pi+ pi+' , 'pi- pi+' /)
       dsigmadp_tot_noColl%total%ydesc(7)= 'pi-/+ pi0'
       dsigmadp_tot_noColl%total%xdesc='p_tot(pi pi) in GeV'
       dsigmadp_tot%total%ydesc(1:6)=(/ 'pi0 pi0', 'pi0 pi+','pi0 pi-','pi+ pi-','pi+ pi+' , 'pi- pi+' /)
       dsigmadp_tot%total%ydesc(7)= 'pi-/+ pi0'
       dsigmadp_tot%total%xdesc='p_tot(pi pi) in GeV'
    end if

    call startRunHist_mr_mc(dsigmadm)
    call startRunHist_mr_mc(dsigmadp_tot)
    call startRunHist_mr_mc(dsigmadp_tot_noColl)




    ensembleLoop :do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       indexLoop : do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)
          flag=.false.
          if (Parts(iEns,iPart)%ID.eq.pion) then
             !Search second Pion
             if (fullEnsemble) then
                endEnsemble=ubound(Parts,dim=1)
             else
                endEnsemble=iEns
             end if

             ensembleLoop_2: do iEns2=iEns,endEnsemble
                if (iEns.eq.iEns2) then
                   if (iPart.eq.ubound(Parts,dim=2)) then
                      cycle ensembleLoop_2
                   else
                      starter=iPart+1
                   end if
                else
                   starter=1
                end if
                indexLoop_2: do iPart2=starter,ubound(Parts,dim=2)
                   if ((Parts(iEns,iPart)%firstEvent.eq.Parts(iEns2,iPart2)%firstEvent)&
                        &.and.(Parts(iEns2,iPart2)%ID.eq.pion)) then
                      if (Parts(iEns,iPart)%perweight-Parts(iEns2,iPart2)%perweight.gt.0.01) then
                         write(*,*) 'Problems twoPion_output'
                         write(*,*) Parts(iEns,iPart)%firstEvent,Parts(iEns2,iPart2)%firstEvent
                         write(*,*) Parts(iEns,iPart)%perweight, Parts(iEns2,iPart2)%perweight
                         stop
                      end if
                      flag=.true. !=found second pion
                      exit ensembleLoop_2
                   end if
                end do indexLoop_2
             end do ensembleLoop_2
          end if

          if (flag) then
             ! Kanäle belegen :
             ! 1 PiNull PiNull
             ! 2 PiNull PiPlus
             ! 3 piNull piMinus
             ! 4 piPlus piMinus
             ! 5 PiPlus PiPlus
             ! 6 piMinus piMinus
             ! 7 piMinus/piPlus + piZero
             charge= Parts(iEns,iPart)%charge + Parts(iEns2,iPart2)%charge

             select case (charge)
             case (-2)   !piMinus piMinus
                channel=6
             case (-1)
                channel=3
             case (0)
                if (Parts(iEns,iPart)%charge.ne.0) then
                   channel=4
                else
                   channel=1
                end if
             case (1)
                channel=2
             case (2)
                channel=5
             end select
             mom1= Parts(iEns,iPart)%mom(1:3)
             mom2= Parts(iEns2,iPart2)%mom(1:3)

             ! ######### dsigma/dm
             m= massInv_twoPi(mom1,mom2 , Parts(iEns,iPart)%charge, Parts(iEns2,iPart2)%charge)

             call AddHist_mr_mc(dsigmadm,m,channel,Parts(iEns,iPart)%perweight)
             if (channel.eq.3.or.channel.eq.2) call AddHist_mr_mc(dsigmadm,m,7,Parts(iEns,iPart)%perweight)

             ! ######### dsigma/dp_tot
             ptot=AbsVec(mom1+mom2)

             ! Check for collisions of the particle using lastCollisionTime=0.
             if (abs(Parts(iEns,iPart)%lastCollTime+Parts(iEns2,iPart2)%lastCollTime).lt.0.0001) then
                call AddHist_mr_mc(dsigmadp_tot_noColl,ptot,channel,Parts(iEns,iPart)%perweight)
                if (channel.eq.3.or.channel.eq.2) call AddHist_mr_mc(dsigmadp_tot_noColl,ptot,7,Parts(iEns,iPart)%perweight)
             end if

             call AddHist_mr_mc(dsigmadp_tot,ptot,channel,Parts(iEns,iPart)%perweight)
             if (channel.eq.3.or.channel.eq.2) call AddHist_mr_mc(dsigmadp_tot,ptot,7,Parts(iEns,iPart)%perweight)

          end if

       end do indexLoop
    end do ensembleLoop


    call endRunHist_mr_mc(dsigmadm)
    call endRunHist_mr_mc(dsigmadp_tot)
    call endRunHist_mr_mc(dsigmadp_tot_noColl)


    ! Write histograms to file
    call WriteHist_mr_mc(dsigmadm,978,&
         & file='twoPi_dsigma_dm_'//realToChar4(ANINT(energy_gamma*1000.))//'_MeV.dat')
    call WriteHist_mr_mc(dsigmadp_tot_noColl,979,&
         & file='twoPi_dsigma_dp_tot_noColl'//realToChar4(ANINT(energy_gamma*1000.))//'_MeV.dat')
    call WriteHist_mr_mc(dsigmadp_tot,977,&
         & file='twoPi_dsigma_dp_tot'//realToChar4(ANINT(energy_gamma*1000.))//'_MeV.dat')
    if (finalFlag) then
       call removeHist_mr_mc(dsigmadm)
       call removeHist_mr_mc(dsigmadp_tot)
       call removeHist_mr_mc(dsigmadp_tot_noColl)
       initFlag=.true.
    end if

  end subroutine twoPi_output_hist


  !****************************************************************************
  !****s* LowPhotonAnalysis/massInv_TwoPi
  ! NAME
  ! function massInv_TwoPi(p,q,chargep,chargeq) Result(m)
  ! INPUTS
  ! * real, dimension(1:3)    ::  p
  ! * real, dimension(1:3)    ::  q
  ! * integer ::  chargep
  ! * integer ::  chargeq
  ! PURPOSE
  ! Evaluates invariant mass of two pions with momentum p,q  and charge
  ! chargep,chargeq
  ! OUTPUT
  ! real  :: m  -- m=SQRT((p+q)*(p+q)) where p and q are the 4-vectors
  ! NOTES
  ! real, parameter :: massPiZero=0.138
  ! real, parameter :: massPiCharged=0.138
  !****************************************************************************
  function massInv_TwoPi(p,q,chargep,chargeq) Result(m)
    use constants, only: mPi

    !Evaluates invariant mass of two pions with momentum pp,qq  and charge
    !chargep,chargeq

    real, INTENT(IN),dimension(1:3)    ::  p
    real, INTENT(IN),dimension(1:3)    ::  q
    integer, INTENT(IN) ::  chargep
    integer, INTENT(IN) ::  chargeq

    real  :: m !m=SQRT((p+q)*(p+q)) where p and q are the 4-vectors

    !    real, parameter :: massPiZero=0.1349764D0
    !    real, parameter :: massPiCharged=0.13956995D0

    ! Use GiBUU masses here, too
    real, parameter :: massPiZero=mPi
    real, parameter :: massPiCharged=mPi

    real ::  massq,massp,Ep,Eq,momentumSquare


    if (chargep.eq.0) then
       massp=massPiZero
    else
       massp=massPiCharged
    end if

    if (chargeq.eq.0) then
       massq=massPiZero
    else
       massq=massPiCharged
    end if

    Ep=SQRT(p(1)**2+p(2)**2+p(3)**2+massp**2)
    Eq=SQRT(q(1)**2+q(2)**2+q(3)**2+massq**2)
    momentumSquare=(p(1)+q(1))**2+(p(2)+q(2))**2+(p(3)+q(3))**2

    m=SQRT((Ep+Eq)**2-momentumSquare)

    if (m.lt.(2*massPiZero)) then
       write(*,*) "Impulse=",p,q
       write(*,*) "Massen=",massP,massQ
       write(*,*) "Invariante Masse",m
       stop
    end if
    return

  end function massInv_TwoPi



end module LowPhotonAnalysis
