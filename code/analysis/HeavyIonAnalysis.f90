!******************************************************************************
!****m* /HeavyIonAnalysis
! NAME
! module HeavyIonAnalysis
!
! PURPOSE
! Contains output routines for heavy ion collisions.
!******************************************************************************
module HeavyIonAnalysis

  use histMP

  implicit none
  private

  !****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputReal
  ! PURPOSE
  ! If .true., then the output of the real particle vector
  ! will be written to the file 'DoHIA.dat'.
  ! SOURCE
  !
  logical, save :: flag_outputReal = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputPert
  ! PURPOSE
  ! If .false., then the output of the perturbative particle vector
  ! will be written to the file 'DoHIA_pert.dat'.
  ! SOURCE
  !
  logical, save :: flag_outputPert = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/flag_outputDetailed
  ! PURPOSE
  ! Print out more detailed information at each time step
  ! from subroutine HeavyIon_evol:
  ! * rhorad_*.dat
  ! * rhoz_*.dat
  ! * rhozx_*.dat
  ! * Fields_*.dat
  ! * pauli_*.dat
  ! * dens_max.dat
  ! SOURCE
  !
  logical, save :: flag_outputDetailed = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* HeavyIonAnalysis/pionAnalysis
  ! PURPOSE
  ! This flag generates various pion spectra (p_T, m_T, y, etc).
  ! The analysis operates under the assumption of a fixed target,
  ! and expects the collision to be performed in the CMS system
  ! (cf. cmsFlag in namelist /heavyIon/).
  ! The analysis matches the one applied to the HADES data in
  ! Agakishiev et al., Eur.Phys.J. A40 (2009) 45-49.
  ! SOURCE
  !
  logical, save :: pionAnalysis = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/etaAnalysis
  ! PURPOSE
  ! This flag generates various eta spectra and eta-related analyses.
  ! SOURCE
  !
  logical, save :: etaAnalysis = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/rapBinningMode
  ! PURPOSE
  ! Select the variable the 'rapBinning' is given for:
  ! * 1: variable is y0 = y/y_cms (the normalised rapidity)
  ! * 2: variable is y
  !
  ! SOURCE
  !
  integer, save :: rapBinningMode = 1
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/rapBinning
  ! PURPOSE
  ! Rapidity binning for the pion and eta analysis
  ! (only used if pionAnalysis = .true. or etaAnalysis = .true. ).
  ! The numbers represent the binning borders in y (or y0, see rapBinningMode).
  ! For each of the bins, a separate pT and/or mT spectrum will be generated.
  !
  ! Only bins, where the upper bound is larger than the lower one are considered
  ! SOURCE
  !
  real, dimension(0:13), save :: rapBinning = (/ -0.75, -0.45, -0.15, 0.15, 0.45, 0.75, 1.05, 1.35, -99.9, -99.9, -99.9, -99.9, -99.9, -99.9 /)
  !
  ! NOTES
  ! * for the Hades AuAu analysis, you should set the bins to
  !   -0.65,-0.55,...0.75
  !****************************************************************************

  integer, save :: nRapBinning ! stores the maximal used bin number

  !****************************************************************************
  !****g* HeavyIonAnalysis/KaonAnalysis
  ! PURPOSE
  ! This flag generates various Kaon spectra and Kaon-related analyses.
  ! SOURCE
  !
  logical, save :: KaonAnalysis = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/nPartAnalysis
  ! PURPOSE
  ! This flag generates output about impact parameter and N_part
  ! SOURCE
  !
  logical, save :: nPartAnalysis = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/DensityPlot
  ! PURPOSE
  ! This flag select printing the density for several time steps
  ! SOURCE
  !
  logical, save :: DensityPlot = .false.
  !****************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/NucleonMassPlot
  ! PURPOSE
  ! This flag select printing the (invariant) mass of the nucleons for several
  ! time steps
  ! SOURCE
  !
  logical, save :: NucleonMassPlot = .false.
  !****************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/do_Tmunu
  ! SOURCE
  logical,save :: do_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output.
  !***************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/do_QRvector
  ! SOURCE
  logical,save :: do_QRvector=.false.
  ! PURPOSE
  ! Switch for QRvector output.
  !***************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/do_Glauber
  ! SOURCE
  logical,save :: do_Glauber=.false.
  ! PURPOSE
  ! Switch for Glauber-MC analysis at timestep 0
  !***************************************************************************


  !***************************************************************************
  !****g* HeavyIonAnalysis/BarMes_Tmunu
  ! SOURCE
  logical,save :: BarMes_Tmunu=.false.
  ! PURPOSE
  ! If .true., then Tmunu is calculated for baryons and mesons separately.
  !***************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/rotateZ_Tmunu
  ! SOURCE
  logical,save :: rotateZ_Tmunu=.false.
  ! PURPOSE
  ! select, whether the particles are first rotated to be aligned to the
  ! z-axis
  !***************************************************************************

  !***************************************************************************
  !****g* HeavyIonAnalysis/correctPot_Tmunu
  ! SOURCE
  integer,save :: correctPot_Tmunu = 0
  ! PURPOSE
  ! select, whether the energy is corrected for the potential or not:
  ! * 0: no correction
  ! * 1: full potential added to p0
  ! * 2: only U_b/2+U_r added to p0
  ! * 3: U_b/2+U_r added to p0 in the LRF
  !***************************************************************************

  !****************************************************************************
  !****g* HeavyIonAnalysis/selectTmunuFormat
  ! SOURCE
  integer,save :: selectTmunuFormat = 2
  ! PURPOSE
  ! select output format of Tmunu (binary encoded):
  ! * 1: ASCII
  ! * 2: Binary
  ! * 3: ASCII + Binary
  !****************************************************************************

  integer, parameter :: nSet = 5
  !****************************************************************************
  !****g* HeavyIonAnalysis/useSet
  ! SOURCE
  logical, dimension(nSet), save :: useSet = (/ .false., .true., .false., .true., .true. /)
  ! PURPOSE
  ! Array to indicate, which particle set will be used for output
  !****************************************************************************
  type(histogramMP), dimension(nSet), save :: &
       hMP_ESet, hMP_mTSet, hMP_pTSet, hMP_ySet

  type tArray
     real, dimension(:), allocatable :: v
  end type tArray
  type(tArray), dimension(nSet), save :: arrMultSet

  logical, save :: initFlag=.true.

  ! string constants may be broken over multiple continuation lines:
  character(*), parameter :: Form5 = &
       "(i4,1x,i2,1x,i8,4(1x,i4),2x,f6.3,1x,i1,3(1x,f8.3),3(1x,f11.3)&
       &1x,e13.6,1x,i5,1x,i4,1x,f6.3)"
  character(*), parameter :: Form6 = &
       "(i4,1x,i2,1x,i6,1x,i10,3(1x,i4),2x,f8.5,1x,i1,4(1x,f10.5)&
       &1x,e13.6,1x,i5,1x,i4,1x,f6.3)"

  public :: DoHeavyIonAnalysis, DoHeavyIonAnalysisTime, HeavyIon_evol
  public :: getRapBinning
  public :: calcGlauber

contains


  !****************************************************************************
  !****s* HeavyIonAnalysis/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads in the namelist "HICanalysis_Input"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !****************************************************************************
  subroutine init
    use output, only: Write_ReadingInput
    use CallStack, only: TRACEBACK
    use RMF, only: getRMF_flag

    integer :: ios, i

    !**************************************************************************
    !****n* HeavyIonAnalysis/HICanalysis_Input
    ! NAME
    ! NAMELIST /HICanalysis_Input/
    ! PURPOSE
    ! Includes the switches:
    ! * flag_outputReal
    ! * flag_outputPert
    ! * flag_outputDetailed
    ! * pionAnalysis
    ! * etaAnalysis
    ! * rapBinningMode
    ! * rapBinning
    ! * KaonAnalysis
    ! * DensityPlot
    ! * NucleonMassPlot
    ! * do_QRvector
    ! * do_Glauber
    ! * do_Tmunu
    ! * BarMes_Tmunu
    ! * rotateZ_Tmunu
    ! * correctPot_Tmunu
    ! * selectTmunuFormat
    ! * useSet
    ! * nPartAnalysis
    !**************************************************************************
    NAMELIST /HICanalysis_Input/ &
         flag_outputReal, flag_outputPert, flag_outputDetailed, &
         pionAnalysis, etaAnalysis, rapBinning, rapBinningMode, KaonAnalysis, &
         DensityPlot, NucleonMassPlot, &
         do_QRvector, do_Glauber, &
         do_Tmunu, &
         BarMes_Tmunu, rotateZ_Tmunu, correctPot_Tmunu, selectTmunuFormat, &
         useSet, nPartAnalysis

    call Write_ReadingInput('HICanalysis_Input',0)
    rewind(5)
    read(5,nml=HICanalysis_Input,iostat=ios)
    call Write_ReadingInput('HICanalysis_Input',0,ios)


    do i=ubound(rapBinning,dim=1),1,-1
       if (rapBinning(i-1) < rapBinning(i)) then
          nRapBinning = i
          exit
       end if
    end do

    write(*,*) 'flag_outputReal    : ', flag_outputReal
    write(*,*) 'flag_outputPert    : ', flag_outputPert
    write(*,*) 'flag_outputDetailed: ', flag_outputDetailed
    write(*,*) 'pionAnalysis       : ', pionAnalysis
    write(*,*) 'etaAnalysis        : ', etaAnalysis
    if (pionAnalysis.or.etaAnalysis) then
       select case (rapBinningMode)
       case (1)
          write(*,*) 'rapBinningMode     : ',rapBinningMode," =y_0"
       case (2)
          write(*,*) 'rapBinningMode     : ',rapBinningMode," =y"
       case default
          call TRACEBACK("wrong rapBinningMode")
       end select
       write(*,'(A,99F7.3)') ' rapBinning         : ', &
            rapBinning(0:nRapBinning)
    end if
    write(*,*) 'KaonAnalysis       : ', KaonAnalysis
    write(*,*) 'nPartAnalysis      : ', nPartAnalysis
    write(*,*) 'DensityPlot        : ', DensityPlot
    write(*,*) 'NucleonMassPlot    : ', NucleonMassPlot
    write(*,*) 'do Q-/R-vector     : ', do_QRvector
    write(*,*) 'do Glauber-MC      : ', do_Glauber
    write(*,*) 'do Tmunu           : ', do_Tmunu
    if (do_Tmunu) then
       write(*,*) '  Tmunu: BarMes    : ', barMes_Tmunu
       write(*,*) '  Tmunu: rotateZ   : ', rotateZ_Tmunu
       if (getRMF_flag()) then
          if (correctPot_Tmunu>0) then
             write(*,*) 'correctPot_Tmunu>0 in RMF mode not possible. Reset.'
             correctPot_Tmunu = 0
          end if
       end if
       write(*,*) '  Tmunu: correctPot: ', correctPot_Tmunu
       write(*,*) '  Tmunu format     : ', selectTmunuFormat

       if (selectTmunuFormat==0) then
          call Traceback("selectTmunuFormat must not == 0")
       end if
    end if
    write(*,*) 'use MP set         : ',useSet

    call Write_ReadingInput('HICanalysis_Input',1)

    initFlag=.false.

  end subroutine init


  !****************************************************************************
  !****s* HeavyIonAnalysis/DoHeavyIonAnalysis
  ! NAME
  ! subroutine DoHeavyIonAnalysis(realParts,pertParts,finalFlag)
  !
  ! PURPOSE
  ! Does some analysis at the end of the run.
  !
  ! INPUTS
  ! * type(particle), dimension(:,:), :: realParts -- real particle vector
  ! * type(particle), dimension(:,:), :: pertParts -- perturb. particle vector
  ! * logical, intent(in) :: finalFlag -- if .true., it is the last run
  !   for one specific energy, therefore final output must be done.
  !
  ! RESULT
  ! some files written to disk
  !****************************************************************************
  subroutine DoHeavyIonAnalysis(realParts, pertParts, finalFlag)

    use IdTable, only: nucleon, pion, EOV, NOP, isMeson
    use particleDefinition
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b, particleId, antiparticle, perturbative
    use initElementary, only: b_ele => impactParameter
    use twoBodyStatistics, only: sqrts_distribution, rate
    use inputGeneral, only: eventtype
    use eventtypes, only: HeavyIon, hadron, elementary, HiLepton, HiPion
    use history, only: history_getParents,history_getGeneration

    type(particle), dimension(:,:), intent(in), target :: realParts
    type(particle), dimension(:,:), intent(in), target :: pertParts
    logical,                        intent(in)         :: finalFlag

    integer :: nEns,nPart
    real :: stossParameter
    integer, save :: isu=0     ! counter of subsequent runs
    type(particle), dimension(1:2) :: dummy

    type(particle), POINTER :: pPart

    if (initFlag) call init

    isu = isu + 1

    select case (eventtype)
    case (HeavyIon)
       stossParameter = b_HI
    case (hadron)
       stossParameter = b_had
    case (elementary)
       stossParameter = b_ele
    case default
       write(*,*) ' Problem in DoHeavyIonAnalysis: ', &
            'impact parameter is not defined'
       stossParameter=0.
    end select

    nEns = size(realParts,dim=1)
    nPart = size(realParts,dim=2)

    if (flag_outputReal) call write_DoHIA
    if (flag_outputPert) call write_DoHIApert

    if (finalFlag) then
       call sqrts_distribution(dummy,0,.true.)
       call rate(dummy,dummy,0.,.true.)
    end if

    if (eventtype.ne.HiLepton) then
       call Spectra
       if (eventtype.ne.HiPion) then
          call countParts(realParts, -99.9, "final_")
       else
          call countParts(pertParts, -99.9, "final_")
       end if
    end if


    if (eventtype==hadron &
         .and. particleId==nucleon &
         .and. antiparticle &
         .and. .not.perturbative ) then
       call histo1
    end if

    if (pionAnalysis) call analyze_Pions
    if (etaAnalysis) call analyze_Etas
    if (KaonAnalysis) call analyze_Kaons

    if (nPartAnalysis) call analyze_Npart

  contains

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/write_DoHIA
    !**************************************************************************
    subroutine write_DoHIA

      integer :: i,j,indFree,parents(1:3),generation,factor

      !***********************************************************************
      !****o* HeavyIonAnalysis/DoHIA.dat
      ! NAME
      ! file DoHIA.dat
      ! PURPOSE
      ! Contains the full dump of the real particle vector at the end of the
      ! simulation, including particle ID, charge, position, momentum, etc.
      !***********************************************************************
      open(30,file='DoHIA.dat',position='Append')
      do i = 1,nEns
         do j = 1,nPart
            pPart => realParts(i,j)
            if (pPart%ID <= 0) cycle

            factor = merge(-1, 1, pPart%anti)
            indFree = merge(1, 0, isFree(realParts,i,j,pPart))

            generation=history_getGeneration(pPart%history)
            parents = history_getParents(pPart%history)
            where (parents>200)
               parents = 200-parents
            end where

            write(30,Form5) &
                 factor*pPart%ID,pPart%charge,&
                 pPart%event(1),parents(1:3),generation,pPart%mass,&
                 indFree,pPart%pos(1:3),&
                 pPart%mom(1:3),&
                 pPart%perweight,i,isu,stossParameter
         end do
      end do
      close(30)
    end subroutine write_DoHIA

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/write_DoHIApert
    !**************************************************************************
    subroutine write_DoHIApert

      integer :: i,j,k, indFree,parents(1:3),factor

      !***********************************************************************
      !****o* HeavyIonAnalysis/DoHIA_pert.dat
      ! NAME
      ! file DoHIA_pert.dat
      ! PURPOSE
      ! Contains the full dump of the perturbative particle vector at the end
      ! of the simulation, including particle ID, charge, position, momentum,
      ! etc.
      !***********************************************************************
      open(30,file='DoHIA_pert.dat',position='Append')
      do i = 1,nEns
         do j = 1,size(pertParts,dim=2) ! nPart for pert vector

            pPart => pertParts(i,j)
            if (pPart%ID == EOV) exit
            if (pPart%ID == NOP) cycle

            factor = merge(-1, 1, pPart%anti)
            indFree = merge(1, 0, isFree(realParts,i,j,pPart))

            parents = history_getParents(pPart%history)
            where (parents>200)
               parents = 200-parents
            end where

            write(30,Form6) &
                 factor*pPart%ID,pPart%charge,pPart%firstEvent,&
                 pPart%event(1),parents(1:3),pPart%mass,&
                 indFree, (pPart%mom(k), k=0,3),&
                 pPart%perWeight,i,isu,stossParameter

         end do
      end do
      close(30)
    end subroutine write_DoHIApert

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/histo1
    !**************************************************************************
    subroutine histo1

      integer, parameter :: N_max=200  ! Max number of pions produced in event
      real, save :: P_Npion(0:N_max)   ! Pion multiplicity distribution
      integer, parameter :: Nmom=100   ! Number of momentum bins
      real, parameter :: dmom=0.02     ! Momentum bin (GeV/c)
      real, save, dimension(1:Nmom,0:10) :: dNpiondMom  ! Pion momentum distr.
      real, save :: pion_events

      integer :: i,j,numPions,ibin,factor
      logical :: flag_not_only_pions
      real :: fnorm,pinumAv,momentumAbs
      type(particle), POINTER :: pPart


      if (isu.eq.1) then
         pion_events=0.
         P_Npion=0.
         dNpiondMom=0.
         open(36,file='FewPionEvents.dat')
      end if

      Ensemble_loop : do i = 1,nEns

         numPions=0
         flag_not_only_pions=.false.

         Particle_loop1 : do j = 1,size(realParts,dim=2)

            if (realParts(i,j)%ID <= 0) cycle Particle_loop1

            if (realParts(i,j)%ID.eq.pion) then
               numPions=numPions+1
            else if (isMeson(realParts(i,j)%ID)) then
               flag_not_only_pions=.true.
            end if

         end do Particle_loop1


         if (.not.flag_not_only_pions) then

            pion_events=pion_events+1.
            if (numPions.le.N_max)  P_Npion(numPions)=P_Npion(numPions)+1.

            Particle_loop2 : do j = 1,size(realParts,dim=2)

               pPart => realParts(i,j)
               if (pPart%ID <= 0) cycle Particle_loop2

               if (numPions.le.1) then
                  factor = merge( -1., 1., pPart%anti)

                  write(36,5) int(factor)*pPart%ID,pPart%charge,&
                       pPart%mass,&
                       pPart%pos(1:3),&
                       pPart%mom(1:3),&
                       i,isu
5                 format(i4,1x,i2,1x,f6.3,3(1x,f8.3),3(1x,f8.3),1x,i5,1x,i4)
               end if

               if (pPart%ID .ne. pion) cycle Particle_loop2


               momentumAbs=sqrt(dot_product(pPart%mom(1:3),pPart%mom(1:3)))
               ibin=nint((momentumAbs-dmom/2.)/dmom)+1
               if (ibin.ge.1 .and. ibin.le.Nmom .and. pPart%charge.ne.0) then
                  dNpiondMom(ibin,0)=dNpiondMom(ibin,0)+1.
                  if (numPions.ge.1 .and. numPions.le.10) &
                       & dNpiondMom(ibin,numPions)=dNpiondMom(ibin,numPions)+1.
               end if

            end do Particle_loop2

         end if


      end do Ensemble_loop

      if (finalFlag) then

         ! Normas:
         !fnorm=1./float(nEns)/float(isu)
         if (pion_events.gt.0.) then
            fnorm=1./pion_events
         else
            fnorm=0.
         end if

         dNpiondMom(:,:)=dNpiondMom(:,:)*fnorm/dmom
         P_Npion(:)=P_Npion(:)*fnorm


         open(36,file='DoHIA1.dat')
         write(36,*)'# Distribution of events in pion number'
         write(36,*)'# Number of events: ', nEns*isu
         write(36,*)'# Npion:    P_Npion:'
         pinumAv=0.
         do j=0,N_max
            pinumAv=pinumAv+float(j)*P_Npion(j)
            write(36,*) j, P_Npion(j)
         end do
         write(36,*)'# Norma: ', sum(P_Npion(:))
         write(36,*)'# Average pion number:', pinumAv
         close(36)

         open(36,file='DoHIA2.dat')
         write(36,*)'# Charged pion momentum distribution'
         write(36,*)'# Number of events: ', nEns*isu
         write(36,*)'# momentum, GeV/c:    dNpiondMom, c/GeV:'
         do j=1,Nmom
            write(36,'(f6.3,8(1x,e10.3))')  (float(j)-0.5)*dmom, dNpiondMom(j,0), &
                 & dNpiondMom(j,1:6), sum(dNpiondMom(j,7:10))
         end do
         write(36,*)'# Norma: ', sum(dNpiondMom(:,0))*dmom
         close(36)

      end if

    end subroutine histo1

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/analyze_Pions
    !**************************************************************************
    subroutine analyze_Pions
      use histMC
      use IDTable, only: EOV
      use minkowski, only: abs4, abs3
      use inputGeneral, only: num_Runs_SameEnergy, num_energies
      use initHeavyIon, only: ekin_lab_Projectile
      use constants, only: mN
      use nucleusDefinition
      use nucleus, only: getProjectile
      use histMC_avg
      use output, only: intToChar2

      type(histogramMC), save :: hist_pT, hist_y, hist_y0, hist_cost
      type(histogramMC), dimension(0:13), save :: hist_mT ! 0=total; 1:..=different rap. bins
      type(histogramMC_avg), save :: hist_v1_y, hist_v2_y, hist_v1_u, hist_v2_u
      type(histogramMC_avg), save :: hist_cos2t
      logical, save :: first = .true.
      real, save :: w, y_cms, u_proj
      integer :: i, j, ch, rb
      real :: mom(0:3)
      real :: beta(1:3)
      real :: pT, mT, m, y, y0, cost, E, p, pabs, v1, v2, ut0, beta_proj, gamma
      real :: cos2t, yTmp
      character(40) :: str
      type(tnucleus), pointer :: proj
      integer, save :: n = 0
      real :: fac

      if (first) then
         call CreateHistMC(hist_pT, 'Pion transverse momentum spectra: dN/dpT',&
              0.,1.,0.01,3)
         call CreateHistMC(hist_y,  'Pion rapidity spectra: dN/dy', &
              -3.,3.,0.1 ,3)
         call CreateHistMC(hist_y0, 'Pion rapidity spectra: dN/dy0', &
              -3.,3.,0.1 ,3)
         call CreateHistMC(hist_mT, 'Pion transverse mass spectra: mT^(-2)*dN/dmT', &
              0.,1.,0.01,3)
         call CreateHistMC(hist_cost, 'Pion polar angle spectra: dN/dcos(theta) with 0.2<p<0.8', &
              -1.,1.,0.05,3)
         call CreateHistMC_avg(hist_v1_y, 'Directed flow of pions: v1(y0)', &
              -2.,2.,0.2 ,3)
         call CreateHistMC_avg(hist_v2_y, 'Elliptic flow of pions: v2(y0)', &
              -2.,2.,0.2 ,3)
         call CreateHistMC_avg(hist_v1_u, 'Directed flow of pions: v1(ut0)', &
              0.,4.,0.2 ,3)
         call CreateHistMC_avg(hist_v2_u, 'Elliptic flow of pions: v2(ur0)', &
              0.,4.,0.2 ,3)
         call CreateHistMC_avg(hist_cos2t, '< cos^2(theta) >(p_cm)', &
              0.,1.,0.01,3)
         do i=1,nRapBinning
            write(str,'(A,i1,A,f5.2,A,f5.2,A)') &
                 " (rap. bin #",i,": y0 = ",rapBinning(i-1),&
                 " ... ",rapBinning(i),")"
            hist_mT(i)%name = trim(hist_mT(i)%name) // str
         end do
         hist_pT%xDesc = 'p_T [GeV]'
         hist_pT%yDesc = (/ "pi-", "pi0", "pi+" /)
         call CopyDesc(hist_y, hist_pT)
         call CopyDesc(hist_y0, hist_pT)
         call CopyDesc(hist_cost, hist_pT)
         call CopyDesc(hist_mT, hist_pT)
         hist_y%xDesc    = 'y'
         hist_y0%xDesc   = 'y0'
         hist_cost%xDesc = 'cos(theta_cm)'
         hist_mT%xDesc   = 'm_T - m [GeV]'

         hist_v1_y%xDesc   = 'y0'
         hist_v1_y%yDesc = (/ "pi-", "pi0", "pi+" /)
         call CopyDesc_avg(hist_v2_y, hist_v1_y)
         call CopyDesc_avg(hist_v1_u, hist_v1_y)
         call CopyDesc_avg(hist_v2_u, hist_v1_y)
         hist_v1_u%xDesc   = 'ut0'
         hist_v2_u%xDesc   = 'ut0'
         call CopyDesc_avg(hist_cos2t, hist_v1_y)
         hist_cos2t%xDesc  = 'p_cm'

         w = 1. / (nEns*num_Runs_SameEnergy*num_energies)   ! weight factor
         E = 2*mN + ekin_lab_Projectile
         p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
         y_cms = 0.5*log((E+p)/(E-p))
         write(*,*) "analyze_Pions: y_cms = ", y_cms
         proj => getProjectile()
         beta_proj = sqrt(sum(proj%vel**2))
         u_proj = beta_proj / sqrt(1. - beta_proj**2)
         write(*,*) "analyze_Pions: u_proj = ", u_proj
         first = .false.
      end if

      do i=1,nEns
         do j=1,nPart
            if (realParts(i,j)%ID == EOV) exit
            if (realParts(i,j)%ID /= pion) cycle

            ch = realParts(i,j)%charge+2
            mom = realParts(i,j)%mom
            beta = mom(1:3) / mom(0)
            gamma = 1./ sqrt( 1. - dot_product(beta,beta) )
            m = abs4(mom)                                 ! mass
            pabs = abs3(mom)                              ! absolute momentum
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! 'normalized' rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass

            cost = mom(3)/pabs
            v1 = mom(1)/pt                       ! v1 = < px/pT >
            v2 = (mom(1)**2 - mom(2)**2)/pt**2   ! v2 = < (px**2-py**2)/pT**2 >
            ut0 = gamma * sqrt(beta(1)**2+beta(2)**2) / u_proj  ! transverse comp. of scaled four-velocity

            cos2t= 2*cost**2 - 1.0 ! = cos(2 theta)

            call AddHistMC(hist_pT,    pt,   ch, w)
            call AddHistMC(hist_y,     y,    ch, w)
            call AddHistMC(hist_y0,    y0,   ch, w)
            call AddHistMC(hist_mT(0), mt-m, ch, w/mt**2)
            if (pabs>0.2 .and. pabs<0.8) &
                 call AddHistMC(hist_cost,  cost, ch, w)
            call AddHistMC_avg(hist_cos2t, pabs, ch, cos2t)
            if (ut0>1.0 .and. ut0<4.2) then
              call AddHistMC_avg(hist_v1_y, y0, ch, v1)
              call AddHistMC_avg(hist_v2_y, y0, ch, v2)
            end if
            if (y0>-1.8 .and. y0<0.) then
              call AddHistMC_avg(hist_v1_u, ut0, ch, v1)
              call AddHistMC_avg(hist_v2_u, ut0, ch, v2)
            end if
            ! determine rap. bin
            select case (rapBinningMode)
            case (1)
               yTmp = y0
            case (2)
               yTmp = y
            end select
            rb = 0
            do while (rb < nRapBinning+1)
               if (yTmp < rapBinning(rb)) exit
               rb = rb + 1
            end do
            if (rb>0 .and. rb< nRapBinning+1) &
                 call AddHistMC(hist_mT(rb), mt-m, ch, w/mt**2)
         end do
      end do

      n = n + 1                                        ! count runs
      fac = num_Runs_SameEnergy*num_energies/float(n)  ! multiplication factor for writing histograms

      call WriteHistMC(hist_pT,  'PionPt.dat',   mul=fac)
      call WriteHistMC(hist_y,   'PionY.dat',    mul=fac)
      call WriteHistMC(hist_y0,  'PionY0.dat',   mul=fac)
      call WriteHistMC(hist_cost,'PionCost.dat', mul=fac)
      do i=0,nRapBinning
         str = ""
         if (i>0) str = "_rapBin"//intToChar2(i)
         call WriteHistMC(hist_mT(i), 'PionMt'//trim(str)//'.dat', mul=fac)
      end do
      call WriteHistMC_avg(hist_v1_y,'PionV1_y.dat')
      call WriteHistMC_avg(hist_v2_y,'PionV2_y.dat')
      call WriteHistMC_avg(hist_v1_u,'PionV1_u.dat')
      call WriteHistMC_avg(hist_v2_u,'PionV2_u.dat')
      call WriteHistMC_avg(hist_cos2t,'PionCos2t.dat', dump=.true.)
    end subroutine analyze_Pions

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/analyze_Etas
    !**************************************************************************
    subroutine analyze_Etas

      use histMC
      use IDTable, only: EOV, pion, eta
      use minkowski, only: abs4, abs3
      use inputGeneral, only: num_Runs_SameEnergy, num_energies
      use initHeavyIon, only: ekin_lab_Projectile
      use constants, only: mN
      use nucleusDefinition
      use nucleus, only: getProjectile
      use output, only: intToChar2

      logical, save :: first = .true.
      real, save :: w, y_cms

      type(histogramMC), save :: hist_y0, hist_cost

      ! 0=total; 1:..=different rap. bins
      type(histogramMC), dimension(0:13), save :: hist_mT, hist_pT

      integer,save :: n = 0
      integer :: i, j, ch, rb
      real :: E, p, m, y, y0, pt, mt, pabs, cost, fac, yTmp
      real :: mom(0:3)
      character(40) :: str

      write(*,*) 'subroutine analyze_Etas'

      if (first) then
         call CreateHistMC(hist_y0, 'rapidity spectra: dN/dy0', &
              -3.,3.,0.1,4)
         call CreateHistMC(hist_cost, 'polar angle spectra: dN/dcos(theta)', &
              -1.,1.,0.05,4)
         call CreateHistMC(hist_pT, 'transverse momentum spectra: dN/dpT',&
              0.,1.,0.01,4)
         call CreateHistMC(hist_mT, 'transverse mass spectra: mT^(-2)*dN/dmT', &
              0.,1.,0.01,4)

         do i=1,nRapBinning
            write(str,'(A,i1,A,f5.2,A,f5.2,A)') &
                 " (rap. bin #",i,": y0 = ",rapBinning(i-1),&
                 " ... ",rapBinning(i),")"
            hist_mT(i)%name = trim(hist_mT(i)%name) // str
            hist_pT(i)%name = trim(hist_pT(i)%name) // str
         end do
         hist_y0%xDesc = 'y0'
         hist_y0%yDesc = (/ "pi-", "pi0", "pi+", "eta" /)
         call CopyDesc(hist_cost, hist_y0)
         call CopyDesc(hist_pT, hist_y0)
         call CopyDesc(hist_mT, hist_y0)
         hist_cost%xDesc = 'cos(theta)'
         hist_pT%xDesc   = 'p_T [GeV]'
         hist_mT%xDesc   = 'm_T - m [GeV]'

         w = 1. / (nEns*num_Runs_SameEnergy*num_energies)   ! weight factor
         E = 2*mN + ekin_lab_Projectile
         p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
         y_cms = 0.5*log((E+p)/(E-p))
         write(*,*) "analyze_Etas: y_cms = ", y_cms
         first = .false.
      end if

      do i=1,nEns
         do j=1,nPart
            select case(realParts(i,j)%ID)
            case(EOV)
               exit
            case(pion)
               ch = realParts(i,j)%charge+2
            case(eta)
               ch = 4
            case default
               cycle
            end select

            mom = realParts(i,j)%mom
            m = abs4(mom)
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! normalized rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass
            pabs = abs3(mom)                              ! absolute momentum
            cost = mom(3)/pabs

            call AddHistMC(hist_y0,    y0,   ch, w)
            call AddHistMC(hist_cost,  cost, ch, w)
            call AddHistMC(hist_pT(0), pt,   ch, w)
            call AddHistMC(hist_mT(0), mt-m, ch, w/mt**2)

            select case (rapBinningMode)
            case (1)
               yTmp = y0
            case (2)
               yTmp = y
            end select
            rb = 0
            do while (rb < nRapBinning+1)
               if (yTmp < rapBinning(rb)) exit
               rb = rb + 1
            end do
            if (rb>0 .and. rb< nRapBinning+1) then
               call AddHistMC(hist_pT(rb), pt,   ch, w)
               call AddHistMC(hist_mT(rb), mt-m, ch, w/mt**2)
            end if

         end do
      end do

      n = n + 1                                        ! count runs
      fac = num_Runs_SameEnergy*num_energies/float(n)  ! multiplication factor for writing histograms

      call WriteHistMC(hist_y0,   'EtaY0.dat',   mul=fac, dump=.true.)
      call WriteHistMC(hist_cost, 'EtaCost.dat', mul=fac, dump=.true.)
      do i=0,nRapBinning
         str = ""
         if (i>0) str = "_rapBin"//intToChar2(i)
         call WriteHistMC(hist_mT(i), 'EtaMt'//trim(str)//'.dat', mul=fac, &
              dump=.true.)
         call WriteHistMC(hist_pT(i), 'EtaPt'//trim(str)//'.dat', mul=fac, &
              dump=.true.)
      end do

    end subroutine analyze_Etas

    !**************************************************************************
    !***is* DoHeavyIonAnalysis/analyze_Kaons
    !**************************************************************************
    subroutine analyze_Kaons
      use histMC
      use initHeavyIon, only: ekin_lab_Projectile
      use idTable, only: Kaon, Kaonbar, isBaryon
      use constants, only: mN
      use minkowski, only: abs4, abs3
      use inputGeneral, only: num_Runs_SameEnergy, num_energies

      type(histogramMC), save :: hist_p, hist_pT, hist_y, hist_y0, &
           hist_cost, hist_mt, hist_parents
      logical, save :: first = .true.
      real, save :: w, y_cms
      integer :: i, j, ch, parents(1:3), channel
      real :: mom(0:3), pT, mT, m, y, y0, cost, E, p, pabs

      if (first) then
         call CreateHistMC(hist_p,  'Kaon momentum spectra: dN/dp', &
              0.,1.5,0.05,4)
         call CreateHistMC(hist_pT, 'Kaon transverse momentum spectra: dN/dpT',&
              0.,1.,0.01,4)
         call CreateHistMC(hist_y,  'Kaon rapidity spectra: dN/dy', &
              -3.,3.,0.1 ,4)
         call CreateHistMC(hist_y0, 'Kaon rapidity spectra: dN/dy0', &
              -3.,3.,0.1 ,4)
         call CreateHistMC(hist_mT, 'Kaon transverse mass spectra: mT^(-2)*dN/dmT', &
              0.,1.,0.01,4)
         call CreateHistMC(hist_cost, 'Kaon polar angle spectra: dN/dcos(theta)', &
              -1.,1.,0.05,4)
         call CreateHistMC(hist_parents, 'Kaon parent sources', &
              0.5,5.5,1.,4)
         hist_p%xDesc = 'p [GeV]'
         hist_p%yDesc = (/ "K0   ", "K+   ", "K-   ", "K0bar" /)
         call CopyDesc(hist_pT,      hist_p)
         call CopyDesc(hist_y,       hist_p)
         call CopyDesc(hist_y0,      hist_p)
         call CopyDesc(hist_cost,    hist_p)
         call CopyDesc(hist_mT,      hist_p)
         call CopyDesc(hist_parents, hist_p)
         hist_pT%xDesc   = 'p_T [GeV]'
         hist_y%xDesc    = 'y'
         hist_y0%xDesc   = 'y0'
         hist_cost%xDesc = 'cos(theta_cm)'
         hist_mT%xDesc   = 'm_T - m [GeV]'
         hist_parents%xDesc = "source channel"
         w = 1. / (nEns*num_Runs_SameEnergy*num_energies)   ! weight factor
         E = 2*mN + ekin_lab_Projectile
         p = sqrt(ekin_lab_Projectile**2 + 2.*mN*ekin_lab_Projectile)  ! assumption: fixed target
         y_cms = 0.5*log((E+p)/(E-p))
         write(*,*) "analyze_Kaons: y_cms = ",y_cms
         first = .false.
      end if

      do i=1,nEns
         do j=1,nPart
            if (realParts(i,j)%ID == EOV) exit
            if (realParts(i,j)%ID == Kaon) then
               ch = realParts(i,j)%charge+1
            else if (realParts(i,j)%ID == Kaonbar) then
               ch = realParts(i,j)%charge+4
            else
               cycle
            end if

            mom = realParts(i,j)%mom
            m = abs4(mom)                                 ! mass
            pabs = abs3(mom)                              ! absolute momentum
            y = 0.5*log((mom(0)+mom(3))/(mom(0)-mom(3)))  ! rapidity
            y0 = y/y_cms                                  ! 'normalized' rapidity
            pt = sqrt(mom(1)**2+mom(2)**2)                ! transverse momentum
            mt = sqrt(m**2+pt**2)                         ! transverse mass
            cost = mom(3)/pabs
            call AddHistMC(hist_p,    pabs, ch, w)
            call AddHistMC(hist_pT,     pt, ch, w)
            call AddHistMC(hist_y,       y, ch, w)
            call AddHistMC(hist_y0,     y0, ch, w)
            call AddHistMC(hist_mT,   mt-m, ch, w/mt**2)
            call AddHistMC(hist_cost, cost, ch, w)
            ! determine where it came from
            parents = history_getParents (realParts(i,j)%history)
            if (parents(2) > 0) then
              ! 2-body collisions
              if (isBaryon(parents(1)) .and. isBaryon(parents(2))) then
                channel = 1  ! channel 1: BB
              else if (isMeson(parents(1)) .and. isMeson(parents(2))) then
                channel = 2  ! channel 2: mm
              else
                channel = 3  ! channel 3: mB
              end if
!               print *,"analyze_Kaons (2-body):", ch, parents(1:3)
            else if (isMeson(parents(1))) then
              channel = 4   ! channel 4: meson decay
!               print *,"analyze_Kaons (meson dec.):", ch, parents(1:3)
            else if (isBaryon(parents(1))) then
              channel = 5   ! channel 5: baryon decay
!               print *,"analyze_Kaons (baryon dec.):", ch, parents(1:3)
            else
              write(*,*) "analyze_Kaons: unknown parents!", ch, parents(1:3), realParts(i,j)%history
              stop
            end if
            call AddHistMC(hist_parents, float(channel), ch, w)
         end do
      end do
      call WriteHistMC(hist_p,       'KaonPlab.dat')
      call WriteHistMC(hist_pT,      'KaonPt.dat')
      call WriteHistMC(hist_y,       'KaonY.dat')
      call WriteHistMC(hist_y0,      'KaonY0.dat')
      call WriteHistMC(hist_mT,      'KaonMt.dat')
      call WriteHistMC(hist_cost,    'KaonCost.dat')
      call WriteHistMC(hist_parents, 'KaonParents.dat')
    end subroutine analyze_Kaons


    !**************************************************************************
    !***is* DoHeavyIonAnalysis/Spectra
    !**************************************************************************
    subroutine Spectra

      use histMP
      use IDtable, only: isHadron

      integer :: i,j, iSet
      real :: mulFak, mom, mom0, mT, pT, y
      type(particle), POINTER :: pPart
      real, parameter :: dy = 0.05
      logical, save :: first = .true.

      if (first) then
         do iSet=1,nSet
            if (.not.useSet(iSet)) cycle
            call CreateHistMP(hMP_ESet(iSet), "dN/pE dE", 0.0, 2.5, 0.02, iSet)
            call CreateHistMP(hMP_mTSet(iSet), "dN/mT2 dmT dy", 0.0, 2.5, 0.02, iSet)
            call CreateHistMP(hMP_pTSet(iSet), "dN/dpT dy", 0.0, 2.5, 0.02, iSet)
            call CreateHistMP(hMP_ySet(iSet), "dN/dy", -3.0, 3.0, 0.02, iSet)
         end do
         first = .false.
      end if

      nEns  = size(realParts,dim=1)
      nPart = size(realParts,dim=2)
      mulfak = 1.0/(nEns)

      do i=1,nEns
         do j=1,nPart
            pPart => realParts(i,j)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            mom = absMom(pPart)
            mom0 = pPart%mom(0)
            mT = sqrt(max(1e-15,pPart%mom(0)**2-pPart%mom(3)**2))
            y = rapidity(pPart)
            pT = sqrt(pPart%mom(1)**2+pPart%mom(2)**2)

            do iSet=1,nSet
               if (.not.useSet(iSet)) cycle
               call AddHistMP(hMP_ESet(iSet), pPart, mom0, 1.0/(mom0*mom), 1.0)
               call AddHistMP(hMP_ySet(iSet), pPart, y, 1.0, 1.0)

               if (abs(y) < dy) then
                  call AddHistMP(hMP_pTSet(iSet), pPart, pT, 1.0, 1.0)
                  call AddHistMP(hMP_mTSet(iSet), pPart, mT, 1.0/(2*dy*mT**2), 1.0)
               end if
            end do

         end do
      end do

      if (finalFlag) then
         do iSet=1,nSet
            if (.not.useSet(iSet)) cycle

            call WriteHistMP(hMP_ESet(iSet), &
                 file='E_Set'//achar(48+iSet)//'_final.dat', &
                 add=1e-20, mul=mulfak, iColumn=1, dump=.true.)
            call WriteHistMP(hMP_mTSet(iSet), &
                 file='mT_Set'//achar(48+iSet)//'_final.dat', &
                 add=1e-20, mul=mulfak, iColumn=1, dump=.true.)
            call WriteHistMP(hMP_pTSet(iSet), &
                 file='pT_Set'//achar(48+iSet)//'_final.dat', &
                 add=1e-20, mul=mulfak, iColumn=1, dump=.true.)
            call WriteHistMP(hMP_ySet(iSet), &
                 file='y_Set'//achar(48+iSet)//'_final.dat', &
                 add=1e-20, mul=mulfak, iColumn=1, dump=.true.)

            call ClearHistMP(hMP_ESet(iSet))
            call ClearHistMP(hMP_mTSet(iSet))
            call ClearHistMP(hMP_pTSet(iSet))
            call ClearHistMP(hMP_ySet(iSet))

         end do
      end if

    end subroutine Spectra


    !**************************************************************************
    !***is* DoHeavyIonAnalysis/analyze_Npart
    !**************************************************************************
    subroutine analyze_Npart

      integer :: i,j,nUnwounded
      type(particle), POINTER :: pPart


      nUnwounded = 0
      do i = 1,nEns
         do j = 1,size(realParts,dim=2)

            pPart => realParts(i,j)
            if (pPart%ID == EOV) exit
!            if (pPart%ID == NOP) cycle
            if (pPart%ID /= 1) cycle
            if (pPart%anti) cycle

            if (pPart%firstEvent /= 0) cycle

            nUnwounded = nUnwounded+1
         end do
      end do

      open(123,file="nUnwounded.txt",status="unknown")
      write(123,*) stossParameter,float(nUnwounded)/nEns
      close(123)

    end subroutine analyze_Npart


  end subroutine DoHeavyIonAnalysis


  !****************************************************************************
  !****s* HeavyIonAnalysis/DoHeavyIonAnalysisTime
  ! NAME
  ! subroutine DoHeavyIonAnalysisTime(realParts, time)
  !
  ! PURPOSE
  ! Do some analysis for every time step
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParts --- real particle vector
  ! * real :: time --- actual time
  !****************************************************************************
  subroutine DoHeavyIonAnalysisTime(realParts, time)
    use particleDefinition
    use inputGeneral, only: time_max, delta_T, timeForOutput, timeSequence, &
         num_runs_sameEnergy
    use output, only: intTochar,intToChar4
    use IdTable, only: EOV, NOP
    use Dilepton_Analysis, only: Dilep_write_CS_time

    type(particle), dimension(:,:), intent(in) :: realParts
    real, intent(in) :: time

    integer, save    :: isut=0      ! number of subsequent runs
    integer, save    :: itime=0

    if (time==0.) itime=0

    if (time>= timeForOutput) then

      if (mod(itime*delta_T,timeSequence)<1E-4) then
        if (flag_outputReal) call writeParticleVectorHI
        call RapiditySpectra
        call InvariantMassSpectra
        call Dilep_write_CS_time(time)
      end if

      itime = itime + 1

    end if

    if (abs(time-time_max) < 1E-4) isut = isut + 1

  contains

    !**************************************************************************
    !****is* DoHeavyIonAnalysisTime/writeParticleVectorHI
    ! NAME
    ! subroutine writeParticleVectorHI()
    !
    ! PURPOSE
    ! write out the particle vectors in the format used for HI analysis
    !**************************************************************************
    subroutine writeParticleVectorHI
      use minkowski, only: abs4
      use inputGeneral, only: eventtype
      use eventtypes, only: HeavyIon, hadron
      use initHeavyIon, only: b_HI => b
      use initHadron, only: b_had => b

      integer :: i, j, fact, indFree
      real    :: stossParameter

      !************************************************************************
      !****o* HeavyIonAnalysis/DoHIATime___.dat
      ! NAME
      ! file DoHIATime___.dat
      ! PURPOSE
      ! Contains the full dump of the real particle vector at a certain time
      ! step.
      ! The time is given in the file name and in the first line of the file.
      ! The columns have the following meaning:
      ! * 1    = particle ID
      ! * 2    = charge
      ! * 3    = vac. mass in GeV
      ! * 4    = (eff. mass)**2 in GeV^2
      ! * 5    = isFree (particle is bound or free?)
      ! * 6-8  = position (x,y,z) in fm
      ! * 9-11 = momentum (px,py,pz) in GeV
      ! * 12   = ensemble no.
      ! * 13   = run no.
      ! * 14   = impact parameter in fm
      !************************************************************************
      open(103,file='DoHIATime'//intToChar(nint(time/delta_T))//'.dat',position='Append')
      if (isut==0) then
         write(103,*) '# time = ',time,' fm/c'
         write(103,*)
      end if

      if (eventtype==HeavyIon) then
         stossParameter = b_HI
      else if (eventtype==hadron) then
         stossParameter = b_had
      else
         write(*,*) ' Problem in DoHeavyIonAnalysisTime: impact parameter is not defined'
         stossParameter=0.
      end if

      Loop_over_ensembles: do i = 1,size(realParts,dim=1)
         Loop_over_particles : do j = 1,size(realParts,dim=2)
            if (realParts(i,j)%id == NOP) cycle Loop_over_particles
            if (realParts(i,j)%id == EOV) exit Loop_over_particles
            if (realParts(i,j)%ID > 0) then
               if (realParts(i,j)%anti) then
                  fact = -1
               else
                  fact = 1
               end if
               if (isFree(realParts,i,j,realParts(i,j))) then
                  indFree=1
               else
                  indFree=0
               end if
               write(103,50) fact*realParts(i,j)%ID, realParts(i,j)%charge, &
                    realParts(i,j)%mass, abs4(realParts(i,j)%mom)**2, &
                    indFree, realParts(i,j)%pos(1:3), &
                    realParts(i,j)%mom(1:3), &
                    i, isut+1, stossParameter
            end if
         end do Loop_over_particles
      end do Loop_over_ensembles
      close(103)

50    format(2i4,2f9.4,i2,6f9.3,2i6,f6.2)

    end subroutine writeParticleVectorHI


    !**************************************************************************
    !****is* DoHeavyIonAnalysisTime/RapiditySpectra
    ! NAME
    ! subroutine RapiditySpectra()
    !
    ! PURPOSE
    ! generate and write out some rapidity spectra
    !**************************************************************************
    subroutine RapiditySpectra
      use idTable, only: nucleon
      use hist

      integer, parameter :: nrap = 60
      real, parameter :: ystart=-6.0, dy=0.1

      type(histogram), allocatable, save :: dNdy(:)
      logical, save :: RapidityInitFLAG=.true.

      integer :: i, j, k, N, nEns
      real :: yb, p0, pz
      character(len=100) :: title

      if (RapidityInitFLAG) then
        N = int((time_max - timeForOutput)/timeSequence)
        ! print *,"RapiditySpectra: allocating histograms: ", N, time_max, timeForOutput, timeSequence
        allocate(dNdy(0:N))
        do k=0,N
          write(title,'(A,f7.2,A)') "nucleon rapidity distribution dN/dy at time = ", timeForOutput + k*timeSequence, " fm/c"
          call createHist(dNdy(k),  title, ystart, ystart+2*nrap*dy, dy)
        end do
        RapidityInitFLAG = .false.
      end if

      k = nint((time - timeForOutput) / timeSequence)

      nEns = size(realParts,dim=1)

      Loop_over_ensembles: do i = 1,nEns
         Loop_over_particles : do j = 1,size(realParts,dim=2)
            if (realParts(i,j)%id == NOP) cycle Loop_over_particles
            if (realParts(i,j)%id == EOV) exit Loop_over_particles

            if (realParts(i,j)%ID == nucleon) then

               p0 = realParts(i,j)%mom(0)
               pz = realParts(i,j)%mom(3)
               yb = 0.5*log( (p0+pz)/(p0-pz) )

               call addHist(dNdy(k), yb, 1.0)
            end if
         end do Loop_over_particles
      end do Loop_over_ensembles

      call writeHist(dNdy(k),  file='RapidityDistributions'//intTochar4(nint(time/delta_T))//'.dat', mul = 1./nEns)

    end subroutine RapiditySpectra


    !**************************************************************************
    !****is* DoHeavyIonAnalysisTime/InvariantMassSpectra
    ! NAME
    ! subroutine InvariantMassSpectra()
    !
    ! PURPOSE
    ! generate and write out invariant mass spectra of some particles
    !**************************************************************************
    subroutine InvariantMassSpectra
      use idTable, only: nucleon,S11_1535,rho,omegaMeson
      use hist

      type(histogram), allocatable, save :: dNdm_nuc(:),dNdm_1535(:), &
           dNdm_rho(:),dNdm_omega(:)
      logical, save :: InitFLAG=.true.

      integer :: i, j, k, N, nEns
      real :: mstar
      character(len=100) :: title

      if (InitFLAG) then
         N = int((time_max - timeForOutput)/timeSequence)
         allocate(dNdm_nuc(0:N),dNdm_1535(0:N),dNdm_rho(0:N),dNdm_omega(0:N))
         do k=0,N
            write(title,'(A,f7.2,A)') &
                 "nucleon invariant mass distribution dN/dm at time = ", &
                 timeForOutput + k*timeSequence, " fm/c"
            call createHist(dNdm_nuc(k),  title, 0.4, 0.4+80*0.01, 0.01)
            write(title,'(A,f7.2,A)') &
                 "N*(1535) invariant mass distribution dN/dm at time = ", &
                 timeForOutput + k*timeSequence, " fm/c"
            call createHist(dNdm_1535(k),  title, 0., 0.+80*0.05, 0.05)
            write(title,'(A,f7.2,A)') &
                 "rho-meson invariant mass distribution dN/dm at time = ", &
                 timeForOutput + k*timeSequence, " fm/c"
                 call createHist(dNdm_rho(k),  title, 0., 0.+50*0.02, 0.02)
            write(title,'(A,f7.2,A)') &
                 "omega-meson invariant mass distribution dN/dm at time = ", &
                 timeForOutput + k*timeSequence, " fm/c"
            call createHist(dNdm_omega(k),  title, 0., 0.+200*0.01, 0.01)
         end do
         InitFLAG = .false.
      end if

      k = nint((time - timeForOutput) / timeSequence)

      nEns = size(realParts,dim=1)

      Loop_over_ensembles: do i = 1,nEns
         Loop_over_particles : do j = 1,size(realParts,dim=2)
            if (realParts(i,j)%id == NOP) cycle Loop_over_particles
            if (realParts(i,j)%id == EOV) exit Loop_over_particles
            select case(realParts(i,j)%ID)
            case (nucleon)
               mstar = sqrtS(realParts(i,j),'InvariantMassSpectra, mstar_nuc')
               call addHist(dNdm_nuc(k),mstar,1.0)
            case(S11_1535)
               mstar = sqrtS(realParts(i,j),'InvariantMassSpectra, mstar_1535')
               call addHist(dNdm_1535(k),mstar,1.0)
            case(rho)
               call addHist(dNdm_rho(k),realParts(i,j)%mass,1.0)
            case(omegaMeson)
               call addHist(dNdm_omega(k),realParts(i,j)%mass,1.0)
            end select
         end do Loop_over_particles
      end do Loop_over_ensembles

      call writeHist(dNdm_nuc(k), &
           file='InvariantMassDistributions_nuc_'//intTochar4(nint(time/delta_T))//'.dat',&
           mul = 1./(nEns*num_runs_sameEnergy))
      call writeHist(dNdm_1535(k), &
           file='InvariantMassDistributions_1535_'//intTochar4(nint(time/delta_T))//'.dat',&
           mul = 1./(nEns*num_runs_sameEnergy))
      call writeHist(dNdm_rho(k), &
           file='InvariantMassDistributions_rho_'//intTochar4(nint(time/delta_T))//'.dat',&
           mul = 1./(nEns*num_runs_sameEnergy))
      call writeHist(dNdm_omega(k), &
           file='InvariantMassDistributions_omega_'//intTochar4(nint(time/delta_T))//'.dat',&
           mul = 1./(nEns*num_runs_sameEnergy))

    end subroutine InvariantMassSpectra


  end subroutine DoHeavyIonAnalysisTime

  !****************************************************************************
  !****s* HeavyIonAnalysis/HeavyIon_evol
  ! NAME
  ! subroutine HeavyIon_evol(realParts, time, timestep)
  !
  ! PURPOSE
  ! time evolution of some global observables
  !
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realParts --- real particle vector
  ! * real :: time --- actual time
  ! * integer :: timestep --- number of time step
  !****************************************************************************
  subroutine HeavyIon_evol(realParts, time, timestep)

    use IdTable
    use particleDefinition
    use densitymodule
    use RMF, only: getRMF_flag, Tens_flag, g_omega, g_rho, g_sigma, &
         ModificationFactor
    use particleProperties, only: hadron
    use output, only: realTochar,intToChar,intToChar4
    use inputGeneral, only: delta_T, eventtype
    use constants, only: pi, hbarc, rhoNull, mN
    use coulomb, only: emfoca
    use PauliBlockingModule, only: pauliBlocking
    use initHadron, only: b,z,p_lab,E_bind
    use collisionNumbering, only: GetCountedEvents
    use eventtypes, only: ET_hadron => hadron, HiLepton
    use thermoDyn, only: temperatureAt, muAt
    use hist, only: histogram, createHist, addHist, writeHist

    type(particle), dimension(:,:), intent(in), target  :: realParts
    real, intent(in) :: time
    integer, intent(in) :: timestep

    type(particle), POINTER :: pPart

    real, dimension(-121:121,-2:2) :: parnum, parnum_free   ! particle numbers
    integer :: nEns,i,j,k,id,charge,Index1,Index2,Index3,npart
    integer, save :: icall=0     ! counter of calls
    integer, save :: isut=0      ! number of subsequent run
    real :: rhobar_points,Qzz,rhorad,rhorad_n,rhorad_p,r1,r2,r, &
         velrad1,velrad2, &
         rhoz,rhoz_bar,rhoz_antibar,rho_bar_max,rho_antibar_max,rholrf, &
         endens,rhobar_local,m_inv
    real :: rhoz_Bar_points,rhoz_AntiBar_points,rhoz_BoundBar_points, &
         rhoz_TargetBar_points,&
         &energy,fnorm
    real :: mstar, pf, sigma_rad, factor
    real, dimension(1:3) :: p2_aver

    ! Check conservation laws:
    real, dimension(0:3) :: p ! total 4-momentum of all particles
    real :: baryon_number,charge_number,strangeness

    ! Properties of the bound system:
    real, parameter :: rho_min_bound=0.01*rhoNull ! minimal density
    real :: B_bound, N_bound, Z_bound ! baryon, neutron and charge numbers
    real :: rho_bound                ! mass averaged density
    real :: E_kinColl_bound          ! collective kinetic energy
    real, dimension(0:3) :: P_bound  ! total 4-momentum (w/o Coulomb contribution in P_bound(0))
    real :: E_Coul_bound             ! Coulomb energy

    real :: cpot
    real, dimension(1:3) :: place,impuls,velColl
    real, dimension(0:3) :: momentum

    logical :: flag_output, blockFlag

    real :: rhoz_n, rhoz_p, Ef_p, Ef_n, pf_p, pf_n, temp, temp_max, mub
    real :: Ecoul   ! EcoulNorm
    real :: pion_ratio

    real, dimension(1:3)     :: number_of_baryons, Qzz_aver, r2_aver, r_aver_sq
    real, dimension(1:3)     :: number_of_neutrons,number_of_protons,np_ratio
    real, dimension(1:3,1:3) :: r_aver
    real, dimension(1:6) :: probability

    real :: mulfak

    type(histogram) :: histMnucleon
    integer :: nNucleons,npoints

    real :: sigmaV
    real, dimension(0:3) :: omegaV, rhoV

    if (initFlag) call init

    if (icall==0) then
       open(32,file='evol.dat')
       write(32,*)'# Column No., quantity:'
       write(32,*)'# 1  time'
       write(32,*)'# 2  nucleon multiplicity'
       write(32,*)'# 3  delta'
       write(32,*)'# 4  S11_1535'
       write(32,*)'# 5  D13_1520'
       write(32,*)'# 6  Other Nonstrange Baryonic Resonances'
       write(32,*)'# 7  pi'
       write(32,*)'# 8  eta'
       write(32,*)'# 9  rho'
       write(32,*)'# 10 K'
       write(32,*)'# 11 Kbar'
       write(32,*)'# 12 Lambda'
       write(32,*)'# 13 Lambda free'
       write(32,*)'# 14 Sigma'
       write(32,*)'# 15 Sigma free'
       write(32,*)'# 16 Xi'
       write(32,*)'# 17 Xi free'
       write(32,*)'# 18 XiStar'
       write(32,*)'# 19 XiStar free'
       write(32,*)'# 20 K^*'
       write(32,*)'# 21 Kbar^*'
       write(32,*)'# 22 Y^* (S=-1 only)'
       write(32,*)'# 23 LambdaBar'
       write(32,*)'# 24 SigmaBar'
       write(32,*)'# 25 Ybar^* (S=1 only)'
       write(32,*)'# 26 XiBar'
       write(32,*)'# 27 XiBar^*'
       write(32,*)'# 28 J/Psi'
       write(32,*)'# 29 D'
       write(32,*)'# 30 Dbar'
       write(32,*)'# 31 D^*'
       write(32,*)'# 32 Dbar^*'

       open(33,file='dens.dat',status='unknown')
       write(33,'(A)') '# time rhocen_baryons rhocen_antibaryons ptot '// &
            'baryon_number charge strangeness '// &
            'temperature(central) mu_B(central) Qzz(central)'

       open(40,file='conserv.dat')
       write(40,*)'# time: energy:  px:  py:  pz:  baryon number:  '// &
            'charge: strangeness:'

       if (getRMF_flag() .and. Tens_flag ) then
          open(42,file='Tensor.dat')
          write(42,*)'# time:   T^00:    T^11:   T^22:   T^33: (GeV/fm^3)'
          open(43,file='BoundSystem.dat')
          write(43,*)'# Properties of the bound system:'
          write(43,*)'# time:   B:   N:   Z:  P(0:3), AGeV:   '// &
               'E_kinColl, AGeV:  E_Coul, AGeV:  <rho>, fm ^-3:'
       end if

       open(44,file='collisions.dat')
       write(44,*)'# time N_2body(integrated) N_2body(timestep)'

       open(50,file='isospinRatios.dat')

    end if

    icall = icall + 1
    if ( abs(time-delta_T) < 1.e-06 ) isut=isut+1
    flag_output=.false.

    parnum= 0.
    parnum_free= 0.
    rhobar_points= 0.
    endens= 0.
    p2_aver= 0.
    r2_aver=0.
    r_aver=0.0
    Qzz_aver=0.
    Qzz=0.
    p(:)=0.
    baryon_number=0.
    charge_number=0.
    strangeness=0.
    number_of_baryons=0.
    number_of_neutrons=0.
    number_of_protons=0.
    E_Coul_bound=0.

    call createHist(histMnucleon, 'dN/dM nucleon', 0.8,0.95,0.002)
    nNucleons = 0

    nEns = size(realParts,dim=1)
    ensemble_loop : do i = 1,nEns
       particle_loop : do j = 1,size(realParts,dim=2)

          pPart => realParts(i,j)

          if (pPart%ID == EOV) exit particle_loop
          if (pPart%ID == NOP) cycle particle_loop

          if (NucleonMassPlot) then
             if (pPart%ID == 1 .and. .not.pPart%anti) then
                call addHist(histMnucleon, sqrtS(pPart), 1.0)
                nNucleons = nNucleons+1
             end if
          end if

          factor = merge(-1., 1., pPart%anti)
          id = pPart%ID*int(factor)
          charge = pPart%charge
          if (abs(id).le.121) then
             parnum(id,charge) = parnum(id,charge) + 1.
             if (id.eq.Lambda .or. id.eq.SigmaResonance .or. id.eq.Xi .or. id.eq.XiStar) then
                if (IsFree(realParts,i,j,pPart)) parnum_free(id,charge) = parnum_free(id,charge) + 1.
             end if
          end if
          p(:)=p(:)+pPart%mom(:)
          if (isBaryon(pPart%ID)) then
             baryon_number=baryon_number+factor
             strangeness=strangeness+real(hadron(pPart%ID)%strangeness)*factor
          else if (isMeson(pPart%ID)) then
             strangeness=strangeness+real(hadron(pPart%ID)%strangeness)
          end if
          charge_number=charge_number+real(pPart%charge)

          if (all(abs(pPart%pos) < 0.5)) then
             !       Central baryon density and Qzz in momentum space:
             if (isBaryon(pPart%ID) .and. .not.pPart%anti) then
                rhobar_points=  rhobar_points + 1.
                Qzz=Qzz+2.*pPart%mom(3)**2-pPart%mom(1)**2-pPart%mom(2)**2
             end if
             !       Central energy density:
             endens= endens + pPart%mom(0)
          end if

          if (isBaryon(pPart%ID)) then

             !       Average squares of momentum components :
             p2_aver(1:3)= p2_aver(1:3) + pPart%mom(1:3)**2

             if (.not.pPart%anti) then

                !         Indexes in large grid:
                Index1=NINT(pPart%pos(1)/gridSpacing(1))
                Index2=NINT(pPart%pos(2)/gridSpacing(2))
                Index3=NINT(pPart%pos(3)/gridSpacing(3))

                if (        abs(Index1).le.gridPoints(1) &
                     & .and. abs(Index2).le.gridPoints(2) &
                     & .and. abs(Index3).le.gridPoints(3)  ) then

                   if (getRMF_flag()) then
                      factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants
                      rhobar_local=(   densField(Index1,Index2,Index3)%baryon(0)  &
                           & + factor*totalDens(Index1,Index2,Index3) )/(1.+factor)  ! density of the baryons
                   else
                      ! w/o RMF "totalDens" is actually (baryon-antibaryon) density and
                      ! densField(0,0,0)%baryon(0) is (baryon+antibaryon) density (see density.f90: addTodensField)
                      rhobar_local=(   densField(Index1,Index2,Index3)%baryon(0)  &
                           & + totalDens(Index1,Index2,Index3) )/2.                 ! density of the baryons
                   end if

                   !cut on density (0.1*rho_sat)
                   if (rhobar_local.gt.0.1*rhoNull) then
                      number_of_baryons(1)=number_of_baryons(1)+1.
                      r2_aver(1)= r2_aver(1) &
                           & + dot_product(pPart%pos,pPart%pos)
                      r_aver(1,:)= r_aver(1,:) + pPart%pos(:)
                      Qzz_aver(1)=Qzz_aver(1) + 2.*pPart%pos(3)**2 &
                           &- pPart%pos(1)**2 - pPart%pos(2)**2
                      if (pPart%ID==1 .and. pPart%charge==0) &
                           & number_of_neutrons(1)=number_of_neutrons(1)+1.
                      if (pPart%ID==1 .and. pPart%charge==1) &
                           & number_of_protons(1)=number_of_protons(1)+1.
                   end if

                   !cut on density (0.01*rho_sat)
                   if (rhobar_local.gt.0.01*rhoNull) then
                      number_of_baryons(2)=number_of_baryons(2)+1.
                      r2_aver(2) = r2_aver(2) &
                           & + dot_product(pPart%pos,pPart%pos)
                      r_aver(2,:)= r_aver(2,:) + pPart%pos(:)
                      Qzz_aver(2)=Qzz_aver(2) + 2.*pPart%pos(3)**2 &
                           &- pPart%pos(1)**2 - pPart%pos(2)**2
                      if (pPart%ID==1 .and. pPart%charge==0) &
                           & number_of_neutrons(2)=number_of_neutrons(2)+1.
                      if (pPart%ID==1 .and. pPart%charge==1) &
                           & number_of_protons(2)=number_of_protons(2)+1.
                   end if

                   !cut on binding energy (E<0)
                   if (getRMF_flag()) then
                      momentum = Particle4MomentumRMF(pPart)
                      energy=momentum(0)
                   else
                      energy=pPart%mom(0)
                   end if
                   place(1:3)= pPart%pos(1:3)
                   impuls(1:3)= pPart%mom(1:3)
                   energy = energy + emfoca(place,impuls,pPart%charge,pPart%ID)
                   if ( energy - pPart%mass < 0.) then
                      number_of_baryons(3)=number_of_baryons(3)+1.
                      r2_aver(3) = r2_aver(3) &
                           & + dot_product(pPart%pos,pPart%pos)
                      r_aver(3,:)= r_aver(3,:) + pPart%pos(:)
                      if (pPart%ID==1 .and. pPart%charge==0) &
                           & number_of_neutrons(3)=number_of_neutrons(3)+1.
                      if (pPart%ID==1 .and. pPart%charge==1) &
                           & number_of_protons(3)=number_of_protons(3)+1.
                   end if

                   if (rhobar_local.gt.rho_min_bound) then
                      place(1:3)= pPart%pos(1:3)
                      impuls(1:3)= pPart%mom(1:3)
                      E_Coul_bound = E_Coul_bound + 0.5*emfoca(place,impuls,pPart%charge,pPart%ID) !see notes in evaluateTotal4Momentum_RMF
                   end if

                end if

             end if

          end if

       end do particle_loop
    end do ensemble_loop

    mulfak = 1.0/float(nEns)

    if (DensityPlot) then
       if (mod(timestep,5)==0) then
          call writeDensityPlane('density_YZ_'//intTochar4(timestep)//'.dat',1)
          call writeDensityPlane('density_XZ_'//intTochar4(timestep)//'.dat',2)
          call writeDensityPlane('density_XY_'//intTochar4(timestep)//'.dat',3)
       end if
    end if

    if (NucleonMassPlot) then
       if (mod(timestep,5)==0) then
          if (nNucleons>0) then
             call WriteHist(histMnucleon, &
                  file='Mnucleon_'//intTochar(timestep)//'.dat', &
                  mul=1.0/nNucleons)
          end if
       end if
    end if

    if (do_QRvector) then
       if (mod(timestep,5)==0) then
          call doQRvector(realParts,timestep)
       end if
    end if

    if (do_Tmunu) then
       if (mod(timestep,5)==0) then
          call doTmunu(realParts,timestep)
       end if
    end if

    if (eventtype.ne.HiLepton) then
       call countParts(realParts, time, "")
    end if

    parnum= parnum*mulfak
    parnum_free= parnum_free*mulfak
    if(rhobar_points.gt.0.) Qzz=Qzz/rhobar_points
    rhobar_points= rhobar_points*mulfak/1.**3
    endens= endens*mulfak/1.**3
    p2_aver(:)= p2_aver(:)*mulfak
    where (number_of_baryons /=0)
       r2_aver(:)= r2_aver(:)/number_of_baryons(:)
       Qzz_aver(:)=Qzz_aver(:)/number_of_baryons(:)
    end where
    do i=1,3
       if (number_of_baryons(i)/=0) &
            r_aver(i,:)  = r_aver(i,:)/number_of_baryons(i)
       r_aver_sq(i) = sqrt(dot_product(r_aver(i,:),r_aver(i,:)))
    end do
    p(:)=p(:)*mulfak
    baryon_number=baryon_number*mulfak
    charge_number=charge_number*mulfak
    strangeness=strangeness*mulfak
    number_of_baryons(:)=number_of_baryons(:)*mulfak
    number_of_protons(:)=number_of_protons(:)*mulfak
    number_of_neutrons(:)=number_of_neutrons(:)*mulfak
    E_Coul_bound=E_Coul_bound*mulfak

    !n/p-ratio:
    do i=1,3
       if ( number_of_protons(i) .ne. 0.0 ) then
          np_ratio(i) = number_of_neutrons(i) / number_of_protons(i)
       else
          np_ratio(i) = 0.0
       end if
    end do

    !pion-ratio pi^{-}/pi^{+}:
    if ( parnum(pion,+1).ne.0.0 ) then
       pion_ratio = parnum(pion,-1)/parnum(pion,+1)
    else
       pion_ratio = 0.0
    end if

    write(32,5) time,sum(parnum(nucleon,:)),sum(parnum(delta,:)),&
         sum(parnum(S11_1535,:)),sum(parnum(D13_1520,:)),&
         sum(parnum(P11_1440,:))+sum(parnum(S11_1650,:))+sum(parnum(S11_2090,:))+sum(parnum(D13_1700:F37_1950,:)),&
         sum(parnum(pion,:)),sum(parnum(eta,:)),sum(parnum(rho,:)),sum(parnum(kaon,:)),sum(parnum(kaonBar,:)),&
         sum(parnum(Lambda,:)),sum(parnum_free(Lambda,:)),&
         sum(parnum(SigmaResonance,:)),sum(parnum_free(SigmaResonance,:)),&
         sum(parnum(Xi,:)),sum(parnum_free(Xi,:)),&
         sum(parnum(XiStar,:)),sum(parnum_free(XiStar,:)),&
         sum(parnum(kaonStar,:)),sum(parnum(kaonStarBar,:)),&
         sum(parnum(Sigma_1385:Sigma_1915,:)),&
         sum(parnum(-Lambda,:)),sum(parnum(-SigmaResonance,:)),&
         sum(parnum(-Sigma_1915:-Sigma_1385,:)),&
         sum(parnum(-Xi,:)),sum(parnum(-XiStar,:)),&
         sum(parnum(JPsi,:)),sum(parnum(dMeson,:)),sum(parnum(dBar,:)),&
         sum(parnum(dStar,:)),sum(parnum(dStarBar,:))

    write(50,5) time,(np_ratio(i),i=1,3),pion_ratio

5   format(1P,50(1x,e13.6))

    rhoz_bar=0.
    rhoz_antibar=0.
    temp=0.
    npoints=0
!    npoints=2
    do k = -npoints,npoints
       do j = -npoints,npoints
          do i = -npoints,npoints
             if (getRMF_flag()) then
                factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants
                rhoz_bar = rhoz_bar + (densField(i,j,k)%baryon(0)+factor*totalDens(i,j,k))/(1.+factor)  ! density of the baryons
                rhoz_antibar = rhoz_antibar + (totalDens(i,j,k)-densField(i,j,k)%baryon(0))/(1.+factor)  ! density of the antibaryons
             else
                ! w/o RMF "totalDens" is actually (baryon-antibaryon) density and
                ! densField(i,j,k)%baryon(0) is (baryon+antibaryon) density (see density.f90: addTodensField)
                rhoz_bar = rhoz_bar + (densField(i,j,k)%baryon(0)+totalDens(i,j,k))/2.
                rhoz_antibar =rhoz_antibar + (densField(i,j,k)%baryon(0)-totalDens(i,j,k))/2.
             end if
             temp = temp + temperatureAt((/float(i)*gridSpacing(1),float(j)*gridSpacing(2),float(k)*gridSpacing(3)/))
          end do
       end do
    end do
    rhoz_bar= rhoz_bar/(2*npoints+1)**3
    rhoz_antibar= rhoz_antibar/(2*npoints+1)**3
    temp=temp/(2*npoints+1)**3

    mub = muAt (rhoz_bar, temp)
    write(33,5) time, rhoz_bar, rhoz_antibar, p(:), baryon_number, charge_number, strangeness, temp, mub, Qzz

    write(40,5) time,(r_aver_sq(i), sqrt(r2_aver(i)), number_of_baryons(i), i=1,3)

    if (getRMF_flag() .and. Tens_flag) then
       write(42,5) time, (Tens(0,0,0,i,i),i=0,3)
       P_bound(:)=0.
       E_kinColl_bound=0.
       B_bound=0.
       N_bound=0.
       Z_bound=0.
       rho_bound=0.
       do Index1=-gridpoints(1),gridpoints(1)
          do Index2=-gridPoints(2),gridPoints(2)
             do Index3=-gridPoints(3),gridPoints(3)
                rhobar_local=densField(Index1,Index2,Index3)%baryon(0)
                if (rhobar_local.gt.rho_min_bound) then
                   P_bound(0:3)=P_bound(0:3)+mom4Dens(Index1,Index2,Index3,0:3)
                   m_inv=mom4Dens(Index1,Index2,Index3,0)**2 &
                        &-dot_product(mom4Dens(Index1,Index2,Index3,1:3),&
                        &mom4Dens(Index1,Index2,Index3,1:3))
                   m_inv=sqrt(max(0.,m_inv))
                   E_kinColl_bound=E_kinColl_bound+mom4Dens(Index1,Index2,Index3,0)&
                        &-m_inv
                   B_bound=B_bound+rhobar_local
                   N_bound=N_bound+densField(Index1,Index2,Index3)%neutron(0)
                   Z_bound=Z_bound+densField(Index1,Index2,Index3)%proton(0)
                   rho_bound=rho_bound+rhobar_local**2
                end if
             end do
          end do
       end do
       P_bound(0:3)=P_bound(0:3)*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       E_kinColl_bound=E_kinColl_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       B_bound=B_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       N_bound=N_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       Z_bound=Z_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
       rho_bound=rho_bound*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)/B_bound

       write(43,5) time, B_bound, N_bound, Z_bound, P_bound(:)/B_bound, E_kinColl_bound/B_bound,&
            &E_Coul_bound/B_bound, rho_bound
    end if

    write(44,5) time, float(getCountedEvents(0,2,1))/float(nEns), float(getCountedEvents(1,2,1))/float(nEns)

    ! flush all files, so that output is actually written to disk
    flush(32)
    flush(33)
    flush(40)
    flush(44)
    flush(50)

    !--------------------------------------------------------------------------
    ! more detailed output, if flag_outputDetailed=.true.
    !--------------------------------------------------------------------------
    if ( .not. flag_outputDetailed ) return
    !--------------------------------------------------------------------------

    if ( abs(time-delta_T) < 1.e-06 .or. abs(time-nint(time)) < 0.5*delta_T ) then
       flag_output=.true.
       open(34,file='rhorad_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(34,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(34,'(A,f10.4)')'# time, fm/c: ', time
       write(34,'(A,i3)')'# run number: ', isut
       open(35,file='rhoz_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(35,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(35,'(A,f10.4)')'# time, fm/c: ', time
       write(35,'(A,i3)')'# run number: ', isut
       open(38,file='rhozx_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(38,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(38,'(A,f10.4)')'# time, fm/c: ', time
       write(38,'(A,i3)')'# run number: ', isut
       if (getRMF_flag()) then
          open(36,file='Fields_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
          write(36,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
          write(36,'(A,f10.4)')'# time, fm/c: ', time
          write(36,'(A,i3)')'# run number: ', isut
       end if
       open(37,file='pauli_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(37,'(A,f10.4)')'# time, fm/c: ', time
       write(37,'(A,i3)')'# run number: ', isut
       write(37,*)'# Momentum, GeV/c:    PauliBlock:'
       open(41,file='temp_'//realTochar(time+1.e-06)//'_'//intTochar(isut)//'.dat',position='append')
       write(41,'(A,3(f10.4,2x))')'# gridSpacing, fm:', gridSpacing
       write(41,'(A,f10.4)')'# time, fm/c: ', time
       write(41,'(A,i3)')'# run number: ', isut
    end if

    if (icall == 1) then
       open(39,file='dens_max.dat',status='unknown')
       if (eventtype == ET_hadron) &
            write(39,'(A,4(e15.8,1x))')" # hadron's b, z, p_lab, E_bind: ", b, z, p_lab, E_bind
       write(39,*)'# time  rho_bar_max  rho_antibar_max'
    end if


    Final_Output : if (flag_output) then

       write(34,*) '# r, fm:   rhorad:  rhorad_n:   rhorad_p:  velrad1/c:   velrad2/c:  PauliBlock(1-6):'
       do k= 0,min(gridpoints(1),gridpoints(2),gridpoints(3))

          r1=max(0.,(float(k)-0.5)*gridSpacing(1))
          r2=(float(k)+0.5)*gridSpacing(1)

          npart=0
          rhorad=0.
          rhorad_n=0.
          rhorad_p=0.
          velrad1=0.
          velrad2=0.
          sigma_rad=0.

          do i = 1,nEns
             particle_loop1 : do j = 1,size(realParts,dim=2)

                pPart => realParts(i,j)

                if (pPart%ID == EOV) exit particle_loop1
                if (pPart%ID == NOP) cycle particle_loop1

                place(1:3)= pPart%pos(1:3)
                r=sqrt(dot_product(place(1:3),place(1:3)))

                if ( r.gt.r1 .and. r.le.r2 ) then

                   npart=npart+1

                   !             Indexes in large grid:
                   Index1=NINT(place(1)/gridSpacing(1))
                   Index2=NINT(place(2)/gridSpacing(2))
                   Index3=NINT(place(3)/gridSpacing(3))

                   if (        abs(Index1).le.gridPoints(1) &
                        & .and. abs(Index2).le.gridPoints(2) &
                        & .and. abs(Index3).le.gridPoints(3)  ) then

                      rhobar_local=densField(Index1,Index2,Index3)%baryon(0)
                      rhorad=rhorad+rhobar_local

                      rhorad_n=rhorad_n+densField(Index1,Index2,Index3)%neutron(0)
                      rhorad_p=rhorad_p+densField(Index1,Index2,Index3)%proton(0)

                      call getFieldRMF(Index1,Index2,Index3, sigmaV)
                      sigma_rad=sigma_rad+sigmaV

                      if (rhobar_local.gt.1.e-06) then
                         velColl(1:3)=densField(Index1,Index2,Index3)%baryon(1:3) &
                              &/rhobar_local
                         if (r.gt.1.e-06) velrad1=velrad1+dot_product(velColl(1:3),place(1:3))/r
                      end if

                      if (allocated(mom4Dens)) then
                         if (mom4Dens(Index1,Index2,Index3,0).gt.1.e-06) then
                            velColl(1:3)=mom4Dens(Index1,Index2,Index3,1:3) &
                                 &/mom4Dens(Index1,Index2,Index3,0)
                            if (r.gt.1.e-06) velrad2=velrad2+dot_product(velColl(1:3),place(1:3))/r
                         end if
                      end if

                   end if

                end if


             end do particle_loop1
          end do

          if (npart.gt.0) then
             rhorad=rhorad/float(npart)
             rhorad_n=rhorad_n/float(npart)
             rhorad_p=rhorad_p/float(npart)
             sigma_rad=sigma_rad/float(npart)
             velrad1=velrad1/float(npart)
             velrad2=velrad2/float(npart)
          end if

          mstar = mN + g_sigma*sigma_rad

          !********** Pauli blocking check, radial dependence: ****************
          do i=1,6
             momentum(1:3)=0.
             momentum(0)=sqrt(mN**2+dot_product(momentum(1:3),momentum(1:3)))
             place(1:3)=0.
             if (i.eq.1) place(1)=float(k)*gridSpacing(1)
             if (i.eq.2) place(1)=-float(k)*gridSpacing(1)
             if (i.eq.3) place(2)=float(k)*gridSpacing(2)
             if (i.eq.4) place(2)=-float(k)*gridSpacing(2)
             if (i.eq.5) place(3)=float(k)*gridSpacing(3)
             if (i.eq.6) place(3)=-float(k)*gridSpacing(3)
             blockFlag=pauliBlocking(momentum,place,1,realParts,probability(i))
          end do
          !********************************************************************

          write(34,5) float(k)*gridSpacing(1), rhorad, rhorad_n, rhorad_p, &
               velrad1, velrad2, probability   !, mstar

       end do


       !********** Pauli blocking check, momentum dependence: *****************
       place(1:3)=0.
       do k=1,201
          do i=1,6
             momentum(1:3)=0.
             if (i.eq.1) momentum(1)=0.00240*float(k-1)
             if (i.eq.2) momentum(1)=-0.00240*float(k-1)
             if (i.eq.3) momentum(2)=0.00240*float(k-1)
             if (i.eq.4) momentum(2)=-0.00240*float(k-1)
             if (i.eq.5) momentum(3)=0.00240*float(k-1)
             if (i.eq.6) momentum(3)=-0.00240*float(k-1)
             momentum(0)=sqrt(mN**2+dot_product(momentum(1:3),momentum(1:3)))
             blockFlag=pauliBlocking(momentum,place,1,realParts,probability(i))
          end do
          write(37,5) 0.00240*float(k-1),probability
       end do
       !***********************************************************************

       if (getRMF_flag()) then
          write(35,'(2A)') &
               & '# z, fm:     rhoz_bar:   rhoz_antibar:',&
               &'rhoz_Bar_points:  rhoz_AntiBar_points:  rhoz_BoundBar_points: rhoz_TargetBar_points:'
       else
          write(35,'(2A)') &
               & '# z, fm:     rhoz_bar:   rhoz_antibar:   rholrf:',&
               &'rhoz_Bar_points:  rhoz_AntiBar_points:  rhoz_BoundBar_points: rhoz_TargetBar_points:'
       end if

       write(41,'(A)') '# z, fm:    temperature, GeV:'

       Loop_over_zGrid1 : do k= -gridpoints(3),gridpoints(3)

          rhoz= densField(0,0,k)%baryon(0)

          !proton & neutron densities
          rhoz_n=densField(0,0,k)%neutron(0)
          rhoz_p=densField(0,0,k)%proton(0)

          if (getRMF_flag()) then

             factor=ModificationFactor(nucleon,.true.)    ! Modification factor of the antinucleon coupling constants

             rhoz_bar=(rhoz+factor*totalDens(0,0,k))/(1.+factor)  ! density of the baryons

             rhoz_antibar=(totalDens(0,0,k)-rhoz)/(1.+factor)     ! density of the antibaryons

             pf=(1.5*pi**2*max(0.,rhoz_bar))**0.333333*hbarc

             pf_p=(3.*pi**2*max(0.,rhoz_p))**0.333333*hbarc
             pf_n=(3.*pi**2*max(0.,rhoz_n))**0.333333*hbarc

             mstar=DiracMass(0,0,k,mN,nucleon,1,.false.)

             call getFieldRMF(0,0,k, omega=omegaV)
             Ef_p=sqrt(pf_p**2+mstar**2)+g_omega*omegaV(0)
             Ef_n=sqrt(pf_n**2+mstar**2)+g_omega*omegaV(0)

          else

             ! Remember: w/o RMF totalDens is (baryon-antibaryon) density
             rhoz_bar=(rhoz+totalDens(0,0,k))/2.         ! density of the baryons
             rhoz_antibar=(rhoz-totalDens(0,0,k))/2.     ! density of the antibaryons

             rholrf= sqrt(max( densField(0,0,k)%baryon(0)**2 &
                  -densField(0,0,k)%baryon(1)**2 &
                  -densField(0,0,k)%baryon(2)**2 &
                  -densField(0,0,k)%baryon(3)**2, 0. ))

          end if

          rhoz_Bar_points=0.
          rhoz_AntiBar_points=0.
          rhoz_BoundBar_points=0.
          rhoz_TargetBar_points=0.
          !        ECoul = 0.0
          !        EcoulNorm = 0.0

          do i = 1,nEns
             particle_loop2 : do j = 1,size(realParts,dim=2)

                pPart => realParts(i,j)

                if (pPart%ID == EOV) exit particle_loop2
                if (pPart%ID == NOP) cycle particle_loop2

                if (pPart%ID >= pion) cycle particle_loop2   ! Count only (anti)baryons

                place(1:3)= pPart%pos(1:3)

                if ( abs(place(1)) <= 0.5*gridSpacing(1) .and. &
                     &abs(place(2)) <= 0.5*gridSpacing(2) .and. &
                     &abs(place(3)-float(k)*gridSpacing(3)) <= 0.5*gridSpacing(3) ) then

                   if ( .not.pPart%anti ) then
                      rhoz_Bar_points=rhoz_bar_points+1.
                      if (getRMF_flag()) then
                         momentum = Particle4MomentumRMF(pPart)
                         energy=momentum(0)
                      else
                         energy=pPart%mom(0)
                      end if
                      impuls(1:3)= pPart%mom(1:3)
                      cpot = emfoca(place,impuls,pPart%charge,pPart%ID)
                      energy=energy+cpot
                      !for testing E_F=E_F(r)
                      !                   if (pPart%charge==1) then
                      !                      Ecoul = Ecoul + cpot
                      !                      EcoulNorm = EcoulNorm + 1.
                      !                   end if
                      ! write(*,*)'Id, energy:', pPart%ID, energy
                      if (energy-pPart%mass < 0) rhoz_BoundBar_points=rhoz_BoundBar_points+1.
                      if (pPart%event(1)==1) &
                           & rhoz_TargetBar_points=rhoz_TargetBar_points+1.
                   else
                      rhoz_antibar_points=rhoz_antibar_points+1.
                   end if

                end if

             end do particle_loop2
          end do

          if (getRMF_flag()) then
             !Coulomb contribution to the Fermi-Energy:
             place=(/0.,0.,float(k)*gridSpacing(3)/)
             impuls=(/0.,0.,pf/)
             Ecoul = emfoca(place,impuls,1,1)
          end if

          fnorm=1./(float(nEns)*gridSpacing(1)*gridSpacing(2)*gridSpacing(3))
          rhoz_Bar_points=rhoz_Bar_points*fnorm
          rhoz_AntiBar_points=rhoz_AntiBar_points*fnorm
          rhoz_BoundBar_points=rhoz_BoundBar_points*fnorm
          rhoz_TargetBar_points=rhoz_TargetBar_points*fnorm

          if (getRMF_flag()) then

             write(35,5)  float(k)*gridSpacing(3), &
                  & rhoz_bar, rhoz_antibar, &
                  & rhoz_Bar_points,rhoz_AntiBar_points,&
                  & rhoz_BoundBar_points,rhoz_TargetBar_points

             call getFieldRMF(0,0,k, sigmaV, omegaV, rhoV)

             write(36,5) float(k)*gridSpacing(3), &
                  & pf_p, pf_n, rhoz_p, rhoz_n, &
                  & sigmaV, mstar, &
                  & g_omega*omegaV(0), &
                  & g_rho*rhoV(0),ECoul, &
                  & (Ef_p+Ecoul+g_rho*rhoV(0)-mN), &
                  & (Ef_n-g_rho*rhoV(0)-mN)
          else

             write(35,5)  float(k)*gridSpacing(3),rhoz_bar,rhoz_antibar,rholrf,rhoz_Bar_points,&
                  &rhoz_AntiBar_points,rhoz_BoundBar_points,rhoz_TargetBar_points

          end if

          write(41,5)  float(k)*gridSpacing(3),temperatureAt((/0.,0.,float(k)*gridSpacing(3)/))

       end do Loop_over_zGrid1

    end if Final_Output

    if (flag_output) write(38,'(A)') '# z, fm:  x, fm:     rhoz_bar, fm^-3:   rhoz_antibar, fm^-3:    temp, GeV:'
    if (getRMF_flag()) factor=ModificationFactor(nucleon,.true.)
    rho_bar_max=-0.1
    rho_antibar_max=-0.1
    temp_max=-0.1

    Loop_over_yGrid2 : do j= -gridpoints(2),gridpoints(2)
       Loop_over_zGrid2 : do k= -gridpoints(3),gridpoints(3)
          Loop_over_xGrid2 : do i= -gridpoints(1),gridpoints(1)

             rhoz= densField(i,j,k)%baryon(0)

             if (getRMF_flag()) then
                rhoz_bar=(rhoz+factor*totalDens(i,j,k))/(1.+factor)  ! density of the baryons
                rhoz_antibar=(totalDens(i,j,k)-rhoz)/(1.+factor)     ! density of the antibaryons
             else
                rhoz_bar=(rhoz+totalDens(i,j,k))/2.         ! density of the baryons
                rhoz_antibar=(rhoz-totalDens(i,j,k))/2.     ! density of the antibaryons
             end if

             temp=temperatureAt((/float(i)*gridSpacing(1),float(j)*gridSpacing(2),float(k)*gridSpacing(3)/))   ! temperature

             if (rhoz_bar.gt.rho_bar_max) rho_bar_max=rhoz_bar
             if (rhoz_antibar.gt.rho_antibar_max) rho_antibar_max=rhoz_antibar
             if (temp.gt.temp_max) temp_max=temp

             if (flag_output .and. j.eq.0) then
                write(38,5)  float(k)*gridSpacing(3),float(i)*gridSpacing(1),&
                     & rhoz_bar, rhoz_antibar, temp
                if (i.eq.gridpoints(1))  write(38,*)
             end if

          end do Loop_over_xGrid2
       end do Loop_over_zGrid2
    end do Loop_over_yGrid2

    write(39,5) time, rho_bar_max, rho_antibar_max, temp_max

  end subroutine HeavyIon_evol


  !****************************************************************************
  !****if* HeavyIonAnalysis/IsFree
  ! NAME
  ! logical function IsFree(realParts,i,j,teilchen)
  !
  ! PURPOSE
  ! Determine whether j-th particle of i-th parallel ensemble is free or not
  !****************************************************************************
  logical function IsFree(realParts,i,j,teilchen)

    use particleDefinition
    use densitymodule, only: Particle4MomentumRMF
    use RMF, only: getRMF_flag
    use coulomb, only: emfoca

    type(particle), dimension(:,:), intent(in)  :: realParts
    integer, intent (in) :: i, j
    type(particle), intent(in)  ::  teilchen  ! particle of interest

    integer, parameter :: imode=1   ! 1 - criterion according to dstmin
    ! 2 - criterion according to binding energy
    real, parameter :: dstmin = 3.  ! minimum distance for free particle (fm)
    real :: tmp,dist2,energy
    integer :: j1
    real, dimension(0:3) :: momentum
    real, dimension(1:3) :: place

    place=teilchen%pos

    if (imode.eq.1) then

       tmp=100.
       do j1=1,size(realParts,dim=2)
          if(teilchen%pert) then
              dist2=(realParts(i,j1)%pos(1)-place(1))**2 &
                   &+(realParts(i,j1)%pos(2)-place(2))**2 &
                   &+(realParts(i,j1)%pos(3)-place(3))**2
              if (dist2.lt.tmp) tmp=dist2
           else if (j.ne.j1) then
              if(teilchen%event(1).lt.1000000 .or. teilchen%event(1).ne.realParts(i,j1)%event(1)) then
                 dist2=(realParts(i,j1)%pos(1)-place(1))**2 &
                     &+(realParts(i,j1)%pos(2)-place(2))**2 &
                     &+(realParts(i,j1)%pos(3)-place(3))**2
                 if (dist2.lt.tmp) tmp=dist2
              end if
          end if
       end do
       IsFree = (tmp.gt.dstmin**2)

    else if (imode.eq.2) then

       if (getRMF_flag()) then
          momentum = Particle4MomentumRMF(teilchen)
          energy=momentum(0)
       else
          energy=teilchen%mom(0)
       end if
       energy=energy+emfoca(place,(/0.,0.,0./),realParts(i,j)%charge,realParts(i,j)%ID)
       IsFree = (energy.gt.teilchen%mass)


    end if

  end function IsFree

  !****************************************************************************
  !****is* HeavyIonAnalysis/countParts
  ! NAME
  ! subroutine countParts(Parts,time,prefix)
  !
  ! PURPOSE
  ! This counts and prints the number of particles at a specific time step or
  ! at the end (final).
  !
  ! NOTES
  ! This routine does not work with multiple events per energy, but only
  ! prints the results of the latest run. So you may average (at least for the
  ! final output) over all lines of the file.
  !****************************************************************************
  subroutine countParts(Parts,time,prefix)

    use histMP, only: Map2HistMP, Map2HistMP_getN, WriteHistMP_Names
    use particleDefinition, only: particle

    type(particle),dimension(:,:),intent(in), target :: Parts
    real, intent(in) :: time
    character*(*),intent(in) :: prefix
!    logical, intent(in) :: finalFlag


    integer :: nHist, nEns, nPart, i, j, iID, iSet
    type(particle), POINTER :: pPart
    real :: mulfak, w
    logical, save :: first = .true.

    !=== allocate arrays at first call ===

    if (first) then
       do iSet=1,nSet
          if (.not.useSet(iSet)) cycle
          nHist = Map2HistMP_getN(iSet)
          allocate( arrMultSet(iSet)%v(0:nHist) )

          ! open all files as 'new' and write the header:

          open(123,file='Mult_Set'//achar(48+iSet)//'.dat', status="unknown")
          call WriteHistMP_Names(iSet,123)
          close(123)

          open(123,file='final_Mult_Set'//achar(48+iSet)//'.dat', status="unknown")
          call WriteHistMP_Names(iSet,123)
          close(123)

       end do
       first = .false.
    end if

    !=== reset all values ===

    do iSet=1,nSet
       if (.not.useSet(iSet)) cycle
       arrMultSet(iSet)%v = 0.0
    end do

    nEns  = size(Parts,dim=1)
    nPart = size(Parts,dim=2)
    mulfak = 1.0/nEns

    !=== count all particles ===

    do i=1,nEns
       do j=1,nPart
          pPart => Parts(i,j)
          if(pPart%Id <  0) exit
          if(pPart%Id == 0) cycle

          w = 1.0
          if (pPart%pert) w = pPart%perweight

          do iSet=1,nSet
             if (.not.useSet(iSet)) cycle
             arrMultSet(iSet)%v(0) = arrMultSet(iSet)%v(0) + w

             iID = Map2HistMP(pPart, iSet)
             if (iID>0) then
                arrMultSet(iSet)%v(iID) = arrMultSet(iSet)%v(iID) + w
             endif
          end do

       end do
    end do

    !=== produce output ===

    ! all files are opened in append mode...
    ! this may be difficult for multiple runs per energy, since these are
    ! listed just behind each other (column 1 is a zikzak)

    do iSet=1,nSet
       if (.not.useSet(iSet)) cycle

       open(123,file=prefix//'Mult_Set'//achar(48+iSet)//'.dat', &
            status="old",position='append')
       write(123,'(f11.5,1P,100E12.4,0P)') time, &
            arrMultSet(iSet)%v(1:)*mulfak, &
            arrMultSet(iSet)%v(0)*mulfak
       close(123)

    end do

  end subroutine countParts

  !****************************************************************************
  !****is* HeavyIonAnalysis/doTmunu
  ! NAME
  ! subroutine doTmunu(realParts,timestep)
  !
  ! PURPOSE
  ! do the Tmunu analysis
  !****************************************************************************
  subroutine doTmunu(realParts,timestep)

    use IdTable
    use particleDefinition
    use TmunuDefinition
    use output, only: intTochar4
    use rotation, only: rotateTo, rotateFrom
    use constants, only: pi
    use densityModule, only: boostToLRF
    use potentialMain, only: potential_LRF, trueEnergy

    type(particle), dimension(:,:), intent(in), target  :: realParts
    integer, intent(in) :: timestep

    logical, save :: initFlagTmunu=.true.

    integer :: iEns,iPart
    type(particle) :: part
    type(tTmunuNmu), dimension(:,:), allocatable, save :: ArrX,ArrY,ArrZ


    real, parameter :: dX = 0.2
    integer :: iBin
    integer, parameter :: nBin = 100
    real, save :: mulFak = 1.0

    type(tTmunuNmu), save :: tTmunuNmu0 ! used to reset the array !
    real :: w
    integer :: iArr, nArr

    if (.not.BarMes_Tmunu) then
       nArr = 1
    else
       nArr = 2
    end if
    iArr = 1


    if (initFlagTmunu) then

       allocate(ArrX(nArr,nBin))
       allocate(ArrY(nArr,nBin))
       allocate(ArrZ(nArr,nBin))

       mulFak = 1.0/(size(realParts,dim=1)*dX**3) ! = 1/(nEns*V)
       initFlagTmunu = .false.
    end if

    ArrX = tTmunuNmu0
    ArrY = tTmunuNmu0
    ArrZ = tTmunuNmu0


    do iEns=1,size(realParts,dim=1)
       do iPart=1,size(realParts,dim=2)

          if (realParts(iEns,iPart)%ID == EOV) exit
          if (realParts(iEns,iPart)%ID == NOP) cycle

          part = realParts(iEns,iPart) ! create local copy
          w = 1.0

          select case (correctPot_Tmunu)

          case (1)
             call boostToLRF(part,1)  ! boost from calculation frame to LRF
             part%mom(0) = part%mom(0) - potential_LRF(part)
             call boostToLRF(part,2)  ! boost from LRF to calculation frame

          case (2)

             part%mom(0) = trueEnergy(part, .true.)

          case (3)

             part%mom(0) = trueEnergy(realParts(iEns,iPart), .true., part)


          end select

          if (rotateZ_Tmunu) then

             part%mom(1:3) = rotateFrom( part%pos(1:3), &
                  part%mom(1:3) )
             part%pos(1:3) = rotateFrom( part%pos(1:3), &
                  part%pos(1:3) )

             w = w * dX**2/(4*pi* part%pos(3)**2)
          end if

          if (BarMes_Tmunu) then
             if (isBaryon(part%ID)) then
                iArr = 1
             else if (isMeson(part%ID)) then
                iArr = 2
             else
                cycle
             end if
          end if

          if ( part%pos(1)>0 &
               .and. abs(part%pos(2))<dX/2 &
               .and. abs(part%pos(3))<dX/2 ) then

             iBin = int( part%pos(1) / dX )+1
             if (iBin <= nBin) call fillTmunu(ArrX(iArr,iBin), part)
          end if

          if ( part%pos(2)>0 &
               .and. abs(part%pos(1))<dX/2 &
               .and. abs(part%pos(3))<dX/2 ) then

             iBin = int( part%pos(2) / dX )+1
             if (iBin <= nBin) call fillTmunu(ArrY(iArr,iBin), part)
          end if

          if ( part%pos(3)>0 &
               .and. abs(part%pos(1))<dX/2 &
               .and. abs(part%pos(2))<dX/2 ) then

             iBin = int( part%pos(3) / dX )+1
             if (iBin <= nBin) call fillTmunu(ArrZ(iArr,iBin), part, w)
          end if

       end do
    end do

    if (iand(selectTmunuFormat,1)==1) then ! ASCII
       do iArr=1,nArr
          open(123,file='Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat', status="unknown")
          write(123,'(A)') headTmunu

          do iBin=1,nBin
             write(123,'(f11.4,1P,100E14.6,0P)') iBin*dX - dX/2, &
                  & ArrX(iArr,iBin)%Tmunu(:)*mulfak, &
                  & ArrX(iArr,iBin)%Nmu(:)*mulfak, &
                  & ArrX(iArr,iBin)%Jmu(:)*mulfak, &
                  & ArrX(iArr,iBin)%B*mulfak, ArrX(iArr,iBin)%S*mulfak
          end do
          write(123,*)
          write(123,*)

          do iBin=1,nBin
             write(123,'(f11.4,1P,100E14.6,0P)') iBin*dX - dX/2, &
                  & ArrY(iArr,iBin)%Tmunu(:)*mulfak, &
                  & ArrY(iArr,iBin)%Nmu(:)*mulfak, &
                  & ArrY(iArr,iBin)%Jmu(:)*mulfak, &
                  & ArrY(iArr,iBin)%B*mulfak, ArrY(iArr,iBin)%S*mulfak
          end do
          write(123,*)
          write(123,*)

          do iBin=1,nBin
             write(123,'(f11.4,1P,100E14.6,0P)') iBin*dX - dX/2, &
                  & ArrZ(iArr,iBin)%Tmunu(:)*mulfak, &
                  & ArrZ(iArr,iBin)%Nmu(:)*mulfak, &
                  & ArrZ(iArr,iBin)%Jmu(:)*mulfak, &
                  & ArrZ(iArr,iBin)%B*mulfak, ArrZ(iArr,iBin)%S*mulfak
          end do
          write(123,*)
          write(123,*)

          close(123)
       end do
    end if

    if (iand(selectTmunuFormat,2)==2) then ! BINARY
       do iArr=1,nArr

          open(123,file='Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat.bin', &
               status="unknown", form="unformatted")
          rewind(123)


          do iBin=1,nBin
             write(123) iBin*dX - dX/2, &
                  & ArrX(iArr,iBin)%Tmunu(:)*mulfak, &
                  & ArrX(iArr,iBin)%Nmu(:)*mulfak, &
                  & ArrX(iArr,iBin)%Jmu(:)*mulfak, &
                  & ArrX(iArr,iBin)%B*mulfak, ArrX(iArr,iBin)%S*mulfak
          end do

          ! how to separate data blocks in binary files????

          do iBin=1,nBin
             write(123) iBin*dX - dX/2, &
                  & ArrY(iArr,iBin)%Tmunu(:)*mulfak, &
                  & ArrY(iArr,iBin)%Nmu(:)*mulfak, &
                  & ArrY(iArr,iBin)%Jmu(:)*mulfak, &
                  & ArrY(iArr,iBin)%B*mulfak, ArrY(iArr,iBin)%S*mulfak
          end do

          do iBin=1,nBin
             write(123) iBin*dX - dX/2, &
                  & ArrZ(iArr,iBin)%Tmunu(:)*mulfak, &
                  & ArrZ(iArr,iBin)%Nmu(:)*mulfak, &
                  & ArrZ(iArr,iBin)%Jmu(:)*mulfak, &
                  & ArrZ(iArr,iBin)%B*mulfak, ArrZ(iArr,iBin)%S*mulfak
          end do

          close(123)
       end do
    end if

  end subroutine doTmunu

  !****************************************************************************
  !****s* HeavyIonAnalysis/doQRvector
  ! NAME
  ! subroutine doQRvector
  !
  ! PURPOSE
  ! calculate the Q- and R-vectors
  !****************************************************************************
  subroutine doQRvector(realParts,timestep)

    use particleDefinition
    use IdTable, only: EOV,NOP
    use EccAndFlow

    type(particle), dimension(:,:), intent(in), target  :: realParts
    integer, intent(in) :: timestep

    integer :: iEns,iPart
    type(particle) :: part

    type(tQRvector), dimension(1:6) :: arrQ
    type(tQRvector), dimension(1:6) :: arrR
    integer :: i
    real :: y, r,phi

    logical, save :: isFirst = .true.

    if (isFirst) then
       open(765,file='QR_vector.dat', status="unknown")
       isFirst = .false.
    end if

    do i=1,6
       call QRvectorInit(arrR(i), i, i)
       call QRvectorInit(arrQ(i), i, 0)
    end do
    call QRvectorInit(arrR(1), 1, 3) ! this is different!

    do iEns=1,size(realParts,dim=1)
       do iPart=1,size(realParts,dim=2)

          if (realParts(iEns,iPart)%ID == EOV) exit
          if (realParts(iEns,iPart)%ID == NOP) cycle

          part = realParts(iEns,iPart) ! create local copy

          if (abs(part%pos(3)) < 0.01) then
             r = sqrt(part%pos(1)**2+part%pos(2)**2)
             phi = atan2(part%pos(1),part%pos(2))
             do i=1,6
                call QRvectorAdd(arrR(i), r, phi)
             end do
          end if

          y = rapidity(part)
          if (abs(y) < 0.05) then
             r = sqrt(part%mom(1)**2+part%mom(2)**2)
             phi = atan2(part%mom(1),part%mom(2))
             do i=1,6
                call QRvectorAdd(arrQ(i), r, phi)
             end do
          end if

       end do
    end do

    write(765,*) timestep, (QRvectorW(arrR(i)),i=1,3),(QRvectorVal(arrR(i)),i=1,6), (QRvectorVal(arrQ(i)),i=1,6)
    flush(765)

  end subroutine doQRvector


  !****************************************************************************
  !****s* HeavyIonAnalysis/getRapBinning
  ! NAME
  ! subroutine getRapBinning(nBins, Bins)
  !
  ! PURPOSE
  ! return the array RapBinning and its size
  !****************************************************************************
  subroutine getRapBinning(nBins, Bins)
    use CallStack, only: TRACEBACK

    integer, intent(out) :: nBins
    real, dimension(:), intent(out) :: Bins

    if (initFlag) call init

    if (ubound(Bins,dim=1) < nRapBinning) &
         call Traceback("input array too small.")
    Bins(lbound(Bins,dim=1):lbound(Bins,dim=1)+nRapBinning) &
         = rapBinning(0:nRapBinning)
    nBins = nRapBinning

  end subroutine getRapBinning

  !****************************************************************************
  !****s* HeavyIonAnalysis/calcGlauber
  ! NAME
  ! subroutine calcGlauber(realParts)
  !
  ! PURPOSE
  ! This routine uses the initialized particles to do a Glauber-MC calculation
  ! of nPart.
  !
  ! This routine should be called directly after the initialization in
  ! timestep 0
  !
  ! NOTES
  ! We are doing a very bad trick in abusing the flag %anti of the particles
  ! in order to indicate, whether the particle is wounded.
  ! After the calculation, all flags are reset to .false.
  !****************************************************************************
  subroutine calcGlauber(realParts,nucA,nucB)
    use particleDefinition
    use nucleusDefinition, only: tNucleus
    use constants, only: pi
    use initHeavyIon, only: b_HI => b

    type(particle), dimension(:,:), intent(inOut), target :: realParts
    type(tNucleus), pointer :: nucA,nucB

    integer :: iEns1, iEns2, iA,iB
    integer :: nEns, nA, nB
    integer :: nPart
    type(particle), pointer :: pA, pB
    real :: s2, sigma
    real, dimension(0:2) :: SS
    integer :: iUnit
    real :: R1=-99.9, R2=-99.9

    real, parameter :: sigma0=23.8 ! in mb, dummy value, to be refined!!!!

    logical, parameter :: doFull = .false.

    if (initFlag) call init
    if (.not.do_Glauber) return

    nEns = size(realParts,dim=1)
    nA = nucA%Mass
    nB = nucB%Mass


    !===== Version 1: parallel ensemble =====

    call setWoundedParallel
    call countWounded
    R1 = SS(1)/SS(0)

    !===== Version 2: full ensemble =====
    ! this stuff is numerically very expensive!!!

    if (doFull) then
       call setWoundedFull
       call countWounded
       R2 = SS(1)/SS(0)
    end if

    write(*,*) 'Glauber: ',b_HI,R1,R2
    open(newunit=iUnit,file='Glauber.dat',status='unknown')
    write(iUnit,*) b_HI,R1,R2
    close(iUnit)


  contains

    subroutine setWoundedParallel
      do iEns1=1,nEns
         do iA=1,nA
            pA => realParts(iEns1,iA)
            do iB=nA+1,nA+nB
               pB => realParts(iEns1,iB)

               if (pA%anti.and.pB%anti) cycle ! not to check anymore

               ! transversal distance squared:
               s2 = (pA%pos(1)-pB%pos(1))**2+(pA%pos(2)-pB%pos(2))**2

               ! cross section:
               sigma = sigma0

               ! check s^2 < sigma/pi: (get fm^2 and mb right!)
               if (10*pi*s2 < sigma) then
                  pA%anti = .true.
                  pB%anti = .true.
               end if

            end do
         end do
      end do
    end subroutine setWoundedParallel

    subroutine setWoundedFull
      do iEns1=1,nEns
         do iA=1,nA
            pA => realParts(iEns1,iA)
            do iEns2=1,nEns
               do iB=nA+1,nA+nB
                  pB => realParts(iEns2,iB)

                  if (pA%anti.and.pB%anti) cycle ! not to check anymore

                  ! transversal distance squared:
                  s2 = (pA%pos(1)-pB%pos(1))**2+(pA%pos(2)-pB%pos(2))**2

                  ! cross section:
                  sigma = sigma0

                  ! check s^2 < sigma/pi: (get fm^2 and mb right!)
                  if (10*pi*s2*nEns < sigma) then
                     pA%anti = .true.
                     pB%anti = .true.
                  end if

               end do
            end do
         end do
      end do
    end subroutine setWoundedFull

    subroutine countWounded
      SS = 0
      ! probably all 'vectorized' ideas do not work...
      do iEns1=1,nEns
         nPart = 0
         do iA=1,nA+nB
            if (realParts(iEns1,iA)%anti) nPart=nPart+1
         end do
         !       write(*,*) iEns1,nPart
         SS = SS + (/1., 1.*nPart, 1.*nPart**2/)
      end do
      !    write(*,*) SS
      !    write(*,*) 'nPart = ',SS(1)/SS(0),'+-',sqrt(SS(1)**2-SS(2))/(SS(0)*(SS(0)-1))
      !    write(*,*) 'nPart = ',SS(1)/SS(0)


      !---- here we also reset the %anti-flag:
      realParts%anti = .false.

    end subroutine countWounded

  end subroutine calcGlauber

end module HeavyIonAnalysis
