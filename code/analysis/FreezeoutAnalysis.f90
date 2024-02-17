!******************************************************************************
!****m* /FreezeoutAnalysis
! NAME
! module FreezeouAnalysis
!
! PURPOSE
! This module provides routines to do a 'freeze out' analysis, i.e. it yields
! access to the position of particles at their last interaction. Particles
! being subject to some potentials have some additional freeze out condition,
! e.g. when the baryon density drops to the value 0.2/fm^3.
!******************************************************************************
module FreezeoutAnalysis

  use CallStack, only: traceback

  implicit none
  private

  public :: getFreezeoutAnalysis_Pert
  public :: getFreezeoutAnalysis_Real
  public :: DoFreezeoutAnalysisPerTime
  public :: DoFreezeoutAnalysisFinalize

  !****************************************************************************
  !****g* FreezeoutAnalysis/FreezeoutAnalysis_Pert
  ! SOURCE
  !
  logical, save :: FreezeoutAnalysis_Pert = .false.
  !
  ! PURPOSE
  ! Flag to do freeze out analysis for perturbative particles
  !****************************************************************************

  !****************************************************************************
  !****g* FreezeoutAnalysis/FreezeoutAnalysis_Real
  ! SOURCE
  !
  logical, save :: FreezeoutAnalysis_Real = .false.
  !
  ! PURPOSE
  ! Flag to do freeze out analysis for real particles
  !****************************************************************************

  !****************************************************************************
  !****g* FreezeoutAnalysis/potThreshold
  ! SOURCE
  !
  real, save :: potThreshold = 0.005
  !
  ! PURPOSE
  ! threshold value in GeV. If the absolute value of the potential is below
  ! this value, the particle is considered to be 'free', e.g. it 'escaped'
  !****************************************************************************

  logical, save :: init = .true.

  real, save, dimension(0:13) :: RapBinning
  integer, save :: nRapBinning

contains

  !****************************************************************************
  !****f* FreezeoutAnalysis/getFreezeoutAnalysis_Pert
  ! NAME
  ! logical function getFreezeoutAnalysis_Pert()
  ! PURPOSE
  ! return the value of FreezeoutAnalysis_Pert
  !****************************************************************************
  logical function getFreezeoutAnalysis_Pert()
    if (init) call initInput
    getFreezeoutAnalysis_Pert = FreezeoutAnalysis_Pert
  end function getFreezeoutAnalysis_Pert

  !****************************************************************************
  !****f* FreezeoutAnalysis/getFreezeoutAnalysis_Real
  ! NAME
  ! logical function getFreezeoutAnalysis_Real()
  ! PURPOSE
  ! return the value of FreezeoutAnalysis_Real
  !****************************************************************************
  logical function getFreezeoutAnalysis_Real()
    if (init) call initInput
    getFreezeoutAnalysis_Real = FreezeoutAnalysis_Real
  end function getFreezeoutAnalysis_Real

  !****************************************************************************
  !****s* FreezeoutAnalysis/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Read namelist 'Freezeout' from jobcard.
  !****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput
    use HeavyIonAnalysis, only: getRapBinning

    !**************************************************************************
    !****n* FreezeoutAnalysis/Freezeout
    ! NAME
    ! namelist /Freezeout/
    ! PURPOSE
    ! Namelist for FreezeoutAnalysis includes:
    ! * FreezeoutAnalysis_Pert
    ! * FreezeoutAnalysis_Real
    ! * potThreshold
    !**************************************************************************
    namelist /Freezeout/ FreezeoutAnalysis_Pert, FreezeoutAnalysis_Real, &
         potThreshold

    integer :: ios !, i,j
    !character(*), parameter :: cPi(-1:1) = (/ "-", "0", "+" /)

    call Write_ReadingInput('Freezeout',0)
    rewind(5)
    read(5,nml=Freezeout,iostat=ios)

    call Write_ReadingInput('Freezeout',0,ios)

    write(*,*) 'Freezeout analysis (real,pert): ', &
               FreezeoutAnalysis_Real, FreezeoutAnalysis_Pert
    write(*,*) 'threshold potential: ',potThreshold,' GeV'
    call Write_ReadingInput('Freezeout',1)

    call getRapBinning(nRapBinning,RapBinning)

!!$    do i=1,nRapBinning
!!$       write(*,*) 'y-Bin ',i,': ',rapBinning(i-1),'...',rapBinning(i)
!!$    end do

!!$    do i=-1,1
!!$       do j=1,nRapBinning
!!$          write(*,'(" (rap. bin #",i2,": y = ",f5.2," ... ",f5.2,") (pi",A,")")') &
!!$               j,rapBinning(j-1),rapBinning(j),cPi(i)
!!$       end do
!!$    end do
!!$
!!$    stop

    init = .false.

  end subroutine initInput

  !****************************************************************************
  !****s* FreezeoutAnalysis/DoFreezeoutAnalysisPerTime
  ! NAME
  ! subroutine DoFreezeoutAnalysisPerTime(iTime, Time, realPart, pertPart)
  ! PURPOSE
  ! Do the analysis after each time step
  !****************************************************************************
  subroutine DoFreezeoutAnalysisPerTime(iTime, Time, realPart, pertPart)
    use particleDefinition
    use particleProperties
    use PIL_freezeout, only: PIL_freezeout_PUT, PIL_freezeout_GET
    use potentialMain, only: potential_LRF
    use densitymodule, only: densityAt
    use dichteDefinition

    integer, intent(in)       :: iTime
    real,    intent(in)       :: Time
    type(particle), intent(in), dimension(:,:),target :: realPart, pertPart

    integer :: i,j
    type(particle), POINTER :: pPart

    real, dimension(0:3) :: pos
    integer :: hist, histOld
    logical :: escaped, escapedOld
    logical :: changed
    real :: pot
    type(dichte) :: dens

    integer :: nNew, nChanged

    if (init) call initInput

    if (FreezeoutAnalysis_Real) then

       nNew = 0
       nChanged = 0

       do i=1,Size(realPart,dim=1)
          do j=1,Size(realPart,dim=2)
             pPart => realPart(i,j)
             if (pPart%Id <  0) exit
             if (pPart%Id <= 0) cycle

             if (pPart%Id /= 101) cycle ! only for pions

             pos(0) = Time
             pos(1:3) = pPart%pos
             hist = pPart%history

             dens=densityAt(pos(1:3))
             pot = potential_LRF(pPart, .true.)
             escaped = .true. ! dummy !!!
!             escaped = (abs(pot)<potThreshold) ! ????

             changed = .not.PIL_freezeout_GET(pPart%number,histOld,escapedOld)

             nNew = nNew+1

             if (.not.changed) then ! i.e. it is already in the list
!                write(*,*) 'got: ',pPart%number,posOld,histOld
                if (hist .ne. histOld) changed = .true.
                if (escaped .neqv. escapedOld) changed = .true.
                if (.not.escaped) changed = .true.

                if (changed) nChanged = nChanged+1
             end if

             if (changed) then
                call PIL_freezeout_PUT(pPart%number,hist,escaped,&
                     pos,pPart%mom(0:3),dens%baryon,pot)
             end if

          end do
       end do

!!$       write(*,*) 'Freezout: new:',nNew,' changed:',nChanged
!!$       write(642,*) time, nNew, nChanged
!!$       flush(642)
    end if

    if (FreezeoutAnalysis_Pert) then

       do i=1,Size(pertPart,dim=1)
          do j=1,Size(pertPart,dim=2)
             pPart => pertPart(i,j)
             if (pPart%Id <  0) exit
             if (pPart%Id <= 0) cycle

             if (hadron(pPart%ID)%stability .ne. 0) cycle

             pos(0) = Time
             pos(1:3) = pPart%pos
             hist = pPart%history

             dens=densityAt(pos)
             pot = potential_LRF(pPart)
!             escaped = .true. ! dummy !!!
             escaped = (abs(pot)<potThreshold) ! ????

             changed = .not.PIL_freezeout_GET(pPart%number,histOld,escapedOld)

             if (.not.changed) then ! i.e. it is already in the list
!                write(*,*) 'got: ',pPart%number,posOld,histOld
                if (hist .ne. histOld) changed = .true.
                if (escaped .neqv. escapedOld) changed = .true.
                if (.not.escaped) changed = .true.
             end if

             if (changed) then
!                write(*,*) 'set: ',pPart%number,pos,history
                call PIL_freezeout_PUT(pPart%number,hist,escaped,pos,&
                     pPart%mom(0:3), dens%baryon,pot)
             end if

          end do
       end do
    end if

  end subroutine DoFreezeoutAnalysisPerTime

  !****************************************************************************
  !****s* FreezeoutAnalysis/DoFreezeoutAnalysisFinalize
  ! NAME
  ! subroutine DoFreezeoutAnalysisFinalize(realPart, pertPart)
  ! PURPOSE
  ! Do the analysis after each time step
  !****************************************************************************
  subroutine DoFreezeoutAnalysisFinalize(realPart, pertPart)
    use particleDefinition

    type(particle), intent(in), dimension(:,:),target :: realPart, pertPart

    if (init) call initInput
    if (FreezeoutAnalysis_Real) call doReal
    if (FreezeoutAnalysis_Pert) call doPert

  contains
    subroutine doReal

      use particleProperties
      use potentialMain, only: potential_LRF
      use PIL_freezeout, only: PIL_freezeout_PUT, PIL_freezeout_GET
      use history, only: history_getParents
      use output, only: intToChar2
      use histMC
      use histMC_avg
      use inputGeneral, only: numTimeSteps, delta_T

      logical, save :: first = .true.

      integer :: nEns,i,j
      type(particle), POINTER :: pPart

      real, dimension(0:3) :: pos, mom, bDens
      integer :: hist
      logical :: escaped
      logical :: inList
      real :: pot, rap
      integer :: iCh,iHist,iBin
      integer, dimension(1:3) :: parents
      integer :: nNew, nChanged

      character(100) :: str
      real :: mulfak, timeH, momH
      type(histogramMC), dimension(-1:1), save :: hMC_time

      type(histogramMC_avg), dimension(-1:1), save :: hMCa_time,hMCa_dens,hMCa_pot
      type(histogramMC_avg), dimension(1:nRapBinning,-1:1) :: hMCa_potY

      character(*), parameter :: sPi(-1:1) = (/ "pim", "pi0", "pip" /)
      character(*), parameter :: cPi(-1:1) = (/ "-", "0", "+" /)

      if (first) then
         first = .false.
      else
         call Traceback('multiple runs nyi.')
      end if

      call CreateHistMC(hMC_time, "dN/dt_freeze", 0., numTimeSteps*delta_T, 0.2, 4)

      call CreateHistMC_avg(hMCa_time, "<t_freeze> vs p", 0., 1., 0.01, 4)
      call CreateHistMC_avg(hMCa_dens, "<rho_freeze> vs p", 0., 1., 0.01, 4)
      call CreateHistMC_avg(hMCa_pot, "<pot_freeze> vs p", 0., 1., 0.01, 4)
      call CreateHistMC_avg(hMCa_potY, "<pot_freeze> vs p", 0., 1., 0.01, 4)

      hMC_time%xDesc = "t_freeze [fm]"
      hMCa_time%xDesc = "p_cm [GeV]"

!!$      write(*,*) 'nRapBinning: ',nRapBinning
!!$      write(*,*) 'lbound:',lbound(hMCa_potY)
!!$      write(*,*) 'ubound:',ubound(hMCa_potY)

      do i=-1,1
         do j=1,nRapBinning
            write(str,'(" (rap. bin #",i2,": y = ",f5.2," ... ",f5.2,") (pi",A,")")') &
                 j,rapBinning(j-1),rapBinning(j),cPi(i)
            hMCa_potY(j,i)%name = trim(hMCa_potY(j,i)%name) // str
         end do
         write(str,'(A,i3,A)') " (charge: ",i,")"
         hMC_time(i)%name = trim(hMC_time(i)%name) // str
         hMCa_time(i)%name = trim(hMCa_time(i)%name) // str
         hMCa_dens(i)%name = trim(hMCa_dens(i)%name) // str
         hMCa_pot(i)%name = trim(hMCa_pot(i)%name) // str

         hMC_time(i)%yDesc = (/ &
           "decay: Delta", "decay: Res  ", &
           "decay: other", "other       " /)

         hMCa_time(i)%yDesc = (/ &
           "decay: Delta", "decay: Res  ", &
           "decay: other", "other       " /)
      end do

      call copyDesc_avg(hMCa_dens, hMCa_time)
      call copyDesc_avg(hMCa_pot, hMCa_time)
      do j=1,nRapBinning
         call copyDesc_avg(hMCa_potY(j,:), hMCa_time)
      end do

      nNew = 0
      nChanged = 0

      nEns  = size(realPart,dim=1)
      mulfak = 1.0/(nEns)

      do i=1,nEns
         do j=1,Size(realPart,dim=2)
            pPart => realPart(i,j)
            if (pPart%Id <  0) exit
            if (pPart%Id <= 0) cycle

            if (pPart%Id /= 101) cycle ! only for pions

            if ((pPart%charge < -1).or.(pPart%charge > 1)) &
                 call Traceback("wrong pion charge")
            iCh = pPart%charge

            inList = PIL_freezeout_GET(pPart%number,hist,escaped,pos,mom,bDens,pot)

            if (.not.inList) then
               write(*,*) 'Oops, particle not in list.'
               cycle
            end if

            parents = history_getParents(hist)
            iHist = 0
            if (parents(2) == 0) then
               select case(parents(1))
               case (2)
                  iHist = 1
               case (3:99)
                  iHist = 2
               case default
                  iHist = 3
               end select
            else
               iHist = 4
            end if

            if (iHist == 0) call Traceback("wrong hist.")

            timeH = pos(0) - 1e-3 ! for avoiding problems at bin boundary
            momH = absMom(pPart)

            call addHistMC(hMC_time(iCh), timeH, iHist, 1.)
            call addHistMC_avg(hMCa_time(iCh), momH, iHist, timeH)
            call addHistMC_avg(hMCa_dens(iCh), momH, iHist, bDens(0))
            call addHistMC_avg(hMCa_pot(iCh), momH, iHist, pot)

            rap = rapidity(pPart)
            iBin = 0
            do while (iBin < nRapBinning+1)
               if (rap < rapBinning(iBin)) exit
               iBin = iBin + 1
            end do
            if (iBin>0 .and. iBin<nRapBinning+1) &
                 call addHistMC_avg(hMCa_potY(iBin,iCh), momH, iHist, pot)


            write(643,*) iCh, iHist, hist, iBin, &
                 pos, mom, pPart%mom, bDens, pot


         end do
      end do

      do i=-1,1
         call writeHistMC(hMC_time(i), file="hMC_time."//sPi(i)//".dat",mul=mulfak,dump=.true.)
         call writeHistMC_avg(hMCa_time(i), file="hMCa_time."//sPi(i)//".dat",dump=.true.)
         call writeHistMC_avg(hMCa_dens(i), file="hMCa_dens."//sPi(i)//".dat",dump=.true.)
         call writeHistMC_avg(hMCa_pot(i), file="hMCa_pot."//sPi(i)//".dat",dump=.true.)
         do j=1,nRapBinning
            call writeHistMC_avg(hMCa_potY(j,i), file="hMCa_potY."//intToChar2(j)//"."//sPi(i)//".dat",dump=.true.)
         end do
      end do

    end subroutine doReal



    subroutine doPert

      call Traceback('DoFreezeoutAnalysisFinalize for pertPart nyi.')

    end subroutine doPert

  end subroutine DoFreezeoutAnalysisFinalize

end module FreezeoutAnalysis
