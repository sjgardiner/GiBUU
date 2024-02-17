!******************************************************************************
!****m* /residue
! NAME
! module residue
!
! PURPOSE
! Calculate a target residue in elementary induced reactions on a nucleus.
! (in perturbative mode)
!
! Using the prescription of particle-hole excitations, 'particles' are
! subtracted from the original nucleus.
!
! NOTES
! * Works both for non-relativistic Skyrme-mode and RMF-mode.
!******************************************************************************
module residue

  implicit none

  private

  !****************************************************************************
  !****g* residue/DetermineResidue
  ! PURPOSE
  ! If .true., then the determination of target residue properties for
  ! every event will be done.
  !
  ! Their output in file 'TargetResidue.dat' at the end of time evolution
  ! is called elsewhere. If nothing is stored, no output is generated.
  ! SOURCE
  logical, save :: DetermineResidue = .true.
  !****************************************************************************

  !****************************************************************************
  !****g* residue/mode
  ! PURPOSE
  ! select the mode, how the residue energy is determined (field res%mom(0)):
  ! * 1: the sum of hole excitation energies
  ! * 2: the sum of energies of the removed particles (with minus sign)
  !
  ! SOURCE
  integer, save :: mode = 1
  !****************************************************************************

  !****************************************************************************
  !****g* residue/switchOutput
  ! PURPOSE
  ! select the output
  ! * 1: write out TargetResidue.dat
  ! * 2: write out TargetResidue.Plot.dat
  ! * 3: write out both files
  !
  ! SOURCE
  integer, save :: switchOutput = 0
  !****************************************************************************


  !****************************************************************************
  !****t* residue/remnant
  ! NAME
  ! type remnant
  ! PURPOSE
  ! Stores properties of nuclear target remnant.
  ! SOURCE
  type remnant
     real, dimension(0:3) :: mom=0.      ! 4-momentum (GeV)
     integer              :: A=0         ! mass number
     integer              :: Z=0         ! charge number
     real                 :: weight=0.   ! perweight
  end type remnant
  !****************************************************************************

  !****************************************************************************
  !****ig* residue/ArrResidue
  ! SOURCE
  type(remnant), allocatable, dimension(:), TARGET, SAVE :: ArrResidue
  ! PURPOSE
  ! Stores mass number, charge number, excitation energy and momentum of
  ! all target residues.
  !****************************************************************************

  logical, save :: initFlag=.true.

  integer, save :: numEnsembles_ = 0
  integer, save :: numParticles_ = 0
  integer, save :: A_orig = 0
  integer, save :: Z_orig = 0

  type(remnant),save :: Remnant0 ! hold the init values (=zero)

  public :: InitResidue
  public :: ResidueAddPH
  public :: ResidueSetWeight
  public :: OutputResidue
  public :: ResidueGetScale

contains
  !****************************************************************************
  !****is* residue/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Read the namelist
  !****************************************************************************
  subroutine readInput

    use output, only: Write_ReadingInput
    use inputGeneral, only: eventType
    use eventtypes

    !**************************************************************************
    !****n* residue/residue_Input
    ! NAME
    ! NAMELIST /residue_Input/
    ! PURPOSE
    ! Includes the switches:
    ! * DetermineResidue
    ! * mode
    ! * switchOutput
    !**************************************************************************
    NAMELIST /residue_Input/ DetermineResidue, mode, switchOutput

    integer :: ios

    if (.not.initFlag) return

    call Write_ReadingInput('residue_Input',0)
    rewind(5)
    read(5,nml=residue_Input,iostat=ios)
    call Write_ReadingInput('residue_Input',0,ios)

    if (DetermineResidue) then
       select case (eventType)
       case (HiPion,HiLepton,Neutrino,hadron,ExternalSource,RealPhoton)
          ! nothing to do

       case default
          write(*,*) "residue not possible for eventtype =", eventType
          write(*,*) "switching it off!"
          DetermineResidue = .false.

       end select
    end if

    write(*,*) 'DetermineResidue : ',DetermineResidue
    write(*,*) 'mode         =',mode
    write(*,*) 'switchOutput =',switchOutput


    ! create empty file:
    if (iand(switchOutput,1)>0) then
       open(171,file='TargetResidue.dat',status='UNKNOWN')
       select case (eventType)
       case (HiLepton,hadron,ExternalSource)

          write(171,'(A)')'#iens: inuc:  A:     Z:    E^* (GeV):       p(1:3) (GeV/c):                                weight:'
       case DEFAULT
          write(171,'(A)')'#iEvent:  A:    Z:    E^* (GeV):       p(1:3) (GeV/c):                                weight:'
       end select
       close(171)
    end if

    call Write_ReadingInput('residue_Input',1)

    initFlag=.false.

  end subroutine readInput



  !****************************************************************************
  !****s* residue/InitResidue
  ! NAME
  ! subroutine InitResidue(numEnsembles,numParticles,A,Z)
  ! PURPOSE
  ! Initialize target residue.
  ! INPUTS
  ! * integer, intent(in) :: numEnsembles,numParticles
  ! * integer, intent(in) :: A,Z -- target nucleus mass and charge numbers
  ! NOTE
  ! Intention is to use this routine in elementary-induced reactions,
  ! where the perturbative event is initialized on every nucleon of every
  ! ensemble. Thus there are two indices identifying the residue for
  ! respective perturbative event.
  !****************************************************************************
  subroutine InitResidue(numEnsembles,numParticles,A,Z)
    use callstack, only: traceback
    use output, only: Write_InitStatus

    integer, intent(in) :: numEnsembles,numParticles
    integer, intent(in) :: A,Z

    call Write_InitStatus("Residue",0)

    if (initFlag) call readInput

    if (.not.DetermineResidue) return

    if (numEnsembles_ > 0) then
       if (numEnsembles_*numParticles_ .ne. numEnsembles*numParticles) then
          call traceback("sizes do not match")
       end if
    end if

    write(*,*)' numEnsembles, numParticles : ', numEnsembles, numParticles
    write(*,*)' A, Z : ', A, Z

    numEnsembles_ = numEnsembles
    numParticles_ = numParticles
    A_orig = A
    Z_orig = Z

    if (.not.allocated(ArrResidue)) &
         allocate(ArrResidue(1:numEnsembles*numParticles))

    Remnant0%A = A
    Remnant0%Z = Z
    ArrResidue = Remnant0 ! initialize the whole array

    call Write_InitStatus("Residue",1)

  end subroutine InitResidue


  !****************************************************************************
  !****s* residue/ResidueAddPH
  ! NAME
  ! subroutine ResidueAddPH(firstEvent,teilchen)
  ! PURPOSE
  ! Add a particle-hole excitation to the residue.
  ! INPUTS
  ! * integer, intent(in) :: firstEvent -- Index of colliding particle
  ! * type(particle), intent(in) :: teilchen -- Struck nucleon
  !****************************************************************************
  subroutine ResidueAddPH(firstEvent,teilchen)

    use particleDefinition
    use densitymodule, only: FermiMomAt, energyDeterminationRMF
    use energyCalc, only: energyDetermination
    use RMF, only: getRMF_flag, g_rho
    use baryonPotentialMain, only: getsymmetryPotFlag_baryon


    integer, intent(in) :: firstEvent
    type(particle), intent(in) :: teilchen

    type(particle) :: Teilchen_Fermi
    real :: pF
    integer :: iEvent
    type(remnant), pointer :: res

    if (initFlag) call readInput

    if(.not.DetermineResidue) return

    iEvent = ResidueCalcIndex(firstEvent)

    if(iEvent.lt.1 .or. iEvent.gt.size(ArrResidue)) then
         write(*,*) ' Wrong firstEvent in the input to ResidueAddPH '
         write(*,*) ' firstEvent : ', firstEvent
         write(*,*) ' iEvent : ', iEvent
         stop
    end if

    res => ArrResidue(iEvent)

    res%A = res%A-1
    res%Z = res%Z-teilchen%charge

    select case (Mode)
    case (1)

       if (getRMF_flag()) then

          if (g_rho/=0.) then
             pF=FermiMomAt(teilchen%pos,teilchen%charge)
          else
             pF=FermiMomAt(teilchen%pos)
          end if

       else

          if (getsymmetryPotFlag_baryon()) then
             pF=FermiMomAt(teilchen%pos,teilchen%charge)
          else
             pF=FermiMomAt(teilchen%pos)
          end if

       end if

       Teilchen_Fermi=teilchen
       Teilchen_Fermi%mom(1:3)=(/pF,0.,0./)

       if(getRMF_flag()) then
          call energyDeterminationRMF(Teilchen_Fermi)
       else
          call energyDetermination(Teilchen_Fermi)
       end if

       res%mom(0)=res%mom(0)+Teilchen_Fermi%mom(0)-teilchen%mom(0)

!           res%mom(0)=res%mom(0)+sqrt(pf**2+0.938**2)&
!                &-sqrt(teilchen%mom(1)**2+teilchen%mom(2)**2+teilchen%mom(3)**2+0.938**2)

       res%mom(1:3)=res%mom(1:3)-teilchen%mom(1:3)

    case (2)
       res%mom(0:3)=res%mom(0:3)-teilchen%mom(0:3)

    end select

  end subroutine ResidueAddPH

  !****************************************************************************
  !****************************************************************************
  real function ResidueGetScale(firstEvent,Q)
    integer, intent(in) :: firstEvent
    integer, intent(in) :: Q

    integer :: iEvent
    type(remnant), pointer :: res
    logical, parameter :: withIsospin = .true.

    ResidueGetScale = 1.0

    if (initFlag) call readInput
    if (.not.DetermineResidue) return

    if (firstEvent==0) return

    iEvent = ResidueCalcIndex(firstEvent)
    res => ArrResidue(iEvent)

    if (withIsospin) then
       select case (Q)
       case (0)
          ResidueGetScale = max(float(res%A-res%Z)/(A_orig-Z_orig), 0.)
       case (1)
          ResidueGetScale = max(float(res%Z)/(Z_orig), 0.)
       end select
    else
       ResidueGetScale = max(float(res%A)/(A_orig), 0.)
    end if

  end function ResidueGetScale

  !****************************************************************************
  !****************************************************************************
  subroutine ResidueSetWeight(firstEvent,weight)

    integer, intent(in) :: firstEvent
    real, intent(in) :: weight

    if(.not.DetermineResidue) return

    ArrResidue(ResidueCalcIndex(firstEvent))%weight = weight

  end subroutine ResidueSetWeight


  !****************************************************************************
  !****************************************************************************
  integer function ResidueCalcIndex(firstEvent)

    use inputGeneral, only: eventType
    use eventtypes
    use callstack, only: traceback

    integer, intent(in) :: firstEvent
    integer :: iEns,iPart,iEvent

    select case (eventType)
    case (HiLepton,hadron,ExternalSource)
       iEns = firstEvent/1000
       iPart= min(mod(firstEvent,1000),numParticles_)
       iEvent = (iEns-1)*numParticles_ + iPart

    case default
       iEvent = firstEvent

    end select

    if (iEvent > size(ArrResidue)) then
       write(*,*) 'iEvent: ',iEvent,"   max:",size(ArrResidue)
       call traceback("out of bounds")
    end if

    ResidueCalcIndex = iEvent

  end function ResidueCalcIndex

  !****************************************************************************
  !****s* residue/OutputResidue
  ! NAME
  ! subroutine OutputResidue
  ! PURPOSE
  ! Output of all residues to the file 'TargetResidue.dat'
  !****************************************************************************
  subroutine OutputResidue

    use inputGeneral, only: eventType
    use eventtypes
    use hist

    integer :: iEvent,iens,inuc
    type(remnant), pointer :: res

    type(histogram), save :: hNinteract, hMom
    real :: fak
    integer, save :: nCalls = 1

    if (initFlag) call readInput
    if(.not.DetermineResidue) return

    !==========================================================================

    if (iand(switchOutput,1)>0) then
       open(171,file='TargetResidue.dat',position='append')
       do iEvent=1,size(ArrResidue)
          res => ArrResidue(iEvent)
          if (eventType==HiLepton .or. eventType==hadron .or. eventType==ExternalSource) then
             iEns = (iEvent-1)/numParticles_+1
             iNuc = mod((iEvent-1),numParticles_)+1
             write(171,'(1x,i4,3x,i3,4x,i3,3x,i3,3x,1P,e14.7,3x,3(e14.7,1x),e13.4)') &
                  iEns,iNuc,res%A,res%Z,res%mom(0:3),res%weight
          else
             write(171,'(1x,i7,4x,i3,3x,i3,3x,1P,e14.7,3x,3(e14.7,1x),e13.4)') &
                  iEvent,res%A,res%Z,res%mom(0:3),res%weight
          end if
       end do
       close(171)
    end if

    !==========================================================================

    if (iand(switchOutput,2)>0) then
       if (nCalls==1) then
          call CreateHist(hNinteract, "number of interactions", 0.5, 20.5, 1.0)
          call CreateHist(hMom, "momentum", 0.0, 5.0, 0.01)
       end if

       do iEvent=1,size(ArrResidue)
          res => ArrResidue(iEvent)

          call addHist(hNinteract, real(Remnant0%A-res%A), res%weight)
          call addHist(hMom, sqrt(sum(res%mom(1:3)**2)), res%weight)
       end do

       fak = 1.0/(nCalls*numEnsembles_)

       open(171,file='TargetResidue.Plot.dat',status='unknown')
       call writeHist(hNinteract, 171, add=1e-20, mul=fak)
       write(171,*)
       write(171,*)
       call writeHist(hMom, 171, add=1e-20, mul=fak)

       close(171)
    end if

!    deallocate(ArrResidue)

    nCalls = nCalls+1

  end subroutine OutputResidue

end module residue


!!$ How 'firstevent' is set in the different inits:
!!$
!!$ LoPion:
!!$ =======
!!$
!!$ incoming Pions:
!!$ pertPart(iEns,iPart)%firstevent = i + size(pertPart,dim=1)*iEns
!!$
!!$ if pertPart was empty at the beginning, then i == iPart; otherwise it
!!$ points to a first epmty place in the vector (unpredictable!)
!!$
!!$ RealPhoton:
!!$ ============
!!$
!!$ pertRun: %firstevent = pert_numbering(realParticles(iEns,iPart))
!!$ realRun: %firstevent = real_numbering()
!!$
!!$ LoLepton:
!!$ ==========
!!$
!!$ private numbering in the init-module:
!!$
!!$ %firstevent = getFirstEvent()
!!$ with
!!$   getFirstEvent(): first = first+1
!!$
!!$ Therefore events are just numbered /1,first/
!!$ (first is reset at every call of init_...)
!!$
!!$ Neutrino:
!!$ ==========
!!$
!!$ equivalent to LoLepton !
!!$
!!$ HiPion:
!!$ ========
!!$
!!$ from 'collisionnumbering':
!!$ %firstEvent=pert_firstnumbering(pair(1),pair(2))
!!$
!!$ (pair(1) = realParts(...), pair(2) = pertPart(iEns,iPart) )
!!$
!!$ pert_firstnumbering() is special for HiPion:
!!$ it is just iCount=iCount+1
!!$ with a reset at ???
!!$
!!$ HiLepton:
!!$ ==========
!!$
!!$ both for realRun and pertRun:
!!$
!!$ %firstevent = iEns*1000+iPart
