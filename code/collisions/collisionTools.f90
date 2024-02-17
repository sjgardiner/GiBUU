!******************************************************************************
!****m* /collisionTools
! NAME
! module collisionTools
! PURPOSE
! Some helper routines for collisions
!******************************************************************************
module collisionTools

  implicit none
  private

  public:: finalCheck
  public:: setEnergyCheck

  logical, save :: initFlag=.true.
  real,    save :: energyCheck=0.01 ! copy of collisionTerm/energyCheck

contains

  !****************************************************************************
  !****s* collisionTools/setEnergyCheck
  ! NAME
  ! subroutine setEnergyCheck(va)
  ! PURPOSE
  ! set the internal variable energyCheck
  !****************************************************************************
  subroutine setEnergyCheck(val)
    real, intent(in) :: val
    energyCheck = val
    initFlag=.false.
  end subroutine setEnergyCheck


  !****************************************************************************
  !****f* collisionTools/finalCheck
  ! NAME
  ! function finalCheck(partIn, partOut, HiEnergyFlagge, woher) result(flag)
  ! PURPOSE
  ! Checks the final state of a collision for the conservation of all
  ! quantum numbers,  also including momentum and energy conservation.
  !
  ! For HiEnergy events we do not check charge and momentum conservation,
  ! this MUST be done separately.
  ! The reason for this is, that some particles which are produced by
  ! Pythia/Fritiof can not be propagated by BUU and therefore do not show up in
  ! the final state vector "partOut"
  ! ("unknown particles wont be propagated").
  !
  ! INPUTS
  ! * type(particle),dimension(:)  :: partIn  -- Incoming particles
  ! * type(particle),dimension(:)  :: partOut -- Outgoing particles
  ! * logical, optional            :: HiEnergyFlag --
  !   .true. if it was a HiEnergy event.
  !   if .true. then energy conservation is not checked
  !   and code does not stop if charge conservation is violated
  ! * character(*), optional      :: woher -- name of calling routine
  !
  ! OUTPUT
  ! * logical :: flag -- .true. if quantum numbers are conserved
  !
  ! NOTES
  ! you may change the behaviour of this routine by changing the internal
  ! parameter 'errcode' at compile time:
  ! * '-1': the checks just produce a error message, but does not stop.
  !   This is just for debugging purposes!
  ! * any other value >0: the code stops
  !
  !****************************************************************************
  function finalCheck(partIn, partOut, HiEnergyFlag, woher) result(flag)

    use particleDefinition
    use IdTable, only: isBaryon, isHadron
    use particleProperties, only: hadron, validCharge
    use callStack, only: traceBack
    use output, only: WriteParticle

    type(particle), intent(in), dimension(:) :: partIn, partOut
    logical, intent(in), optional :: HiEnergyFlag
    character(*), intent(in), optional :: woher
    logical :: flag

    integer :: iQ_In, iQ_Out, iB_In, iB_Out, iS_In, iS_Out, k
    real, dimension(0:3) :: mom_In, mom_Out
    integer, parameter :: errcode = 2 ! -1: code continues, >0: code stops

    flag=.false.

    if (initFlag) call Traceback('finalCheck: not initialized!')

    ! Check charge of each particle
    do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (.not.validCharge(partOut(k))) then
          if (present(woher)) write(*,*) 'Problem in ',trim(woher),' !!!!!!'
          write(*,*) 'Charge of particle is not valid in finalCheck'
          call WriteParticle(6,1,k,partOut(k))
          write(*,*) 'Charge=',partOut(k)%charge
          write(*,*) 'Id=',partOut(k)%id
          write(*,*) 'Antiparticle=',partOut(k)%anti
          write(*,*) 'Momentum=',partOut(k)%mom
          write(*,*)
          if (present(HiEnergyFlag)) then
             if (HiEnergyFlag) then
                write(*,*) 'It was a HiEnergy event!'
                call PYLIST(2)
             else
                write(*,*) 'It was a low-energy event!'
             end if
          else
             write(*,*) 'It was a decay-event!'
          end if
          write(*,*)
          write(*,'(79("*"))')
          call printPart(partIn, .true.)
          call printPart(partOut, .false.)
          call writeParticle(6,1,partIn)
          call writeParticle(6,2,partOut)
          call Traceback('Severe problem. Code stops!!!!')
       end if
    end do

    if (present(HiEnergyFlag)) then
       if (HiEnergyFlag) then
          flag = .true.
          return
       end if
    end if

    ! Check baryon number conservation:
    iB_In = 0
    do k=lbound(partIn,dim=1),ubound(partIn,dim=1)
       if (isBaryon(partIn(k)%Id)) then
          if (.not.partIn(k)%anti) then
             iB_In = iB_In + 1
          else
             iB_In = iB_In - 1
          end if
       end if
    end do
    iB_Out = 0
    do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (isBaryon(partOut(k)%Id)) then
          if (.not.partOut(k)%anti) then
             iB_Out = iB_Out + 1
          else
             iB_Out = iB_Out - 1
          end if
       end if
    end do
    if (iB_In /= iB_Out) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)&
               //' : Baryon number conservation'
       else
          write(*,*) 'Problems in collisionTerm: Baryon number conservation'
       end if
       write(*,*) 'Baryon number in: ', iB_In
       write(*,*) 'Baryon number out:', iB_Out
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       call Traceback('stop in finalCheck, baryon number check', errcode)
       if (errcode < 0) return ! --> failure
    end if

    ! Check strangeness conservation:
    iS_In = 0
    do k=lbound(partIn,dim=1),ubound(partIn,dim=1)
       if (.not. isHadron(partIn(k)%Id)) cycle
       if (.not. partIn(k)%anti) then
          iS_In = iS_In + hadron(partIn(k)%Id)%strangeness
       else
          iS_In = iS_In - hadron(partIn(k)%Id)%strangeness
       end if
    end do
    iS_Out = 0
    do k=lbound(partOut,dim=1),ubound(partOut,dim=1)
       if (.not. isHadron(partOut(k)%Id)) cycle
       if (.not. partOut(k)%anti) then
          iS_Out = iS_Out + hadron(partOut(k)%Id)%strangeness
       else
          iS_Out = iS_Out - hadron(partOut(k)%Id)%strangeness
       end if
    end do
    if (iS_In /= iS_Out) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)//' : Strangeness conservation'
       else
          write(*,*) 'Problems in collisionTerm: Strangeness conservation'
       end if
       write(*,*) 'Strangeness in: ', iS_In
       write(*,*) 'Strangeness out:', iS_Out
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       !  Check conservation of srts
       do k=0,3
          mom_In(k) =Sum(partIn%mom(k))
          mom_Out(k)=Sum(partOut%mom(k))
       end do
       write(*,*) ' (Missing mass)^2 :', &
            (mom_In(0)-mom_Out(0))**2 &
            - sum((mom_In(1:3)-mom_Out(1:3))**2)
       call Traceback('stop in finalCheck, strangeness check', errcode)
       if (errcode < 0) return ! --> failure
    end if


    !  Check total charge:
    iQ_In  = Sum(partIn%charge)
    iQ_Out = Sum(partOut%charge)
    if (iQ_In /= iQ_Out) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)//' : Charge conservation'
       else
          write(*,*) 'Problems in collisionTerm: Charge conservation'
       end if
       write(*,*) 'Charge in: ', iQ_In
       write(*,*) 'Charge out:', iQ_Out
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       call Traceback('stop in finalCheck, charge check', errcode)
       if (errcode < 0) return ! --> failure
    end if

    !  if (debug) Print *, 'CollisionTerm: SqrtS=',sqrtS(partIn(1),partIn(2))

    !  Check conservation of 4-momentum:
    do k=0,3
       mom_In(k) =Sum(partIn%mom(k))
       mom_Out(k)=Sum(partOut%mom(k))
    end do

    if (any(abs(mom_In-mom_Out).gt.energyCheck)) then
       if (present(woher)) then
          write(*,*) 'Problems in '//trim(woher)&
               //' : Energy/momentum conservation'
       else
          write(*,*) 'Problems in collisionTerm: Energy/momentum conservation'
       end if
       write(*,'(A,4G13.5)')'Momentum in: ', mom_In
       write(*,'(A,4G13.5)')'Momentum out:', mom_Out
       write(*,*)
       write(*,'(79("*"))')
       call printPart(partIn, .true.)
       call printPart(partOut, .false.)
       call writeParticle(6,1,partIn)
       call writeParticle(6,2,partOut)
       call Traceback('stop in finalCheck, momentum check', errcode)

       if (errcode < 0) call print899
       if (errcode < 0) return ! --> failure
    end if

    flag=.true.

  contains

    !**************************************************************************
    ! for debugging purposes...
    !**************************************************************************
    subroutine printPart(part, in)
      use mediumModule, only: mediumAt
      use dichteDefinition, only: dichte
      use densityModule, only: densityAt
      use output, only: writeParticle_debug
      type(particle), intent(in), dimension(:) :: part
      logical, intent(in) :: in
      integer :: i
      type(dichte) :: density
      do i=lbound(part,dim=1) , ubound(part,dim=1)
         if (part(i)%ID==0) cycle
         if (in) then
            write(*,*) 'Incoming Particle number #', i
         else
            write(*,*) 'Outgoing Particle number #', i
         end if
         density = densityAt(part(i)%pos)
         call writeParticle_debug(part(i),mediumAt(part(i)%pos),&
              density%baryon)
      end do
      write(*,'(79("*"))')

    end subroutine printPart

    !**************************************************************************
    ! for debugging purposes...
    !**************************************************************************
    subroutine print899

      use dichteDefinition
      use densitymodule, only: densityAt
      use lorentzTrafo, only: lorentz
      use potentialMain, only: potential_LRF
      use coulomb, only: emfoca

      integer :: k
      type(dichte) :: dens
      real, dimension(1:3) :: beta
      type(particle) :: partDummy
      real :: pot, potC

      write(899,'(A,4G13.5)')'Momentum in: ', mom_In
      write(899,'(A,4G13.5)')'Momentum out:', mom_Out
      call writeParticle(899,1,partIn)
      call writeParticle(899,2,partOut)

      do k=1,size(partIn)
         if (partIn(k)%ID==0) cycle
         dens = densityAt(partIn(k)%pos)
         if (dens%baryon(0)>1E-8) then
            beta = dens%baryon(1:3)/dens%baryon(0)
            partDummy = partIn(k)
            call lorentz( beta,partDummy%mom)
            pot = potential_LRF(partDummy, dens, addCoulomb=.false.)
         else
            pot = 0.
         end if
         potC = emfoca(partIn(k)%pos, partIn(k)%mom(1:3), &
              partIn(k)%charge, partIn(k)%ID)
         write(899,'(i3,1P,6e13.5)') k,dens%baryon,pot, potC
      end do

      do k=1,size(partOut)
         if (partOut(k)%ID==0) cycle
         dens = densityAt(partOut(k)%pos)
         if (dens%baryon(0)>1E-8) then
            beta = dens%baryon(1:3)/dens%baryon(0)
            partDummy = partOut(k)
            call lorentz( beta,partDummy%mom)
            pot = potential_LRF(partDummy, dens, addCoulomb=.false.)
         else
            pot = 0.
         end if
         potC = emfoca(partOut(k)%pos, partOut(k)%mom(1:3), &
              partOut(k)%charge, partOut(k)%ID)
         write(899,'(i3,1P,6e13.5)') k,dens%baryon,pot,potC
      end do


      flush(899)

    end subroutine print899

  end function finalCheck

end module collisionTools
