!******************************************************************************
!****m* /PIL_freezeout
! NAME
! module PIL_freezeout
! PURPOSE
! Provide some storage method for the freeze out position of some particles
!
! This is closely connected to the module PILIndex
!
! INPUTS
! (none)
!
! NOTES
! * "PIL" stands for "PartInfoList"
!******************************************************************************
module PIL_freezeout

  use PILIndex

  implicit none
  private

  !****************************************************************************
  !****t* PIL_freezeout/tValueEntry
  ! NAME
  ! type tValueEntry
  ! PURPOSE
  ! a container of the information to be stored
  !
  ! NOTES
  ! * The time of the last interaction we get directly from the particle
  !   (lastCollisionTime), but also from th e0-th component of position
  ! * The history field should be identical with that of the particle. It is
  !   used to decide, whether the particle has suffered an elastic interaction
  !   (is this true??? or does an elastic scatter also changes the number?)
  !
  ! SOURCE
  !
  type tValueEntry
     real, dimension(0:3) :: pos = 0
     real, dimension(0:3) :: mom = 0
     real, dimension(0:3) :: dens = 0
     integer :: history = 0
     logical :: escaped = .false.
     real :: pot = 0.
  end type
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_freezeout/IndexList
  ! PURPOSE
  ! The list, were the particle numbers are connected with the (physical)
  ! storage position
  !
  ! SOURCE
  !
  type(tIndexList), save :: IndexList
  !****************************************************************************


  !****************************************************************************
  !****g* PIL_freezeout/ValueList
  ! PURPOSE
  ! The list, were the information is stored
  !
  ! SOURCE
  !
  type(tValueEntry), allocatable, save :: ValueList(:)
  !****************************************************************************


  public :: PIL_freezeout_PUT, PIL_freezeout_GET
  public :: PIL_freezeout_DeAllocate
  public :: PIL_freezeout_ZERO
!   PUBLIC :: PIL_freezeout_Print

contains

  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_DeAlloc
  ! NAME
  ! subroutine PIL_freezeout_DeAlloc
  ! PURPOSE
  ! Deallocate the memory for this list and the corresponding index list.
  !****************************************************************************
  subroutine PIL_freezeout_DeAllocate()
    call PILIndex_DeAllocate(IndexList)
    if (allocated(ValueList)) deallocate(ValueList)
  end subroutine PIL_freezeout_DeAllocate


  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_ZERO
  ! NAME
  ! subroutine PIL_freezeout_ZERO()
  ! PURPOSE
  ! Reset the list by setting the counter of stored information to 0.
  ! No allocation or deallocation of memory happens.
  !****************************************************************************
  subroutine PIL_freezeout_ZERO()
    IndexList%nEntry = 0
  end subroutine PIL_freezeout_ZERO


  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_PUT
  ! NAME
  ! subroutine PIL_freezeout_PUT(number,hist,escaped,r,mom,dens,pot)
  ! PURPOSE
  ! Store the information "r" connected with particle "number" in the list.
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! * real    :: pos    -- the 4-position to store
  ! * real    :: mom    -- the 4-momentum to store
  ! * real    :: dens   -- the 4-density (baryon) to store
  ! * integer :: hist   -- the history of the particle
  ! * logical :: escaped-- flag to indicate, whether particle has 'escaped'
  ! * logical :: pot    -- value of the potential
  ! OUTPUT
  ! none
  !****************************************************************************
  subroutine PIL_freezeout_PUT(number,hist,escaped,pos,mom,dens,pot)
    integer, intent(IN) :: number
    integer, intent(IN) :: hist
    logical, intent(IN) :: escaped
    real, dimension(0:3), intent(IN) :: pos
    real, dimension(0:3), intent(IN) :: mom
    real, dimension(0:3), intent(IN) :: dens
    real, intent(IN)    :: pot

    integer :: iEntry

    iEntry = PILIndex_PUT(IndexList, number,"freezeout")

!    write(*,*) '###PUT: ',number,r,iEntry

    if (iEntry>0) then      ! everything is ok
       ValueList(iEntry)%pos = pos
       ValueList(iEntry)%mom = mom
       ValueList(iEntry)%dens = dens
       ValueList(iEntry)%history = hist
       ValueList(iEntry)%escaped = escaped
       ValueList(iEntry)%pot = pot
    else
       call PIL_freezeout_Allocate() ! do (re)allocate
       ValueList(-iEntry)%pos = pos
       ValueList(-iEntry)%mom = mom
       ValueList(-iEntry)%dens = dens
       ValueList(-iEntry)%history = hist
       ValueList(-iEntry)%escaped = escaped
       ValueList(-iEntry)%pot = pot
    end if
  end subroutine PIL_freezeout_PUT


  !****************************************************************************
  !****f* PIL_freezeout/PIL_freezeout_GET
  ! NAME
  ! logical function PIL_freezeout_GET(number,hist,escaped,pos,mom,dens,pot)
  ! PURPOSE
  ! Get the stored information of particle "number"
  ! INPUTS
  ! * integer :: number -- the (unique) particle number
  ! OUTPUT
  ! * real    :: pos    -- the 4-position to store
  ! * real    :: mom    -- the 4-momentum to store
  ! * real    :: dens   -- the 4-density (baryon) to store
  ! * real    :: time   -- the time the particle was stored
  ! * integer :: hist   -- the history of the particle
  ! * logical :: escaped-- flag to indicate, whether particle has 'escaped'
  ! * logical :: pot    -- value of the potential
  ! * the (logical) return value signals, whether information about
  !   this particle was found in the list or not.
  !****************************************************************************
  logical function PIL_freezeout_GET(number,hist,escaped,pos,mom,dens,pot)
    integer, intent(IN)  :: number
    integer, intent(OUT) :: hist
    logical, intent(OUT) :: escaped
    real, dimension(0:3), intent(OUT), optional :: pos
    real, dimension(0:3), intent(OUT), optional :: mom
    real, dimension(0:3), intent(OUT), optional :: dens
    real, intent(OUT), optional :: pot

    integer :: iEntry

    iEntry = PILIndex_FIND(IndexList,number)

    ! ATTENTION: iEntry is the line of information in the IndexList.
    ! The information connected to the particle is stored in the line
    ! IndexList%Entry(iEntry) in the array ValueList !!!!!

    if (iEntry > 0) then
       if (present(pos)) &
            pos = ValueList(IndexList%Entry(iEntry))%pos
       if (present(mom)) &
            mom = ValueList(IndexList%Entry(iEntry))%mom
       if (present(dens)) &
            dens = ValueList(IndexList%Entry(iEntry))%dens
       hist = ValueList(IndexList%Entry(iEntry))%history
       escaped = ValueList(IndexList%Entry(iEntry))%escaped
       if (present(pot)) &
            pot = ValueList(IndexList%Entry(iEntry))%pot
       PIL_freezeout_GET = .TRUE.
    else
       if (present(pos)) pos = 0.0
       if (present(mom)) mom = 0.0
       if (present(dens)) dens = 0.0
       hist = 0
       escaped = .false.
       if (present(pot)) pot = 0.0
       PIL_freezeout_GET = .FALSE.
    end if

  end function PIL_freezeout_GET


  !****************************************************************************
  !****is* PIL_freezeout/PIL_freezeout_Allocate
  ! NAME
  ! subroutine PIL_freezeout_Allocate
  ! PURPOSE
  ! Do the allocation and reallocation of the value vector.
  ! The new size is taken from the size of the IndexList vectors.
  ! NOTES
  ! For security one should insert here checks, whether the memory allocations
  ! failed and stop execution in these cases.
  !
  ! If the compiler has MOVE_ALLOC(), one could do it as
  !   call allocate(L0(n1))
  !   L0(1:n0) = L(:)
  !   call move_alloc(L0, ValueList)
  ! or
  !   call move_alloc(ValueList, L0)
  !   call allocate(ValueList(n1))
  !   ValueList(1:n0) = L0(:)
  !****************************************************************************
  subroutine PIL_freezeout_Allocate
    integer :: n0, n1,i
    type(tValueEntry), allocatable :: L0(:)

    n1 = size(IndexList%PartNumber) ! new size

    if (.not.allocated(ValueList)) then
       allocate(ValueList(n1))
       return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    do i=1,n0
       L0(i)%pos = ValueList(i)%pos
       L0(i)%mom = ValueList(i)%mom
       L0(i)%dens = ValueList(i)%dens
       L0(i)%history = ValueList(i)%history
       L0(i)%escaped = ValueList(i)%escaped
       L0(i)%pot = ValueList(i)%pot
    end do
    deallocate(ValueList)
    allocate(ValueList(n1))
    do i=1,n0
       ValueList(i)%pos = L0(i)%pos
       ValueList(i)%mom = L0(i)%mom
       ValueList(i)%dens = L0(i)%dens
       ValueList(i)%history = L0(i)%history
       ValueList(i)%escaped = L0(i)%escaped
       ValueList(i)%pot = L0(i)%pot
    end do
    deallocate(L0)

  end subroutine PIL_freezeout_Allocate


  !****************************************************************************
  !****s* PIL_freezeout/PIL_freezeout_Print
  ! NAME
  ! subroutine PIL_freezeout_Print(file)
  ! PURPOSE
  ! Print the list to file
  !****************************************************************************
!   subroutine PIL_freezeout_Print(file)
!     integer, intent(IN) :: file
!
!     integer :: i
!
!     write(file,*) '****** PIL_freezeout:',IndexList%nEntry
!     do i=1,IndexList%nEntry
!        write(file,'(i8.0,i8.0,3g12.5)') i,IndexList%PartNumber(i),ValueList(IndexList%Entry(i))%pos
!     end do
!
!   end subroutine PIL_freezeout_Print



end module PIL_freezeout
