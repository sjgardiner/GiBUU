!******************************************************************************
!****m* /RMFflag
! NAME
! module RMFflag
! PURPOSE
! Offer the possibility to check, whether RMF mode and Skyrme mode are
! initialized simultanously. This is a state which should be considered
! to be a bug.
!
! It is necessary to put this feature in a specialised module in order to
! circumvent cyclic dependencies
!******************************************************************************
module RMFflag

  use CallStack, only: TRACEBACK


  implicit none
  private

  integer,save :: RMFflag_status = -1
  integer,save :: EQS_status = -1

  logical, parameter :: dieHard = .true.


  public :: setRMF
  public :: setSkyrme


contains

  !****************************************************************************
  !****s* RMFflag/setRMF
  ! NAME
  ! subroutine setRMF(flag)
  ! PURPOSE
  ! store the value of the RMF_flag.
  ! Throws an error, if Skyrme already initialized.
  !****************************************************************************
  subroutine setRMF(flag)

    logical, intent(in) :: flag

    if (flag) then

       if (dieHard) then
          if (EQS_status .gt. -1) then
             call TRACEBACK("RMF and Skyrme not allowed")
          end if
       else
          if (EQS_status .gt. 0) then
             call TRACEBACK("RMF and Skyrm EQS>0 not allowed")
          end if
       end if

       RMFflag_status = 1

    else
       RMFflag_status = 0
    end if

  end subroutine setRMF

  !****************************************************************************
  !****s* RMFflag/setSkyrme
  ! NAME
  ! subroutine setSkyrme(EQS)
  ! PURPOSE
  ! store the value of EQS_type in Skyrme mode.
  ! Throws an error, if RMF already initialized.
  !****************************************************************************
  subroutine setSkyrme(EQS)

    integer, intent(in) :: EQS

    if (RMFflag_status.eq.1) then
       if (dieHard) then
          call TRACEBACK("RMF and Skyrme not allowed")
       else
          if (EQS>0) then
             call TRACEBACK("RMF and Skyrm EQS>0 not allowed")
          end if
       end if
    end if

    EQS_status = EQS

  end subroutine setSkyrme

end module RMFflag
