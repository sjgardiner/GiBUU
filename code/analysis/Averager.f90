!******************************************************************************
!****m* /Averager
! NAME
! module Averager
! PURPOSE
! Encapsulate all routines and datas for a statistical 'averager'
!
! Calculates mean, deviation etc of a statistical variable
!
! INPUTS
! ...(This module needs no input)
!
!******************************************************************************
module Averager

  implicit none
  private

  !****************************************************************************
  !****t* Averager/tAverager
  ! NAME
  ! type Averager
  ! PURPOSE
  ! Type definition to store all information of an averager
  ! SOURCE
  !
  type, public :: tAverager
     integer :: count = 0 ! number of entries counted
     real,dimension(0:4) :: sum = 0.0 ! sum of entries to different powers
  end type tAverager
  !****************************************************************************

  public :: AveragerClear
  public :: AveragerAdd
  public :: AveragerCount
  public :: AveragerNorm
  public :: AveragerMean
  public :: AveragerVariance
  public :: AveragerStdDev
  public :: AveragerStdErr
  public :: AveragerCentralMom
  public :: AveragerWrite


contains

  !****************************************************************************
  !****s* Averager/AveragerClear
  ! elemental subroutine AveragerClear(A)
  ! PURPOSE
  ! reset the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to clear
  ! OUTPUT
  ! * type(tAverager) :: A --- the Averager to clear
  !****************************************************************************
  elemental subroutine AveragerClear(A)

    type(tAverager), intent(inout) :: A

    A%count = 0
    A%sum = 0.0

  end subroutine AveragerClear

  !****************************************************************************
  !****s* Averager/AveragerAdd
  ! elemental subroutine AveragerAdd(A, val, weight)
  ! PURPOSE
  ! add a value with some weight to the averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  ! * real :: val --- the value to add
  ! * real, OPTIONAL :: weight --- te weight of the value
  ! OUTPUT
  ! * type(tAverager) :: A --- the Averager to use
  !****************************************************************************
  elemental subroutine AveragerAdd(A, val, weight)

    type(tAverager), intent(inout) :: A
    real, intent(in) :: val
    real, intent(in), optional :: weight

    real :: w

    w = 1.0
    if (present(weight)) w = weight

    A%count = A%count + 1
    A%sum = A%sum + w * (/ 1.0, val, val**2, val**3, val**4 /)

  end subroutine AveragerAdd

  !****************************************************************************
  !****f* Averager/AveragerCount
  ! elemental integer function AveragerCount(A)
  ! PURPOSE
  ! return the count of the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  !****************************************************************************
  elemental integer function AveragerCount(A)
    type(tAverager), intent(in) :: A
    AveragerCount = A%count
  end function AveragerCount

  !****************************************************************************
  !****f* Averager/AveragerNorm
  ! elemental real function AveragerNorm(A)
  ! PURPOSE
  ! return the normalization (i.e. sum of weights) of the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  !****************************************************************************
  elemental real function AveragerNorm(A)
    type(tAverager), intent(in) :: A
    AveragerNorm = A%sum(0)
  end function AveragerNorm

  !****************************************************************************
  !****f* Averager/AveragerMean
  ! elemental real function AveragerMean(A)
  ! PURPOSE
  ! return the average value of the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  !****************************************************************************
  elemental real function AveragerMean(A)
    type(tAverager), intent(in) :: A

    if (A%Sum(0) > 0) then
       AveragerMean = A%Sum(1)/A%Sum(0)
    else
       AveragerMean = 0
    end if

  end function AveragerMean

  !****************************************************************************
  !****f* Averager/AveragerVariance
  ! elemental real function AveragerVariance(A)
  ! PURPOSE
  ! return the variance (i.e. the second central moment) of the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  ! NOTES
  ! * this is just a shortcut to AveragerCentralMom(A, 2)
  !****************************************************************************
  elemental real function AveragerVariance(A)
    type(tAverager), intent(in) :: A
    AveragerVariance = AveragerCentralMom(A, 2)
  end function AveragerVariance

  !****************************************************************************
  !****f* Averager/AveragerStdDev
  ! elemental real function AveragerStdDev(A)
  ! PURPOSE
  ! return the standard deviation (i.e. the square root of the variance) of
  ! the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  ! NOTES
  ! * this is just a shortcut to sqrt(AveragerVariance)
  !****************************************************************************
  elemental real function AveragerStdDev(A)
    type(tAverager), intent(in) :: A
    AveragerStdDev = sqrt(AveragerVariance(A))
  end function AveragerStdDev

  !****************************************************************************
  !****f* Averager/AveragerStdErr
  ! elemental real function AveragerStdErr(A)
  ! PURPOSE
  ! return the standard error of the given Averager
  !
  ! The standard error of the arithmetic average is given by the standard
  ! deviation, divided by the square root of the nomber of entries
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  !****************************************************************************
  elemental real function AveragerStdErr(A)
    type(tAverager), intent(in) :: A

    if (A%Sum(0) > 0) then
       AveragerStdErr = sqrt(AveragerVariance(A)/A%Sum(0))
    else
       AveragerStdErr = 0
    end if
  end function AveragerStdErr

  !****************************************************************************
  !****f* Averager/AveragerCentralMom
  ! elemental real function AveragerCentralMom(A, i)
  ! PURPOSE
  ! return the ith central moment of the given Averager
  ! INPUTS
  ! * type(tAverager) :: A --- the Averager to use
  ! * integer :: i --- the moment to give
  !****************************************************************************
  elemental real function AveragerCentralMom(A, i)
    type(tAverager), intent(in) :: A
    integer, intent(in) :: i

    real, dimension(0:4) :: h

    if (A%Sum(0) <= 0) then
       AveragerCentralMom = 0.0
       return
    end if

    h = A%Sum/A%Sum(0) ! the moments around origin

    select case (i)
    case (0,1)
       AveragerCentralMom = 0.0
    case (2)
       AveragerCentralMom = h(2) - h(1)**2
    case (3)
       AveragerCentralMom = h(3) - 3*h(1)*h(2) + 2*h(1)**3
    case (4)
       AveragerCentralMom = h(4) - 4*h(1)*h(3) + 6*h(1)**2*h(2) - 3*h(1)**4
    end select

  end function AveragerCentralMom

  !****************************************************************************
  !****s* Averager/AveragerWrite
  ! subroutine AveragerWrite(iFile, num, A)
  ! PURPOSE
  ! write the given Averager to some file
  ! INPUTS
  ! * integer :: iFile --- file number for output
  ! * real :: num --- some number to prepend to the output
  ! * type(tAverager) :: A --- the Averager to clear
  !****************************************************************************
  subroutine AveragerWrite(iFile, num, A)

    integer, intent(in) :: ifile
    real, intent(in) :: num
    type(tAverager), intent(in) :: A

    write(iFile,'(f9.3,i7,1P,5e13.4)') num, &
         AveragerCount(A), AveragerNorm(A), AveragerMean(A), &
         AveragerCentralMom(A,2), &
         AveragerCentralMom(A,3), &
         AveragerCentralMom(A,4)
  end subroutine AveragerWrite

end module Averager
