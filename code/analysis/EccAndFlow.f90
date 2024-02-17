!******************************************************************************
!****m* /EccAndFlow
! NAME
! module EccAndFlow
! PURPOSE
! Encapsulate all routines and datas for calculation off eccentricity
! and flow according the Q- and R-vector method
!
! INPUTS
! ...(This module needs no input)
!
!******************************************************************************
module EccAndFlow

  implicit none
  private

  !****************************************************************************
  !****t* EccAndFlow/tQRvector
  ! NAME
  ! type tQRvector
  ! PURPOSE
  ! Type definition to store all information of a Q- or R-vector
  !
  ! R = r^m exp(i n phi_r)
  ! Q = p^m exp(i n phi_p)
  !
  ! SOURCE
  !
  type, public :: tQRvector
     integer :: n = 1
     integer :: m = 1
     real, dimension(1:2) :: Q = 0  ! the vector
     real :: sumW = 0.      ! sum of weights
     real :: sumN = 0.      ! sum of entries = sum of perweights
  end type tQRvector
  !****************************************************************************

  public :: QRvectorInit
  public :: QRvectorAdd
  public :: QRvectorVal
  public :: QRvectorW



contains

  subroutine QRvectorInit(QR, n, m)
    type(tQRvector) :: QR
    integer, intent(in) :: n
    integer, intent(in) :: m

    QR%n = n
    QR%m = m
    QR%Q = 0
    QR%sumW = 0.
    QR%sumN = 0.
  end subroutine QRvectorInit

  subroutine QRvectorAdd(QR, r, phi, w )
    type(tQRvector) :: QR
    real, intent(in) :: r
    real, intent(in) :: phi
    real, intent(in), OPTIONAL :: w

    real :: w_

    w_ = 1.0
    if (present(w)) w_ = w

    QR%Q = QR%Q + r**QR%m * (/cos(QR%n*phi),sin(QR%n*phi)/) * w_
    QR%sumW = QR%sumW + r**QR%m * w_
    QR%sumN = QR%sumN + w_

  end subroutine QRvectorAdd

  pure function QRvectorVal(QR) result(val)
    type(tQRvector), intent(in) :: QR
    real, dimension(1:2) :: val

    if (abs(QR%sumW)>1e-20) then
       val = QR%Q/QR%sumW
    else
       val = 0.
    end if

  end function QRvectorVal

  pure real function QRvectorW(QR)
    type(tQRvector), intent(in) :: QR

    if (abs(QR%sumN)>1e-20) then
       QRvectorW = QR%sumW/QR%sumN
    else
       QRvectorW = 0.
    end if
  end function QRvectorW

end module EccAndFlow
