!******************************************************************************
!****m* /createAgAgSmear
! NAME
! module createAgAgSmear
! PURPOSE
! This module creates a random generator for smearing the momenta according
! a 2D array
!******************************************************************************
module createAgAgSmear

  implicit none

  private

  real(4), dimension(0:200,0:300) ::  Arr

  public :: fillArr
  public :: integrateArr
  public :: printLine
  public :: bisectLine
  public :: dumpArr


contains

  subroutine setBinContent(pIn,ratio,val)
    real, intent(in) :: pIn,ratio
    integer, intent(in) :: val

    integer :: ix,iy
    real(4) :: val4


    ix = int(pIn/5.0)+1
    iy = int(ratio/0.01)+1
    val4 = val*1.0

    if ((ix>300).or.(iy>200)) then
       write(*,*) "Error: ",pIn,ratio,ix,iy
       stop
    end if

!    write(*,*) pIn,ratio,ix,iy,val4
    Arr(iy,ix) = val4

  end subroutine setBinContent

  subroutine fillArr
    integer :: iLine
    integer,parameter :: nLine = 60000

    real :: pIn,ratio
    integer :: val

    open(10,file="Mom_Smearing_Data.txt",status="OLD")

    Arr = 0.
    do iLine=1,nLine
       read(10,*) pIn,ratio,val
       call setBinContent(pIn,ratio,val)
    end do

    close(10)

  end subroutine fillArr

  subroutine integrateArr
    integer :: iLine,i
    real(4) :: s

    do iLine=1,300
       s = 0.
       do i=0,200
          s = s + Arr(i,iLine)
          Arr(i,iLine) = s
       end do
       Arr(:,iLine) = Arr(:,iLine) * (1./s)
    end do
  end subroutine integrateArr

  subroutine printLine(iLine)
    integer, intent(in) :: iLine

    integer :: i

    do i=0,200
       write(*,*) iLine*5.0,i*0.01,Arr(i,iLine)
    end do

  end subroutine printLine

  integer function bisectLine(iLine, val8)
    integer, intent(in) :: iLine
    real, intent(in) :: val8

    integer :: i1,iM,i2
    real(4) :: val

    val=val8
    i1 = 0
    i2 = 200

    if (val < Arr(i1,iLine)) then
       write(*,*) 'Error 1'
       stop
    end if
    if (val > Arr(i2,iLine)) then
       write(*,*) 'Error 2'
       stop
    end if

    do while(i2-i1 > 1)
       iM = (i1+i2)/2
       if (val > Arr(iM,iLine)) then
          i1 = iM
       else
          i2 = iM
       end if
    end do

    bisectLine = i1

  end function bisectLine

  subroutine dumpArr

    open(10,file="Mom_Smearing_Data.bin",status="UNKNOWN",form="UNFORMATTED")
    write(10) Arr
    close(10)

  end subroutine dumpArr


end module createAgAgSmear


program useAgAgSmear

  use createAgAgSmear, only: fillArr, integrateArr, printLine, bisectLine, &
       dumpArr

  implicit none

  call fillArr
!  call printLine(300)
  call integrateArr
!  call printLine(200)

  call dumpArr


  ! some test examples:
  write(*,*) bisectLine(200, 0.1)
  write(*,*) bisectLine(200, 0.2)
  write(*,*) bisectLine(200, 0.5)
  write(*,*) bisectLine(200, 0.7)
  write(*,*) bisectLine(200, 0.9)



  ! How to generate a random number:
!!$  r = ran()
!!$  randomIndex = bisectLine(iLine, r)


end program useAgAgSmear
