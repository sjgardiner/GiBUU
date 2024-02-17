!******************************************************************************
!****m* /HadesAgAgSmear
! NAME
! module HadesAgAgSmear
! PURPOSE
! This module implements a random generator for smearing the momenta according
! a 2D array.
! This array may be generated with testRun/createAgAgSmear.f90 (Ag+Ag@1.58GeV)
! or testRun/createPPSmear.f90 (p+p@1.58 GeV)
!
! Input files are generated on basis of private communications of:
! * Ag+Ag: Jan-Hendrik Otto
! * p+p: Karina Scharmann
!******************************************************************************
module HadesAgAgSmear

  use callstack, only: Traceback

  implicit none

  private

  integer, parameter :: nRatio = 200
  integer, parameter :: npIn = 300

  real, parameter :: dRatio = 0.01
  real, parameter :: dpIn = 0.005 ! in GeV

  real(4), dimension(0:nRatio,0:npIn), save ::  Arr
  logical, save :: initFlag = .true.

  logical, parameter :: doExtrapolate = .true.

  character(1000), save :: filename
  ! possibilities are :
  ! * Mom_Smearing_Data.bin
  ! * Mom_Smearing_Data_PP.bin

  public :: doSmear
  public :: setFileName

contains

  subroutine setFileName(s)
    use inputGeneral, only: path_To_Input, ExpandPath

    character*(*), intent(in) :: s

    filename = trim(s)

    if (len_trim(filename)>0) then
       if (index(filename,"/")>0) then
          call ExpandPath(filename)
          filename = trim(filename)
       else
          filename = trim(path_to_Input)//'/hades/' &
               //trim(filename)
       end if
       write(*,*) 'hadesSmearFile : ', trim(filename)
    end if


  end subroutine setFileName

  subroutine fetchArr

    use inputGeneral, only: path_To_Input
    use output, only: Write_ReadingInput

    integer :: ios

    call Write_ReadingInput(trim(filename),0)
    open(1013,file=filename,status="OLD",form="UNFORMATTED",iostat=ios)
    if (ios.ne.0) call Traceback("File does not exist.")
    read(1013,iostat=ios) Arr
    if (ios.ne.0) call Traceback("Error while reading file.")
    close(1013)
    call Write_ReadingInput(trim(filename),1)

    initFlag = .false.
  end subroutine fetchArr

  integer function bisectLine(iLine, val8)
    integer, intent(in) :: iLine
    real, intent(in) :: val8

    integer :: i1,iM,i2
    real(4) :: val

    val=val8
    i1 = 0
    i2 = nRatio

    if (val < Arr(i1,iLine)) call Traceback("Error 1")
    if (val > Arr(i2,iLine)) call Traceback("Error 2")

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

  real function getRandomRatio(pIn)

    use random, only: rn

    real, intent(in) :: pIn

    integer :: ipIn, iRatio
    real :: r, h

    if (initFlag) call fetchArr

    ipIn = int(pIn/dpIn)+1
    if (ipIn > npIn) then
       if (doExtrapolate) then
          ipIn = npIn
       else
          call Traceback("pIn too large")
       end if
    end if

    r = rn()
    iRatio = bisectLine(ipIn, r)

    if (iRatio == nRatio) then
       getRandomRatio = nRatio*dRatio
       return
    end if

    ! now do a linear interpolation:
    h = (r-Arr(iRatio,ipIn))/(Arr(iRatio+1,ipIn)-Arr(iRatio,ipIn))
    getRandomRatio = (iRatio+h)*dRatio

  end function getRandomRatio

  subroutine doSmear(p)
    use minkowski, only: abs3

    real, dimension(0:3), intent(inout) :: p

    real :: r, pIn

!    write(*,*) 'in:  ',p

    pIn = abs3(p)
    r = getRandomRatio(pIn)
    !m2 = p(0)**2-pIn**2
    !E2new = m2 + (r*pIn)**2
    p(1:3) = p(1:3)*r
    p(0) = sqrt(p(0)**2-(1.-r**2)*pIn**2)

!    write(*,*) 'out: ',p,r

  end subroutine doSmear

end module HadesAgAgSmear
