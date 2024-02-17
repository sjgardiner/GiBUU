program writeHist2Dintegrals

  use hist
  use hist2D

  implicit none

  type(histogram2D) :: H2D
  type(histogram) :: HX, HY

  character*(*),parameter :: fName = "initHiLep.nuQ2.000.dat"
  real :: add, mul
  logical :: flagOK

  call FetchHist2D(H2D, fName//".bin", add=add, mul=mul, flagOK=flagOK)
  if (.not.flagOK) then
     write(*,*) 'read not okay.'
     stop
  end if

  call IntegrateHist2D(H2D, HX, 1)
  call IntegrateHist2D(H2D, HY, 2)

  call writeHist(HX,add=add,mul=mul,file=fname//".X.dat")
  call writeHist(HY,add=add,mul=mul,file=fname//".Y.dat")


end program writeHist2Dintegrals
