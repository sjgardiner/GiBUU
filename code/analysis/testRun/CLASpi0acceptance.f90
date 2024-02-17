program CLASpi0acceptance

  use inputGeneral
  use particleProperties
  use clebschGordan
  use hist
  use hist2D

  implicit none



  call readInputGeneral
  call initParticleProperties

  call doCalcMC2

contains

  subroutine doCalcMC()

    use random, only: rnFlat
    use constants, only: pi

    use minkowski, only: abs4

    real, parameter :: mPi = 0.135 ! GeV, maybe not the default value (0.138)
    real, parameter :: pmax=4.0 ! GeV
    real, dimension(0:3) :: mom1, mom2, mompi

    real :: cost1, cost2, sint1, sint2, phi, p1, p2, b, pmin
    real :: costpi, ppi, w0,w1

    integer, parameter :: nMC = 1000
    integer :: iMC

    real, parameter :: costMax = 0.98901586336191682981 ! 8.5°
    real, parameter :: costMin = 0.7071067811865475244 ! 45°


    do iMC=1,nMC
       cost1 = rnFlat(costMin,costMax)
       cost2 = rnFlat(costMin,costMax)
       phi   = rnFlat(0.0,pi)

       sint1 = sqrt(1.0-cost1**2)
       sint2 = sqrt(1.0-cost2**2)

       b = 2*(1.0 - sint1*sint2*cos(phi) - cost1*cost2)
!       write(*,*) b

       if (abs(b) .lt. 1e-5) then
          ! ...
       end if

       if (b .lt. 0) then
          write(*,*) b
       end if

       pmin = mPi**2/(b * pmax)
       if (pmin.gt.pmax) cycle

       p1 = rnFlat(pmin,pmax)
       p2 = mPi**2/(b * p1)

       mom1 = p1 * (/ 1.0, sint1, 0.0, cost1 /)
       mom2 = p2 * (/ 1.0, sint2*cos(phi), sint2*sin(phi), cost2 /)

       mompi = mom1+mom2

!       write(*,*) abs4(mompi)


       ppi = sqrt(Dot_Product(mompi(1:3),mompi(1:3)))
       costpi = mompi(3)/ppi

       w0 = 1.0/(b**2 * p1)
       w1 = w0 * AccGamma(p1,cost1) * AccGamma(p2,cost2)

!       write(*,*) p1,p2,ppi, w0,w1
       write(*,*) ppi,costpi, w0,w1

    end do

  end subroutine doCalcMC



  real function AccGamma(p,cost)

    use constants, only: pi

    real, intent(in) :: p,cost
    real :: theta ! theta in radian

    theta = acos(cost)

!    write(*,*) cost,theta,theta*180./pi

    AccGamma = 0.0
    if (p .lt. 0.3) return
    if (p .gt. 4.0) return
    if (theta .lt. 8.5/180.0*pi) return
    if (theta .gt. 45.0/180.0*pi) return
    AccGamma = 1.0

  end function AccGamma

  subroutine doCalcMC2()

    use random, only: rnFlat
    use constants, only: pi

    use minkowski, only: abs4

    real, parameter :: mPi = 0.135 ! GeV, maybe not the default value (0.138)
    real, dimension(0:3) :: mom1, mom2, mompi

    real :: cost1, cost2, sint1, sint2, phi, p1, p2, a,b
    real :: costpi, ppi, w0,w1

    integer, parameter :: nMC = 10000000
    integer :: iMC

!!$    real, parameter :: costMax = 0.98901586336191682981 ! 8.5°
!!$    real, parameter :: costMin = 0.7071067811865475244 ! 45°
!!$
!!$    real, parameter :: pmin=0.3 ! GeV
!!$    real, parameter :: pmax=4.0 ! GeV

    real, parameter :: costMax = 1.0
    real, parameter :: costMin = -1.0

    real, parameter :: pmin=0.0 ! GeV
    real, parameter :: pmax=5.0 ! GeV

    type(histogram) :: hP, hTh
    type(histogram2D) :: hPTh

    real :: addFak, mulFak

    call createHist(hP,"P",0.0,5.0,0.1)
    call createHist(hTh,"theta",0.0,90.0,1.0)
    call createHist2D(hPTh, "p vs. theta", (/0.0,0.0/), (/5.0,90.0/), (/0.1,1.0/))

    iMC = 0
!    do iMC=1,nMC
    do
       cost1 = rnFlat(costMin,costMax)
       cost2 = rnFlat(costMin,costMax)

       sint1 = sqrt(1.0-cost1**2)
       sint2 = sqrt(1.0-cost2**2)

       p1 = rnFlat(pmin,pmax)
       p2 = rnFlat(pmin,pmax)

       a = mpi**2 - 2*p1*p2*(1.-cost1*cost2)
       b = 2*p1*p2*sint1*sint2


       if (b<-a) cycle
!       write(*,*) a,b,b>-a

       if (abs(a).ge.abs(b)) cycle

       phi = acos(-a/b)

       mom1 = p1 * (/ 1.0, sint1, 0.0, cost1 /)
       mom2 = p2 * (/ 1.0, sint2*cos(phi), sint2*sin(phi), cost2 /)

       mompi = mom1+mom2

!       write(*,*) abs4(mompi),phi,-a/b

       ppi = sqrt(Dot_Product(mompi(1:3),mompi(1:3)))
       costpi = mompi(3)/ppi

       w0 = p1*p2/sqrt(b**2-a**2)
       w1 = w0 * AccGamma(p1,cost1) * AccGamma(p2,cost2)

!       write(*,*) p1,p2,ppi, w0,w1
!       write(*,*) ppi,costpi, w0,w1
!       write(*,*) ppi,acos(costpi)*180./pi, w0,w1

       call addHist(hP,ppi,w0,w1)
       call addHist(hTh,acos(costpi)*180./pi,w0,w1)
       call addHist2D(hPTh,(/ppi,acos(costpi)*180./pi/),w0,w1)

       if (w1>0) then
          iMC = iMC+1
          if (iMC>nMC) exit
       end if
    end do

    addFak = 1e-20
    !    mulFak = 1.0/nMC
    mulFak = 1.0

    call writeHist(hP, add=addFak,mul=mulFak, doAve=.true., file="hP.dat")
    call writeHist(hTh, add=addFak,mul=mulFak, doAve=.true., file="hTh.dat")
    call writeHist2D_Gnuplot(hPTh, add=addFak,mul=mulFak, doAve=.true., file="hPTh.dat")

  end subroutine doCalcMC2

end program CLASpi0acceptance
