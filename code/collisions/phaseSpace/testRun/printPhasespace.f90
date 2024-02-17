program printPhaseSpace

  use inputGeneral
  use particleDefinition
  use particleProperties
  use threeBodyPhaseSpace, only: Integrate_3BodyPS
  use version, only: printVersion

  implicit none

  real :: srts0, srts, ps
!  real :: m1=0.5478 ! eta
  real :: m1=0.7826 ! omega
  real :: m2=0.9380 ! p
  real :: m3=0.9380 ! p
  real :: h,h1,h2,deltah
  integer :: i, n=100

  call printVersion
  call readInputGeneral
  call initParticleProperties

  h1 = log(m1+m2+m3)
  h2 = log(m1+m2+m3+1.0)
  deltah = (h2-h1)/n

  do i=1,n
     h = exp(h1+i*deltah)
     ps = Integrate_3BodyPS(h,m1,m2,m3)
     write(101,'(f7.4,1P,e13.4)') h,ps

  end do


end program printPhaseSpace
