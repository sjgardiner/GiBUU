program plot_E_A

  use skyrme, only: U, E_div_A
  use particleProperties, only: initParticleProperties
  use inputGeneral
  use constants, only : pi, hbarc

  implicit none

  !****************************************************************************
  ! Parameters of potential:
  ! Units: alpha [MeV], beta [MeV], c [MeV], lambda [1/fm], eta [fm^4]
  ! 'c' and 'lambda' determine the momemtum dependence.
  ! 'eta' concerns the surface term (cf. SurfacePotFlag)
  real, parameter, dimension(1:14) :: &
        alpha  = (/ -108.619,  -9.972, -286.99993 , -124.2611 , -29.253, 0., 0., 0.,  -22.9 ,    5.72 ,   56.7 , -209.2, -135.2 , -124.8  /), &
        beta   = (/  136.779,  38.032,  233.6517  ,   71.03007,  57.248, 0., 0., 0.,   33.9 ,   13.6  ,    4.75,  156.4,  176.8 ,  204.8 /), &
        tau    = (/    1.259,   2.404,    1.225707,    2.00108,   1.76 , 0., 0., 0.,    1.82,    2.61 ,    3.69,   1.35,    1.22,    1.2/), &
        lambda = (/    2.13 ,   2.126,    0.      ,    0.     ,   2.13 , 0., 0., 0.,    0.96,    0.807,    1.18,    0. ,    2.92,    4.0/), &
        c      = (/  -63.601, -63.601,    0.      ,    0.     , -63.516, 0., 0., 0., -102.  , -144.   , -145.  ,    0. ,  -61.23,  -77.29 /), &
        eta    = (/    0.4  ,   0.3  ,    0.4     ,    0.27   ,   0.36 , 0., 0., 0.,    0.  ,    0.   ,    0.  ,    0. , 0. , 0. /)
  !****************************************************************************


  integer, parameter :: iSet = 5
  integer :: iRho
  real, parameter :: rho0 = 0.168*hbarc**3
  real, parameter :: mN = 0.938
  real :: rho, pF, U0, U1, mStar0, mStar1, EA

  do iRho=0,150
     rho = iRho*0.01*rho0
     pF=(3./2.*pi**2*rho)**(1./3.)

     U0 = U(rho,0.,rho0, alpha(iSet)/1e3,beta(iSet)/1e3,c(iSet)/1e3,tau(iSet),lambda(iSet)*hbarc)
     U1 = U(rho,pF,rho0, alpha(iSet)/1e3,beta(iSet)/1e3,c(iSet)/1e3,tau(iSet),lambda(iSet)*hbarc)
     EA = E_div_A(rho,rho0, alpha(iSet)/1e3,beta(iSet)/1e3,c(iSet)/1e3,tau(iSet),lambda(iSet)*hbarc)

     mStar0 = U0 + mN
     mStar1 = sqrt( (U1+sqrt(mN**2+pF**2))**2-pF**2 )

     write(*,*) rho/hbarc**3, EA, U0, U1, mStar0,mStar1, sqrt( (U1+sqrt(1.535**2+pF**2))**2-pF**2 )

  end do






end program plot_E_A
