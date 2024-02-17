program calcBranching

  use inputGeneral
  use IDTable, only: D13_1520,S11_1535
  use constants, only: pi
  use particleProperties, only: InitParticleProperties,nDecays,hadron
  use decayChannels
  use baryonWidthVacuum, only: vacuumWidth
  use CallStack, only: Traceback

  implicit none

  real :: gammaTotal
  real, dimension(nDecays) :: ratio,rho_AB_Pole,partialG,sPartialG
  real :: sSpectral

  integer :: ID, iM
  real :: M,mass0,spectral,spectral2


  call readInputGeneral
  call InitParticleProperties

!  ID = D13_1520
  ID = S11_1535
  mass0 = hadron(id)%mass

  sPartialG = 0
  sSpectral = 0

  do iM=100,300
     M = iM*0.01

     gammaTotal = vacuumWidth(M,ID,ratio,rho_AB_Pole)
     partialG = gammaTotal*ratio

     spectral = 2./pi * M**2 * gammaTotal / ((M**2-mass0**2)**2+M**2*gammaTotal**2)
     spectral2 = 2./pi * M**2 * gammaTotal / ((M**2-(mass0-0.25)**2)**2+M**2*gammaTotal**2)
     sPartialG = sPartialG + partialG*spectral*0.01
     sSpectral = sSpectral + spectral*0.01
     write(1001,*) M,ratio,gammaTotal,spectral,spectral2
     write(1002,*) M,sPartialG/sSpectral,sSpectral
     write(1003,*) M,sPartialG/sum(sPartialG),sum(sPartialG)
  end do


end program calcBranching
