program test_XsectionRatios

use inputGeneral, only : readInputGeneral
use XsectionRatios, only : Ratio
use RMF, only: getRMF_parSet, walecka
use constants, only: mn,mpi
use particleProperties, only: initParticleProperties, hadron
use IdTable

implicit none

real :: Elab, k_infty2, k2, mass_seed, rhobar, shift, VectorPotential, U
real :: mnstar, mDeltaStar, Estar, E, tmp, srtsstar
real, dimension(1:2) :: XsRatio
integer :: i

character(len=50) :: f

call readInputGeneral
call initParticleProperties

write(f,'(A,i2,A)') 'XsectionRatio_2gev_RMF_set',getRMF_parSet(),'.dat'

Elab = 2. !GeV

k_infty2 = 2.*Elab*mn + Elab**2

open(1,file=trim(f),status='unknown')
write(1,*)'# The ratio of the in-medium and vacuum NN -> NDelta and NN -> NNpi cross sections'
write(1,*)'# Elab, GeV:', Elab
write(1,*)'# k_infty2, GeV^2:', k_infty2
write(1,*)'# rhobar, fm^-3:  sqrt(s^*), GeV:   m_N^*, GeV:       m_Delta^*, GeV:   Ratio_NNND:       Ratio_NNNNpi:'

mass_seed = mn

do i=0,2000
   rhobar = i*0.001
    
   call walecka (rhobar, shift, em0=mass_seed, V=VectorPotential, potential=U)


   k2 = k_infty2 - 2.*mn*(U+VectorPotential/mn*Elab)

   mnstar = mn-shift
   mDeltaStar = hadron(Delta)%mass-shift

   Estar = sqrt(k2+mnstar**2)
   E = Estar + VectorPotential

   tmp = abs(E-Elab-mn)
   if(tmp .gt. 1.e-05) then
      write(*,*) 'Energy nonconservation :', tmp
      stop
   end if

   srtsstar = sqrt((Estar+mnstar)**2 - k2)

   XsRatio(1) = Ratio(srtsstar,(/mnstar,mnstar,mnstar,mDeltaStar/),(/mn,mn,mn,hadron(Delta)%mass/))

   XsRatio(2) = Ratio(srtsstar,(/mnstar,mnstar,mnstar,mnstar,mpi/),(/mn,mn,mn,mn,mpi/))
   
   write(1,'(9(e13.6,5x))') rhobar, srtsstar, mnstar, mDeltaStar, XsRatio

   mass_seed = mn - shift

end do

end program test_XsectionRatios
