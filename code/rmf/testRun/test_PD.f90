program test_PD

  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties
  use RMF, only: getRMF_parSet, PD, m_nucleon, m_minus, PD_rhoB, PD_2
  use constants, only : pi, hbarc, f_pi

  implicit none

  character(len=50) :: f
  integer :: numSet

  call readInputGeneral
  call initParticleProperties

  numSet=getRMF_parSet()

  if(numSet.le.30) then
     write(*,*)' wrong RMF set for the parity doublet model : ',  numSet
     stop
  end if

  write(f,'(A,i2,A)') 'RMF_set',numSet,'.dat'
  open(1,file=trim(f),status='unknown')
  write(1,*)'# Parity doublet model, RMF set : ', numSet
  write(1,*)'# m_nucleon, GeV : ', m_nucleon
  write(1,*)'# m_minus, GeV : ', m_minus
  write(1,*)'#################################'

!  call run1
!  call run2
  call run3

contains

  subroutine run1

    real :: mubStar, mubar, sigma_seed, rhoNuc, rhoPart, sigma, shift(2), rhos, energyDensity, press
    real :: ScalarPotential(2), VectorPotential(2), U(2)
    integer :: i

    write(1,*)'# Col. No.:    quantity:'
    write(1,*)'# 1            mubStar, GeV'
    write(1,*)'# 2            mubar, GeV'
    write(1,*)'# 3            rhoNuc, fm^-3'
    write(1,*)'# 4            rhoPart, fm^-3'
    write(1,*)'# 5            rhos, fm^-3'
    write(1,*)'# 6            sigma, GeV'
    write(1,*)'# 7            m^*_nuc, GeV'
    write(1,*)'# 8            m^*_part, GeV'
    write(1,*)'# 9            endens-endens_vac, GeV/fm^3'
    write(1,*)'# 10           P-P_vac, GeV/fm^3'
    write(1,*)'# 11           S_nuc, GeV'
    write(1,*)'# 12           V_nuc, GeV'
    write(1,*)'# 13           U_nuc, GeV'
    write(1,*)'# 14           S_part, GeV'
    write(1,*)'# 15           V_part, GeV'
    write(1,*)'# 16           U_part, GeV'

    sigma_seed = 0.12 !0.01  !f_pi  !0.  !0.8



!    muBStar=0.935
    muBStar=0.900
    call PD(mubStar,sigma,shift,flagPlot=.true.,sigma_inp=sigma_seed,mub=mubar, &
         rhoPlus=rhoNuc,rhoMinus=rhoPart,rhoscalar=rhos,&
         endens=energyDensity,pressure=press,S=ScalarPotential,V=VectorPotential,potential=U)
!!$    stop

    do i=0,400
       mubStar = 0.9086 +  i*0.001
       !   mubStar = 0.924 - i*0.0001
       !   muBStar = 0.825 +  i*0.001
       !   mubStar = m_nucleon - i*0.001

       mubStar = 0.935 -  i*0.001

       call PD(mubStar,sigma,shift,flagPlot=.false.,sigma_inp=sigma_seed,mub=mubar, &
            rhoPlus=rhoNuc,rhoMinus=rhoPart,rhoscalar=rhos,&
            endens=energyDensity,pressure=press,S=ScalarPotential,V=VectorPotential,potential=U)

       sigma_seed = sigma

       write(1,'(1P,16(e13.6,1x))') mubStar, mubar, rhoNuc, rhoPart, rhos, sigma, &
            m_nucleon-shift(1), m_minus-shift(2),&
            energyDensity, press,&
            ScalarPotential(1), VectorPotential(1), U(1),&
            ScalarPotential(2), VectorPotential(2), U(2)
       flush(1)

    end do

    close(1)

  end subroutine run1

  subroutine run2

    integer :: i
    real :: rhoB

    do i=0,168
       rhoB = i*0.001
       call PD_rhoB(rhoB)
    end do



  end subroutine run2

  subroutine run3

    integer :: i
    real :: mubStar

    do i=0,400
    !    do i=0,20
       !   mubStar = 0.9086 +  i*0.001
       !   mubStar = 0.924 - i*0.0001
       !   muBStar = 0.825 +  i*0.001
       !   mubStar = m_nucleon - i*0.001

       mubStar = 0.945 -  i*0.001
!       mubStar = 0.935 -  i*0.001
!       mubStar = 0.900 -  i*0.001

       call PD_2(mubStar)
    end do
  end subroutine run3





end program test_PD
