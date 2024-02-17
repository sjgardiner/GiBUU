program test
use inputGeneral
use particleProperties, only: initParticleProperties

call readInputGeneral
call initParticleProperties

!call tester_nuc
!call tester_bar
call tester_mes
!call writePhotoProduction

end program test

subroutine tester_nuc
  use baryonWidthMedium_tables, only : dumpTable

  call dumpTable(113,1, -1,1,1,1, 1)
  write(113,*)
  call dumpTable(113,1, 10,-1,3,3, 2)

  call dumpTable(114,1, -1,-1,3,3, 12)




end subroutine tester_nuc


subroutine tester_bar
  use output
  use baryonWidthMedium_tables, only : get_inMediumWidth
!  use IdTable, only : nres
  implicit none
  integer :: i,j,k,particleID
  real,parameter :: rho=0.16
!  real,parameter :: rho=0.001
  real :: rhoN,rhoP,pabs,mass
  real, dimension(0:3) :: momentum

  particleID=1
  do i=1,4
     Open(57,File='WidthBar_'//inttochar(particleID)//'_rho_'//intToChar(i)//'.dat')
     write(*,*) "i",i
     rhoN=rho/4.*float(i)/2.
     rhoP=rho/4.*float(i)/2.
     write(57,*) '#rhoN', rhoN
     write(57,*) '#rhoP', rhoP

     write(57,*) '#rho=',rhoN+rhoP
     write(57,*) '#ID=',particleID
     mass=0.938
     do k=0,100
        pabs=float(k)*0.01
        momentum=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
        write(57,'(5E18.6)') pabs,mass, get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,1)&
             , get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,2), get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,3)
     end do
     close(57)
  end do


  do particleID=4,4
     do i=0,4
        Open(57,File='WidthBar_'//IntToChar(particleID)//'_rho_'//intToChar(i)//'.dat')
        rhoN=rho/4.*float(i)/2.
        rhoP=rho/4.*float(i)/2.
        write(57,*) '#rho=',rhoN+rhoP
        write(57,*) '#ID=',particleID
        do j=0,100
!        write(*,*) "j",j
           mass=0.8+float(j)*0.02
           do k=0,34
              pabs=float(k)*0.03
              momentum=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
              write(57,'(5G18.6)') pabs,mass, get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,1)&
                   & , get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,2),&
                   & get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,3)

           end do
           write(57,*)
        end do
     end do
     close(57)
  end do


end subroutine tester_bar


subroutine tester_mes
  use output
  use mesonWidthMedium_tables, only : get_inMediumWidth,numSteps_absP,numSteps_mass
  use mediumDefinition
  implicit none
  integer :: i,j,k,particleID
  real,parameter :: rho=0.20
  real,    save :: max_absP_Mes = 3.0 ! max. momentum for tabulation
  real, parameter :: min_mass=0.01
  real, parameter :: max_mass=3.0
  real :: rhoN,rhoP,pabs,mass,delta_absP,delta_mass
  type(medium) :: med

  write(*,*)' numSteps_absP : ', numSteps_absP
  write(*,*)' numSteps_mass : ', numSteps_mass

  delta_absP=max_absP_Mes/float(numSteps_absP)
  delta_mass=(max_mass-min_mass)/float(numSteps_mass)

!  particleID=103    ! rho
!  particleID=101   ! pion
  particleID=102   ! eta
!  particleID=104   ! sigmaMeson
!  particleID=105   ! omegaMeson
  do i=4,4 !1,4
     !     Open(57,File='WidthMes_'//inttochar(particleID)//'_rho_'//intToChar(i)//'.dat')
     !     Open(57,File='WidthMes_'//inttochar(particleID)//'_pion_'//intToChar(i)//'.dat')
          Open(57,File='WidthMes_'//inttochar(particleID)//'_eta_'//intToChar(i)//'.dat')
     !Open(57,File='WidthMes_'//inttochar(particleID)//'_sigma_'//intToChar(i)//'.dat')
     !Open(57,File='WidthMes_'//inttochar(particleID)//'_omega_'//intToChar(i)//'.dat')
     write(*,*) "i",i
     rhoN=rho/4.*float(i)/2.
     rhoP=rho/4.*float(i)/2.
     med%densityNeutron = rhoN
     med%densityProton = rhoP

     write(57,*) '#rhoN', rhoN
     write(57,*) '#rhoP', rhoP

     write(57,*) '#rho=',rhoN+rhoP
     write(57,*) '#ID=',particleID

     do j=0,numSteps_mass
         mass=min_mass+float(j)*delta_mass
         do k=0,numSteps_absP
             pabs=float(k)*delta_absP
             write(57,'(3E18.6)') pabs,mass,get_inMediumWidth(particleID,pabs,mass,med)

         end do
     end do
  end do
  close(57)

end subroutine tester_mes

subroutine writePhotoProduction
  use output
  use baryonWidthMedium_tables, only : get_inMediumWidth
!  use IdTable, only : nres
  implicit none
  integer :: i,j,k,particleID
  real,parameter :: rho=0.16
!  real,parameter :: rho=0.001
  real :: rhoN,rhoP,pabs,mass
  real, dimension(0:3) :: momentum
  real :: mN

  do particleID=4,4
     do i=0,4
        Open(57,File='PhotoProd_'//IntToChar(particleID)//'_rho_'//intToChar(i)//'.dat')
        rhoN=rho/4.*float(i)/2.
        rhoP=rho/4.*float(i)/2.
!!$        write(57,*) '#rho=',rhoN+rhoP
!!$        write(57,*) '#ID=',particleID

        mN = 0.938

        do j=100,300
           !        write(*,*) "j",j
           !           mass=0.8+float(j)*0.02
           mass = j*0.01

           if (mass<mN) cycle
           pabs = (mass**2-mN**2)/(2*mN)
           momentum=(/sqrt(mass**2+pabs**2),pabs,0.,0./)
           write(57,'(5G18.6)') pabs,mass, &
                get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,1),&
                get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,2),&
                get_inMediumWidth(particleID,momentum,mass,rhoN,rhoP,3)

        end do
     end do
     close(57)
  end do

end subroutine writePhotoProduction
