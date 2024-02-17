program test
use inputGeneral
use particleProperties, only: initParticleProperties

call readInputGeneral
call initParticleProperties

call tester_rho

end program test


subroutine tester_rho

  use rhoWidthMedium_tables, only : get_GammaColl_rho, readTable, writeTable
  use inMediumWidth_rho, only : GammaColl
  use mesonWidthVacuum, only: vacuumWidth
  use mediumDefinition
  use constants
  use ParticleProperties, only : hadron
  use output, only: intToChar

  implicit none
  
  integer, parameter :: n_mass = 100   ! number of mass points
  integer, parameter :: n_p= 50        ! number of momentum points
  integer, parameter :: n_dens = 30    ! number of density points
  integer, parameter :: n_temp = 10    ! number of temperature points

  real,    parameter :: min_mass=2.*mElec  ! min. mass (GeV)
  real,    parameter :: max_mass=3.0   ! max. mass (GeV)  
  real,    parameter :: p_max = 3.0    ! max. momentum for tabulation (GeV/c)
  real,    parameter :: dens_max = 0.48 ! max. baryon density for tabulation (fm^-3)
  real,    parameter :: temp_max = 0.1 ! max. temperature for tabulation (GeV)

  integer, save :: imode=2  ! 1 -- read partial files GammaColl.103.dat.bz2 produced by tabulate_GammaColl_rho
                            !      and combine them in one file GammaColl.103.dat.bz2
                            ! 2 -- plot collisional width of rho-meson for selected temperatures and densities  

  real :: delta_mass,delta_p,delta_dens,delta_temp,mass,p,dens,temp
  real :: m0,GammaColl1,GammaColl2,GammaVac,width,SpFun
  integer :: idir,index_mass,index_p,index_dens,index_temp  
  type(medium) :: med
  character(200) :: fileName  ! name of bz2-file
                                                      
  select case(imode)
  case(1)
     
     do idir=1,200
        fileName='./'//trim(intToChar(idir))//'/GammaColl.103.dat.bz2'
        call readTable(fileName)
     end do
   
     call writeTable
     
  case(2)

     delta_mass=(max_mass-min_mass)/float(n_mass)
     delta_p=p_max/float(n_p)    
     delta_dens=dens_max/float(n_dens)
     delta_temp=temp_max/float(n_temp)

     m0=hadron(103)%mass
     
     temp = 0.    ! GeV     
     dens = 0.16    ! fm^-3
     
     med%temperature = temp     
     med%density = dens

     Open(57,File='GammaCollRho_T0gev_dens0.16fm-3_p0gev_res_mod_rhoN.dat',status='unknown')
     write(57,*) '# temperature, GeV : ', med%temperature
     write(57,*) '# density, fm^-3 : ', med%density
     write(57,'(A)')'#   p, GeV/c:        mass, GeV:       Gamma_coll (table,calc), GeV:     Gamma_vac, GeV:  SpFun, GeV^-2:'
     open(58,file='GammaCollRes_T0gev_dens0.16fm-3_p0gev_mod_rhoN.dat',status='unknown')
     write(58,*) '# temperature, GeV : ', med%temperature
     write(58,*) '# density, fm^-3 : ', med%density 
     write(58,'(A)')'# p, GeV/c:          mass, GeV:       Gamma_coll_Res (id=2-31), GeV:'     
     do index_p=0,0 !20,10  !n_p
         p=float(index_p)*delta_p
         do index_mass=0,n_mass
             mass=min_mass+float(index_mass)*delta_mass
             GammaColl1=get_GammaColl_rho(mass,p,med)
             GammaColl2=GammaColl(103,mass,p,dens,temp)
             GammaVac=vacuumWidth(mass,103)
             width=GammaColl2+GammaVac
             SpFun=mass*width/pi/((mass**2-m0**2)**2+(mass*width)**2)
             write(57,'(6(4x,e13.7))') p,mass,GammaColl1,GammaColl2,GammaVac,SpFun
         end do
     end do
     close(57)
     close(58)

  end select
  
end subroutine tester_rho




