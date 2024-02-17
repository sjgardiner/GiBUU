program plotInMediumWidth_rho

  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  use version, only: PrintVersion
  use inMediumWidth
  use inMediumWidth_rho
  use rhoWidthMedium_tables

  implicit none

  call PrintVersion

  call readinputGeneral
  call initParticleProperties

  call loop1

contains

  subroutine loop1

    use mediumDefinition

    integer :: ID
    type(medium) :: med
    real :: mom, mass, width
    integer :: iM

    ID = 103

    med%densityNeutron = 0.2 * 0.168 ! in fm^-3
    med%densityProton  = 0.2 * 0.168
    med%density = med%densityNeutron+med%densityProton
    med%temperature = 0.100

    mom = 0.5

    do iM=10,1000,10
       mass = iM * 1e-3

       width = get_GammaColl_rho(mass,mom,med)

       write(*,*) mom,mass, med%densityNeutron,med%densityProton, width
    end do


  end subroutine loop1

end program plotInMediumWidth_rho
