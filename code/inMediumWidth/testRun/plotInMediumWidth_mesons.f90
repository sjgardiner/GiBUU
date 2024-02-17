program plotInMediumWidth_mesons

  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  use version, only: PrintVersion
  use inMediumWidth
  use inMediumWidth_rho
  use mesonWidthMedium_tables

  implicit none

  call PrintVersion

  call readinputGeneral
  call initParticleProperties

  !    ID = 103
  !    ID = 105

  !    fakDens= 0.2
  !    fakDens= 1.0

  !    mom = 0.0
  !    mom = 0.2
  !    mom = 0.5
  !    mom = 1.0

  !  call loop1(103, 0.2, 0.0)

  call loop2(103, 0.2, 0.0)
  call loop2(103, 0.2, 0.2)
  call loop2(103, 0.2, 0.5)
  call loop2(103, 0.2, 1.0)


contains

  subroutine loop1(ID,fakDens,mom)

    use mediumDefinition
    use mesonWidth, only: fullWidthMeson

    integer,intent(in) :: ID
    real, intent(in) :: fakDens,mom

    type(medium) :: med, med1
    real :: mass, width, widthVac
    real :: h, dMom, dMass, dMassVac, dMed
    integer :: iM


    med%densityNeutron = fakDens * 0.168 ! in fm^-3
    med%densityProton  = fakDens * 0.168

    med1%densityNeutron = fakDens * 0.168 *(1.+1e-4) ! in fm^-3
    med1%densityProton  = fakDens * 0.168 *(1.+1e-4)

    !    do iM=10,1000,10
    do iM=10,1000,1
       mass = iM * 1e-3

       width = get_inMediumWidth(ID,mom,mass,med)
       widthVac = fullWidthMeson(ID,mass)

       h = get_inMediumWidth(ID,mom+1e-4,mass,med)
       dMom = (h-width)/(1e-4)

       h = get_inMediumWidth(ID,mom,mass+1e-4,med)
       dMass = (h-width)/(1e-4)

       h = fullWidthMeson(ID,mass+1e-4)
       dMassVac = (h-widthVac)/(1e-4)

       h = get_inMediumWidth(ID,mom,mass,med1)
       dMed = (h-width)/(fakDens*0.168*1e-4)

       write(197,*) mom,mass, med%densityNeutron,med%densityProton, width,widthVac, dMom,dMass,dMed, dMassVac
    end do


  end subroutine loop1


  subroutine loop2(ID,fakDens,mom)

    use mediumDefinition
    use mesonWidth, only: fullWidthMeson
    use particleProperties, only: hadron

    integer,intent(in) :: ID
    real, intent(in) :: fakDens,mom

    type(medium) :: med
    real, dimension(2) :: width
    real :: mass, mPole
    real :: muH, momH, mom0H
    real, dimension(2) :: widthH
    real :: mom0, mu, chi
    real, parameter :: eps = 1e-4
    integer :: iM, iH
    real, dimension(2) :: dMom0, dMom
    real, dimension(0:5) :: ham

    med%densityNeutron = fakDens * 0.168 ! in fm^-3
    med%densityProton  = fakDens * 0.168

    mPole = hadron(ID)%mass

    !    do iM=10,1000,10
    do iM=10,1000,2
       mass = iM * 1e-3

       mu = mass
       width(1) = get_inMediumWidth(ID,mom,mu,med)
       width(2) = fullWidthMeson(ID,mu)

       ! Hamilton:

       do iH=0,5
          chi = iH*1.0
          if (mass < mPole) chi = -chi

          ham(iH) = sqrt(mPole**2 + mom**2 + chi*mu*(width(1)+width(2)))
       end do

       ! dMom0:

       mom0 = sqrt(mass**2+mom**2)
       mom0H = mom0+eps
       muH = sqrt(mom0H**2-mom**2)

       widthH(1) = get_inMediumWidth(ID,mom,muH,med)
       widthH(2) = fullWidthMeson(ID,muH)

       dMom0 = (muH*widthH - mu*width)/eps

       ! dMom:

       momH = mom+eps
       if (momH>mom0) then
          muH = -99.9
       else
          muH = sqrt(mom0**2-momH**2)
       end if

       widthH(1) = get_inMediumWidth(ID,mom,muH,med)
       widthH(2) = fullWidthMeson(ID,muH)

       dMom = (muH*widthH - mu*width)/eps


       write(197,*) mom,mass, med%densityNeutron,med%densityProton, &
            width,dMom0,dMom, ham
    end do

    write(197,*)
    write(197,*)

  end subroutine loop2



end program plotInMediumWidth_mesons
