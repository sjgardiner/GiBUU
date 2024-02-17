program test
  use inputGeneral

  implicit none

  write(*,*) '**************************************************'
  write(*,*) 'Initializing database'
  call readinputGeneral
  call init_Database
  write(*,*) '**************************************************'

  write(*,*) '**************************************************'
  write(*,*) 'Testing Rho Delta -> '
  write(*,*) '**************************************************'

  call testRhoDelta


  contains


  subroutine testRhoDelta
    use IDTABLE
    use mediumDefinition
    use particleDefinition
    use particleProperties
    use rhoDelta_resonance
    use preEventDefinition
    implicit none
    real                                          :: srts                  ! sqrt(s) in the process
    type(particle),dimension(1:2)   :: teilchenIn        ! colliding particles
    type(medium)                      :: mediumATcollision    ! Medium informations at the position of the collision

    logical                         :: plotFlag          ! Switch on plotting of the  Xsections
    real,dimension(0:3)       :: momentumLRF        ! Total Momentum in LRF
    type(preEvent),dimension(1:3) :: teilchenOut     ! colliding particles
    real                            :: sigmaTot         ! total Xsection
    real                            :: sigmaElast      ! elastic Xsecti
    integer :: i,chargeMeson,chargeNuk, k
    real :: plab,ekin
    real :: piN,rhoN, pipiN,resonanz
    real :: dens
    integer :: numTries=10000

    NAMELIST /initTest/ chargeMeson, chargeNuk, dens
    write (*,*)
    write (*,*) '**Initializing testing parameter'
    write(*,*)  ' Reading input ....'
    rewind(5)
    read(5,nml=initTest)
    write(*,*) ' Set charge to ', chargeMeson,'.'
    write(*,*) ' Set nuk charge to ',  chargeNuk,'.'
    write(*,*) ' Set density to ',  dens,'.'

    mediumAtCollision%useMedium=.true.
    mediumAtCollision%densityProton=dens/2.
    mediumAtCollision%densityNeutron=dens/2.

    teilchenIN(1)%Id=rho
    teilchenIN(2)%Id=Delta
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(2)%charge=chargeNuk
    teilchenIN(1)%mass=meson(rho)%mass-0.1
    teilchenIN(2)%mass=baryon(delta)%mass
    teilchenIN(2)%mom=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%vel=(/0.,0.,0./)

    Open(301,file='RhoDeltaXsections.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,250
       plab=i*0.01
!    do i=0,300
!       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )
       teilchenIN(1)%mom=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%vel=(/teilchenIN(1)%mom(1)/teilchenIN(1)%mom(0),0.,0./)

       srts=SQRT((SQRT(teilchenIN(1)%mass**2+plab**2)+baryon(Delta)%mass)**2-plab**2)
       momentumLRF=(/SQRT(teilchenIN(1)%mass**2+plab**2)+baryon(Delta)%mass ,plab,0.,0./)
       ekin=SQRT(teilchenIN(1)%mass**2+plab**2)-teilchenIN(1)%mass
       call rhoDelta(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.2,.true.)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

       Do k=lBOund(teilchenOut,dim=1),uBound(teilchenOut,dim=1)
          If ( (teilchenOut(k) %ID.eq.Delta).and.(teilchenOut(k)%charge.lt.0)) then
             Print *, teilchenOut%charge
             Print *, teilchenOut%ID
             stop 'nukcharge'
          end if
       End do
    end do
    close(301)

    piN=0.
    rhoN=0.
    pipiN=0.
    Do i=1,numTries
       plab=1.0
       srts=SQRT((SQRT(teilchenIN(1)%mass**2+plab**2)+baryon(Delta)%mass)**2-plab**2)
       teilchenIN(1)%mom=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%vel=(/teilchenIN(1)%mom(1)/teilchenIN(1)%mom(0),0.,0./)

       momentumLRF=(/SQRT(teilchenIN(1)%mass**2+plab**2)+baryon(Delta)%mass ,plab,0.,0./)
       call rhoDelta(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)
       If ((teilchenOut(1)%ID.eq.rho).and.(teilchenOut(2)%ID.eq.Delta)) then
          rhoN=rhoN+sigmatot
       else If (teilchenOut(1)%ID.eq.S11_1535)then
          resonanz=Resonanz+sigmatot
       end if
    end do
    write(*,*) 'Xsections at plab= 1.0 GeV'
    write(*,*) 'Simga for pi P11(1440). :', rhoN/float(numTries)
    write(*,*) 'Simga for S11_1535 prod. :', resonanz/float(numTries)


    teilchenIN(1)%Id=rho
    teilchenIN(2)%Id=Delta
    teilchenIN(1)%charge=chargeMeson
    teilchenIN(2)%charge=-chargeNuk
    teilchenIN(2)%anti=.true.
    teilchenIN(1)%mass=meson(rho)%mass
    teilchenIN(2)%mass=baryon(delta)%mass
    teilchenIN(2)%mom=(/teilchenIN(2)%mass,0.,0.,0./)
    teilchenIN(2)%vel=(/0.,0.,0./)

    Open(301,file='RhoDeltaXsections_anti.dat',status='unknown')
    write(301,*) '# nukCharge=',chargeNuk,'    mesonCharge=', chargeMeson
    Do i=1,250
       plab=i*0.01
!    do i=0,300
!       plab = 0.5*exp( (log(500.)-log(0.5))*i/300. )
       teilchenIN(1)%mom=(/sqrt(teilchenIN(1)%mass**2+plab**2),plab,0.,0./)
       teilchenIN(1)%vel=(/teilchenIN(1)%mom(1)/teilchenIN(1)%mom(0),0.,0./)

       srts=SQRT((SQRT(meson(rho)%mass**2+plab**2)+baryon(Delta)%mass)**2-plab**2)
       momentumLRF=(/SQRT(meson(rho)%mass**2+plab**2)+baryon(Delta)%mass ,plab,0.,0./)
       ekin=SQRT(meson(rho)%mass**2+plab**2)-meson(rho)%mass
       call rhoDelta(srts,teilchenIN,mediumATcollision,momentumLRF,teilchenOUT,sigmaTot,sigmaElast,.true.,2.3,.false.)

       write(301,'(5F8.3)') plab,sigmaTot,sigmaElast

       Do k=lBOund(teilchenOut,dim=1),uBound(teilchenOut,dim=1)
          If ( (teilchenOut(k) %ID.eq.Delta).and.(teilchenOut(k)%charge.gt.0)) then
             Print *, teilchenOut%charge
             Print *, teilchenOut%ID
             stop 'nukcharge'
          end if
       End do
    end do
    close(301)



  end subroutine testRhoDelta

end program test
