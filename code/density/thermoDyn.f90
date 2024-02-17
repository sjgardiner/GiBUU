!******************************************************************************
!****m* /thermoDyn
! NAME
! module thermoDyn
! PURPOSE
! Includes routines for the determination of temperature T and chemical potential mu.
!******************************************************************************
module thermoDyn
  implicit none
  private

  real,dimension(:,:,:),Allocatable :: pSquared
  real,dimension(:,:,:),Allocatable :: temperature


  !****************************************************************************
  !****g*  thermoDyn/temperatureSwitch
  ! SOURCE
  !
  integer,save :: temperatureSwitch = 1
  !
  ! PURPOSE
  ! * 1=groundstate calculations (T=0,mu=E_F)
  ! * 2=the full procedure
  !****************************************************************************


  !****************************************************************************
  !****g*  thermoDyn/linearExtrapolation
  ! SOURCE
  !
  logical,save :: linearExtrapolation=.true.
  !
  ! PURPOSE
  ! * .true.= Use linear extrapolation for temperature between gridPoints
  ! * .false.= Do not use it
  !****************************************************************************


  logical, save :: initFlag=.true.
  logical,parameter :: debugFlag=.false.

  public :: updateTemperature, temperatureAt, muAt, cleanUp

Contains

  subroutine init

    use output, only: Write_ReadingInput
    use densityModule, only: getGridPoints

    integer :: ios
    integer, dimension (1:3) :: gridPoints

    !**************************************************************************
    !****n*  thermoDyn/initThermoDynamics
    ! NAME
    ! NAMELIST initThermoDynamics
    ! PURPOSE
    ! Includes the input switches:
    ! * temperatureSwitch
    ! * linearExtrapolation
    !**************************************************************************
    NAMELIST /initThermoDynamics/ temperatureSwitch, linearExtrapolation

    call Write_ReadingInput('initThermoDynamics',0)
    rewind(5)
    read(5,nml=initThermoDynamics,iostat=ios)
    call Write_ReadingInput('initThermoDynamics',0,ios)

    write(*,*) 'Set temperature switch to ',temperatureSwitch,'.'
    write(*,*) 'Set linear extrapolation switch to ',linearExtrapolation,'.'
    call Write_ReadingInput('initThermoDynamics',1)

    if (temperatureSwitch == 2) then
      gridPoints = getGridPoints()
      allocate(pSquared(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),-gridPoints(3):gridPoints(3)))
      allocate(temperature(-gridPoints(1):gridPoints(1),-gridPoints(2):gridPoints(2),-gridPoints(3):gridPoints(3)))
    end if

    initFlag = .false.

  end subroutine init


  subroutine cleanUp
    if (temperatureSwitch == 2) then
       if (allocated(psquared))    deallocate(pSquared)
       if (allocated(temperature)) deallocate(temperature)
    end if
  end subroutine


  !****************************************************************************
  !****s*  thermoDyn/upDateTemperature
  ! PURPOSE
  ! Calculation of local 'temperature'.
  !****************************************************************************
  subroutine upDateTemperature (teilchen)

    use particleDefinition
    use output, only: DoPR

    type(particle),dimension(:,:),intent(in) :: teilchen

    if (initFlag) call init

    select case (temperatureSwitch)
    case (1)
       return
    case (2)
       if (DoPr(2)) write(*,*) 'Updating Temperature'
       ! determine pSquared(x,y,z) according to particle vector
       call evaluatePSquared
       ! now determine temperature from <p**2>
       call evaluateTemperature
    case default
       write(*,*) 'This temperatureSwitch is not valid:', temperatureSwitch
       Stop
    end select

  contains

    !**************************************************************************
    !****s*  upDateTemperature/evaluatePSquared
    ! PURPOSE
    ! Calculation of pSquared(x,y,z).
    ! This is stored in the field psquared(-gridPoints(1:3):gridPoints(1:3)).
    ! Usage of the smearing weigths to smear particles' momentum over
    ! gridpoints in their neighborhood.
    !**************************************************************************
    subroutine evaluatePSquared

      use IdTable, only: isBaryon
      use densitymodule, only: gridSpacing, gridPoints, nLargePoints, nSmallPoints, smearingWeights, boostToLRF

      type(particle) :: baryon

      integer :: i,j !loop variables over particle vector
      integer :: small,large,Index1,Index2,Index3
      integer, dimension(1:3) :: posOrig, indexSmall
      real :: momSquare

      pSquared(:,:,:)=0.

      do i=1,Size(Teilchen,dim=2)    ! Loop over all Particles in Ensemble
         do j=1,Size(Teilchen,dim=1) ! Loop over all Ensembles

            ! Count only baryons, but no antibaryons :
            if ((.not. isBaryon(teilchen(j,i)%ID)) .or. teilchen(j,i)%anti) cycle

            !position in large grid:
            posOrig=NINT(Teilchen(j,i)%pos(1:3)/gridSpacing)
            !position in small grid:
            indexSmall=NINT((Teilchen(j,i)%pos(1:3)/gridSpacing-posOrig)*(2.*nSmallPoints+1.))

            ! Test for errors:
            if ((abs(indexSmall(1)).gt.nSmallPoints)&
                 &.or.(abs(indexSmall(2)).gt.nSmallPoints)&
                 &.or.(abs(indexSmall(3)).gt.nSmallPoints)) then
               write(*,*) 'Problem in updateTemperature, module density.f90'
               write(*,*) IndexSmall, 'too big, choose different grid'
               stop
            end if

            small=1+(nSmallPoints+indexSmall(3))                               &
                 &  +(nSmallPoints+indexSmall(2))*(2*nSmallPoints+1)       &
                 &  +(nSmallPoints+indexSmall(1))*(2*nSmallPoints+1)**2
            large=0

            !Smearing particle over points in Neighborhood:
            do Index1=posOrig(1)-nLargePoints,posOrig(1)+nLargePoints
               do Index2=posOrig(2)-nLargePoints,posOrig(2)+nLargePoints
                  do Index3=posOrig(3)-nLargePoints,posOrig(3)+nLargePoints
                     large=large+1
                     if ((abs(Index1).le.gridPoints(1)) &
                          &.and.(abs(Index2).le.gridPoints(2)) &
                          &.and.(abs(Index3).le.gridPoints(3)) &
                          &.and.(smearingweights(small,large).gt.0.0)) then  !Point is inside Grid
                        !Evaluate pSquared in local restframe
                        baryon=teilchen(j,i)
                        call boostToLRF(baryon,1)
                        momSquare=Dot_Product(baryon%mom(1:3),baryon%mom(1:3))
                        pSquared(Index1,Index2,Index3)=pSquared(Index1,Index2,Index3) &
                             &   +momSquare*smearingWeights(small,large)
                     end if
                  end do !loop over Index1
               end do !loop over Index2
            end do !!loop over Index3
         end do !"End do" over all particles in one ensemble
      end do !"End do" over ensembles

    end subroutine evaluatePSquared


    !**************************************************************************
    !****s*  upDateTemperature/evaluateTemperature
    ! PURPOSE
    ! Calculation of temperature(x,y,z).
    ! This is stored in the field temperature(-gridPoints(1:3):gridPoints(1:3)).
    ! Use pSquared(x,y,z) and rho(x,y,z) to determine the temperature.
    !**************************************************************************
    subroutine evaluateTemperature

      use densitymodule, only: gridspacing, gridpoints, densField
      use minkowski, only: abs4

      integer :: ix,iy,iz
      real :: rhoLRF, pSquaredAverage  ! density in LRF and <p**2>

      temperature(:,:,:)=0.

      do ix=-gridPoints(1),gridPoints(1)
         do iy=-gridPoints(2),gridPoints(2)
            do iz=-gridPoints(3),gridPoints(3)
               if (pSquared(ix,iy,iz)>0.) then
                  if (densField(ix,iy,iz)%baryon(0)>1.e-02) then
                     pSquaredAverage=pSquared(ix,iy,iz)/densField(ix,iy,iz)%baryon(0)
                     rhoLRF=abs4(densField(ix,iy,iz)%baryon)
                     temperature(ix,iy,iz)=tempe(pSquaredAverage,rhoLRF)
                  end if
               end if
            end do
         end do
      end do


      if (debugFlag) then
         open(10,file='temperaturesReal.dat')
         open(20,File='TemperatureZAxis.dat')
         do ix=-gridPoints(1),gridPoints(1)
            do iy=-gridPoints(2),gridPoints(2)
               do iz=-gridPoints(3),gridPoints(3)
                  write(10,'(4F12.4)') ix*Gridspacing(1),iy*Gridspacing(2),iz*Gridspacing(3),temperature(ix,iy,iz)
                  if (ix==0 .and. iy==0) write(20,'(2F12.4)') iz*Gridspacing(3),temperature(ix,iy,iz)
               end do
            end do
         end do
         close(10)
         close(20)
         write(*,*) 'Stopping after evaluateTemperature due to debugFlag=.true.!'
         stop
      end if

    end subroutine evaluateTemperature

  end subroutine updateTemperature


  !****************************************************************************
  !****f*  thermoDyn/temperatureAt
  ! PURPOSE
  ! Evaluates temperature at some space point r. Therefore it uses the values
  ! which are stored in the field temperature(-gridPoints(1:3):gridPoints(1:3)).
  !****************************************************************************
  real function temperatureAt (r)
    use densitymodule, only: gridSpacing, gridSize, gridPoints

    real, intent(in), dimension(1:3)    :: r         !position where density should be calculated
    integer,dimension(1:3) :: indizes

    if (initFlag) call init

    select case (temperatureSwitch)
    case (1) ! Assume ground state
       temperatureAt = 0.
    case (2) !Dynamic temperature according to test-particle density
      if ((abs(r(1)).gt.gridSize(1)).or.(abs(r(2)).gt.gridSize(2)).or.(abs(r(3)).gt.gridSize(3))) then
        !outside grid
        temperatureAt = 0.
      else if (linearExtrapolation) then
        temperatureAt = extrapolate(r)
      else !No extrapolation
        indizes = Int(r/gridSpacing)
        temperatureAt = temperature(indizes(1),indizes(2),indizes(3))
      end if
    case default
      write(*,*) 'This temperature switch is not valid:',temperatureSwitch
      write(*,*) 'Stop BUU'
      stop
    end select

  contains

    real function extrapolate(r)
      !For linear extrapolation between gridpoints
      real, dimension(1:3),intent(in) :: r !position where temperature should be calculated

      integer, dimension(1:3) :: lowIndex
      integer :: i,j,k
      integer,dimension(0:7,1:3) :: grid  !field with all corners of 3D-box where r is situated in
      real,dimension(1:3) :: gridPos      !position of gridpoints
      real :: factor

      !The point "r" where the temperature should be calculated
      !is sitting in a 3D-box with lowest corner "LowIndex" on the
      !density grid. We first construct this point. Then all other corners of the 3D-box
      !are constructed.
      !In the end a simple linear extrapolation is used to make the density smooth
      !inside the box.

      ! (1.) Construct Lowest lying point (most negative point!)
      LowIndex=Int(r/gridSpacing)!Always chooses integers which are closer to the origin
      do i=1,3 !Convert them to integers which are more negative!
         if (r(i)/gridsize(i).lt.0) then
            LowIndex(i)=LowIndex(i)-1
         end if
      end do

      ! (2.) Define GridPoints on the corners of the 3D box in which the point "r" is situated
      do i=0,1
         do j=0,1
            do k=0,1
               grid(i+2*j+4*k,1:3)=LowIndex+(/i,j,k/)
            end do
         end do
      end do

      ! (3.) Do linear extrapolation
      extrapolate=0.
      do i=0,7
         gridPos(1:3)=grid(i,1:3)*gridSpacing(1:3) !position of grid point
         factor=1.  !evaluate weight for linear extrapolation for each grid point
         do j=1,3
            factor=factor*Abs(gridSpacing(j)-abs(r(j)-gridPos(j)))/gridspacing(j)
         end do
         if ((abs(grid(i,1)).le.gridPoints(1)).and.(Abs(grid(i,2)).le.gridPoints(2))&
              & .and.(abs(grid(i,3)).le.gridPoints(3))) then
            extrapolate=extrapolate+(temperature(grid(i,1),grid(i,2),grid(i,3))*factor)
         end if
      end do
    end function extrapolate

  end function temperatureAt


  !****************************************************************************
  !****f*  thermoDyn/tempe
  ! NAME
  ! real function tempe (p2av, rho)
  ! INPUTS
  ! * real, intent(in) :: p2av     ! average momentum squared <p**2> in GeV**2
  ! * real, intent(in) :: rho      ! baryon density in fm^-3
  ! OUTPUT
  ! * temperature in GeV
  ! PURPOSE
  ! Evaluates temperature as a function of rho(0) and <p**2>. At first call it
  ! initializes a field "temSave" which holds this information. Therefore no
  ! further calculation is necessary after initializing this field. It is
  ! generated in the subroutine "initTempe".
  !****************************************************************************
  real function tempe (p2av, rho)
    use constants, only: rhoNull

    real, intent(in) :: p2av, rho

    integer :: index_rho,index_p2av
    real, parameter :: delta_rho  = 0.1      ! stepsize for rho in units of rho0
    real, parameter :: delta_p2av = 0.002     ! stepsize for <p**2>
    integer, parameter :: maxIndex_rho  = 60
    integer, parameter :: maxIndex_p2av = 25000
    real, save, allocatable :: temSave(:,:)  ! Field where results are stored
    logical, save :: init=.true.   ! Checks wether temSave is initialized by "initTempe"

    if (init) call initTempe

    index_rho=nint(rho/delta_rho/rhoNull)
    index_p2av=nint(p2av/delta_p2av)
    if (index_rho>maxIndex_rho) then
       write(*,*) 'Problem in tempe: density higher than expected:', rho,  maxIndex_rho*delta_rho*rhoNull
       index_rho=maxIndex_rho
    end if
    if (index_p2av>maxIndex_p2av) then
       write(*,*) 'Problem in tempe: <p**2> higher than expected:',p2av,maxIndex_p2av*delta_p2av
       index_p2av=maxIndex_p2av
    end if

    if (index_rho==0) then
       tempe=0.
    else
       tempe=temSave(index_rho,index_p2av)
    end if

  contains

    subroutine initTempe
      use constants, only: hbarc

      integer :: i,j
      real, parameter :: delta_temperature = 0.0002  ! exactness of temperature in GeV
      real :: p2average ! <p**2>
      real :: temp      ! temperature
      real :: mu        ! chemical potential
      real :: rho       ! density
      real :: p2        ! <p**2>

      write(*,*) '**In initTempe'

      allocate(temSave(1:maxIndex_rho,0:maxIndex_p2av))

      do i=1,maxIndex_rho  ! loop over density
         rho=float(i)*delta_rho*rhoNull*hbarc**3
         temp=0.
         do j=0,maxIndex_p2av  ! loop over <p**2>
            p2average=float(j)*delta_p2av
            do  ! find temperature
               mu=muAt(rho,temp)
               p2=integral(temp,mu,2)/rho
               if (p2>p2average) then
                  temSave(i,j)=temp
                  exit
               else  ! raise temperature
                  temp=temp+delta_temperature
               end if
            end do
         end do
      end do

      if (debugFlag) then
         open(10,File='Temperature.dat')
         write(10,*) '# rho[fm^-3],<p^2>/rho,temperature[GeV]'
         do i=1,maxIndex_Rho
            do j=0,maxIndex_p2av
               write(10,*) i*delta_Rho*rhoNull,j*delta_p2av,temsave(i,j)
            end do
         end do
         close(10)
      end if

      init=.false.

    end subroutine initTempe

  end function tempe


  !****************************************************************************
  !****f*  thermoDyn/muAt
  ! PURPOSE
  ! Determine chemical potential as function of temperature and rho.
  ! Everything in GeV!!!
  ! INPUTS
  ! real :: rho ,  density in GEV**3
  ! real :: temp ,  temperature in GeV
  !****************************************************************************
  real function muAt (rho, temp)

    use constants, only: pi,mn

    real, intent(in) :: rho, temp

    real  :: pf,mum(2),R,f(3),mu,aux
    real, parameter :: eps=1e-04 !Exactness of result for mu
    logical :: flag
    integer :: iter
    integer, parameter :: imode=2   ! 1 - bisection method , 2 - Newton method

    if(rho.lt.1.e-06) then
       muAt=0.
       return
    end if

    pf=(1.5*pi**2*rho)**(1./3.)
    mum(1)=sqrt(pf**2+mn**2)

    if(temp.lt.1.e-03) then
       muAt=mum(1)
       return
    end if

    R=1./rho**0.333333   ! interparticle distance (GeV^-1)

    if(R .gt. 2.*pi/sqrt(3.*mn*temp)) then  ! Boltzmann limit
       aux=integral(temp,mn,4)
       if(aux.gt.1.e-06) then
          mum(1)=temp*log(rho/aux)+mn
       else
          write(*,*)'In muAt -- rho, temp, aux : ', rho, temp, aux
          stop
       end if
!       write(*,*)' Boltzmann limit, mu : ', mum(1)
    end if

    select case(imode)

    case(1)

       !Starting values:
       mum(1)=-5.
       mum(2)=5.
       f(1)=integral(temp,mum(1),1)-rho
       f(2)=integral(temp,mum(2),1)-rho
       if (f(1)*f(2).gt.0) then
          write(*,*) 'problems with f(1),f(2)',f(1),f(2),rho,temp
       end if

       flag=.true.
       do while(flag)
          mu=(mum(1)+mum(2))/2.
          f(3)=integral(temp,mu,1)-rho

          if (f(3)*f(1).gt.0) then
             mum(1)=mu
             f(1)=f(3)
          else
             mum(2)=mu
             f(2)=f(3)
          end if
          if (mum(2)-mum(1).lt.eps) flag=.false.
       end do
       muAt=(mum(2)+mum(1))/2.

    case(2)

       iter=0
       do
          iter=iter+1
          f(1)=integral(temp,mum(1),1)-rho

!          write(*,*)' iter, mum(1), f(1) : ', iter, mum(1), f(1)

          if(abs(f(1)).lt.1.e-07) exit

          f(3)=integral(temp,mum(1),3)

          if(abs(f(3)).gt.1.e-06) mum(2)=mum(1)-f(1)/f(3)

          mum(1)=mum(2)

          if(iter.eq.20) then
             write(*,*)' exit after 20 iterations, mum(1), f(1) : ', mum(1), integral(temp,mum(1),1)-rho
             write(*,*)' rho, temp : ', rho, temp
             exit
          end if

       end do

       muAt=mum(1)

    end select


  end function muAt


  !****************************************************************************
  !****f*  thermoDyn/integral
  ! NAME
  ! real function integral(temperature,mu,Switch)
  ! PURPOSE
  ! Evaluates for a gas of degenerate neutrons and protons of given "temperature" and
  ! given chemical potential "mu" the following:
  ! Switch=1 : rho=<1>=4*Integral 1/(1+exp((E(p)-mu)/T)) (dp)**3)/(2pi)**3   over the full p-space.
  ! Switch=2 : <p**2>=4*Integral p**2/(1+exp((E(p)-mu)/T)) (dp)**3)/(2pi)**3  over the full p-space.
  ! Switch=3 : \partial\rho/\partial\mu = (4/T)*Integral exp((E(p)-mu)/T)/(1+exp((E(p)-mu)/T))**2 (dp)**3)/(2pi)**3   over the full p-space.
  ! Switch=4 : rho=4*Integral exp((-E(p)+mu)/T) (dp)**3)/(2pi)**3  over the full p-space. Boltzmann limit.
  ! The first Factor of 4 is due to spin&isospin degeneracy.
  !****************************************************************************
  real function integral (temperature, mu, Switch)

    use constants, only: mN, pi

    real,    intent(in) :: temperature, mu
    integer, intent(in) :: Switch

    real, parameter ::  dp=0.01 ! approximate stepsize for integral
    real :: dp2              ! real stepsize
    real :: pmax             ! cut-off for integral in p-space
    real :: eps              ! energy

    integer np !number of points for integral
    integer i
    real :: p, a, b, value

    ! Cut-off for integral
    pmax=max(sqrt(max((10*temperature+max(0.,mu))**2-mN**2,0.)),0.1)

!    write(*,*)' pmax : ', pmax

    !Number of points for Riemann-Integral
    np=int(pmax/dp)+1

    !Delta(p) for Integral
    dp2=pmax/float(np)

    integral=0
    do i=1,np
       p=(float(i)-0.5)*dp2

       eps=sqrt(p**2+mN**2) !Energy

       if (temperature.gt.0.) then
          a=(eps-mu)/temperature
          if (a.lt.80) then !Cut off for Integral, other contribution can be neglected
             b=exp(a)
             select case(Switch)
                case(1)
                   value=p**2/(b+1.)
                case(2)
                   value=p**4/(b+1.)
                case(3)
                   value=p**2*b/(b+1.)**2/temperature
                case(4)
                   value=p**2/b
             end select
          else
             value=0.
          end if

       else if (temperature.eq.0.) then
          if (eps.gt.mu) then !Above fermi-surface
             value=0.
          else               !Below fermi-surface
             select case(Switch)
                case(1)
                   value=p**2
                case(2)
                   value=p**4
                case(3)
                   value=0.
             end select
          end if

       else
          write(*,*) 'temperature negative in subroutine "integral":', temperature
          write(*,*) 'Critical Error. Stop !'
       end if
       integral=integral+value*dp2
    end do
    ! Multiply by 4*4pi/(2pi**3)=2./pi**2 : isospin*spin*raumwinkel/(2pi)**3
    integral=integral*2./pi**2

  end function integral


end module thermoDyn
