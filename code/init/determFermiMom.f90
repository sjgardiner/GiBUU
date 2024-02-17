!******************************************************************************
!****m* /determFermiMom
! NAME
! module determFermiMom
!
! PURPOSE
! determine the fermi momentum
!
! "fmbybE" stands for "fermi momentum by binding energy"
!
!******************************************************************************
module determFermiMom
  use particleDefinition, only: particle
  implicit none

  private

  public :: determine_fermimom
  public :: determine_fermiNucDLDA


  type(particle), save :: teilchen_fermi ! dummy particle

  real, save :: baryondensity ! global variable for root finding

contains

  !****************************************************************************
  !****f* determFermiMom/determine_fermimom
  ! NAME
  ! real function determine_fermimom(rhoB)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! * real :: rhoB -- ???
  ! OUTPUT
  ! * real :: mom -- fermi momentum
  !
  !****************************************************************************
  function determine_fermimom(rhoB) result (mom)

    use constants, only: mN
    use particleDefinition
    use cernlib, only: dzerox

    real, intent(in) :: rhoB
    real :: mom

    logical, save::initflag=.true.

    if (initflag) then

       call setToDefault(teilchen_fermi)
       teilchen_fermi%mass=mN
       ! finalstate%charge =charge_out
       ! teilchen_fermi%mom=p_out
       teilchen_fermi%pos=(/0.,0.,0./)
       teilchen_fermi%anti=.false.
       teilchen_fermi%id=1
       teilchen_fermi%pert=.false.
       ! finalState%prodTime=0.
       ! finalState%lastCollTime=0.
       ! finalState%formTime=0.
       ! finalState%scaleCS=1.
       ! finalState%inF=.false.
       write(*,*) 'determin_fermimomentum initialized'
       initflag=.false.
    end if

    baryondensity=rhoB ! set the external parameter for root finding
    mom = DZEROX(0d0,1d0,1d-9,1000,func,1)

    if (mom .lt. 0.01) then
       write(*,*) 'mom =',mom,' rho=',rhoB
    end if

  end function determine_fermimom

  !****************************************************************************
  !****f* determFermiMom/determine_fermiNucDLDA
  ! NAME
  ! real function determine_fermiNucDLDA(pos,rhocent)
  ! PURPOSE
  ! ...
  ! INPUTS
  ! * real, dimension(1:3) :: pos -- ???
  ! * real                 :: rhocent -- ???
  ! OUTPUT
  ! * real :: mom -- fermi momentum
  !
  !****************************************************************************
  function determine_fermiNucDLDA(pos,rhocent) result (mom)

    use constants, only: hbarc, pi, mN
    use NucDLDA, only: getEParticleLaplace
    use densitymodule, only: gridSpacing
    use baryonPotentialMain, only: rhoLaplace

    real ::  mom
    real, intent(in) :: rhocent
    real, dimension(1:3),intent(in) :: pos

    real, dimension(1:4):: Etemp
    real :: rhoptemp, pFermi2temp, difftemp
    integer, save :: notBound=0
    logical, save :: firsttime=.true.

    rhoptemp=rhoLaplace(pos,gridSpacing)
    pFermi2temp=(1.5*pi**2*rhocent)**(2./3.)*hbarc**2
    call getEParticleLaplace(Etemp,rhocent,rhoptemp,pFermi2temp)
    mom=0.
    if (Etemp(1).ge.0.) then
       difftemp=pFermi2temp-(2.*mN)/(1000.)*Etemp(1)
       if (difftemp.ge.0.) then
          mom=SQRT(difftemp)
          if (firsttime) then
             open(97,file='NotBound.dat')
             firsttime=.false.
          else
             open(97,file='NotBound.dat',position='Append')
          end if
          notBound=notBound+1
          write(97,fmt='(I5,1X,F10.4,1X,F10.4,1X,F10.4)') notBound,pFermi2temp/(2.*mN)*1000. &
               &, difftemp/(2.*mN)*1000.,Etemp(1)
          close(97)
       else
          mom=0.
       end if
    else
       mom=SQRT(pFermi2temp)
    end if

  end function determine_fermiNucDLDA


  !****************************************************************************
  !****if* determFermiMom/func
  ! NAME
  ! double precision function func
  ! PURPOSE
  ! function used to find a zero in determine_fermimomentum
  !****************************************************************************
  double precision function func(p)
    use baryonPotentialMain, only: BaryonPotential
    use mediumDefinition
    real, intent(in):: p
    type(medium)    :: med
    teilchen_fermi%mom(1:3)=(/p,0.,0./)
    med%temperature    =0.
    med%useMedium      =.true.
    med%density        = baryondensity
    med%densityProton  = baryondensity/2
    med%densityNeutron = baryondensity/2
    func=sqrt((teilchen_fermi%mass+BaryonPotential(teilchen_fermi,med,.true.))**2 +p**2)-teilchen_fermi%mass+0.016
    !non relativistic energy  -Binding energy
  end function func

end module determFermiMom
