!
! compile this with: make ARGS="-fopenmp"
!
program ThermalModel2

  use, intrinsic :: iso_c_binding
  use inputGeneral
  use particleProperties
  use IdTable, only: pion, nmes
  use CallStack, only: Traceback
  use omp_lib

  implicit none

  real, parameter :: massMax = 3.0 ! upper boundary of integration

  real, parameter :: dMass = 0.2
  integer, parameter :: nMass = (3.3-1.0)/dMass     ! +0.3 for security
  real, dimension(2,0:nMass) :: arrTabBoltzmannPot,arrTabBoltzmannPotSign


  ! final values (after root finding):
  real :: out_T = 0.0       ! temperature
  real :: out_muB = 0.0     ! baryo chemical potential
  real :: out_lambda = 0.0  ! overall fugacity

  real :: out_eps = 0.0    ! energy density
  real :: out_rhoB = 0.0   ! (total) baryon density
  real :: out_rho  = 0.0   ! particle density

  ! input values
  real :: in_T = 0.200     ! temperature
  real :: in_muB = 0.01    ! baryo chemical potential
  real :: in_lambda = 1.0  ! overall fugacity

  real :: in_eps = 0.519  ! energy density
  real :: in_rhoB = 0.0   ! (total) baryon density
  real :: in_rho  = 0.1   ! particle density

  integer :: modus = 0 ! 0: do not run solver, 1: run solver2, 2: run solver3

  real :: eps_part = 0.

  !****************************************************************************
  !****g* ThermalModel2/FileNameHydro
  ! PURPOSE
  ! The absolute filename of the file containing hydro info
  !
  ! possible values:
  ! * if not set, default is '[path_To_Input]/???.dat'
  ! * if given, but does not contain '/':
  !   default is '[path_To_Input]/[FileNameHydro]'
  ! * otherwise: filename is absolute, including path
  !
  ! NOTE
  ! if you want to use the file 'XXX.dat' in the actual directory,
  ! give it as './XXX.dat'
  !
  ! SOURCE
  !
  character(300), save :: FileNameHydro = ''
  !****************************************************************************

  integer, save :: colEps = 23
  integer, save :: colRho = 24
  integer, save :: colN0 = 12
  integer, save :: colB0 = 21
  integer, save :: colMax = 35


  call readInputGeneral
  call initParticleProperties
  call readInput

!$OMP PARALLEL
  !$ write(*,*) 'Thread: ',OMP_GET_THREAD_NUM()
!$OMP END PARALLEL


!!$  out_T=1.0
!!$  do
!!$     write(*,*) out_T, log(out_T)
!!$     out_T = out_T/10
!!$  end do

  select case(modus)
  case (0)
     call calcENB(in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho, eps_part)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho, eps_part


  case (1)
     call runSolver2a
     call calcENB(out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho, eps_part)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, in_eps,in_rhoB,in_rho
     write(*,*) 'OUT:',out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho, eps_part

  case (2)
     call runSolver2b
     call calcENB(out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho, eps_part)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, in_eps,in_rhoB,in_rho
     write(*,*) 'OUT:',out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho, eps_part


  case (3)
     call runSolver3
     call calcENB(out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho, eps_part)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, in_eps,in_rhoB,in_rho
     write(*,*) 'OUT:',out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho, eps_part

  case (10)
     call runMode10

  case (11)
     call Traceback("nyi")

  case (12)
     call runMode12

  case (99)
     call runMode99

  end select

contains

  !************************************************************************
  !****s* ThermalModel2/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "Thermalmodel2".
  !************************************************************************
  subroutine readInput
    use output, only: Write_ReadingInput
    use inputGeneral, only: path_To_Input, ExpandPath

    NAMELIST /ThermalModel2/ modus, &
         & in_T, in_muB, in_lambda, &
         & in_eps, in_rhoB, in_rho, &
         & FileNameHydro

    integer :: ios

    call Write_ReadingInput('ThermalModel2',0)
    rewind(5)
    read(5,nml=ThermalModel2,IOSTAT=ios)
    call Write_ReadingInput('ThermalModel2',0,ios)




    write(*,*) "modus: ",modus
    select case(modus)
    case (0)
       write(*,*) "--- do not run solver"
       write(*,*) "T      = ",in_T
       write(*,*) "muB    = ",in_muB
       write(*,*) "lambda = ",in_lambda
    case (1)
       write(*,*) "--- run solver2a"
       write(*,*) "eps    = ", in_eps
       write(*,*) "rhoB   = ", in_rhoB
       !       write(*,*) "rho    = ", in_rho
    case (2)
       write(*,*) "--- run solver2b"
       write(*,*) "eps    = ", in_eps
       !       write(*,*) "rhoB   = ", in_rhoB
       write(*,*) "rho    = ", in_rho
    case (3)
       write(*,*) "--- run solver3"
       write(*,*) "eps    = ", in_eps
       write(*,*) "rhoB   = ", in_rhoB
       write(*,*) "rho    = ", in_rho

    case (10)
       write(*,*) "--- write eps,rhoB,rho from file"

    case (11)
       write(*,*) "--- fit eps,rhoB from file"

    case (12)
       write(*,*) "--- fit eps,rho from file"

    case (99)
       write(*,*) "--- write table to file"

    case default
       call Traceback("wrong modus")
    end select

    if (modus > 9) then
       if (len_trim(FileNameHydro)>0) then
          if (index(FileNameHydro,"/")>0) then
             call ExpandPath(FileNameHydro)
             FileNameHydro = trim(FileNameHydro)
          else
             FileNameHydro = trim(path_to_Input)//'/'//trim(FileNameHydro)
          end if
          write(*,*) 'Hydro.dat : ', trim(FileNameHydro)
       else
          FileNameHydro = trim(path_to_Input)//'/Hydro.dat'
          write(*,*) 'Hydro.dat : -- hardwired --'
       end if
    end if

    call Write_ReadingInput('ThermalModel2',1)

  end subroutine readInput

  !************************************************************************
  !****s* ThermalModel2/calcENB
  ! NAME
  ! subroutine calcENB
  ! PURPOSE
  ! calculates energy density, particle density and baryon density
  !
  ! we first have to calculate rhoB for the calculation of the potential
  ! part
  !************************************************************************
  subroutine calcENB(T,muB,lambda, eps,rhoB,rho, eps_part)

    use particleProperties, only: hadron

    real, intent(in) :: T
    real, intent(in) :: muB
    real, intent(in) :: lambda

    real, intent(out) :: eps
    real, intent(out) :: rhoB
    real, intent(out) :: rho
    real, intent(out) :: eps_part ! the particle contribution only

    integer :: ID

    real:: tmpN, tmpE, fak, epsI,rhoI

    logical, parameter :: doPot = .false.
    logical, save :: isFirst = .true.
    logical :: isnan

    eps = 0.
    rho = 0.
    rhoB = 0.
    eps_part = 0.

    write(*,*) 'calcENB:',T,muB,lambda

    isnan = (T/=T)
    if (isnan) call Traceback("T is NaN")

    !===== BARYONS =====

    !----- Nucleons -----

    fak = lambda*(hadron(1)%Spin*2+1)*(hadron(1)%isoSpinTimes2+1)

    call numFermiEN(hadron(1)%mass, 1., T, muB, tmpN,tmpE)
    tmpN = tmpN*fak
    tmpE = tmpE*fak

    eps = eps + tmpE
    rho = rho + tmpN

    call numFermiEN(hadron(1)%mass,-1., T, muB, tmpN,tmpE)
    tmpN = tmpN*fak
    tmpE = tmpE*fak

    eps = eps + tmpE
    rho = rho + tmpN

    !----- Resonances -----

    epsI = 0
    rhoI = 0
!$OMP PARALLEL DO  &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(T,muB,lambda) &
!$OMP PRIVATE(ID,tmpN,tmpE) &
!$OMP REDUCTION(+:epsI,rhoI) &
!$OMP IF (.not.isFirst)
    do ID=2,61
       call integrationBaryon(ID, T,muB, tmpN,tmpE)
       tmpN = tmpN*lambda
       tmpE = tmpE*lambda

       epsI = epsI + tmpE
       rhoI = rhoI + tmpN
    end do
    eps = eps + epsI
    rho = rho + rhoI

    rhoB = rho

    !===== MESONS =====

    epsI = 0
    rhoI = 0
!$OMP PARALLEL DO  &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(T,muB,lambda) &
!$OMP PRIVATE(ID,tmpN,tmpE) &
!$OMP REDUCTION(+:epsI,rhoI) &
!$OMP IF (.not.isFirst)
    do ID=101,122
       call integrationMeson(ID, T,muB, tmpN,tmpE)
       tmpN = tmpN*lambda
       tmpE = tmpE*lambda

       epsI = epsI + tmpE
       rhoI = rhoI + tmpN
    end do
    eps = eps + epsI
    rho = rho + rhoI

    eps_part = eps

    !===== POTENTIAL =====

    if (doPot) then
       call prepareTabBoltzmannPot(T,muB,rhoB)

       !----- Nucleons -----

       fak = lambda*(hadron(1)%Spin*2+1)*(hadron(1)%isoSpinTimes2+1)
       tmpE= numFermiPot(hadron(1)%mass, 1., T, muB, rhoB)+numFermiPot(hadron(1)%mass,-1., T, muB, rhoB)
       !    write(*,*) 1,tmpE*fak

       eps = eps + tmpE*fak

       !----- Resonances -----

       do ID=2,61
          call integrationBaryonPot(ID, T,muB,rhoB, tmpE)
          !       write(*,*) ID,tmpE*fak
          eps = eps + tmpE*lambda
       end do

    end if

    isFirst = .false.

  end subroutine calcENB


  !************************************************************************
  !****f* ThermalModel/BoltzmannN
  ! NAME
  ! real function BoltzmannN(mass,nB, T,muB)
  ! PURPOSE
  ! calculate the particle density for given parameters
  ! INPUTS
  ! * real :: mass
  ! * real :: nB  -- baryon number
  ! * real :: T   -- temperature
  ! * real :: muB -- baryo chemical potential
  ! OUTPUT
  ! * function value
  !************************************************************************
  real function BoltzmannN(mass,nB, T,muB)

    use constants, only: pi
    use besselK, only: BesselK2,BesselKexp2

    real, intent(in) :: mass,nB, T,muB

    real :: betamu, betam, res
    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    if (T<0.0005) then
       call Traceback("T too small")
    end if

    if (nB*muB < 400*T) then

       betamu = (nB*muB)/T
       betam  = mass/T
       res = exp(betamu) * BesselK2(betam)
       BoltzmannN = fak*res * mass**2 * T

    else

       betamu = (nB*muB)/T
       betam  = mass/T
       res = exp(betamu-betam) * BesselKexp2(betam)
       BoltzmannN = fak*res * mass**2 * T

    end if

    return
  end function BoltzmannN

  !************************************************************************
  !****f* ThermalModel/BoltzmannE
  ! NAME
  ! real function BoltzmannE(mass,nB, T,muB)
  ! PURPOSE
  ! calculate the energy density for given parameters
  ! INPUTS
  ! * real :: mass
  ! * real :: nB  -- baryon number
  ! * real :: T   -- temperature
  ! * real :: muB -- baryo chemical potential
  ! OUTPUT
  ! * function value
  !************************************************************************
  real function BoltzmannE(mass,nB, T,muB)

    use constants, only: pi
    use besselK, only: BesselK1,BesselK2,BesselKexp1,BesselKexp2

    real, intent(in) :: mass,nB, T,muB

    real :: betamu, betam, res
    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    if (T<0.0005) then
       call Traceback("T too small")
    end if

    if (nB*muB < 400*T) then

       betamu = (nB*muB)/T
       betam  = mass/T
       res = exp(betamu) * (BesselK2(betam)+betam/3.*BesselK1(betam))
       BoltzmannE = fak * res * 3 * mass**2 * T**2

    else

       betamu = (nB*muB)/T
       betam  = mass/T
       res = exp(betamu-betam) * (BesselKexp2(betam)+betam/3.*BesselKexp1(betam))
       BoltzmannE = fak * res * 3 * mass**2 * T**2

    end if

    return

  end function BoltzmannE


!!$  !************************************************************************
!!$  !****f* ThermalModel/FermiN
!!$  ! NAME
!!$  ! real function FermiN(mass,nB, T,muB)
!!$  ! PURPOSE
!!$  ! calculate the particle density for given parameters
!!$  !    n = 4\pi/(2\pi)^3 \int \frac{p^2 dp}{\exp(\beta(E-\mu))+1}
!!$  ! INPUTS
!!$  ! * real :: mass
!!$  ! * real :: nB  -- baryon number
!!$  ! * real :: T   -- temperature
!!$  ! * real :: muB -- baryo chemical potential
!!$  ! OUTPUT
!!$  ! * function value
!!$  ! NOTES
!!$  ! Does this by series expansion
!!$  !************************************************************************
!!$  real function FermiN(mass,nB, T,muB)
!!$
!!$    use constants, only: pi
!!$    use besselK, only: BesselK2, BesselKexp2
!!$
!!$    real, intent(in) :: mass,nB, T,muB
!!$
!!$    real :: betamu, betam, res, val, s
!!$    real, parameter :: fak = 4*pi/(2*pi*0.197)**3
!!$    integer :: i
!!$
!!$    if (T<0.001) then
!!$       call Traceback("T too small")
!!$    end if
!!$
!!$    if (nB*muB < 20*T) then
!!$
!!$       write(*,*) 'A'
!!$
!!$       betamu = (nB*muB)/T
!!$       betam  = mass/T
!!$       res = 0.
!!$       s = 1.0
!!$
!!$       do i=1,20
!!$
!!$          val = s/i*exp(i*betamu) * BesselK2(i*betam)
!!$          res = res + val
!!$          write(*,*) i,val,res
!!$          s = -s
!!$          if (abs(val) < 1e-9*res) exit
!!$       end do
!!$       FermiN = fak*res * mass**2 * T
!!$
!!$    else
!!$
!!$       write(*,*) 'B'
!!$
!!$       betamu = (nB*muB)/T
!!$       betam  = mass/T
!!$       res = 0.
!!$       s = 1.0
!!$
!!$       do i=1,20
!!$
!!$
!!$          val = s/i*exp(i*(betamu-betam)) * BesselKexp2(i*betam)
!!$          write(*,*) val
!!$          res = res + val
!!$          s = -s
!!$          if (abs(val) < 1e-9*res) exit
!!$       end do
!!$       FermiN = fak*res * mass**2 * T
!!$
!!$       write(*,*) fak*res * mass**2 * T
!!$
!!$    end if
!!$
!!$    return
!!$  end function FermiN
!!$
!!$  !************************************************************************
!!$  !****f* ThermalModel/FermiE
!!$  ! NAME
!!$  ! real function FermiE(mass,nB, T,muB)
!!$  ! PURPOSE
!!$  ! calculate the energy density for given parameters
!!$  !    e = 4\pi/(2\pi)^3 \int \frac{E p^2 dp}{\exp(\beta(E-\mu))+1}
!!$  ! INPUTS
!!$  ! * real :: mass
!!$  ! * real :: nB  -- baryon number
!!$  ! * real :: T   -- temperature
!!$  ! * real :: muB -- baryo chemical potential
!!$  ! OUTPUT
!!$  ! * function value
!!$  ! NOTES
!!$  ! Does this by series expansion
!!$  !************************************************************************
!!$  real function FermiE(mass,nB, T,muB)
!!$
!!$    use constants, only: pi
!!$    use besselK, only: BesselK1,BesselK2,BesselKexp1,BesselKexp2
!!$
!!$    real, intent(in) :: mass,nB, T,muB
!!$
!!$    real :: betamu, betam, res, val, s
!!$    real, parameter :: fak = 4*pi/(2*pi*0.197)**3
!!$    integer :: i
!!$
!!$    if (T<0.001) then
!!$       call Traceback("T too small")
!!$    end if
!!$
!!$    if (nB*muB < 400*T) then
!!$
!!$       betamu = (nB*muB)/T
!!$       betam  = mass/T
!!$       res = 0.
!!$       s = 1.0
!!$
!!$       do i=1,20
!!$          val = s/(i**2)*exp(i*betamu) * (BesselK2(i*betam)+i*betam/3.*BesselK1(i*betam))
!!$          res = res + val
!!$          s = -s
!!$          if (abs(val) < 1e-9*res) exit
!!$       end do
!!$
!!$       FermiE = fak*res * 3 * mass**2 * T**2
!!$    else
!!$
!!$       betamu = (nB*muB)/T
!!$       betam  = mass/T
!!$       res = 0.
!!$       s = 1.0
!!$
!!$       do i=1,20
!!$          val = s/(i**2)*exp(i*(betamu-betam)) * (BesselKexp2(i*betam)+i*betam/3.*BesselKexp1(i*betam))
!!$          res = res + val
!!$          s = -s
!!$          if (abs(val) < 1e-9*res) exit
!!$       end do
!!$
!!$       FermiE = fak*res * 3 * mass**2 * T**2
!!$
!!$    end if
!!$
!!$    return
!!$  end function FermiE

  !************************************************************************
  !****s* ThermalModel/numFermiEN
  ! NAME
  ! subroutine numFermiEN(mass,nB, T,muB, tmpN,tmpE)
  ! PURPOSE
  ! calculate energy and particle density for fermions by numerical
  ! integration of the p integral (with E = \sqrt{p^2+m^2})
  !    n = 4\pi/(2\pi)^3 \int \frac{p^2 dp}{\exp(\beta(E-\mu))+1}
  !    e = 4\pi/(2\pi)^3 \int \frac{E p^2 dp}{\exp(\beta(E-\mu))+1}
  !
  ! NOTES
  ! * The procedure is inspired by Numerical Recipies, 3rd edition
  ! * For values of \beta\mu < 15, the integration is done by trapezoidal
  !   rule. A variable transformation p = exp(v-exp(-v)) is performed.
  !   The considered v-range is ...
  !************************************************************************
  subroutine numFermiEN(mass,nB, T,muB, tmpN,tmpE)

    use constants, only: pi

    real, intent(in) :: mass,nB, T,muB
    real, intent(out) :: tmpN
    real, intent(out) :: tmpE

    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    real, dimension(2) :: sum, sumOld, hsum
    integer :: it
    real :: beta, betamu

    !--- for trapezoidal:
    real, parameter :: v1 = -4.5
    real, parameter :: v2 =  5.0
    integer :: iN, i
    real :: rN, dv, v
    logical, dimension(2) :: mask
    real, parameter :: epsAbs = 1e-6


    if (T<0.0005) then
       call Traceback("T too small")
    end if

    betamu = nB*muB/T ! "= alpha"
    beta = 1/T

    if (nB*muB < 1500000*T) then
       !===== do it with trapezoidal rule =====

       ! it = 1:
       sum = 0.5*(v2-v1)*(fFermiEN(v1,mass,beta,betamu)+fFermiEN(v2,mass,beta,betamu))
       sumOld = sum

       do it=2,20
          iN = 2**(it-2)
          rN = iN
          dv = (v2-v1)/rN

          v = v1 + 0.5*dv
          hsum = 0.
          do i=1,iN
             hsum = hsum + fFermiEN(v,mass,beta,betamu)
             v = v+dv
          end do

          sum=0.5*(sum+hsum*dv)

          if (it>5) then
             ! some covergence test
!             mask = abs(sum-sumOld) < epsAbs*abs(sumOld))
          end if

!          write(*,*) it, sum*fak

          sumOld = sum
       end do

       sum = sum*fak

    else
       !===== do it with splitted integrals and DE rule =====

       write(*,*) "mass,nB, T,muB:",mass,nB, T,muB
       call Traceback("not yet implemented")
    end if

    tmpN = sum(1)
    tmpE = sum(2)


  end subroutine numFermiEN

  function fFermiEN(v, mass, beta,betamu)
    real, dimension(2) :: fFermiEN
    real, intent(in) :: v, mass, beta,betamu

    real :: p,E

    p = exp(v-exp(-v))
    E = sqrt(p**2+mass**2)
    fFermiEN = p**2/(exp(beta*E-betamu)+1) * p*(1+exp(-v)) * (/1.,E/)

  end function fFermiEN

  !************************************************************************
  !****s* ThermalModel/numFermidEN
  ! NAME
  ! subroutine numFermidEN(mass,nB, T,muB, tmpN,tmpE)
  ! PURPOSE
  ! calculate integrals needed for derivatives of energy and particle
  ! density for fermions by numerical integration of the p integral
  ! (with E = \sqrt{p^2+m^2})
  !    F_n = 4\pi/(2\pi)^3 \int \frac{E^n p^2 dp \exp(\beta(E-\mu))}{(\exp(\beta(E-\mu))+1}
  ! for n=0,1,2
  !
  ! NOTES
  ! * The procedure is inspired by Numerical Recipies, 3rd edition
  ! * For values of \beta\mu < 15, the integration is done by trapezoidal
  !   rule. A variable transformation p = exp(v-exp(-v)) is performed.
  !   The considered v-range is ...
  !************************************************************************
  subroutine numFermidEN(mass,nB, T,muB, res)

    use constants, only: pi

    real, intent(in) :: mass,nB, T,muB
    real, dimension(0:2), intent(out) :: res

    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    real, dimension(0:2) :: sum, sumOld, hsum
    integer :: it
    real :: beta, betamu

    !--- for trapezoidal:
    real, parameter :: v1 = -4.5
    real, parameter :: v2 =  5.0
    integer :: iN, i
    real :: rN, dv, v
    logical, dimension(0:2) :: mask
    real, parameter :: epsAbs = 1e-6


    if (T<0.0005) then
       call Traceback("T too small")
    end if

    betamu = nB*muB/T ! "= alpha"
    beta = 1/T

    if (nB*muB < 1500000*T) then
       !===== do it with trapezoidal rule =====

       ! it = 1:
       sum = 0.5*(v2-v1)*(fFermidEN(v1,mass,beta,betamu)+fFermidEN(v2,mass,beta,betamu))
       sumOld = sum

       do it=2,20
          iN = 2**(it-2)
          rN = iN
          dv = (v2-v1)/rN

          v = v1 + 0.5*dv
          hsum = 0.
          do i=1,iN
             hsum = hsum + fFermidEN(v,mass,beta,betamu)
             v = v+dv
          end do

          sum=0.5*(sum+hsum*dv)

          if (it>5) then
             ! some covergence test
!             mask = abs(sum-sumOld) < epsAbs*abs(sumOld))
          end if

!          write(*,*) it, sum*fak

          sumOld = sum
       end do

       sum = sum*fak

    else
       !===== do it with splitted integrals and DE rule =====

       write(*,*) "mass,nB, T,muB:",mass,nB, T,muB
       call Traceback("not yet implemented")
    end if

    res = sum


  end subroutine numFermidEN

  function fFermidEN(v, mass, beta,betamu)
    real, dimension(3) :: fFermidEN
    real, intent(in) :: v, mass, beta,betamu

    real :: p,E

    p = exp(v-exp(-v))
    E = sqrt(p**2+mass**2)
    fFermidEN = p**2*exp(beta*E-betamu)/(exp(beta*E-betamu)+1)**2 * p*(1+exp(-v)) * (/1.,E,E**2/)

  end function fFermidEN

  !************************************************************************
  !****f* ThermalModel/numFermiPot
  ! NAME
  ! real function numFermiPot(mass,nB, T,muB,rho)
  ! PURPOSE
  ! calculate potential energy density for fermions by numerical
  ! integration of the p integral (with E = \sqrt{p^2+m^2})
  !    U = 4\pi/(2\pi)^3 \int \frac{p^2 dp}{\exp(\beta(E-\mu))+1} U_b(p,rho)
  !************************************************************************
  real function numFermiPot(mass,nB, T,muB,rho)

    use constants, only: pi

    real, intent(in) :: mass,nB, T,muB, rho

    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    real :: sum, sumOld, hsum
    integer :: it
    real :: beta, betamu

    !--- for trapezoidal:
    real, parameter :: v1 = -4.5
    real, parameter :: v2 =  5.0
    integer :: iN, i
    real :: rN, dv, v
    logical :: mask
    real, parameter :: epsAbs = 1e-6


    if (T<0.0005) then
       call Traceback("T too small")
    end if

    betamu = nB*muB/T ! "= alpha"
    beta = 1/T

    if (nB*muB < 1500000*T) then
       !===== do it with trapezoidal rule =====

       ! it = 1:
       sum = 0.5*(v2-v1)*(fFermiPot(v1,mass,beta,betamu,rho)+fFermiPot(v2,mass,beta,betamu,rho))
       sumOld = sum

       do it=2,20
          iN = 2**(it-2)
          rN = iN
          dv = (v2-v1)/rN

          v = v1 + 0.5*dv
          hsum = 0.
          do i=1,iN
             hsum = hsum + fFermiPot(v,mass,beta,betamu,rho)
             v = v+dv
          end do

          sum=0.5*(sum+hsum*dv)

          if (it>5) then
             ! some covergence test
!             mask = abs(sum-sumOld) < epsAbs*abs(sumOld))
          end if

!          write(*,*) it, sum*fak

          sumOld = sum
       end do

       sum = sum*fak

    else
       !===== do it with splitted integrals and DE rule =====

       write(*,*) "mass,nB, T,muB:",mass,nB, T,muB
       call Traceback("not yet implemented")
    end if

    numFermiPot = sum
  end function numFermiPot

  function fFermiPot(v, mass, beta,betamu,rho)
    real :: fFermiPot
    real, intent(in) :: v, mass, beta,betamu, rho

    real :: p,E

    p = exp(v-exp(-v))
    E = sqrt(p**2+mass**2)
    fFermiPot = p**2/(exp(beta*E-betamu)+1) * p*(1+exp(-v)) * Ubare(p,rho)/2

  end function fFermiPot

  !************************************************************************
  !****f* ThermalModel/numBoltzmannPot
  ! NAME
  ! real function numBoltzmannPot(mass,nB, T,muB,rho)
  ! PURPOSE
  ! calculate potential energy density for fermions by numerical
  ! integration of the p integral (with E = \sqrt{p^2+m^2})
  !    U = 4\pi/(2\pi)^3 \int \frac{p^2 dp}{\exp(\beta(E-\mu))+1} U_b(p,rho)/2
  !************************************************************************
  real function numBoltzmannPot(mass,nB, T,muB,rho)

    use constants, only: pi

    real, intent(in) :: mass,nB, T,muB, rho

    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    real :: sum, sumOld, hsum
    integer :: it
    real :: beta, betamu

    !--- for trapezoidal:
    real, parameter :: v1 = -4.5
    real, parameter :: v2 =  5.0
    integer :: iN, i
    real :: rN, dv, v
    logical :: mask
    real, parameter :: epsAbs = 1e-6

    if (T<0.0005) then
       call Traceback("T too small")
    end if

    betamu = nB*muB/T ! "= alpha"
    beta = 1/T

    if (nB*muB < 1500000*T) then
       !===== do it with trapezoidal rule =====

       ! it = 1:
       sum = 0.5*(v2-v1)*(fBoltzmannPot(v1,mass,beta,betamu,rho)+fBoltzmannPot(v2,mass,beta,betamu,rho))
       sumOld = sum

       do it=2,20
          iN = 2**(it-2)
          rN = iN
          dv = (v2-v1)/rN

          v = v1 + 0.5*dv
          hsum = 0.
          do i=1,iN
             hsum = hsum + fBoltzmannPot(v,mass,beta,betamu,rho)
             v = v+dv
          end do

          sum=0.5*(sum+hsum*dv)

          if (it>5) then
             ! some covergence test
!             mask = abs(sum-sumOld) < epsAbs*abs(sumOld))
          end if

!          write(*,*) it, sum*fak

          sumOld = sum
       end do

       sum = sum*fak

    else
       !===== do it with splitted integrals and DE rule =====

       write(*,*) "mass,nB, T,muB:",mass,nB, T,muB
       call Traceback("not yet implemented")
    end if


!!$    if (nB>0) then
!!$       write(917,*) mass,sum,sum*exp(mass/T)*T/mass
!!$    else
!!$       write(918,*) mass,sum,sum*exp(mass/T)*T/mass
!!$    end if

    numBoltzmannPot = sum
  end function numBoltzmannPot

  function fBoltzmannPot(v, mass, beta,betamu,rho)
    real :: fBoltzmannPot
    real, intent(in) :: v, mass, beta,betamu, rho

    real :: p,E

    p = exp(v-exp(-v))
    E = sqrt(p**2+mass**2)
    fBoltzmannPot = p**2*exp(-beta*E+betamu) * p*(1+exp(-v)) * Ubare(p,rho)/2

  end function fBoltzmannPot


  subroutine prepareTabBoltzmannPot(T,muB,rho)
    real, intent(in)  :: T,muB,rho

    real :: m
    integer :: iM
    logical :: flag
    real, dimension(2) :: v

    flag = (3.3/T < 100.) ! check argument of exp()

    write(419,*) '# T,muB,rho=',T,muB,rho

    do iM=0,nMass
       m = 1.0 + iM*dMass
       v = (/ numBoltzmannPot(m, 1., T,muB,rho),numBoltzmannPot(m,-1., T,muB,rho) /)
       if (flag) then
          arrTabBoltzmannPot(:,iM) = v * exp(m/T)*T/m
       else
          arrTabBoltzmannPotSign(:,iM) = sign(1.,v)
          arrTabBoltzmannPot(:,iM) = log( abs(v) )
       end if

       write(419,*) iM,m,v,arrTabBoltzmannPot(:,iM),arrTabBoltzmannPotSign(:,iM)
    end do

  end subroutine prepareTabBoltzmannPot

  real function tabBoltzmannPotS(mass,T)

    real, intent(in)  :: mass,T

    integer :: iM
    real :: m0, w
    real, dimension(2) :: v,v1,v2
    logical :: flag

    flag = (3.3/T < 100.) ! check argument of exp()

    iM = floor((mass - 1.0)/dMass)
    m0 = 1.0+iM*dMass
    w = (mass - m0)/dMass

    if (flag) then
       v = arrTabBoltzmannPot(:,iM)*(1.-w) + arrTabBoltzmannPot(:,iM+1)*w
       tabBoltzmannPotS = (v(1)+v(2))*exp(-mass/T)*mass/T
    else
       ! check sign change in interpolation bin
       v =  (arrTabBoltzmannPotSign(:,iM)*arrTabBoltzmannPotSign(:,iM+1))
       if ((v(1) < 0).or.(v(2) < 0)) then
          write(*,*) "--- Problem sign"
          tabBoltzmannPotS = 0
       else
          v1 = arrTabBoltzmannPot(:,iM)*(1.-w)
          v2 = arrTabBoltzmannPot(:,iM+1)*(w)
          v = arrTabBoltzmannPotSign(:,iM)*exp(v1+v2)
          tabBoltzmannPotS = v(1)+v(2)
       end if
    end if

  end function tabBoltzmannPotS


  !************************************************************************
  !****s* ThermalModel/integrateBaryon
  ! NAME
  ! subroutine integrationBaryon(ID,nB, T,muB, tmpN,tmpE)
  ! PURPOSE
  ! calculate energy and particle density for baryons (Boltzmann)
  !************************************************************************
  subroutine integrationBaryon(ID, T,muB, tmpN,tmpE)

    use constants, only: pi
    use particleProperties, only: hadron
    use baryonWidth, only: fullWidthBaryon

    integer, intent(in) :: ID
    real, intent(in)  :: T,muB
    real, intent(out) :: tmpN
    real, intent(out) :: tmpE

    integer :: im, nm
    real :: mass0, gamma0, nS, nPi
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac, deg
    real, parameter :: dy0 = pi/200.

    tmpN = 0.
    tmpE = 0.

    mass0  = hadron(id)%mass
    gamma0 = hadron(id)%width
    deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

    if (gamma0 < 1e-3) then
       if (massMax > mass0) then
          tmpN = BoltzmannN(mass0, 1., T,muB)+BoltzmannN(mass0,-1., T,muB)
          tmpE = BoltzmannE(mass0, 1., T,muB)+BoltzmannE(mass0,-1., T,muB)
          tmpN = tmpN*deg
          tmpE = tmpE*deg
       end if
       return
    end if

    mmin = hadron(id)%minmass
    mmax = massMax
    if (mmax < mmin) return ! integral is zero

    ymax = 2*atan(2*(mmax-mass0) / gamma0)
    ymin = 2*atan(2*(mmin-mass0) / gamma0)

    nm = max(int((ymax-ymin)/dy0),1)
    dy  = (ymax-ymin)/float(nm)

    do im=1,nm
       y = ymin+(float(im)-0.5)*dy
       m = 0.5*tan(0.5*y)*gamma0 + mass0
       m = min(max(m,mmin),mmax)
       gamma = fullWidthBaryon(id, m)
       spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
       intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
       tmpN = tmpN + spectral/intfac*(BoltzmannN(m, 1., T,muB)+BoltzmannN(m,-1., T,muB))
       tmpE = tmpE + spectral/intfac*(BoltzmannE(m, 1., T,muB)+BoltzmannE(m,-1., T,muB))
    end do
    tmpN = tmpN*deg*dy
    tmpE = tmpE*deg*dy

  end subroutine integrationBaryon

  !************************************************************************
  !****s* ThermalModel/integrateMeson
  ! NAME
  ! subroutine integrationMeson(ID,nB, T,muB, tmpN,tmpE)
  ! PURPOSE
  ! calculate energy and particle density for mesons (Boltzmann)
  !************************************************************************
  subroutine integrationMeson(ID, T,muB, tmpN,tmpE)

    use constants, only: pi
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson

    integer, intent(in) :: ID
    real, intent(in)  :: T,muB
    real, intent(out) :: tmpN
    real, intent(out) :: tmpE

    integer :: im, nm
    real :: mass0, gamma0, nS, nPi
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac, deg
    real, parameter :: dy0 = pi/200.

    tmpN = 0.
    tmpE = 0.

    mass0  = hadron(id)%mass
    gamma0 = hadron(id)%width
    deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

    if (gamma0 < 1e-3) then
       if (massMax > mass0) then
          tmpN = BoltzmannN(mass0, 0., T,muB)
          tmpE = BoltzmannE(mass0, 0., T,muB)
          tmpN = tmpN*deg
          tmpE = tmpE*deg
       end if
       return
    end if

    mmin = hadron(id)%minmass
    mmax = massMax
    if (mmax < mmin) return ! integral is zero

    ymax = 2*atan(2*(mmax-mass0) / gamma0)
    ymin = 2*atan(2*(mmin-mass0) / gamma0)

    nm = max(int((ymax-ymin)/dy0),1)
    dy  = (ymax-ymin)/float(nm)

    do im=1,nm
       y = ymin+(float(im)-0.5)*dy
       m = 0.5*tan(0.5*y)*gamma0 + mass0
       m = min(max(m,mmin),mmax)
       gamma = fullWidthMeson(id, m)
       spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
       intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
       tmpN = tmpN + spectral/intfac*BoltzmannN(m, 0., T,muB)
       tmpE = tmpE + spectral/intfac*BoltzmannE(m, 0., T,muB)
    end do
    tmpN = tmpN*deg*dy
    tmpE = tmpE*deg*dy

  end subroutine integrationMeson

  !****************************************************************************
  !****s* ThermalModel/integrateBaryonPot
  ! NAME
  ! subroutine integrationBaryonPot(ID,nB, T,muB,rho, tmpE)
  ! PURPOSE
  ! calculate potential energy density baryons (Boltzmann)
  !************************************************************************
  subroutine integrationBaryonPot(ID, T,muB,rho, tmpE)

    use constants, only: pi
    use particleProperties, only: hadron
    use baryonWidth, only: fullWidthBaryon

    integer, intent(in) :: ID
    real, intent(in)  :: T,muB,rho
    real, intent(out) :: tmpE

    integer :: im, nm
    real :: mass0, gamma0, nS, nPi
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac, deg
    real, parameter :: dy0 = pi/200.

    tmpE = 0.

    mass0  = hadron(id)%mass
    gamma0 = hadron(id)%width
    deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

    if (gamma0 < 1e-3) then
       if (massMax > mass0) then
          !          tmpE = numBoltzmannPot(mass0, 1., T,muB, rho)+numBoltzmannPot(mass0,-1., T,muB, rho)
          !          tmpE = tmpE*deg
          tmpE = tabBoltzmannPotS(mass0,T)*deg

       end if
       return
    end if

    mmin = hadron(id)%minmass
    mmax = massMax
    if (mmax < mmin) return ! integral is zero

    ymax = 2*atan(2*(mmax-mass0) / gamma0)
    ymin = 2*atan(2*(mmin-mass0) / gamma0)

    nm = max(int((ymax-ymin)/dy0),1)
    dy  = (ymax-ymin)/float(nm)

    do im=1,nm
       y = ymin+(float(im)-0.5)*dy
       m = 0.5*tan(0.5*y)*gamma0 + mass0
       m = min(max(m,mmin),mmax)
       gamma = fullWidthBaryon(id, m)
       spectral = 2./pi * m**2 * gamma / ((m**2-mass0**2)**2+m**2*gamma**2)
       intfac = gamma0 / ((m-mass0)**2+gamma0**2/4.)
       !       tmpE = tmpE + spectral/intfac*(numBoltzmannPot(m, 1., T,muB, rho)+numBoltzmannPot(m,-1., T,muB, rho))
       tmpE = tmpE + spectral/intfac*tabBoltzmannPotS(m,T)
    end do
    tmpE = tmpE*deg*dy

  end subroutine integrationBaryonPot



  !****************************************************************************
  !****if* ThermalModel2/Ubare
  ! NAME
  ! real function Ubare(p,rho)
  ! INPUTS
  ! * real :: p  -- momentum of particle (in LRF)
  ! * real :: rho       -- baryon density in LRF
  ! PURPOSE
  ! This function provides the bare potential U_b alone
  ! (please note: you will need 1/2*U_b)
  !****************************************************************************
  real function Ubare(p,rho)
    use constants, only: rhoNull

    real, intent(in) ::  p,rho
    integer, parameter :: EQS = 5

    !****************************************************************************
    ! Parameters of potential:
    ! Units: alpha [MeV], beta [MeV], c [MeV], lambda [1/fm], eta [fm^4]
    ! 'c' and 'lambda' determine the momemtum dependence.
    ! 'eta' concerns the surface term (cf. SurfacePotFlag)
    real, parameter, dimension(1:14) :: &
         alpha  = (/ -108.619,  -9.972, -286.99993 , -124.2611 , -29.253, 0., 0., 0.,  -22.9 ,    5.72 ,   56.7 , -209.2, -135.2 , -124.8  /), &
         beta   = (/  136.779,  38.032,  233.6517  ,   71.03007,  57.248, 0., 0., 0.,   33.9 ,   13.6  ,    4.75,  156.4,  176.8 ,  204.8 /), &
         tau    = (/    1.259,   2.404,    1.225707,    2.00108,   1.76 , 0., 0., 0.,    1.82,    2.61 ,    3.69,   1.35,    1.22,    1.2/), &
         lambda = (/    2.13 ,   2.126,    0.      ,    0.     ,   2.13 , 0., 0., 0.,    0.96,    0.807,    1.18,    0. ,    2.92,    4.0/), &
         c      = (/  -63.601, -63.601,    0.      ,    0.     , -63.516, 0., 0., 0., -102.  , -144.   , -145.  ,    0. ,  -61.23,  -77.29 /), &
         eta    = (/    0.4  ,   0.3  ,    0.4     ,    0.27   ,   0.36 , 0., 0., 0.,    0.  ,    0.   ,    0.  ,    0. , 0. , 0. /)
  !****************************************************************************

    Ubare = (alpha(EQS)*(rho/rhoNull) &
         + 2./(tau(EQS)+1.)*beta(EQS) *(rho/rhoNull)**tau(EQS))/1000.

    Ubare = Ubare + momentumDependentPart(p,c(EQS),lambda(EQS),rho)

  end function Ubare


  !****************************************************************************
  !****if* ThermalModel2/momentumDependentPart
  ! NAME
  ! real function momentumDependentPart(pin,c,lambda,rho)
  ! PURPOSE
  ! This function provides the analytical momdep. potential.
  ! INPUTS
  ! * real :: pIn       -- absolute momentum in GeV in LRF
  ! * real :: c         -- parameter of potential
  ! * real :: lambda    -- parameter of potential
  ! * real :: rho       -- baryon density in LRF
  ! NOTES
  ! see effenberger, dr.-thesis, pages 14-16
  !****************************************************************************
  real function momentumDependentPart(pin,c,lambda,rho)
    use constants, only: pi, rhoNull, hbarc

    !input
    real, intent(in) ::  c
    real, intent(in) ::  lambda
    real, intent(in) ::  pIn
    real, intent(in) ::  rho

    integer,parameter ::  isum=5

    real      pfermi, p
    real      pot, t1, t2, t3, t4, t5, t6, t7
    real      xtest, temp, zwi
    integer   isu


    p = pin/hbarc !  convert p in units 1/fm
    if (p.lt.1.0e-10) then
       p = 1.0e-10
    end if

    pfermi = (3./2.*pi*pi*rho )**(1./3.)

    !further on everything is calculated in 1 / fm
    !the constant c,lambda converts then the potential in gev

    !determine whether the small momentum expansion has to be used
    xtest = 2.0*pfermi*p/(pfermi**2+lambda**2)

    if (xtest .gt. 1.0e-06) then

       t1 = pi*lambda**3*4.0/(2.0*pi)**3
       t2 = (pfermi**2+lambda**2-p**2)/(2.0*p*lambda)
       t3 = (p+pfermi)**2+lambda**2
       t4 = (p-pfermi)**2+lambda**2
       t5 = 2.0*pfermi/lambda
       t6 = (p+pfermi)/lambda
       t7 = (p-pfermi)/lambda
       pot = t2*log(t3/t4) + t5 - 2.0*(atan(t6)-atan(t7))
       pot = pot*t1
       pot = pot*2.0*c/rhoNull

    else

       t1   = pi*lambda**3*4.0/(2.0*pi)**3
       t2   = (pfermi**2 + lambda**2 - p**2)/lambda
       t3   = 0.0
       temp = 2.0*pfermi/(pfermi**2+lambda**2)
       zwi = temp
       do isu = 0, isum
          if (abs(zwi)>1e-100) then
             t3 = t3 + zwi/float(2*isu+1)
             zwi = zwi * p*p * temp*temp
          else
             zwi = 0.0
          end if
       end do
       t5 = 2.0*pfermi/lambda
       t6 = (p+pfermi)/lambda
       t7 = (p-pfermi)/lambda

       pot = t2*t3+t5-2.0*(atan(t6)-atan(t7))
       pot = pot*t1
       pot = pot*2.0*c/rhoNull

    end if
    momentumDependentPart = pot/1000. !convert to GeV
    return
  end function momentumDependentPart

  !************************************************************************
  !************************************************************************

  real function f_T2v(T)
    real, intent(in) :: T
    f_T2v = log(T)
  end function f_T2v

  real function f_v2T(v)
    real, intent(in) :: v
    f_v2T = exp(v)
  end function f_v2T

  !************************************************************************

  real function f_muB2v_0(muB)
    real, intent(in) :: muB
    f_muB2v_0 = muB
  end function f_muB2v_0

  real function f_v2muB_0(v)
    real, intent(in) :: v
    f_v2muB_0 = v
  end function f_v2muB_0
  !------------------------------------------------------------------------
  real function f_muB2v(muB)
    use constants, only: pi
    real, intent(in) :: muB
    f_muB2v = 1./tan(pi*(1.2-muB)/1.2)
  end function f_muB2v

  real function f_v2muB(v)
    use constants, only: pi
    real, intent(in) :: v
    f_v2muB = 1.2/pi*(atan(v)+pi/2)
  end function f_v2muB

  !************************************************************************

  real function f_lambda2v_0(lambda)
    real, intent(in) :: lambda
    f_lambda2v_0 = lambda
  end function f_lambda2v_0

  real function f_v2lambda_0(v)
    real, intent(in) :: v
    f_v2lambda_0 = v
  end function f_v2lambda_0
  !------------------------------------------------------------------------
  real function f_lambda2v(lambda)
    use constants, only: pi
    real, intent(in) :: lambda
    f_lambda2v = 1./tan(pi*(2.0-lambda)/2.0)
  end function f_lambda2v

  real function f_v2lambda(v)
    use constants, only: pi
    real, intent(in) :: v
    f_v2lambda = 2.0/pi*(atan(v)+pi/2)
  end function f_v2lambda

  !************************************************************************


  !************************************************************************
  !****f* ThermalModel2/funcToMinimize2a
  ! NAME
  ! integer(c_int) function funcToMinimize2a(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver2a: eps,rhoB -> T,muB
  !************************************************************************
  function funcToMinimize2a(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align
    use particleProperties, only: hadron

    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize2a

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda, eps,rhoB,rho,eps_part

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = f_v2T(xv(1))
    muB = f_v2muB_0(xv(2))
    lambda = 1.0

    call calcENB(T,muB,lambda, eps,rhoB,rho, eps_part)

    yv(1) = eps  - par(1)
    yv(2) = rhoB - par(2)

!!$    write(*,*) '>>',xv(1:2)
!!$    write(*,*) '::',yv(1:2)
!!$    write(*,*) '..',eps,par(1)
!!$    write(*,*) '..',rhoB,par(2)

    funcToMinimize2a = fgsl_success
  end function funcToMinimize2a

  !************************************************************************
  !****f* ThermalModel2/funcToMinimize2b
  ! NAME
  ! integer(c_int) function funcToMinimize2b(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver2b: eps,rho -> T,muB
  !************************************************************************
  function funcToMinimize2b(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align
    use particleProperties, only: hadron

    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize2b

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda, eps,rhoB,rho,eps_part
    real, parameter :: c = 2.0
    real :: penalty


    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

!!$    write(*,*) 'xv:',xv

    call calcENB(T,muB,lambda, eps,rhoB,rho, eps_part)

    if (c*rhoB<(rho-rhoB)) then
       penalty = (rho-rhoB)/(c*rhoB)
    else
       penalty = 1.0
    end if

    yv(1) = eps  - par(1)
    yv(2) = rho - par(2)
!    yv(2) = rho*penalty - par(2)




!!$    write(*,*) '>>',xv(1:2)
!!$    write(*,*) '::',yv(1:2)
!!$    write(*,*) '..',eps,par(1)
!!$    write(*,*) '..',rho,par(2)


    funcToMinimize2b = fgsl_success
  end function funcToMinimize2b

  !************************************************************************
  !****f* ThermalModel2/funcToMinimize3
  ! NAME
  ! integer(c_int) function funcToMinimize3(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver3
  !************************************************************************
  function funcToMinimize3(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align
    use particleProperties, only: hadron
    use constants, only: pi

    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize3

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda, eps,rhoB,rho,eps_part
    real, parameter :: c = 2.0
    real :: penalty

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 3 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    write(*,*) 'xv:',xv

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = f_v2lambda(xv(3))

    call calcENB(T,muB,lambda, eps,rhoB,rho, eps_part)

    if (c*rhoB<(rho-rhoB)) then
       penalty = (rho-rhoB)/(c*rhoB)
    else
       penalty = 1.0
    end if

    yv(1) = eps  - par(1)
    yv(2) = rhoB - par(2)
    yv(3) = rho  - par(3)
!    yv(3) = rho*penalty  - par(3)

!!$    write(*,*) '>>',xv
!!$    write(*,*) '::',yv
!!$    write(*,*) '..',eps,par(1)
!!$    write(*,*) '..',rhoB,par(2)
!!$    write(*,*) '..',rho,par(3)

    funcToMinimize3 = fgsl_success
  end function funcToMinimize3

  !************************************************************************
  !****s* ThermalModel2/runSolver2a
  ! NAME
  ! subroutine runSolver2a()
  ! PURPOSE
  ! The routine to run a 2 dimensional root finding algorithm: given
  ! (eps,rhoB), it finds (T,muB)
  !************************************************************************
  subroutine runSolver2a()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 2
    type(fgsl_multiroot_fsolver) :: solver
    type(fgsl_multiroot_function) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 500
    integer :: iIt


    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrids,nParam)

    param = (/ in_eps, in_rhoB /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_init(funcToMinimize2a,nParam,pParam)

    xv = (/ f_T2v(in_T), f_muB2v_0(in_muB) /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)

    ! setup solver:
    status = fgsl_multiroot_fsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)

    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:2),fptr(1:2)
    end do

    call WriteStatus(status)

    out_T = f_v2T(xptr(1))
    out_muB = f_v2muB_0(xptr(2))
    out_lambda = 1.0

  end subroutine runSolver2a

  !************************************************************************
  !****s* ThermalModel2/runSolver2b
  ! NAME
  ! subroutine runSolver2b()
  ! PURPOSE
  ! The routine to run a 2 dimensional root finding algorithm: given
  ! (eps,rho), it finds (T,muB)
  !************************************************************************
  subroutine runSolver2b()

    use fgsl ! too much stuff for 'only'
    use constants, only: pi

    integer(fgsl_size_t),parameter :: nParam = 2
    type(fgsl_multiroot_fsolver) :: solver
    type(fgsl_multiroot_function) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 500
    integer :: iIt



!    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrids,nParam)
    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrid,nParam)
!    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_dnewton,nParam)

    param = (/ in_eps, in_rho /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_init(funcToMinimize2b,nParam,pParam)

    xv = (/ f_T2v(in_T), f_muB2v(in_muB) /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)

    ! setup solver:
    status = fgsl_multiroot_fsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)

    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:2),fptr(1:2)
    end do

    call WriteStatus(status)

    out_T = f_v2T(xptr(1))
    out_muB = f_v2muB(xptr(2))
    out_lambda = 1.0

  end subroutine runSolver2b


  !************************************************************************
  !****s* ThermalModel2/runSolver3
  ! NAME
  ! subroutine runSolver3()
  ! PURPOSE
  ! The routine to run a 3 dimensional root finding algorithm: given
  ! (eps,rhoB,rho), it finds (T,muB,lambda)
  !************************************************************************
  subroutine runSolver3()

    use fgsl ! too much stuff for 'only'
    use constants, only: pi

    integer(fgsl_size_t),parameter :: nParam = 3
    type(fgsl_multiroot_fsolver) :: solver
    type(fgsl_multiroot_function) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 100
    integer :: iIt

!    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrids,nParam)
    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_hybrid,nParam)
!    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_dnewton,nParam)
!    solver = fgsl_multiroot_fsolver_alloc(fgsl_multiroot_fsolver_broyden,nParam)

    param = (/ in_eps, in_rhoB, in_rho /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_init(funcToMinimize3,nParam,pParam)

    xv = (/ f_T2v(in_T), f_muB2v(in_muB), f_lambda2v(in_lambda) /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)

    ! setup solver:
    status = fgsl_multiroot_fsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)

    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:3),fptr(1:3)
    end do

    call WriteStatus(status)

    out_T = f_v2T(xptr(1))
    out_muB = f_v2muB(xptr(2))
    out_lambda = f_v2lambda(xptr(3))

  end subroutine runSolver3

  !************************************************************************
  !****s* ThermalModel2/WriteStatus
  ! NAME
  ! subroutine WriteStatus(status)
  ! PURPOSE
  ! Write a gsl status string to stdout
  !************************************************************************
  subroutine WriteStatus(status)
    use fgsl, only: fgsl_int, fgsl_char, fgsl_strmax
    use fgsl, only: fgsl_strerror

    integer(fgsl_int),intent(in) :: status
    character(kind=fgsl_char, len=fgsl_strmax) :: message
    message = fgsl_strerror(status)
    write(*,*)
    write(*,*) 'Status: ',message
    write(*,*)

  end subroutine WriteStatus

  !************************************************************************
  !************************************************************************
  !************************************************************************



  subroutine runMode10

    integer :: status

    real, dimension(100) :: arr
    real ::  arr1save = 1e22
    logical :: isnan

    real:: start_T, start_muB

    open(13,file=FileNameHydro ,status='old',action='read',iostat=status)
    if (status/=0) then
       write(*,*) 'status =',status
       call Traceback()
    end if

    start_T = in_T
    start_muB = in_muB

    do
       read(13,*,iostat=status) arr(1:2)
       if (status>0) cycle
       if (status<0) exit

       if (arr(1)<arr1save) then
          write(*,*) '----------------'
          write(421,*)
          write(421,*)

          arr1save = arr(1)
       end if
       arr1save = arr(1)


       in_T = arr(1)
       in_muB = arr(2)
       in_lambda = 1.0

       call calcENB(in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho, eps_part)
       write(*,*) 'IN: ',in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho, eps_part
       write(421,*) in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho, eps_part


    end do

  end subroutine runMode10


  subroutine runMode12

    integer :: status

    real, dimension(100) :: arr
    real ::  arr1save = 1e22
    logical :: isnan

    real:: start_T, start_muB

    open(13,file=FileNameHydro ,status='old',action='read',iostat=status)
    if (status/=0) then
       write(*,*) 'status =',status
       call Traceback()
    end if

    start_T = in_T
    start_muB = in_muB

    do
       read(13,*,iostat=status) arr(1:colMax)
       if (status>0) cycle
       if (status<0) exit

       if (arr(1)<arr1save) then
          write(*,*) '----------------'
          write(421,*)
          write(421,*)

          arr1save = arr(1)
       end if
       arr1save = arr(1)

       isnan = (arr(colEps)/=arr(colEps))

       write(*,*) arr(1),arr(colEps),arr(colRho),isnan
       if (isnan) then
          write(421,*) arr(1),arr(colEps),arr(colRho),isnan, 0.,0.,0.
          cycle
       end if

       if (arr(colEps)< 0.05) then
          write(421,*) arr(1),arr(colEps),arr(colRho),isnan, 0.,0.,0.
          cycle
       end if

       write(*,*) 'rhoM/rho=',(arr(colN0)-arr(colB0))/arr(colN0)
       if ((arr(colN0)-arr(colB0))/arr(colN0) < 1e-20) then
          ! meson contribution vanishes -> it must be T=0
          write(421,*) arr(1),arr(colEps),arr(colRho),isnan, 0.,0.,0.
          cycle
       end if


       in_T = start_T
       in_muB = start_muB

       in_eps = arr(colEps)
       in_rho = arr(colRho)


       call runSolver2b

       write(421,*) arr(1),arr(colEps),arr(colRho),isnan, out_T,out_muB,out_lambda
       flush(421)

    end do
  end subroutine runMode12

  subroutine runMode99

    integer :: iMu, iT
    real :: mu,T


    do iT=0,32
       T = iT*0.005
       if (iT==0) T = 0.001

       do iMu=0,50
          mu=iMu*0.02
          if (iMu==0) mu = 0.001

          call calcENB(T,mu,1.0, out_eps,out_rhoB,out_rho, eps_part)
          write(499,*) T,mu,1.0, out_eps,out_rhoB,out_rho, eps_part
       end do
       write(499,*)

    end do

  end subroutine runMode99

end program ThermalModel2
