!
! compile this with: make ARGS="-fopenmp"
!
! compile this with: make MODE=opt0 FPE=0
!
program ThermalModel3

  use, intrinsic :: iso_c_binding
  use inputGeneral
  use particleProperties
  use CallStack, only: Traceback
  use omp_lib

  implicit none

  real, parameter :: massMax = 3.0 ! upper boundary of integration


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

  real, parameter :: fit_T0 = 0.
  real, parameter :: fit_muB1 = 0.
  real, parameter :: fit_muB2 = 1.2
  real, parameter :: fit_lambda1 = 0.
  real, parameter :: fit_lambda2 = 2.0

  real :: h_eps,h_rhoB,h_rho


  !****************************************************************************
  !****g* ThermalModel3/FileNameHydro
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

  real, dimension(3,3) :: hJ

!  real, dimension(2,2) :: try

  call readInputGeneral
  call initParticleProperties
  call readInput

!$OMP PARALLEL
  !$ write(*,*) 'Thread: ',OMP_GET_THREAD_NUM()
!$OMP END PARALLEL


  ! valid range: exp(-708)...exp(709) = 3.3e-308 ... 8.2e+307

!  do out_T=10,1000
!     write(*,*) out_T,exp(-out_T)
!  end do

!!$  out_T=1.0
!!$  do
!!$     write(*,*) out_T, log(out_T)
!!$     out_T = out_T/10
!!$  end do

!!$  try = reshape( source = (/ 1.1,1.2, 2.1,2.2  /), &
!!$       shape = (/ 2, 2 /))
!!$
!!$  write(*,*) try(1,:)
!!$  write(*,*) try(2,:)
!!$  stop

  select case(modus)
  case (0)
     call calcENB(in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho
     write(*,*)
     write(*,*) 'in_eps =',out_eps,'; in_rhoB =',out_rhoB,'; in_rho =',out_rho,'; in_T=',in_T,'; in_muB =',in_muB,'; in_lambda =',in_lambda
     write(*,*)

     call calcJacobian(in_T,in_muB,in_lambda, hJ)
     write(*,*) hJ(1,:)
     write(*,*) hJ(2,:)
     write(*,*) hJ(3,:)

     write(*,*) 'dT:'
     call calcENB(in_T+0.00001,in_muB,in_lambda, h_eps,h_rhoB,h_rho)
     write(*,*) (h_eps-out_eps)/0.00001,(h_rhoB-out_rhoB)/0.00001,(h_rho-out_rho)/0.00001
     write(*,*) 'dmu:'
     call calcENB(in_T,in_muB+0.0001,in_lambda, h_eps,h_rhoB,h_rho)
     write(*,*) (h_eps-out_eps)/0.0001,(h_rhoB-out_rhoB)/0.0001,(h_rho-out_rho)/0.0001
     write(*,*) 'dlambda:'
     call calcENB(in_T,in_muB,in_lambda+0.0001, h_eps,h_rhoB,h_rho)
     write(*,*) (h_eps-out_eps)/0.0001,(h_rhoB-out_rhoB)/0.0001,(h_rho-out_rho)/0.0001

     call invertJ(hJ)

  case (1)
     call runSolver2a
     call calcENB(out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, in_eps,in_rhoB,in_rho
     write(*,*) 'OUT:',out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho

  case (2)
     call runSolver2b
     call calcENB(out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, in_eps,in_rhoB,in_rho
     write(*,*) 'OUT:',out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho


  case (3)
     call runSolver3
     call calcENB(out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho)
     write(*,*) 'IN: ',in_T,in_muB,in_lambda, in_eps,in_rhoB,in_rho
     write(*,*) 'OUT:',out_T,out_muB,out_lambda, out_eps,out_rhoB,out_rho

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
  !****s* ThermalModel3/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "ThermalModel2".
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

    call Write_ReadingInput('ThermalModel3',1)

  end subroutine readInput

  !************************************************************************
  !****s* ThermalModel3/calcENB
  ! NAME
  ! subroutine calcENB
  ! PURPOSE
  ! calculates energy density, particle density and baryon density
  !
  !************************************************************************
  subroutine calcENB(T,muB,lambda, eps,rhoB,rho)

    use particleProperties, only: hadron

    real, intent(in) :: T
    real, intent(in) :: muB
    real, intent(in) :: lambda

    real, intent(out) :: eps
    real, intent(out) :: rhoB
    real, intent(out) :: rho

    integer :: ID

    real:: tmpN, tmpE, fak, epsI,rhoI


    logical, save :: isFirst = .true.
    logical :: isnan

    real, dimension(0:2) :: res, resP,resM

    eps = 0.
    rho = 0.
    rhoB = 0.

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

!!$    rhoB = rho ! for debug
!!$    return     ! for debug

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
       call integrationBaryon(ID, T,muB, resP,resM)
       epsI = epsI + (resP(1)+resM(1))*lambda
       rhoI = rhoI + (resP(0)+resM(0))*lambda
    end do
    eps = eps + epsI
    rho = rho + rhoI

    rhoB = rho

!!$    return     ! for debug

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
       call integrationMeson(ID, T, res)
       epsI = epsI + res(1)*lambda
       rhoI = rhoI + res(0)*lambda
    end do
    eps = eps + epsI
    rho = rho + rhoI

    isFirst = .false.

  end subroutine calcENB

  !************************************************************************
  !****s* ThermalModel3/calcJacobian
  ! NAME
  ! subroutine calcJacobian(T,muB,lambda, J)
  ! PURPOSE
  ! Calculate the derivatives with respect to T,muB,lambda
  !
  ! column/row order adjusted using some test outputs...
  !************************************************************************
  subroutine calcJacobian(T,muB,lambda, J)

    use particleProperties, only: hadron

    real, intent(in) :: T
    real, intent(in) :: muB
    real, intent(in) :: lambda
    real, dimension(1:3,1:3), intent(out) :: J !
    ! index 1: T,mu,lambda
    ! index 2: eps,rhoB,rho,


    real :: tmpN, tmpE, fak
    logical :: isnan
    real, dimension(0:2) :: res, resP,resM
    integer :: ID

    J = 0.

    isnan = (T/=T)
    if (isnan) call Traceback("T is NaN")

    ! we first calculate everything neglecting lambda and multiply
    ! it at the end

    !----- Nucleons -----

    fak = (hadron(1)%Spin*2+1)*(hadron(1)%isoSpinTimes2+1)

    call numFermiEN(hadron(1)%mass, 1., T, muB, tmpN,tmpE)
    res = numFermidEN(hadron(1)%mass, 1., T, muB)

    J(1,1) = J(1,1) + fak*(res(2)-muB*res(1))/T**2      ! = de/dT
    J(2,1) = J(2,1) + fak*res(1)/T                      ! = de/dmu
    J(3,1) = J(3,1) + fak*tmpE                          ! = de/dlambda
    J(1,2) = J(1,2) + fak*(res(1)-muB*res(0))/T**2      ! = drhoB/dT
    J(2,2) = J(2,2) + fak*res(0)/T                      ! = drhoB/dmu
    J(3,2) = J(3,2) + fak*tmpN                          ! = drhoB/dlambda

    call numFermiEN(hadron(1)%mass,-1., T, muB, tmpN,tmpE)
    res = numFermidEN(hadron(1)%mass,-1., T, muB)

    J(1,1) = J(1,1) + fak*(res(2)+muB*res(1))/T**2      ! = de/dT
    J(2,1) = J(2,1) - fak*res(1)/T                      ! = de/dmu
    J(3,1) = J(3,1) + fak*tmpE                          ! = de/dlambda
    J(1,2) = J(1,2) + fak*(res(1)+muB*res(0))/T**2      ! = drhoB/dT
    J(2,2) = J(2,2) - fak*res(0)/T                      ! = drhoB/dmu
    J(3,2) = J(3,2) + fak*tmpN                          ! = drhoB/dlambda

!!$    J(:,3) = J(:,2) ! rho == rhoB ! for debug
!!$    J(1:2,:) = J(1:2,:)*lambda ! for debug
!!$    return                     ! for debug

    !----- Resonances -----

    do ID=2,61
       call integrationBaryon(ID, T,muB, resP,resM)
       J(1,1) = J(1,1) + ((resP(2)+resM(2))-muB*(resP(1)-resM(1)))/T**2 ! = de/dT
       J(2,1) = J(2,1) + (resP(1)-resM(1))/T                            ! = de/dmu
       J(3,1) = J(3,1) + (resP(1)+resM(1))                              ! = de/dlambda
       J(1,2) = J(1,2) + ((resP(1)+resM(1))-muB*(resP(0)-resM(0)))/T**2 ! = drhoB/dT
       J(2,2) = J(2,2) + (resP(0)-resM(0))/T                            ! = drhoB/dmu
       J(3,2) = J(3,2) + (resP(0)+resM(0))                              ! = drhoB/dlambda
    end do

    J(:,3) = J(:,2) ! rho == rhoB

!!$    J(1:2,:) = J(1:2,:)*lambda ! for debug
!!$    return                     ! for debug

    !----- Mesons -----

    do ID=101,122
       call integrationMeson(ID, T, res)
       J(1,1) = J(1,1) + res(2)/T**2        ! = de/dT
       J(3,1) = J(3,1) + res(1)             ! = de/dlambda
       J(1,3) = J(1,3) + res(1)/T**2        ! = drho/dT
       J(3,3) = J(3,3) + res(0)             ! = drho/dlambda
    end do

    J(1:2,:) = J(1:2,:)*lambda ! = d/dT, d/dmu

  end subroutine calcJacobian

  !************************************************************************
  !****s* ThermalModel3/numFermiEN
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
!    logical, dimension(2) :: mask
!    real, parameter :: epsAbs = 1e-6


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

    if (beta*E-betamu < 400) then
       fFermiEN = p**2/(exp(beta*E-betamu)+1) * p*(1+exp(-v)) * (/1.,E/)
    else
       fFermiEN = 0.
    end if

  end function fFermiEN

  !************************************************************************
  !****f* ThermalModel/numFermidEN
  ! NAME
  ! function numFermidEN(mass,nB, T,muB) result(res)
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
  function numFermidEN(mass,nB, T,muB) result(res)

    use constants, only: pi

    real, intent(in) :: mass,nB, T,muB
    real, dimension(0:2) :: res

    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    real, dimension(0:2) :: sum, sumOld, hsum
    integer :: it
    real :: beta, betamu

    !--- for trapezoidal:
    real, parameter :: v1 = -4.5
    real, parameter :: v2 =  5.0
    integer :: iN, i
    real :: rN, dv, v
!    logical, dimension(0:2) :: mask
!    real, parameter :: epsAbs = 1e-6


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


  end function numFermidEN

  function fFermidEN(v, mass, beta,betamu)
    real, dimension(3) :: fFermidEN
    real, intent(in) :: v, mass, beta,betamu

    real :: p,E

    ! please note: we square exp(...), therefore the valid range of the
    ! argument has to be halved

    p = exp(v-exp(-v))
    E = sqrt(p**2+mass**2)
    if (beta*E-betamu < 200) then
       fFermidEN = p**2*exp(beta*E-betamu)/(exp(beta*E-betamu)+1)**2 * p*(1+exp(-v)) * (/1.,E,E**2/)
    else
       fFermidEN = 0.
    end if

  end function fFermidEN

  !************************************************************************
  !****s* ThermalModel/BoltzmandEN
  ! NAME
  ! function BoltzmanndEN(mass,nB, T,muB) result(res)
  ! PURPOSE
  ! calculate integrals needed for derivatives of energy and particle
  ! density for Boltzmanns by numerical integration of the p integral
  ! (with E = \sqrt{p^2+m^2})
  !    B_n = 4\pi/(2\pi)^3 \int \frac{E^n p^2 dp} \exp(\beta(E-\mu))
  ! for n=0,1,2
  !
  ! Thus: n_B = B_0, e_B = B_1
  !
  ! NOTES
  ! * The procedure is inspired by Numerical Recipies, 3rd edition
  ! * For values of \beta\mu < 15, the integration is done by trapezoidal
  !   rule. A variable transformation p = exp(v-exp(-v)) is performed.
  !   The considered v-range is ...
  !************************************************************************
  function BoltzmanndEN(mass,nB, T,muB) result(res)
    use constants, only: pi
    use besselK, only: BesselK0,BesselK1,BesselK2,BesselKexp0,BesselKexp1,BesselKexp2

    real, intent(in) :: mass,nB, T,muB
    real, dimension(0:2) :: res

    real :: betamu, betam
    real, parameter :: fak = 4*pi/(2*pi*0.197)**3

    if (T<0.0005) then
       call Traceback("T too small")
    end if

    res = 0.

    if ((nB*muB < 400*T).and.(mass < 400*T)) then

       betamu = (nB*muB)/T
       betam  = mass/T

       res(0) = BesselK2(betam) * mass**2 * T
       res(1) = (3.*BesselK2(betam)+betam*BesselK1(betam)) * mass**2 * T**2
       res(2) = (betam*(betam**2+12.)*BesselK0(betam)+(5*betam**2+24.)*BesselK1(betam)) * mass * T**4

       res = res * fak * exp(betamu)
    else

       betamu = (nB*muB-mass)/T ! attention: abuse
       betam  = mass/T

       if (betamu>-400) then
          res(0) = BesselKexp2(betam) * mass**2 * T
          res(1) = (3.*BesselKexp2(betam)+betam*BesselKexp1(betam)) * mass**2 * T**2
          res(2) = (betam*(betam**2+12.)*BesselKexp0(betam)+(5*betam**2+24.)*BesselKexp1(betam)) * mass * T**4

          res = res * fak * exp(betamu)
       end if

    end if

  end function BoltzmanndEN

  !************************************************************************
  !****s* ThermalModel/integrateBaryon
  ! NAME
  ! subroutine integrationBaryon(ID,nB, T,muB, resP, resB)
  ! PURPOSE
  ! calculate energy and particle density for baryons (Boltzmann)
  ! (n = res(0), e = res(1))
  ! (resP = baryons, resM = antibaryons)
  !************************************************************************
  subroutine integrationBaryon(ID, T,muB, resP, resM)

    use constants, only: pi
    use particleProperties, only: hadron
    use baryonWidth, only: fullWidthBaryon

    integer, intent(in) :: ID
    real, intent(in)  :: T,muB
    real, dimension(0:2), intent(out) :: resP,resM

    integer :: im, nm
    real :: mass0, gamma0
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac, deg
    real, parameter :: dy0 = pi/200.


    resP = 0.
    resM = 0.

!!$    res = BoltzmanndEN(0.938, 1., T,muB)
!!$    write(*,*) 'res=',res*4
!!$    stop

    mass0  = hadron(id)%mass
    gamma0 = hadron(id)%width
    deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

    if (gamma0 < 1e-3) then
       if (mass0 < massMax) then
          resP = BoltzmanndEN(mass0, 1., T,muB)*deg
          resM = BoltzmanndEN(mass0,-1., T,muB)*deg
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
       resP = resP + spectral/intfac*BoltzmanndEN(m, 1., T,muB)
       resM = resM + spectral/intfac*BoltzmanndEN(m,-1., T,muB)
    end do

    resP = resP*deg*dy
    resM = resM*deg*dy

  end subroutine integrationBaryon

  !************************************************************************
  !****s* ThermalModel/integrateMeson
  ! NAME
  ! subroutine integrationMeson(ID,nB, T, res)
  ! PURPOSE
  ! calculate energy and particle density for mesons (Boltzmann)
  ! (n = res(0), e = res(1))
  !************************************************************************
  subroutine integrationMeson(ID, T, res)

    use constants, only: pi
    use particleProperties, only: hadron
    use mesonWidth, only: fullWidthMeson

    integer, intent(in) :: ID
    real, intent(in)  :: T
    real, dimension(0:2), intent(out) :: res

    integer :: im, nm
    real :: mass0, gamma0
    real :: mmin, mmax, m
    real :: ymin, ymax, dy, y
    real :: gamma, spectral, intfac, deg
    real, parameter :: dy0 = pi/200.

    res = 0.

    mass0  = hadron(id)%mass
    gamma0 = hadron(id)%width
    deg = (hadron(ID)%Spin*2+1)*(hadron(ID)%isoSpinTimes2+1)

    if (gamma0 < 1e-3) then
       if (mass0 < massMax) then
          res =  BoltzmanndEN(mass0, 0., T,0.)*deg
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
       res = res + spectral/intfac*BoltzmanndEN(m, 0., T,0.)
    end do

    res = res*deg*dy

  end subroutine integrationMeson

  !************************************************************************
  !************************************************************************

  real function f_T2v(T)
    real, intent(in) :: T
    !    f_T2v = log(T-fit_T0)
    f_T2v = T
  end function f_T2v

  real function f_v2T(v)
    real, intent(in) :: v
    !    f_v2T = fit_T0+exp(v)
    f_v2T = v
  end function f_v2T

  real function f_d_v2T(v)
    real, intent(in) :: v
    !    f_d_v2T = exp(v)
    f_d_v2T = 1.
  end function f_d_v2T

  !************************************************************************

  real function f_muB2v(muB)
    use constants, only: pi
    real, intent(in) :: muB
!    f_muB2v = 1./tan(pi*(fit_muB1-muB)/(fit_muB2-fit_muB1))
    f_muB2v = muB
  end function f_muB2v

  real function f_v2muB(v)
    use constants, only: pi
    real, intent(in) :: v
!    f_v2muB = (fit_muB2-fit_muB1)/pi*(atan(v)+pi/2)+fit_muB1
    f_v2muB = v
  end function f_v2muB

  real function f_d_v2muB(v)
    use constants, only: pi
    real, intent(in) :: v
!    f_d_v2muB = (fit_muB2-fit_muB1)/(pi*(1+v**2))
    f_d_v2muB = 1.
  end function f_d_v2muB

  !************************************************************************

  real function f_lambda2v(lambda)
    use constants, only: pi
    real, intent(in) :: lambda
!    f_lambda2v = 1./tan(pi*(fit_lambda1-lambda)/(fit_lambda2-fit_lambda1))
    f_lambda2v = lambda
  end function f_lambda2v

  real function f_v2lambda(v)
    use constants, only: pi
    real, intent(in) :: v
!    f_v2lambda = (fit_lambda2-fit_lambda1)/pi*(atan(v)+pi/2)+fit_lambda1
    f_v2lambda = v
  end function f_v2lambda

  real function f_d_v2lambda(v)
    use constants, only: pi
    real, intent(in) :: v
!    f_d_v2lambda = (fit_lambda2-fit_lambda1)/(pi*(1+v**2))
    f_d_v2lambda = 1.
  end function f_d_v2lambda

  !************************************************************************
  !************************************************************************
  !************************************************************************


  !************************************************************************
  !****f* ThermalModel3/funcToMinimize2a
  ! NAME
  ! integer(c_int) function funcToMinimize2a(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver2a: eps,rhoB -> T,muB
  !************************************************************************
  function funcToMinimize2a_f(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align

    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize2a_f

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

    write(*,*) '   _f  :',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yv(1) = J(3,1)*lambda - par(1) ! e
    yv(2) = J(3,2)*lambda - par(2) ! rhoB

!!$    write(*,*) '>>',xv(1:2)
!!$    write(*,*) '::',yv(1:2)
!!$    write(*,*) '..',eps,par(1)
!!$    write(*,*) '..',rhoB,par(2)

    funcToMinimize2a_f = fgsl_success

  end function funcToMinimize2a_f

  function funcToMinimize2a_df(x, params, df) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector, fgsl_matrix
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align, fgsl_matrix_align

    type(c_ptr), value :: x, params, df
    integer(c_int) :: funcToMinimize2a_df

    type(fgsl_vector) :: fx
    type(fgsl_matrix) :: dff
    real(fgsl_double), pointer :: par(:), xv(:), yd(:,:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(dff, df)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_matrix_align(yd, dff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

    write(*,*) '   _df :',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yd(1:2,1:2) = reshape( source = (/ &
         J(1,1)*f_d_v2T(xv(1)),         & ! de/dT
         J(2,1)*f_d_v2muB(xv(2)),       & ! de/dmuB
         J(1,2)*f_d_v2T(xv(1)),         & ! drhoB/dT
         J(2,2)*f_d_v2muB(xv(2)) /),    & ! drhoB/dmuB
         shape = (/ 2, 2 /))

    funcToMinimize2a_df = fgsl_success
  end function funcToMinimize2a_df

  function funcToMinimize2a_fdf(x, params, f, df) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector, fgsl_matrix
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align, fgsl_matrix_align

    type(c_ptr), value :: x, params, f, df
    integer(c_int) :: funcToMinimize2a_fdf

    type(fgsl_vector) :: fx, ff
    type(fgsl_matrix) :: dff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:), yd(:,:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call fgsl_obj_c_ptr(dff, df)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)
    status = fgsl_matrix_align(yd, dff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

    write(*,*) '   _fdf:',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yv(1) = J(3,1)*lambda - par(1) ! e
    yv(2) = J(3,2)*lambda - par(2) ! rhoB

    yd(1:2,1:2) = reshape( source = (/ &
         J(1,1)*f_d_v2T(xv(1)),         & ! de/dT
         J(2,1)*f_d_v2muB(xv(2)),       & ! de/dmuB
         J(1,2)*f_d_v2T(xv(1)),         & ! drhoB/dT
         J(2,2)*f_d_v2muB(xv(2)) /),    & ! drhoB/dmuB
         shape = (/ 2, 2 /))

    funcToMinimize2a_fdf = fgsl_success
  end function funcToMinimize2a_fdf



  !************************************************************************
  !****f* ThermalModel3/funcToMinimize2b
  ! NAME
  ! integer(c_int) function funcToMinimize2b(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver2b: eps,rho -> T,muB
  !************************************************************************
  function funcToMinimize2b_f(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align

    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize2b_f

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

    write(*,*) '   _f  :',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yv(1) = J(3,1)*lambda - par(1) ! e
    yv(2) = J(3,3)*lambda - par(2) ! rho

!!$    write(*,*) '>>',xv(1:2)
!!$    write(*,*) '::',yv(1:2)
!!$    write(*,*) '..',eps,par(1)
!!$    write(*,*) '..',rhoB,par(2)

    funcToMinimize2b_f = fgsl_success

  end function funcToMinimize2b_f

  function funcToMinimize2b_df(x, params, df) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector, fgsl_matrix
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align, fgsl_matrix_align

    type(c_ptr), value :: x, params, df
    integer(c_int) :: funcToMinimize2b_df

    type(fgsl_vector) :: fx
    type(fgsl_matrix) :: dff
    real(fgsl_double), pointer :: par(:), xv(:), yd(:,:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(dff, df)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_matrix_align(yd, dff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

    write(*,*) '   _df :',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yd(1:2,1:2) = reshape( source = (/ &
         J(1,1)*f_d_v2T(xv(1)),         & ! de/dT
         J(2,1)*f_d_v2muB(xv(2)),       & ! de/dmuB
         J(1,3)*f_d_v2T(xv(1)),         & ! drho/dT
         J(2,3)*f_d_v2muB(xv(2)) /),    & ! drho/dmuB
         shape = (/ 2, 2 /))

    funcToMinimize2b_df = fgsl_success
  end function funcToMinimize2b_df

  function funcToMinimize2b_fdf(x, params, f, df) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector, fgsl_matrix
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align, fgsl_matrix_align

    type(c_ptr), value :: x, params, f, df
    integer(c_int) :: funcToMinimize2b_fdf

    type(fgsl_vector) :: fx, ff
    type(fgsl_matrix) :: dff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:), yd(:,:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call fgsl_obj_c_ptr(dff, df)
    call c_f_pointer(params, par, (/ 2 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)
    status = fgsl_matrix_align(yd, dff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = 1.0

    write(*,*) '   _fdf:',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yv(1) = J(3,1)*lambda - par(1) ! e
    yv(2) = J(3,3)*lambda - par(2) ! rho

    yd(1:2,1:2) = reshape( source = (/ &
         J(1,1)*f_d_v2T(xv(1)),         & ! de/dT
         J(2,1)*f_d_v2muB(xv(2)),       & ! de/dmuB
         J(1,3)*f_d_v2T(xv(1)),         & ! drho/dT
         J(2,3)*f_d_v2muB(xv(2)) /),    & ! drho/dmuB
         shape = (/ 2, 2 /))

    funcToMinimize2b_fdf = fgsl_success
  end function funcToMinimize2b_fdf


  !************************************************************************
  !****f* ThermalModel3/funcToMinimize3
  ! NAME
  ! integer(c_int) function funcToMinimize3(x, params, f)
  ! PURPOSE
  ! The function to minimize with solver3
  !************************************************************************
  function funcToMinimize3_f(x, params, f) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align

    type(c_ptr), value :: x, params, f
    integer(c_int) :: funcToMinimize3_f

    type(fgsl_vector) :: fx, ff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call c_f_pointer(params, par, (/ 3 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = f_v2lambda(xv(3))

    write(*,*) '   _f  :',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yv(1) = J(3,1)*lambda - par(1) ! e
    yv(2) = J(3,2)*lambda - par(2) ! rhoB
    yv(3) = J(3,3)*lambda - par(3) ! rho

!!$    write(*,*) '>>',xv(1:2)
!!$    write(*,*) '::',yv(1:2)
!!$    write(*,*) '..',eps,par(1)
!!$    write(*,*) '..',rhoB,par(2)

    funcToMinimize3_f = fgsl_success

  end function funcToMinimize3_f

  function funcToMinimize3_df(x, params, df) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector, fgsl_matrix
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align, fgsl_matrix_align

    type(c_ptr), value :: x, params, df
    integer(c_int) :: funcToMinimize3_df

    type(fgsl_vector) :: fx
    type(fgsl_matrix) :: dff
    real(fgsl_double), pointer :: par(:), xv(:), yd(:,:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(dff, df)
    call c_f_pointer(params, par, (/ 3 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_matrix_align(yd, dff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = f_v2lambda(xv(3))

    write(*,*) '   _df :',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yd(1:3,1:3) = reshape( source = (/ &
         J(1,1)*f_d_v2T(xv(1)),         & ! de/dT
         J(2,1)*f_d_v2muB(xv(2)),       & ! de/dmuB
         J(3,1)*f_d_v2lambda(xv(3)),    & ! de/dlambda
         J(1,2)*f_d_v2T(xv(1)),         & ! drhoB/dT
         J(2,2)*f_d_v2muB(xv(2)),       & ! drhoB/dmuB
         J(3,2)*f_d_v2lambda(xv(3)),    & ! drhoB/dlambda
         J(1,3)*f_d_v2T(xv(1)),         & ! drho/dT
         J(2,3)*f_d_v2muB(xv(2)),       & ! drho/dmuB
         J(3,3)*f_d_v2lambda(xv(3)) /), & ! drho/dlambda
         shape = (/ 3, 3 /))

    funcToMinimize3_df = fgsl_success
  end function funcToMinimize3_df

  function funcToMinimize3_fdf(x, params, f, df) bind(c)
    use fgsl, only: fgsl_success
    use fgsl, only: fgsl_int, fgsl_double, fgsl_vector, fgsl_matrix
    use fgsl, only: fgsl_obj_c_ptr, fgsl_vector_align, fgsl_matrix_align

    type(c_ptr), value :: x, params, f, df
    integer(c_int) :: funcToMinimize3_fdf

    type(fgsl_vector) :: fx, ff
    type(fgsl_matrix) :: dff
    real(fgsl_double), pointer :: par(:), xv(:), yv(:), yd(:,:)
    integer(fgsl_int) :: status

    real :: T,muB,lambda
    real, dimension(1:3,1:3) :: J

    call fgsl_obj_c_ptr(fx, x)
    call fgsl_obj_c_ptr(ff, f)
    call fgsl_obj_c_ptr(dff, df)
    call c_f_pointer(params, par, (/ 3 /))
    status = fgsl_vector_align(xv, fx)
    status = fgsl_vector_align(yv, ff)
    status = fgsl_matrix_align(yd, dff)

    T = f_v2T(xv(1))
    muB = f_v2muB(xv(2))
    lambda = f_v2lambda(xv(3))

    write(*,*) '   _fdf:',T,muB,lambda

    call calcJacobian(T,muB,lambda, J)

    yv(1) = J(3,1)*lambda - par(1) ! e
    yv(2) = J(3,2)*lambda - par(2) ! rhoB
    yv(3) = J(3,3)*lambda - par(3) ! rho

    yd(1:3,1:3) = reshape( source = (/ &
         J(1,1)*f_d_v2T(xv(1)),         & ! de/dT
         J(2,1)*f_d_v2muB(xv(2)),       & ! de/dmuB
         J(3,1)*f_d_v2lambda(xv(3)),    & ! de/dlambda
         J(1,2)*f_d_v2T(xv(1)),         & ! drhoB/dT
         J(2,2)*f_d_v2muB(xv(2)),       & ! drhoB/dmuB
         J(3,2)*f_d_v2lambda(xv(3)),    & ! drhoB/dlambda
         J(1,3)*f_d_v2T(xv(1)),         & ! drho/dT
         J(2,3)*f_d_v2muB(xv(2)),       & ! drho/dmuB
         J(3,3)*f_d_v2lambda(xv(3)) /), & ! drho/dlambda
         shape = (/ 3, 3 /))

    funcToMinimize3_fdf = fgsl_success
  end function funcToMinimize3_fdf



  !************************************************************************
  !****s* ThermalModel3/runSolver2a
  ! NAME
  ! subroutine runSolver2a()
  ! PURPOSE
  ! The routine to run a 2 dimensional root finding algorithm: given
  ! (eps,rhoB), it finds (T,muB)
  !************************************************************************
  subroutine runSolver2a()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 2
    type(fgsl_multiroot_fdfsolver) :: solver
    type(fgsl_multiroot_function_fdf) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 500
    integer :: iIt


    solver = fgsl_multiroot_fdfsolver_alloc(fgsl_multiroot_fdfsolver_gnewton,nParam)

    param = (/ in_eps, in_rhoB /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_fdf_init( &
         funcToMinimize2a_f,funcToMinimize2a_df,funcToMinimize2a_fdf, &
         nParam,pParam)

    xv = (/ f_T2v(in_T), f_muB2v(in_muB) /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)

    ! setup solver:
    status = fgsl_multiroot_fdfsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fdfsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)

    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fdfsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fdfsolver_iterate(solver);
       if (status /= fgsl_success .or. iIt > nIt) exit

       status = fgsl_multiroot_test_residual(fvec, 1.0d-10)
       if (status == fgsl_success) exit

       write(*,*) iIt,xptr(1:2),fptr(1:2)
    end do

    call WriteStatus(status)

    out_T = f_v2T(xptr(1))
    out_muB = f_v2muB(xptr(2))
    out_lambda = 1.0

  end subroutine runSolver2a

  !************************************************************************
  !****s* ThermalModel3/runSolver2b
  ! NAME
  ! subroutine runSolver2b()
  ! PURPOSE
  ! The routine to run a 2 dimensional root finding algorithm: given
  ! (eps,rho), it finds (T,muB)
  !************************************************************************
  subroutine runSolver2b()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 2
    type(fgsl_multiroot_fdfsolver) :: solver
    type(fgsl_multiroot_function_fdf) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 500
    integer :: iIt


    solver = fgsl_multiroot_fdfsolver_alloc(fgsl_multiroot_fdfsolver_gnewton,nParam)

    param = (/ in_eps, in_rho /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_fdf_init( &
         funcToMinimize2b_f,funcToMinimize2b_df,funcToMinimize2b_fdf, &
         nParam,pParam)

    xv = (/ f_T2v(in_T), f_muB2v(in_muB) /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)

    ! setup solver:
    status = fgsl_multiroot_fdfsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fdfsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)

    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fdfsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fdfsolver_iterate(solver);
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
  !****s* ThermalModel3/runSolver3
  ! NAME
  ! subroutine runSolver3()
  ! PURPOSE
  ! The routine to run a 3 dimensional root finding algorithm: given
  ! (eps,rhoB,rho), it finds (T,muB,lambda)
  !************************************************************************
  subroutine runSolver3()

    use fgsl ! too much stuff for 'only'

    integer(fgsl_size_t),parameter :: nParam = 3
    type(fgsl_multiroot_fdfsolver) :: solver
    type(fgsl_multiroot_function_fdf) :: func
    real(fgsl_double), target :: param(nParam), xv(nParam)
    type(fgsl_vector) :: xvec, fvec
    real(fgsl_double), pointer :: fptr(:), xptr(:)
    type(c_ptr) :: pParam
    integer(fgsl_int) :: status
    integer, parameter :: nIt = 500
    integer :: iIt

!    solver = fgsl_multiroot_fdfsolver_alloc(fgsl_multiroot_fdfsolver_hybridsj,nParam)
!    solver = fgsl_multiroot_fdfsolver_alloc(fgsl_multiroot_fdfsolver_newton,nParam)
    solver = fgsl_multiroot_fdfsolver_alloc(fgsl_multiroot_fdfsolver_gnewton,nParam)


    param = (/ in_eps, in_rhoB, in_rho /)
    pParam = c_loc(param)
    func = fgsl_multiroot_function_fdf_init( &
         funcToMinimize3_f,funcToMinimize3_df,funcToMinimize3_fdf, &
         nParam,pParam)

    xv = (/ f_T2v(in_T), f_muB2v(in_muB), f_lambda2v(in_lambda) /) ! start values

    ! create vector of values from array:
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xv,nParam,xvec,nParam,0_fgsl_size_t,1_fgsl_size_t)

    ! setup solver:
    status = fgsl_multiroot_fdfsolver_set(solver, func, xvec)

    fvec = fgsl_multiroot_fdfsolver_f(solver)
    status = fgsl_vector_align(fptr, fvec)

    call fgsl_vector_free(xvec)
    xvec = fgsl_multiroot_fdfsolver_root(solver)
    status = fgsl_vector_align(xptr, xvec)

    iIt = 0
    do
       iIt = iIt+1
       status = fgsl_multiroot_fdfsolver_iterate(solver);
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
  !****s* ThermalModel3/WriteStatus
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
!    logical :: isnan

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

       call calcENB(in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho)
       write(*,*) 'IN: ',in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho
       write(421,*) in_T,in_muB,in_lambda, out_eps,out_rhoB,out_rho


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

          call calcENB(T,mu,1.0, out_eps,out_rhoB,out_rho)
          write(499,*) T,mu,1.0, out_eps,out_rhoB,out_rho
       end do
       write(499,*)

    end do

  end subroutine runMode99


  subroutine invertJ(J)

    use fgsl ! too much stuff for 'only'

    real, dimension(3,3), intent(in) :: J

    real(fgsl_double), target :: hJ(3,3), invf(3,3)
    integer(fgsl_int) :: status, signum
    type(fgsl_matrix) :: a,inv
    type(fgsl_permutation) :: p

    a = fgsl_matrix_init(1.0_fgsl_double)
    status = fgsl_matrix_align(hJ, 3_fgsl_size_t, 3_fgsl_size_t, 3_fgsl_size_t,a)

    inv = fgsl_matrix_init(1.0_fgsl_double)
    status = fgsl_matrix_align(invf, 3_fgsl_size_t, 3_fgsl_size_t, 3_fgsl_size_t,inv)

    p = fgsl_permutation_alloc(3_fgsl_size_t)

    hJ = J

    status = fgsl_linalg_LU_decomp(a, p, signum)
    status = fgsl_linalg_lu_invert(a, p, inv)

    write(*,*) 'J='
    write(*,*) J(:,1)
    write(*,*) J(:,2)
    write(*,*) J(:,3)

    write(*,*) 'J^-1='
    write(*,*) invf(:,1)
    write(*,*) invf(:,2)
    write(*,*) invf(:,3)

  end subroutine invertJ

end program ThermalModel3
