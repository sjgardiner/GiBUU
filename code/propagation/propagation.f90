!******************************************************************************
!****m* /propagation
! NAME
! module propagation
! PURPOSE
! Module which includes the propagation of the test-particles.
!******************************************************************************
module propagation

  use CallStack, only: Traceback
  use hist2D

  implicit none

  private

  !****************************************************************************
  !****g* propagation/Mode
  ! SOURCE
  !
  integer, save :: Mode = 2
  ! PURPOSE
  ! define the type of propagation:
  ! * 0: Cascade
  ! * 1: Euler
  ! * 2: PredictorCorrector
  ! * 3: no propagation, random placement in box
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/delta_E
  ! PURPOSE
  ! Delta energy in derivatives
  ! SOURCE
  !
  real, save :: delta_E=0.01
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/delta_P
  ! PURPOSE
  ! Delta Momentum in derivatives
  ! SOURCE
  !
  real, save :: delta_P=0.01
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/UseCoulombDirectly
  ! PURPOSE
  ! Whether to use coulomb force directly in propagation or not. (If switched
  ! off while coulomb is switched on in module coulomb, the effect of the
  ! coulomb potential comes in via the gradient of the potentials. With this
  ! flag you can not switch on/off coulomb, you just select, how it is
  ! treated.)
  ! SOURCE
  !
  logical,  save :: UseCoulombDirectly=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/UseHadronic
  ! PURPOSE
  ! Whether to use hadronic potentials in propagation
  ! SOURCE
  !
  logical,  save :: UseHadronic=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/FreezeNonint
  ! SOURCE
  !
  logical, save :: FreezeNonint=.false.
  ! PURPOSE
  ! If switched on, the real particles which did not interact will have zero
  ! velocities, i.e. will be "frozen". This is important for stability of the
  ! nuclear ground state in real particle simulations.
  ! Note that this flag influences only when freezeRealParticles=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/dh_dp0_switch
  ! PURPOSE
  ! Switch which decides whether we use dh_dp0.
  ! SOURCE
  !
  logical,  save :: dh_dp0_switch=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/RungeKuttaOrder
  ! PURPOSE
  ! Order of Runge-Kutta in derivatives:
  ! * 1 = first order Runge-Kutta
  ! * 2 = second order Runge-Kuttay
  ! SOURCE
  !
  integer,  save :: RungeKuttaOrder=1
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/offShellInfo
  ! PURPOSE
  ! print out offShellInfo: set to .true. automatically if offShellTransport
  ! is used
  ! SOURCE
  !
  logical,save :: offShellInfo=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/offShellInfoDetail
  ! PURPOSE
  ! print out detailed offShellInfo
  ! SOURCE
  !
  logical,save :: offShellInfoDetail=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation/tachyonDebug
  ! PURPOSE
  ! ...
  ! SOURCE
  !
  logical,save :: tachyonDebug=.false.
  !****************************************************************************

  logical,save :: startFlag = .true.
  logical,parameter :: warnNegativeMass2 = .false.

  !****************************************************************************
  !****s* propagation/updateVelocity
  ! NAME
  ! subroutine updateVelocity(teilchen)
  ! PURPOSE
  ! Updates velocities of a single particle or a field of particles
  ! and checks for v<c.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: teilchen
  ! or:
  ! * type(particle), dimension(:)   :: teilchen
  ! or:
  ! * type(particle) :: teilchen
  !****************************************************************************
  Interface updateVelocity
     Module Procedure updateVelocity_1, updateVelocity_field, &
          updateVelocity_matrix
  End Interface updateVelocity

  ! Tachyon histograms:
  type(histogram2D),save :: hist2D_rho,hist2D_omega,hist2D_phi
  logical, save :: doInitTachyonHists = .true.
  logical, save :: doWriteTachyonHists = .false.

  public :: propagate
  public :: propagate_euler
  public :: propagate_cascade
  public :: updateVelocity
  public :: checkVelo,checkVelos
  public :: gradients

contains


  !****************************************************************************
  !****s* propagation/init
  ! NAME
  ! subroutine init
  ! PURPOSE
  ! Reads out namelist "propagation"
  ! INPUTS
  ! * (none)
  ! OUTPUT
  ! * Initializes global module variables
  !****************************************************************************
  subroutine init
    use output
    use offshellPotential

    integer :: ios

    !**************************************************************************
    !****n* propagation/Propagation
    ! NAME
    ! NAMELIST propagation
    ! PURPOSE
    ! Namelist which includes the input switches:
    ! * delta_P
    ! * delta_E
    ! * UseCoulombDirectly
    ! * UseHadronic
    ! * FreezeNonint
    ! * RungeKuttaOrder
    ! * Mode
    ! * dh_dp0_switch
    ! * offShellInfoDetail
    ! * tachyonDebug
    !**************************************************************************
    NAMELIST /propagation/ delta_P, delta_E, UseCoulombDirectly, UseHadronic, &
         FreezeNonint,RungeKuttaOrder, Mode, &
         dh_dp0_switch, offShellInfoDetail, tachyonDebug

    character(22), dimension(0:4), parameter :: cMode = (/ &
         'Cascade           ',&
         'Euler             ',&
         'PredictorCorrector',&
         'Random Placement  ',&
         'Cascade+hard core ' /)


    call Write_ReadingInput('propagation',0)
    rewind(5)
    read(5,nml=propagation,iostat=ios)
    if (ios>0) then
       write(*,*)
       write(*,*) 'Maybe you use an old jobcard.'
       write(*,*) 'Variabes have been renamed:'
       write(*,*) '   coulomb --> UseCoulombDirectly'
       write(*,*) '   hadronic --> Usehadronic'
       write(*,*) '   predictorCorrector --> Mode'
       write(*,*) '   DerivativeType --> RungeKuttaOrder'
    end if
    call Write_ReadingInput('propagation',0,ios)

    write(*,*) 'delta_P = ',delta_P
    write(*,*) 'delta_E = ',delta_E
    write(*,*) 'use coulomb force directly = ',UseCoulombDirectly
    write(*,*) 'hadronic potential flag    = ',UseHadronic
    write(*,*) 'FreezeNonint               = ',FreezeNonint
    write(*,*) 'RungeKuttaOrder            = ',RungeKuttaOrder
    select case (Mode)
    case (0:3)
       write(*,*) 'Propagation mode           : ',Mode,"= ",&
            cMode(Mode)
    case default
       write(*,*) 'Propagation mode           : ',Mode,"= WRONG!"
       call traceBack('wrong input value')
    end select
    write(*,*) 'dh_dp0_switch: ',dh_dp0_switch

    if (tachyonDebug) then
       write(*,*) 'Do tachyon debug!!'
    end if

    if (get_useOffShellPotentialBaryons().or.get_useOffShellPotentialMesons()) then
       offShellInfo=.true.
       if (tachyonDebug) offShellInfoDetail=.true.
       write(*,*) 'Detailed off-shell info?', offshellInfoDetail
    end if

    call Write_ReadingInput('propagation',1)

    startFlag = .false.
  end subroutine init




  !****************************************************************************
  !****s* propagation/gradients
  ! NAME
  ! subroutine gradients(part,Grad_P,Grad_R,flagOk)
  ! PURPOSE
  ! determine gradients of the Hamiltonfunction
  ! INPUTS
  ! * type(particle)  :: part  --
  !   particle whose gradients we want to calculate
  ! OUTPUT
  ! * real, dimension(1:3) :: Grad_P -- momentum gradient
  ! * real, dimension(0:3) :: Grad_R -- space gradients =(d/dt , d/dr(1:3) )
  !   [Note that there is no minus sign!!]
  ! * logical, optional :: flagOk -- indicates failure
  ! NOTES
  ! This routine is also called in RMF mode for HeavyIon init for rho, omega
  ! and phi mesons. In this case EnergyDeterminationRMF should be called
  !****************************************************************************
  subroutine gradients(partIn,Grad_P,Grad_R,flagOk)

    use particleDefinition
    use coulomb, only: emfoca
    use densityModule
    use energyCalc
    use offShellPotential
    use hist
    use idTable, only: photon
    use derivatives, only: finiteDifference
    use minkowski, only: SP
    use ieee_arithmetic, only: ieee_is_nan
    use RMF, only: getRMF_flag

    type(particle),intent(in)  :: partIN
    real, intent(out), dimension(1:3) :: Grad_P
    real, intent(out), dimension(0:3),optional :: Grad_R
    logical, intent(out), optional :: flagOk

    real :: delta_R(1:3)
    real,dimension(-2:2,0:3) :: energy
    type(particle) :: part
    integer :: i,j
    real :: cpot                     ! Coulomb potential
    real, dimension(1:3) :: emForce  ! electromagnetic Force
    real :: dH_dp0           ! Derivative with respect to p_0: dH/dp_0
    type(histogram),save :: dh_dp0_hist
    type(histogram2D),save :: dh_dp0_hist2D,dh_dp0_hist2Dmu2
    integer, save :: numCalls=0
    logical, save :: DoInit=.true.
    logical :: outOfBounds!,success
    integer :: deriv_scheme
    real :: mu2
    logical :: isOk

    if (present(flagOk)) flagOk = .true.

    numCalls=numCalls+1

    if (startFlag) call init

    if (get_offshell_debug()) then
       if (DoInit) then
          call createHist(dH_dp0_hist,'dH/dp_0',-2.,5.,0.1)
          call createHist2D(dH_dp0_hist2D,'Num events as function of dH/dp_0 and mass',(/0.0,-2./),(/1.3,5./),(/0.02,0.1/))
          call createHist2D(dH_dp0_hist2Dmu2,'Num events as function of dH/dp_0 and mu2',(/0.0,-2./),(/1.3,5./),(/0.02,0.1/))
          DoInit=.false.
       end if
    end if

    ! Check Input
    if (partIn%ID <= 0) then !no particle
       write(*,*) 'ID =',partIn%ID
       call TRACEBACK('Error. Particle-ID<=0  in Gradient.')
    else if (partIn%ID.eq.photon) then
       ! Photons do not undergo any potentials
       if (present(grad_R)) grad_R=0.
       grad_P(1:3)  =  partIn%mom(1:3)/partIn%mom(0)
       if (Dot_Product(grad_P,grad_P).gt.1.0000000001) then
          write(*,*) 'Problem with photon in gradients'
          write(*,*) '|Velo|^2=',Dot_Product(grad_P,grad_P)
          write(*,*) 'Velo=',grad_P
          write(*,*) 'Momentum=',partIn%mom
          return
       end if
    end if

    grad_P=0.

    if (getRMF_flag()) then
       if (UseHadronic .and. (.not.UseCoulombDirectly)) then
          call TRACEBACK('RMF works only with UseCoulombDirectly')
       end if
    end if

    if (present(grad_R)) then
       grad_R=0.
       if (UseCoulombDirectly) then ! special treatment of electromagnetic part
          cpot = emfoca(partIn%pos(1:3),partIn%mom(1:3),partIn%charge,partIn%ID,emForce)
          grad_R(1:3)=grad_R(1:3)-emForce  !Grad(Potential)=-Force
       end if
    end if

    if (UseHadronic) then ! with hadronic potentials

       if (present(grad_R)) then
          delta_R(1:3)=getGridSpacing() ! gridsize in position space

          ! Derivative with respect to r: grad_r H
          part=partIN
          energy=0.
          do i=1,3 !loop over x,y,z
             deriv_scheme=0
             part%pos=partIn%pos
             do j=-RungeKuttaOrder,RungeKuttaOrder
                part%pos(i)=partIN%pos(i)+float(j)*delta_R(i)
                ! determine energy due to hadronic interaction
                if (treatParticleOffShell(part%Id,part%offshellPar)) then
                   energy(j,i)=hamiltonFunc_offshell(part,outOfBounds,&
                        flagOk=isOk)
                   if (.not.isOk) then
                      call treatNotOk('r')
                      return
                   end if
                   if (outOfBounds) then
                      call OutOfBoundsMessage('r')
                      deriv_scheme=-j
                   end if
                else
                   if (getRMF_flag()) then
                      call energyDeterminationRMF(part)
                      ! see error message above for (.not.UseCoulombDirectly)
                   else
                      call energyDetermination(part,&
                           ForbidCoulomb=UseCoulombDirectly,skipTestE0=.true.)
                   end if
                   energy(j,i)=part%mom(0)
                end if
             end do
             grad_R(i)=grad_R(i)+finiteDifference(energy(:,i),delta_R(i),RungeKuttaOrder,deriv_scheme)
          end do
       end if



       ! Derivative with respect to p: grad_p H
       part=partIN
       energy=0.

       do i=1,3 !loop over x,y,z
          deriv_scheme=0
          part%mom=partIn%mom
          do j=-RungeKuttaOrder,RungeKuttaOrder
             part%mom(i)=partIN%mom(i)+float(j)*delta_P
             ! determine energy due to hadronic interaction
             if (treatParticleOffShell(part%Id,part%offshellPar)) then
                energy(j,i)=hamiltonFunc_offshell(part,outOfBounds, &
                     flagOk=isOk)
                if (.not.isOk) then
                   call treatNotOk('p')
                   return
                end if
                if (outOfBounds) then
                   call OutOfBoundsMessage('p')
                   deriv_scheme=-j
                end if
             else
                if (getRMF_flag()) then
                   call energyDeterminationRMF(part)
                   ! see error message above for (.not.UseCoulombDirectly)
                else
                   call energyDetermination(part,&
                        ForbidCoulomb=UseCoulombDirectly,skipTestE0=.true.)
                end if
                energy(j,i)=part%mom(0)
             end if
          end do
          grad_P(i)=grad_P(i)+finiteDifference(energy(:,i),delta_P,RungeKuttaOrder,deriv_scheme)
       end do



       !grad_R(0)=0.

       if (treatParticleOffShell(partIn%Id,partIn%offshellPar).and.dh_dp0_switch) then
          !********************************************************************
          ! Derivative with respect to p_0: dH/dp_0
          part=partIN
          energy=0.
          part%mom=partIn%mom
          deriv_scheme=0
          do j=-RungeKuttaOrder,RungeKuttaOrder
             part%mom(0)=partIN%mom(0)+float(j)*delta_E
             ! determine energy due to hadronic interaction
             energy(j,0)=hamiltonFunc_offshell(part,outofBounds, &
                  flagOk=isOk)
             if (.not.isOk) then
                call treatNotOk('E')
                return
             end if
             if (outOfBounds) then
                call OutOfBoundsMessage('E')
                deriv_scheme=-j
             end if
          end do
          dH_dp0=finiteDifference(energy(:,0),delta_E,RungeKuttaOrder,deriv_scheme)

          if (offShellInfoDetail) then
             write(*,*) 'dH_dp0=', dh_dp0
             if (ieee_is_nan(dh_dp0)) then
                write(*,*) "NaN in dH_dp0   ID=",partIn%Id, &
                     "     offShellParam=",partIn%offshellPar
                write(*,*) "energy=", energy(:,0)
                call Traceback()
             end if
          end if


          ! See Oliver's notes for the case of an energy dependent Hamiltonian:
          if (present(grad_R)) grad_R(1:3)=grad_R(1:3)/(1-dH_dp0)
          grad_P=grad_P/(1-dH_dp0)


          if (get_offshell_debug()) then

             mu2 = SP(partIN%mom,partIN%mom)
             call addHist(dH_dp0_hist,dH_dp0,1.)
             call AddHist2D(dH_dp0_hist2D, (/part%mass,dH_dp0/),1.)
             call AddHist2D(dH_dp0_hist2Dmu2, (/mu2,dH_dp0/),1.)

             if (mod(numCalls,1000).eq.0) then
                call writeHist(dh_dp0_hist,44,mul=dH_dp0_hist%xbin,file='dH_dp0.dat')
                call writeHist2D_Gnuplot(dh_dp0_hist2D,44,mul=dH_dp0_hist2D%xbin(1)*dH_dp0_hist2D%xbin(2),file='dH_dp0_2D.dat')
                call writeHist2D_Gnuplot(dh_dp0_hist2Dmu2,44,mul=dH_dp0_hist2D%xbin(1)*dH_dp0_hist2D%xbin(2),file='dH_dp0_2D_mu2.dat')
             end if
          end if
       end if
    else !no hadronic potentials
       !grad_R=0.
       grad_P=grad_P+partIN%mom(1:3)/FreeEnergy(partIN)
    end if

  contains

    subroutine OutOfBoundsMessage(C)

      character*(*), intent(in) :: C
      if (offShellInfoDetail) &
           write(*,'(A,A,i4)') C,': Hamilton function off shell is out of bounds!!',j
      if (j.eq.0) then ! Out of bounds at central value
         write(*,'(A,A,I4,A,I4)') C,': Problem with OutOfBounds: j=',j
      else if (deriv_scheme*j.gt.0) then ! Out of bounds on both sides
         write(*,'(A,A,I4,A,I4)') C,': Problem with OutOfBounds: deriv_scheme='&
              & ,deriv_scheme," new deriv_scheme=",-j
      end if

    end subroutine OutOfBoundsMessage

    subroutine treatNotOk(C)

      use output, only: DoPr
      use particleProperties, only: partName

      character*(*), intent(in) :: C

      if (present(flagOk)) flagOk = .false.

      if (DoPR(2)) write(*,'(A,A,A,A,i10)') &
           "WARNING: gradients '",C,"' not okay: ",&
           partName(partIn),partIn%number

      Grad_P = 0.
      if (present(Grad_R)) Grad_R = 0.

    end subroutine treatNotOk


  end subroutine gradients

  !****************************************************************************
  !****s* propagation/propagate
  ! NAME
  ! subroutine propagate(realPart, pertPart, delta_T, timeStep)
  ! PURPOSE
  ! This routine propagates the particle vectors.
  ! INPUTS
  ! * type(particle), dimension(:,:)  :: realPart
  !   -- real particle vector which should be propagated
  ! * type(particle), dimension(:,:)  :: pertPart
  !   -- perturbative particle vector which should be propagated
  ! * real    :: delta_T    -- time step size (fm/c)
  ! * integer :: timeStep   -- time step number
  ! OUTPUT
  ! * type(particle), dimension(:,:)  :: realPart
  ! * type(particle), dimension(:,:)  :: pertPart
  ! NOTES
  ! * Uses cascade, Euler or predictor-corrector time stepping.
  !****************************************************************************
  subroutine propagate(realPart, pertPart, delta_T, timeStep)
    use particleDefinition
    use inputGeneral, only: freezeRealParticles
    use densityModule, only: updateDensity
    use yukawa, only: updateYukawa
    use coulomb, only: updateCoulomb
    use offshellPotential, only: get_useOffShellPotentialMesons, get_useOffShellPotentialBaryons

    type(particle),intent(inOut),dimension(:,:) :: realPart
    type(particle),intent(inOut),dimension(:,:) :: pertPart
    real, intent(in) :: delta_T
    integer, intent(in) :: timeStep

    logical :: doReal
    ! Gradients of predictor step
    real, dimension(:,:,:), save, Allocatable :: gradR_real, gradP_real
    real, dimension(:,:,:), save, Allocatable :: gradR_pert, gradP_pert

    if (startFlag) call init

    doReal = (.not.freezeRealParticles) .or. (timeStep<=1)

    select case (Mode)
    case (0) ! === Cascade ===

       if (doReal) call propagate_cascade(realPart,delta_T)
       call propagate_cascade(pertPart,delta_T)

    case (1) ! === Euler ===

       if (doReal) call propagate_euler(realPart,delta_T)
       call propagate_euler(pertPart,delta_T)

    case (2) ! === PredictorCorrector ===

       ! Allocate fields to store the gradients
       if (doReal) then
          if(.not.allocated(gradR_real)) &
               allocate(gradR_real(0:3, &
               lbound(realPart,dim=1):ubound(realPart,dim=1), &
               lbound(realPart,dim=2):ubound(realPart,dim=2)))
          if(.not.allocated(gradP_real)) &
               allocate(gradP_real(1:3, &
               lbound(realPart,dim=1):ubound(realPart,dim=1), &
               lbound(realPart,dim=2):ubound(realPart,dim=2)))
       end if
       if(.not.allocated(gradR_pert)) &
            allocate(gradR_pert(0:3, &
            lbound(pertPart,dim=1):ubound(pertPart,dim=1), &
            lbound(pertPart,dim=2):ubound(pertPart,dim=2)))
       if(.not.allocated(gradP_pert)) &
            allocate(gradP_pert(1:3, &
            lbound(pertPart,dim=1):ubound(pertPart,dim=1), &
            lbound(pertPart,dim=2):ubound(pertPart,dim=2)))

       ! (1) Predictor Step: Define the predicted value of the particle at Delta_T

       if (doReal) call predictorStep(realPart,GradP_real,GradR_real,delta_T)
       call predictorStep(pertPart,GradP_pert,GradR_pert,delta_T)

       ! (2) Update potentials
       if (doReal) then
          call updateDensity(realPart)
          call updateCoulomb
          call updateYukawa(.false.)
       end if

       ! (3) Corrector step: consider also the predicted gradients
       if (doReal) call correctorStep(realPart,GradP_real,GradR_real,delta_T)

       call correctorStep(pertPart,GradP_pert,GradR_pert,delta_T)

    case (3)  ! === Random Placement ===

       if (doReal) call propagate_random(realPart)
       call propagate_random(pertPart)

    end select

    !  ** Set masses properly:
    if (get_useOffShellPotentialBaryons() .or. get_useOffShellPotentialMesons()) then
       if (doReal) call setMass(realPart)
       call setMass(pertPart)
    end if

    if (doReal) call checkVelos(realPart)
    call checkVelos(pertPart)

  end subroutine propagate


  !****************************************************************************
  !****s* propagation/setMass
  ! NAME
  ! subroutine setMass(Parts)
  ! PURPOSE
  ! * Resets %mass according to the offshell parameter and the actual
  !   momentum and position of the particle.
  ! * if particles is NOT meant to be treated offshell, nothing is changed.
  ! INPUTS
  ! * type(particle),dimension(:,:) :: Parts
  !
  ! OUTPUT
  ! * type(particle),dimension(:,:) :: Parts
  !****************************************************************************
  subroutine setMass(Parts)
    use particleDefinition
    use offshellPotential, only: treatParticleOffShell
    use minkowski, only: SP
    use potentialMain, only: massDetermination, scapot

    type(particle),dimension(:,:),intent(inout),target :: Parts

    integer :: iEns,iPart
    logical :: success
    type(particle), pointer :: pPart

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID.le.0) cycle

          if (.not.treatParticleOffShell(pPart%Id,pPart%offshellPar)) cycle

          call massDetermination(pPart,success=success)

          if (.not.success) then
             write(*,*) "Can't set baremass in setMass!!!"
             write(*,*) "Particle is deleted!!!"
             call protocolFailure(pPart,'setMass1')
             pPart%ID=0
             cycle
          end if
          if (SP(pPart%mom,pPart%mom)<=0.01**2) then ! 10 MeV
             write(*,*)
             write(*,*) 'WARNING in setMass: p^mu p_mu too small! ',SP(pPart%mom,pPart%mom)
             write(*,*) 'particle is deleted!!'
             call protocolFailure(pPart,'setMass2')
             pPart%ID = 0 ! delete particle
             cycle
          end if
          if (pPart%mass<0.01) then ! 10MeV
             write(*,*)
             write(*,*) 'WARNING in setMass: baremass too small! ',pPart%mass
             write(*,*) 'particle is deleted!!'
             call protocolFailure(pPart,'setMass3')
             pPart%ID = 0 ! delete particle
             cycle
          else if (pPart%mass+scapot(pPart)>pPart%mom(0)) then ! inv. mass > energy !
             write(*,*)
             write(*,*) 'WARNING in setMass: mass too large! ',pPart%mass
             write(*,*) 'particle is deleted!!'
             call protocolFailure(pPart,'setMass4')
             pPart%ID = 0 ! delete particle
             cycle
          end if

       end do
    end do
  end subroutine setMass


  !****************************************************************************
  !****s* propagation/predictorStep
  ! NAME
  ! subroutine predictorStep(Parts,gradP,gradR,delta_T)
  ! PURPOSE
  ! propagate particle vector using an Euler method and returns the gradients.
  ! INPUTS
  ! * real, intent(in) :: delta_T
  !
  ! OUTPUT
  ! * type(particle),dimension(:,:) :: Parts --
  !   vector at predicted position and momentum
  ! * real, dimension(:,:,:) :: gradR, gradP ! Gradients in R and P-Space
  !
  ! NOTES
  ! * Used as  predictor step in a predictor/corrector scheme
  !****************************************************************************
  subroutine predictorStep(Parts,gradP,gradR,delta_T)
    use offshellPotential, only: get_useOffShellPotentialBaryons,&
         get_useOffShellPotentialMesons,SetOffShellEnergy
    use particleDefinition
    use minkowski, only: SP
    use output
    use IdTable, only: isBaryon

    real, dimension(1:,:,:),intent(out) :: gradP
    real, dimension(0:,:,:),intent(out) :: gradR
    type(particle),intent(inOut),dimension(:,:),target :: Parts

    real, intent(in) :: delta_T
    integer :: iEns,iPart
    type(particle), pointer :: pPart

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          if(FreezeNonint .and. .not.pPart%pert .and.pPart%event(1)<1000000) then
             pPart%vel=0.
             cycle
          end if

          ! Evaluate dp/dt=-dH/dr and dr/dt=dH/dp
          call gradients(pPart,GradP(:,iEns,iPart),GradR(:,iEns,iPart))

          ! Set the perturbative particle vector to the predictor values:
          pPart%mom(1:3)=pPart%mom(1:3)-delta_T*gradR(1:3,iEns,iPart)
          pPart%mom( 0 )=pPart%mom( 0 )+delta_T*gradR(0,iEns,iPart)
          pPart%pos(1:3)=pPart%pos(1:3)+delta_T*gradP(1:3,iEns,iPart)
          pPart%vel(1:3)=gradP(1:3,iEns,iPart)

          call SetOffShellEnergy(pPart)

          if (SP(pPart%mom,pPart%mom) < -epsilon(0.0D0)) then
             if (warnNegativeMass2) then
                write(*,*) 'problems in predictorstep, SP(p,p).lt.0:', &
                     SP(pPart%mom,pPart%mom)
                write(*,*) 'gradR = ', gradR(1:3,iEns,iPart)
                call writeParticle_debug(pPart)
             end if
             if (isBaryon(pPart%ID)) then
                write(*,*) 'Deleting Particle'
                call protocolFailure(pPart,'predictor')
                pPart%ID=0
             end if
          end if

       end do
    end do

    if (get_useOffShellPotentialBaryons().or.get_useOffShellPotentialMesons()) &
         call setMass(Parts)

  end subroutine predictorStep



  !****************************************************************************
  !****s* propagation/correctorStep
  ! NAME
  ! subroutine correctorStep(Parts,gradP,gradR,delta_T)
  !
  ! PURPOSE
  ! Performs a corrector step on a particle vector previously undergoing a
  ! simple Euler method.
  !
  ! INPUTS
  ! * real :: delta_T
  ! * real, dimension(:,:,:) :: gradR, gradP --
  !   gradients in R and P-Space of predictor step
  ! * type(particle),dimension(:,:) :: Parts --
  !   particle vector at predicted values
  !
  ! OUTPUT
  ! * type(particle),dimension(:,:) :: Parts -- particle vector at t+Delta(t)
  !
  ! NOTES
  ! * Used as  corrector step in a predictor/corrector scheme
  !****************************************************************************
  subroutine correctorStep(Parts,gradP,gradR,delta_T)
    use particleDefinition
    use densityModule, only: gridsize
    use inputGeneral, only: continousBoundaries
    use offshellPotential, only: SetOffShellEnergy

    real, dimension(1:,:,:),intent(in) :: gradP
    real, dimension(0:,:,:),intent(in) :: gradR
    type(particle),intent(inOut),dimension(:,:),target :: Parts

    real, intent(in) :: delta_T
    real, dimension(1:3) :: Grad_P_Predictor
    real, dimension(0:3) :: Grad_R_Predictor
    integer :: i
    integer :: iEns,iPart
    logical :: checkVelo_flag
    type(particle), pointer :: pPart

    ! (3) Corrector Step: Evaluate Gradients at predicted value

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          if(FreezeNonint .and. .not.pPart%pert .and.pPart%event(1)<1000000) then
             pPart%vel=0.
             cycle
          end if

          call gradients(pPart,Grad_P_Predictor,Grad_R_Predictor) ! Evaluate dH/dr and dH/dp

          ! (3) Do time stepping considering also the predicted gradients
          !
          !       momentum(1:3,t+delta t)=momentum(1:3,t)-delta_T*(grad_R+grad_R_Predictor)/2.
          !       momentum_predicted(1:3,t+delta t)=momentum(1:3,t)-delta_T*grad_R
          !
          ! => momentum(1:3,t+delta t)=momentum_Predicted(1:3)-delta_T*(-grad_R+grad_R_Predictor)/2.
          ! And be reminded that "Parts" is the predictor (therefore the - sign in front of gradR and gradP).
          !
          ! The same argument also holds for the position.
          !
          pPart%mom(1:3)=pPart%mom(1:3)-delta_T*(-gradR(1:3,iEns,iPart)+grad_R_Predictor(1:3))*0.5
          pPart%mom( 0 )=pPart%mom( 0 )+delta_T*(-gradR( 0, iEns,iPart)+grad_R_Predictor( 0 ))*0.5
          pPart%pos(1:3)=pPart%pos(1:3)+delta_T*(-gradP(1:3,iEns,iPart)+grad_P_Predictor(1:3))*0.5

          ! (4) Save velocity
          pPart%vel(1:3)=(gradP(1:3,iEns,iPart)+grad_P_Predictor(1:3))*0.5

          call SetOffShellEnergy(pPart)

          if (continousBoundaries) then
             do i=1,3
                if (pPart%pos(i).gt.gridsize(i)) then
                   pPart%pos(i)= pPart%pos(i)-2*gridsize(i)
                else if (pPart%pos(i).lt.-gridsize(i)) then
                   pPart%pos(i)= pPart%pos(i)+2*gridsize(i)
                end if
             end do
          end if

          if (tachyonDebug) checkVelo_flag=checkVelo(pPart)
       end do
    end do

  end subroutine correctorStep



  !****************************************************************************
  !****s* propagation/propagate_euler
  ! NAME
  ! subroutine propagate_euler(Parts,delta_T)
  ! PURPOSE
  ! propagate a particle vector
  ! INPUTS
  ! * type(particle),dimension(:,:)  :: Parts
  !   -- particle vector which should be propagated
  ! * real :: delta_T ! time step (fm/c)
  ! OUTPUT
  ! * type(particle),dimension(:,:)  :: Parts
  !   -- particle vector which should be propagated
  ! NOTES
  ! Uses simple Euler time stepping.
  ! See also files in "Documentation_Extra/propagation/".
  !****************************************************************************
  subroutine propagate_euler(Parts,delta_T)

    use particleDefinition
    use nucleusDefinition
    use densityModule, only: gridsize
    use inputGeneral, only: continousBoundaries

    type(particle),intent(inOut),dimension(:,:),target :: Parts
    real, intent(in) :: delta_T

    integer :: iEns,iPart
    type(particle), pointer :: pPart

    real,dimension(0:3) :: grad_R
    real,dimension(1:3) :: grad_P
    integer :: i
    logical :: checkVelo_flag

    if (startFlag) call init

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          if(FreezeNonint .and. .not.pPart%pert .and.pPart%event(1)<1000000) then
             pPart%vel=0.
             cycle
          end if

          call gradients(pPart,Grad_P,Grad_R) ! Evaluate dH/dr and dH/dp
          pPart%mom(1:3)=pPart%mom(1:3)-delta_T*grad_R(1:3)
          pPart%mom( 0 )=pPart%mom( 0 )+delta_T*grad_R( 0 )
          pPart%pos(1:3)=pPart%pos(1:3)+delta_T*grad_P(1:3)
          pPart%vel(1:3)=grad_P(1:3)
          if (tachyonDebug) checkVelo_Flag=checkVelo(pPart)

          if (continousBoundaries) then
             do i=1,3
                if (pPart%pos(i).gt.gridsize(i)) then
                   pPart%pos(i)= pPart%pos(i)-2*gridsize(i)
                else if (pPart%pos(i).lt.-gridsize(i)) then
                   pPart%pos(i)= pPart%pos(i)+2*gridsize(i)
                end if
             end do
          end if

       end do
    end do

  end subroutine propagate_euler

  !****************************************************************************
  !****s* propagation/propagate_cascade
  ! NAME
  ! subroutine propagate_cascade(Parts,delta_T)
  ! PURPOSE
  ! Routine propagates particles in case of cascade mode.
  ! Useful also for doing initial step in the RMF mode.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts
  !   -- particles which should be propagated
  ! * real :: delta_T                            -- time step (fm/c)
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: Parts
  ! NOTES
  ! Straight line trajectories; no momentum change (mean fields switched off).
  !****************************************************************************
  subroutine propagate_cascade(Parts,delta_T)
    use nucleusDefinition
    use particleDefinition
    use eventtypes, only: HeavyIon,InABox,InABox_pion,InABox_delta,Box,elementary
    use inputGeneral, only: eventtype, continousBoundaries
    use nucleus, only: getTarget, getProjectile
    use densityModule, only: gridsize

    type(particle),intent(inOut),dimension(:,:),target :: Parts
    real, intent(in) :: delta_T

    type(tNucleus), pointer :: proj, targ
    type(particle), pointer :: pPart
    integer :: iEns,iPart, i
    real, dimension(1:3) :: t_vel,p_vel
    logical :: defreezeFlag

    if (startFlag) call init

    defreezeFlag = eventtype==InABox .or. eventtype==InABox_pion &
         .or. eventtype==InABox_delta .or. eventtype==Box &
         .or. eventtype==elementary

    if (eventtype==HeavyIon) then

       targ => getTarget()
       proj => getProjectile()

       t_vel=targ%vel
       p_vel=proj%vel

    else

       t_vel=0.
       p_vel=0.

    end if

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          if (.not.defreezeFlag .and. .not.pPart%pert .and. 0<pPart%event(1) .and. pPart%event(1)<1000000) then

             ! The particle is assumed to be "frozen"
             ! until it collides with other particle, i.e.
             ! it propagates with the speed of either target or projectile.

             if (mod(pPart%event(1),2)==1) then
                pPart%vel(1:3) = t_vel(1:3)
             else
                pPart%vel(1:3) = p_vel(1:3)
             end if

          else

             pPart%vel(1:3) = pPart%mom(1:3) / FreeEnergy( pPart )

          end if

          pPart%pos(1:3) = pPart%pos(1:3) + delta_T * pPart%vel(1:3)

          if (continousBoundaries) then
             do i=1,3
                if (pPart%pos(i).gt.gridsize(i)) then
                   pPart%pos(i)= pPart%pos(i)-2*gridsize(i)
                else if (pPart%pos(i).lt.-gridsize(i)) then
                   pPart%pos(i)= pPart%pos(i)+2*gridsize(i)
                end if
             end do
          end if

       end do
    end do

  end subroutine propagate_cascade

  !****************************************************************************
  !****s* propagation/propagate_random
  ! NAME
  ! subroutine propagate_random(Parts)
  ! PURPOSE
  ! Routine places particles randomly inside the box.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts
  !   -- particles which should be propagated
  ! OUTPUT
  ! * type(particle), dimension(:,:) :: Parts
  ! NOTES
  ! Only useful for box calculations!
  !****************************************************************************
  subroutine propagate_random(Parts)
    use particleDefinition
    use inputGeneral, only: eventtype
    use eventtypes, only: Box
    use densityModule, only: gridsize
    use random, only: rnFlat

    type(particle),intent(inOut),dimension(:,:),target :: Parts

    integer :: iEns,iPart, i
    type(particle), pointer :: pPart

    if (startFlag) call init

    if (eventtype.ne.Box) then
       call Traceback("ranom placement only for box calculations!")
    end if

    do iEns=lbound(Parts,dim=1),ubound(Parts,dim=1)
       do iPart=lbound(Parts,dim=2),ubound(Parts,dim=2)

          pPart => Parts(iEns,iPart)
          if (pPart%ID < 0) exit
          if (pPart%ID == 0) cycle

          pPart%vel(1:3) = pPart%mom(1:3) / FreeEnergy( pPart )


          do i=1,3
             pPart%pos(i) = rnFlat(-gridsize(i),gridsize(i))
          end do

      end do
    end do

  end subroutine propagate_random

  !****************************************************************************
  !****f* propagation/checkVelo
  ! NAME
  ! logical function checkVelo(part,verbose)
  ! PURPOSE
  ! Checks the velocity of a given particle
  ! INPUTS
  ! * type(particle),intent(inout) :: part
  ! * logical, optional :: verbose -- flag to produce additional output
  ! NOTES
  ! * Also does some statistical stuff
  !****************************************************************************
  logical function checkVelo(part,verbose)

    use IdTable, only: photon
    use particleDefinition
    use minkowski, only: abs4
    use offshellPotential, only: treatParticleOffShell, get_offshell_debug, &
         get_useOffShellPotentialMesons
    use output
    use mediumDefinition
    use mediumModule, only: mediumAt

    type(particle),intent(inout) :: part
    logical, intent(in), optional :: verbose

    integer, save :: failures=0
    integer, save :: numtry=0

    type(medium) :: med
    logical :: verbose_

    checkVelo=.true.
    if (part%ID.le.0 .or. part%ID==photon) return
    numtry=numtry+1

    verbose_ = .true.
    if (present(verbose)) verbose_ = verbose

    if (startFlag) call init

    if (get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
       call addTachyonHists(part, 0.,1.)
    end if

    if (Dot_Product(part%vel(1:3),part%vel(1:3)) .ge. 1.) then

       failures=failures+1
       checkVelo=.false.

       med = mediumAt(part%pos)

       if (DoPR(2).and.verbose_) then
          write(*,*)
          write(*,'(A,G12.5)') 'Problems in CheckVelo : v = ',&
               sqrt( Dot_Product(part%vel(1:3),part%vel(1:3)))
          write(*,'(A,G12.5,A,G12.5,A,G12.5)') 'm = ',abs4(part%mom),&
               '; pabs = ',absMom(part), &
               '; rho = ',med%densityproton+med%densityneutron
          call WriteParticle(6,0,0,part)
          write(*,*) 'Number of failures :',failures,' (',float(failures)/float(numtry)*100., '%)'
       end if
       call protocolFailure(part,'checkVelo')

       if (get_offshell_debug() .and. get_useOffShellPotentialMesons()) then
          call addTachyonHists(part, 1.,0.)
       end if

       if (treatParticleOffShell(part%Id,part%offshellPar)) then
          !ACCEPTABLE ONLY FOR OFFSHELL TRANSPORT, OTHERWISE SERIOUS ERROR!!!!!!!!!!
          if (DoPR(2)) write(*,*) '--- this particle is now deleted! (1)'
          part%ID=0
       else
          call traceback('Stop!')
       end if
    end if

  end function checkVelo


  !****************************************************************************
  !****s* propagation/checkVelos
  ! NAME
  ! subroutine checkVelos(part)
  ! PURPOSE
  ! Checks the velocities
  ! INPUTS
  ! *type(particle),dimension(:,:),intent(inout) :: part
  !****************************************************************************
  subroutine checkVelos(part)
    use particleDefinition

    type(particle),dimension(:,:),intent(inout) :: part
    integer :: iEns,iPart
    logical :: dummy

    do iEns=lbound(part,dim=1),ubound(part,dim=1)
       do iPart=lbound(part,dim=2),ubound(part,dim=2)
          if (part(iEns,iPart)%ID<=0) cycle
          dummy=checkVelo(part(iEns,iPart))
       end do
    end do

    if (doWriteTachyonHists) call writeTachyonHists

  end subroutine checkVelos

  !****************************************************************************
  subroutine initTachyonHists
    if (doInitTachyonHists) then
       call createHist2D(hist2D_rho,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
       call createHist2D(hist2D_omega,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
       call createHist2D(hist2D_phi,'number of tachyons as function of mass and momentum',(/0.,0./),(/2.,2./),(/0.02,0.02/))
       doInitTachyonHists = .false.
       doWriteTachyonHists = .true.
    end if
  end subroutine initTachyonHists

  !****************************************************************************
  subroutine addTachyonHists(part, v1,v2)
    use particleDefinition
    use IdTable, only: rho,omegaMeson,phi
    use minkowski, only: abs4

    type(particle),intent(in) :: part
    real, intent(in) :: v1,v2
    integer, save :: nCall = 0


    if (doInitTachyonHists) call initTachyonHists

    select case (part%ID)
    case (rho)
       call AddHist2D(hist2D_rho,(/abs4(part%mom),absMom(part)/),v1,v2)
       doWriteTachyonHists = .true.
       nCall = nCall + 1
    case (omegaMeson)
       call AddHist2D(hist2D_omega,(/abs4(part%mom),absMom(part)/),v1,v2)
       doWriteTachyonHists = .true.
       nCall = nCall + 1
    case (phi)
       call AddHist2D(hist2D_phi,(/abs4(part%mom),absMom(part)/),v1,v2)
       doWriteTachyonHists = .true.
       nCall = nCall + 1
    end select

    if (mod(nCall,20)==0) then
       call writeTachyonHists
       nCall = 0
    end if

  end subroutine addTachyonHists

  !****************************************************************************
  subroutine writeTachyonHists
    if (doWriteTachyonHists) then
       call writeHist2D_Gnuplot(hist2D_rho,44,file='Tachyons2D_rho.dat')
       call writeHist2D_Gnuplot(hist2D_omega,45,file='Tachyons2D_omega.dat')
       call writeHist2D_Gnuplot(hist2D_phi,46,file='Tachyons2D_phi.dat')
       doWriteTachyonHists = .false.
    end if
  end subroutine writeTachyonHists


  !****************************************************************************
  ! cf. Interface updateVelocity
  !****************************************************************************
  subroutine updateVelocity_matrix(Parts)
    use particleDefinition
    use output, only: DoPR

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut), dimension(:,:) :: Parts
    integer :: i,j
    logical :: success

    if (startFlag) call init

    if (DoPr(2)) write(*,*) 'Updating particle velocities'

    ensLoop: do i=lbound(Parts, dim=1),ubound(Parts, dim=1)
       do j=lbound(Parts, dim=2),ubound(Parts, dim=2)
          if (Parts(i,j)%ID <= 0) cycle ensLoop
          call gradients(Parts(i,j),Grad_P) ! Evaluate dH/dp
          Parts(i,j)%vel(1:3)=grad_P
          success=checkVelo(Parts(i,j))
       end do
    end do ensLoop
  end subroutine updateVelocity_matrix

  !****************************************************************************
  ! cf. Interface updateVelocity
  !****************************************************************************
  subroutine updateVelocity_field(Parts)
    use particleDefinition

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut), dimension(:) :: Parts
    integer :: i
    logical :: success

    if (startFlag) call init

    do i=lbound(Parts, dim=1),ubound(Parts, dim=1)
       if (Parts(i)%ID <= 0) cycle
       call gradients(Parts(i),Grad_P) ! Evaluate dH/dp
       Parts(i)%vel(1:3)=grad_P
       success=checkVelo(Parts(i))
    end do
  end subroutine updateVelocity_field

  !****************************************************************************
  ! cf. Interface updateVelocity
  !****************************************************************************
  subroutine updateVelocity_1(part,success)
    use particleDefinition
    use output, only: writeParticle_debug

    real,dimension(1:3)  :: grad_P
    type(particle),intent(inOut) :: part
    logical,intent(out),optional :: success
    logical :: c

    if (startFlag) call init

    if (part%ID <= 0) return
    call gradients(part,Grad_P) ! Evaluate dH/dp
    part%vel(1:3)=grad_P

    c=checkVelo(part)

    if (present(success)) then
       success = c
       if (tachyonDebug.and..not.success) then
          call writeParticle_debug(part)
          stop
       end if
    end if
  end subroutine updateVelocity_1


  !****************************************************************************
  subroutine protocolFailure(part,kindOfFailure)
    use particleDefinition

    character(*),   intent(in) :: kindOfFailure
    type(particle), intent(in) :: part

    logical, save :: firstTime=.true.
    character(20) :: form
    integer, save :: numFailures=0
    form='(A40,4G13.5)'

    numFailures=numFailures+1

    if (firstTime) then
       open(22,file='propa_failures.txt')
       write(22,*)
       firstTime=.false.
    else
       open(22,file='propa_failures.txt',position='append')
    end if
    write(22,'(A,A)')'Problems in ', kindOfFailure
    write(22,form)'|v| = '                      , sqrt(Dot_Product(part%vel(1:3),part%vel(1:3)))
    write(22,form)'Particle ID: '               , part%ID
    write(22,form)'Particle mass: '             , part%mass
    write(22,form)'Particle offshellparameter: ', part%offshellPar
    write(22,form)'Particle number: '           , part%number
    write(22,form)'Particle firstevent: '       , part%firstevent
    write(22,form)'Particle velocity: '         , part%vel
    write(22,form)'Particle momentum: '         , part%mom
    write(22,form)'Particle perturbative?: '    , part%pert
    write(22,form)'Particle charge: '           , part%Charge
    write(22,form)'Particle position: '         , part%pos
    write(22,*)   'Number of failures: '        , numFailures
    write(22,*)
    close(22)
  end subroutine protocolFailure


end module propagation
