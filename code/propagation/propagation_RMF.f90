!******************************************************************************
!****m* /propagation_RMF
! NAME
! module propagation_RMF
! PURPOSE
! Module which includes the propagation of the test-particles
! in relativistic mean field. It substitutes the module propagation
! in case if RMF is used.
!******************************************************************************
module propagation_RMF

  implicit none

  private




  !****************************************************************************
  !****g* propagation_RMF/predictorCorrector
  ! PURPOSE
  ! Switch for predictor-corrector method in the propagation.
  ! If .false. then simple Euler method is used
  ! (i.e. only predictor step is done)
  ! SOURCE
  !
  logical, save :: predictorCorrector=.true.
  !****************************************************************************

  !****************************************************************************
  !****g* propagation_RMF/deleteTachyons
  ! PURPOSE
  ! Switch for treatment of particles with velocity > 1.
  !
  ! Possible values:
  ! * if .true., these particles are deleted by setting ID=0
  ! * if .false., these particles are propagated, but with modified velocity
  !
  ! SOURCE
  !
  logical, save :: deleteTachyons = .false.
  !****************************************************************************

  logical, save :: startFlag=.true.
  logical, parameter :: debugFlag=.false.

  public :: propagate_RMF ! subroutine to propagate particles

contains


  !****************************************************************************
  subroutine init
    use output

    integer :: ios

    !**************************************************************************
    !****n* propagation_RMF/propagation_RMF_input
    ! NAME
    ! NAMELIST propagation_RMF_input
    ! PURPOSE
    ! Namelist which includes the input switches:
    ! * predictorCorrector
    ! * deleteTachyons
    !**************************************************************************
    NAMELIST /propagation_RMF_input/ predictorCorrector, deleteTachyons

    call Write_ReadingInput('propagation_RMF_input',0)
    rewind(5)
    read(5,nml=propagation_RMF_input,iostat=ios)
    call Write_ReadingInput('propagation_RMF_input',0,ios)

    write(*,*) ' Predictor-Corrector mode: ',predictorCorrector
    write(*,*) ' Delete Tachyons:          ',deleteTachyons

    call Write_ReadingInput('propagation_RMF_input',1)
  end subroutine init

  !****************************************************************************
  !****s* propagation_RMF/propagate_RMF
  ! NAME
  ! subroutine propagate_RMF(realPart, pertPart, delta_T, TimeStep)
  ! PURPOSE
  ! This routine propagates particles in a relativistic mean field.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: realPart
  ! * type(particle), dimension(:,:) :: pertPart
  !   -- real and perturbative particle arrays which should be propagated
  ! * real :: delta_T -- size of time step (fm/c)
  ! * integer :: TimeStep -- number of time step
  ! NOTES
  ! Particle coordinates and kinetic momentum are propagated.
  ! Uses predictor corrector scheme or simple Euler time stepping.
  !****************************************************************************
  subroutine propagate_RMF(realPart, pertPart, delta_T, TimeStep)

    use particleDefinition
    use inputGeneral, only: freezeRealParticles
    use densitymodule, only: updateRMF, energyDeterminationRMF, storeFieldsRMF
    use coulomb, only: updateCoulomb
    use offShellPotential, only: treatParticleOffShell,SetOffShellEnergy
    use output, only: WriteParticle, DoPR
    use propagation, only: checkVelo, checkVelos
    use minkowski, only: abs4Sq
    use particleProperties, only: partName

    type(particle),intent(inOut),dimension(:,:),target :: realPart
    type(particle),intent(inOut),dimension(:,:),target :: pertPart
    real, intent(in) :: delta_T
    integer, intent(in) :: TimeStep

    ! Working variables:
    real, save, Allocatable, dimension(:,:,:) :: rk_re, pk_re, rk_pe, pk_pe    ! arrays to store the coordinates and momenta
                                                                               ! at the current time
    real, save, Allocatable, dimension(:,:,:) :: pdot_old_re, pdot_old_pe      ! arrays to store dp^*/dt calculated
                                                                               ! at the current time
    real, save, Allocatable, dimension(:,:,:) :: rdot_old_re, rdot_old_pe    ! arrays to store the particle velocities
                                                                                     ! at the current time

    type(particle), pointer :: pPart
    integer, save :: n_ensembles, n_particles_re, n_particles_pe
    integer :: i, j, errCode
    real, dimension(1:3) :: pdot, rdot
    real, save :: delta_T_old  ! Previous time step size
    logical :: flagOk

    if (startFlag) then !Initialisierung im ersten Aufruf
       call init
       n_ensembles = Size(realPart,dim=1)
       n_particles_re = Size(realPart,dim=2)
       n_particles_pe = Size(pertPart,dim=2)
       allocate(rk_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(rk_pe(1:n_ensembles,1:n_particles_pe,1:3))
       allocate(pk_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(pk_pe(1:n_ensembles,1:n_particles_pe,1:3))
       allocate(pdot_old_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(pdot_old_pe(1:n_ensembles,1:n_particles_pe,1:3))
       allocate(rdot_old_re(1:n_ensembles,1:n_particles_re,1:3))
       allocate(rdot_old_pe(1:n_ensembles,1:n_particles_pe,1:3))
       startFlag=.false.
    end if

    if ( TimeStep == 1 ) delta_T_old = delta_T

    !--------------------------------------------------------------------------
    ! 1) Predictor Step:
    !--------------------------------------------------------------------------

    Loop_over_ensembles_1 : do i=1,n_ensembles

       !-----------------------------------------------------------------------
       ! 1a) real particles
       !-----------------------------------------------------------------------

       if (.not.freezeRealParticles) then
          do j=1,n_particles_re

             pPart => realPart(i,j)
             if ( pPart%id == 0 ) cycle
             if ( pPart%id < 0  ) exit

             call SetOffShellEnergy(pPart, errCode)
             if (errCode>0) then
                call reportError(pPart, errCode)
                cycle
             end if

             call rhs(pPart,delta_T_old,pdot,rdot)

             ! if(pPart%id == 101) write(*,*)' pion velo: ', rdot

             if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
                if (DoPR(2)) write(*,'(A,A,i10,G12.3,i3)'), &
                     "RMF-predictor: real tachyon: ", &
                     partName(pPart), pPart%number, &
                     sqrt(dot_product(rdot(1:3),rdot(1:3))),errCode
                if (deleteTachyons) then
                   pPart%vel(1:3)=rdot(1:3)
                   flagOk = checkVelo(pPart,.false.) ! will set ID=0
                   cycle
                else
                   if (DoPR(2)) then
                      call WriteParticle(6,i,j,pPart)
                      write(*,*) '  offshellparameter :', &
                           pPart%offshellPar, abs4Sq(pPart%mom)
                   end if
                   call reportError(pPart, 6)

                   rdot(1:3)=pPart%mom(1:3)/pPart%mom(0)
                   pdot(1:3)=0.
                end if
             end if

             rk_re(i,j,1:3) = pPart%pos(1:3)
             pk_re(i,j,1:3) = pPart%mom(1:3)

             pPart%pos(1:3) = pPart%pos(1:3)+delta_T*rdot(1:3)
             pPart%mom(1:3) = pPart%mom(1:3)+delta_T*pdot(1:3)

             pPart%vel(1:3) = rdot(1:3)

             pdot_old_re(i,j,1:3) = pdot(1:3)
             rdot_old_re(i,j,1:3) = rdot(1:3)

             call SetOffShellEnergy(pPart, errCode)
             if (errCode>0) then
                call reportError(pPart, errCode)
                cycle
             end if

          end do
       end if


       !-----------------------------------------------------------------------
       ! 1b) perturbative particles
       !-----------------------------------------------------------------------

       do j=1,n_particles_pe

          pPart => pertPart(i,j)
          if ( pPart%id == 0 ) cycle
          if ( pPart%id < 0  ) exit

          call SetOffShellEnergy(pPart)
          call rhs(pPart,delta_T_old,pdot,rdot)

          if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
             if (DoPR(2)) then
                write(*,'(A,A,i10,G12.3)'), &
                     "RMF-predictor: perturbative tachyon: ", &
                     partName(pPart), pPart%number, &
                     sqrt(dot_product(rdot(1:3),rdot(1:3)))
                call WriteParticle(6,i,j,pPart)
                write(*,*) '  offshellparameter :', &
                     pPart%offshellPar, abs4Sq(pPart%mom)
             end if

             rdot(1:3)=pPart%mom(1:3)/pPart%mom(0)
             pdot(1:3)=0.
          end if

          rk_pe(i,j,1:3) = pPart%pos(1:3)
          pk_pe(i,j,1:3) = pPart%mom(1:3)

          pPart%pos(1:3) = pPart%pos(1:3)+delta_T*rdot(1:3)
          pPart%mom(1:3) = pPart%mom(1:3)+delta_T*pdot(1:3)

          call SetOffShellEnergy(pPart)
          pPart%vel(1:3)=rdot(1:3)

          pdot_old_pe(i,j,1:3) = pdot(1:3)
          rdot_old_pe(i,j,1:3) = rdot(1:3)

       end do

    end do Loop_over_ensembles_1


    ! Store the space components of the omega- and rho-field:
    if (.not.freezeRealParticles) call storeFieldsRMF()

    !--------------------------------------------------------------------------
    ! 2) Corrector step (optional):
    !--------------------------------------------------------------------------

    if (predictorCorrector) then

       if (.not.freezeRealParticles) then
          call updateRMF(realPart)
          call updateCoulomb()
       end if

       Loop_over_ensembles_2 : do i=1,n_ensembles

          !--------------------------------------------------------------------
          ! 2a) real particles
          !--------------------------------------------------------------------

          if (.not.freezeRealParticles) then
             do j=1,n_particles_re

                pPart => realPart(i,j)
                if ( pPart%id == 0 ) cycle
                if ( pPart%id < 0  ) exit

                call rhs(pPart,delta_T,pdot,rdot)

                rdot(1:3)=0.5*(rdot_old_re(i,j,1:3)+rdot(1:3))
                pdot(1:3)=0.5*(pdot_old_re(i,j,1:3)+pdot(1:3))

                if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then

                   if (DoPR(2)) write(*,'(A,A,i10,G12.3)'), &
                     "RMF-corrector: real tachyon: ", &
                     partName(pPart), pPart%number, &
                     sqrt(dot_product(rdot(1:3),rdot(1:3)))

                   if (deleteTachyons) then
                      pPart%vel(1:3)=rdot(1:3)
                      flagOk = checkVelo(pPart,.false.) ! will set ID=0
                      cycle
                   else
                      if (DoPR(2)) then
                         call WriteParticle(6,i,j,pPart)
                         write(*,*) '  offshellparameter :', &
                              pPart%offshellPar, abs4Sq(pPart%mom)
                      end if
                      call reportError(pPart, 7)

                      rdot(1:3)=pPart%mom(1:3)/pPart%mom(0)
                      pdot(1:3)=0.
                   end if
                end if

                pPart%pos(1:3) = rk_re(i,j,1:3)+delta_T*rdot(1:3)
                pPart%mom(1:3) = pk_re(i,j,1:3)+delta_T*pdot(1:3)

                pPart%vel(1:3)=rdot(1:3)

                call SetOffShellEnergy(pPart, errCode)
                if (errCode>0) then
                   call reportError(pPart, errCode)
                   cycle
                end if

             end do
          end if

          !--------------------------------------------------------------------
          ! 2b) perturbative particles
          !--------------------------------------------------------------------

          do j=1,n_particles_pe

             pPart => pertPart(i,j)
             if ( pPart%id == 0 ) cycle
             if ( pPart%id < 0  ) exit

             if(.not.treatParticleOffShell(pPart%id,pPart%offshellPar)) &
                  call energyDeterminationRMF( pPart )

             call rhs(pPart,delta_T,pdot,rdot)

             rdot(1:3)=0.5*(rdot_old_pe(i,j,1:3)+rdot(1:3))
             pdot(1:3)=0.5*(pdot_old_pe(i,j,1:3)+pdot(1:3))

             if(dot_product(rdot(1:3),rdot(1:3)).gt.1.) then
                if (DoPR(2)) then
                   write(*,'(A,A,i10,G12.3)'), &
                        "RMF-corrector: perturbative tachyon: ", &
                        partName(pPart), pPart%number, &
                        sqrt(dot_product(rdot(1:3),rdot(1:3)))
                   call WriteParticle(6,i,j,pPart)
                   write(*,*) '  offshellparameter :', &
                        pPart%offshellPar, abs4Sq(pPart%mom)
                end if

                rdot(1:3)=pPart%mom(1:3)/pPart%mom(0)
                pdot(1:3)=0.
             end if

             pPart%pos(1:3) = rk_pe(i,j,1:3)+delta_T*rdot(1:3)
             pPart%mom(1:3) = pk_pe(i,j,1:3)+delta_T*pdot(1:3)

             call SetOffShellEnergy(pPart)
             pPart%vel(1:3)=rdot(1:3)

          end do

       end do Loop_over_ensembles_2

    end if

    if (.not.freezeRealParticles) then
       call updateRMF(realPart)
       call updateCoulomb
    end if

    !--------------------------------------------------------------------------
    ! 3) set energies of perturbative on-shell particles:
    !--------------------------------------------------------------------------

    Loop_over_ensembles_3 : do i=1,n_ensembles
       do j=1,n_particles_pe

          pPart => pertPart(i,j)
          if ( pPart%id == 0 ) cycle
          if ( pPart%id < 0  ) exit

          if(treatParticleOffShell(pPart%id,pPart%offshellPar)) cycle

          call energyDeterminationRMF( pPart )
          pPart%vel(1:3) = pPart%mom(1:3)/pPart%mom(0)
       end do
    end do Loop_over_ensembles_3

    delta_T_old = delta_T

    call checkVelos(realPart)
    call checkVelos(pertPart)

  end subroutine propagate_RMF


  !****************************************************************************
  !****************************************************************************
  subroutine rhs(Part,dt,pdot,rdot)

    use particleDefinition
    use IdTable
    use densitymodule
    use RMF
    use coulomb, only: emfoca
    use propagation, only: gradients
    use eventtypes, only: HeavyIon
    use inputGeneral, only: eventtype
    use callstack, only: traceback

    real,                 intent(in) :: dt
    type(particle),       intent(in) :: Part
    real, dimension(1:3), intent(out) :: pdot, rdot

    real, dimension(0:3) :: V             ! vector component of the particle self-energy (vector field) [GeV]

    real, dimension(0:3,1:3) :: dV_dr     ! space derivatives of the vector field V (GeV/fm),
                                          ! 1-st index -> field component,
                                          ! 2-nd index -> space index of the derivative

    real,  dimension(1:3) :: dsigma_dr ! space derivatives of the scalar field (GeV/fm)

    real, dimension(1:3) :: Grad_P     ! auxiliary, needed for off-shell particles
    real, dimension(0:3) :: Grad_R     ! auxiliary, needed for off-shell particles

    integer :: i1, k, k_max,id,charge
    integer :: index1, index2, index3
    real, dimension(1:3) :: pdot_space, pdot_time, place, impuls, emForce
    real :: cpot, factor, mass
    real :: V1, V2, S1, S2
    logical :: anti, AddCoulomb, flagOk


    if (lorentz_flag) then
      k_max=3    ! With space components of the omega-field
    else
      k_max=0    ! W/o space components of the omega-field
    end if

    id = Part%id
    charge = Part%charge
    anti = Part%anti
    mass = Part%mass

    ! Modification factor for the coupling constants:
    factor = ModificationFactor(Id,anti)

    pdot_space(1:3) = 0.
    pdot_time(1:3) = 0.
    rdot = Part%vel
    place(1:3) = Part%pos(1:3)
    impuls(1:3) = Part%mom(1:3)

    AddCoulomb = .true.

    ! nuclear part (only for baryons, kaons and antikaons):
    if ( (isBaryon(id)  .or. id.eq.Kaon .or. id.eq.kaonBar) .and. factor.gt.0. ) then

        index1 = nint( place(1) / gridSpacing(1) )
        index2 = nint( place(2) / gridSpacing(2) )
        index3 = nint( place(3) / gridSpacing(3) )

        if ( abs(index1) <= gridPoints(1)-1 .and. &
             abs(index2) <= gridPoints(2)-1 .and. &
             abs(index3) <= gridPoints(3)-1 ) then     ! Inside grid:

           do k = 0,k_max

              V(k) = Selfenergy_vector(Index1,Index2,Index3,k,id,charge,anti)

              V1 = Selfenergy_vector(Index1+1,Index2,Index3,k,id,charge,anti)
              V2 = Selfenergy_vector(Index1-1,Index2,Index3,k,id,charge,anti)

              dV_dr(k,1) = ( V1 - V2) / (2.*gridSpacing(1))

              V1 = Selfenergy_vector(Index1,Index2+1,Index3,k,id,charge,anti)
              V2 = Selfenergy_vector(Index1,Index2-1,Index3,k,id,charge,anti)

              dV_dr(k,2) = ( V1 - V2) / (2.*gridSpacing(2))

              V1 = Selfenergy_vector(Index1,Index2,Index3+1,k,id,charge,anti)
              V2 = Selfenergy_vector(Index1,Index2,Index3-1,k,id,charge,anti)

              dV_dr(k,3) = ( V1 - V2) / (2.*gridSpacing(3))

           end do

           !isoscalar, scalar contributions (dm*/dr_i):

           S1 = Selfenergy_scalar(Index1+1,Index2,Index3,id,anti)
           S2 = Selfenergy_scalar(Index1-1,Index2,Index3,id,anti)

           dsigma_dr(1) = ( S1 - S2 ) / (2.*gridSpacing(1))

           S1 = Selfenergy_scalar(Index1,Index2+1,Index3,id,anti)
           S2 = Selfenergy_scalar(Index1,Index2-1,Index3,id,anti)

           dsigma_dr(2) = ( S1 - S2 ) / (2.*gridSpacing(2))

           S1 = Selfenergy_scalar(Index1,Index2,Index3+1,id,anti)
           S2 = Selfenergy_scalar(Index1,Index2,Index3-1,id,anti)

           dsigma_dr(3) = ( S1 - S2 ) / (2.*gridSpacing(3))

           if ( kaonpot_flag .and. (Part%id == kaon .or. Part%id == kaonBar) ) then
              do i1=1,3
                 dsigma_dr(i1) = -dsigma_dr(i1) &
                      &          + ( dot_product( V(0:k_max),dV_dr(0:k_max,i1) ) ) / (2.*gridSpacing(i1))
              end do
              dsigma_dr(1:3) = dsigma_dr(1:3) / DiracMass(Index1,Index2,Index3,mass,id,charge,anti)
           end if


           do i1 = 1,3

             !**** Space derivative contributions
             !**** to the rhs of equation dp^*/dt=pdot
             !**** pdot_space(i1) = -domega_dr(0,i1) -drho_dr(0,i1)

             !**** vector field contributions:
             pdot_space(i1) = -dV_dr(0,i1)

             if (lorentz_flag) then
               pdot_space(i1) = pdot_space(i1) &
                    & + dot_product( Part%vel(1:3), (dV_dr(1:3,i1)-dV_dr(i1,1:3)) )
             end if

             ! For antibaryons or antikaons the vector field changes sign:
             if ( anti .or. id.eq.kaonBar ) pdot_space(i1) = -pdot_space(i1)

             !**** Add-up scalar field contributions:
             pdot_space(i1) = pdot_space(i1) -(DiracMass(Index1,Index2,Index3,mass,id,charge,anti)) &
                                           &  / Part%mom(0) *  dsigma_dr(i1)

             !**** Time derivative contribution
             !**** to the rhs of equation dp^*/dt=pdot

             if (lorentz_flag) then
                pdot_time(i1) = ( &
                     &   SelfEnergy_vector_old(Index1,Index2,Index3,i1,id,charge,anti) &
                     & - SelfEnergy_vector(Index1,Index2,Index3,i1,id,charge,anti)  &
                     &           ) / dt
                ! For antibaryons or antikaons the vector field changes sign:
                if ( anti .or. id.eq.kaonBar ) pdot_time(i1) = -pdot_time(i1)
             else
                pdot_time(i1) = 0.
             end if

           end do

        end if

        rdot = Part%vel

        AddCoulomb = .true.

    else if(eventtype==HeavyIon .and. (id.eq.rho .or. id.eq.omegaMeson .or. id.eq.phi)) then

        call gradients(Part,Grad_P,Grad_R, flagOk)

        if (.not.flagOk.and.deleteTachyons) then
           !           call Traceback('*** gradients, flagOk=F')
           Grad_P = (/ 1., 1., 1. /) ! this is too large --> checkvelo failure
        end if

        rdot = Grad_P
        pdot_space(1:3) = -Grad_R(1:3)

        AddCoulomb = .false.

    end if

    if(AddCoulomb) then
       !**** electromagnetic part
       cpot = emfoca(place,impuls,Part%charge,Part%ID,emForce)
       pdot_space(1:3) = pdot_space(1:3) + emForce(1:3)
    end if

    !**** Total derivative dp^*/dt
    pdot = pdot_space + pdot_time

  end subroutine rhs

  !****************************************************************************
  !****is* propagation_RMF/reportError
  ! NAME
  ! subroutine reportError(Part, errCode)
  ! PURPOSE
  ! generate histograms reporting, at which masses errors occured
  !****************************************************************************
  subroutine reportError(Part, errCode)
    use particleDefinition
    use histMC
    use minkowski, only: abs4Sq
    use Callstack, only: Traceback
    use output, only: DoPR
    use offshellPotential, only: get_offshell_debug

    type(particle),intent(inout)  :: part
    integer, intent(in) :: errCode

    logical, save :: doInit = .true.
    integer, save :: nCalls = 0

    type(histogramMC), save :: hErrMass, hErrMu
    real :: mu

    if (errCode == 0) return ! no error

    if (get_offshell_debug()) then
       nCalls=nCalls+1

       if (doInit) then
          call CreateHistMC(hErrMass, "mass", 0.0, 1.5, 0.01, 7)
          call CreateHistMC(hErrMu, "mu", 0.0, 1.5, 0.01, 7)
          doInit = .false.
       end if

       if (mod(nCalls,10).eq.0) then
          call WriteHistMC(hErrMass,file='hErrMass.dat')
          call WriteHistMC(hErrMu,file='hErrMu.dat')
       end if

       if (Part%ID .eq. 103) then ! only rho mesons
          mu = abs4Sq(part%mom)
          if (mu>0) then
             mu = sqrt(mu)
          else
             mu = -sqrt(-mu)
          end if

          call addHistMC(hErrMass, Part%mass, errCode, 1.)
          call addHistMC(hErrMu, mu, errCode, 1.)
       end if

    end if

    if (deleteTachyons) then
       if (DoPR(2)) write(*,*) '--- this particle is now deleted! (2)'
       part%ID=0
    end if


  end subroutine reportError

end module propagation_RMF
