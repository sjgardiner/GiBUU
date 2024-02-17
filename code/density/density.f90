!******************************************************************************
!****m* /densitymodule
! NAME
! module densitymodule
!
! PURPOSE
! Administrates the calculation of the density.
!******************************************************************************
module densitymodule
  use dichteDefinition
  use constants, only: singlePrecision
  use Callstack, only: TRACEBACK
  use hist

  implicit none

  private

  !****************************************************************************
  !****g* densitymodule/densitySwitch
  ! SOURCE
  !
  integer, save :: densitySwitch=1
  ! PURPOSE
  ! This switch decides whether the density is static or dynamic during
  ! the run. ("Static" makes sense only for fixed target scenarios!)
  !
  ! One can use a static density if the nucleus stays roughly in its
  ! ground state during the collision.
  !
  ! possible values:
  ! * 0: Density is set to 0.
  ! * 1: Dynamic density according to test-particle distribution.
  ! * 2: Static density (not for heavy-ion collisions).
  ! * 3: Resting matter: Density is given by the two input parameters
  !   "densityInput_neutron" and "densityInput_proton".
  ! * 4: Dynamic density in a box. Assumes the same density everywhere, but
  !   also calculates the momentum distribution
  !****************************************************************************


  !****************************************************************************
  !****g* densitymodule/linearInterpolation
  ! SOURCE
  !
  logical, save :: linearInterpolation = .true.
  ! PURPOSE
  ! If this switch is 'true', then the dynamic-density mode uses linear
  ! interpolation to determine the density in between the gridpoints.
  !****************************************************************************


  ! Input parameters for the grid where density is calculated upon:

  !****************************************************************************
  !****g* densitymodule/gridSize
  ! SOURCE
  !
  real, dimension(1:3), save, public :: gridSize = (/12.,12.,12./)
  ! PURPOSE
  ! Size of density grid in fm.
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/gridPoints
  ! SOURCE
  !
  integer, dimension(1:3), save, public :: gridPoints = (/30,30,30/)
  ! PURPOSE
  ! Number of gridpoints in each space direction.
  !****************************************************************************


  !****************************************************************************
  !****g* densitymodule/gridSpacing
  ! SOURCE
  !
  real, dimension(1:3), save, public :: gridSpacing = 1.
  ! PURPOSE
  ! Spacing of the grid in each dimension, determined by gridSize/gridPoints.
  !****************************************************************************


  !****************************************************************************
  !****g* densitymodule/densityInput_proton
  ! SOURCE
  !
  real, save :: densityInput_proton=0.084
  ! PURPOSE
  ! Assumed proton density if densitySwitch=3
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/densityInput_neutron
  ! SOURCE
  !
  real, save :: densityInput_neutron=0.084
  ! PURPOSE
  ! Assumed neutron density if densitySwitch=3
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/setnewsmearing
  ! SOURCE
  !
  logical, save :: setnewsmearing=.false.
  ! PURPOSE
  ! Readjust the smearing to a different width if .true.
  !****************************************************************************

  !****************************************************************************
  !****g* densitymodule/newsmearing
  ! SOURCE
  !
  real, save :: newsmearing=1.
  ! PURPOSE
  ! Use a smearing width as in a grid wtih newsmearing times the gridspacing
  !****************************************************************************


  type(dichte),          save, Allocatable, dimension(:,:,:), public :: densField      ! Field where density is saved in
  real(singlePrecision), save, Allocatable, dimension(:,:,:), public :: totalDens      ! density of baryons plus density of antibaryons (fm^-3)

  ! Fields needed when the RMF is used:
  real,                  save, Allocatable, dimension(:,:,:)   :: sigmaField      ! sigma-meson field (GeV)
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: omegaField      ! omega-meson field (GeV)
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: rhoField        ! rho-meson field (GeV)

  ! stored values from previous time step:
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: baryonCurrent_old  ! baryon density (fm^-3)
  real(singlePrecision), save, Allocatable, dimension(:,:,:)   :: sigmaField_old     ! sigma-meson field
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: omegaField_old     ! omega-meson field
  real(singlePrecision), save, Allocatable, dimension(:,:,:,:) :: rhoField_old       ! rho-meson field


  real(singlePrecision), save, Allocatable, dimension(:,:,:) :: scalarDens   ! scalar density (fm^-3)
  real(singlePrecision), save, Allocatable, dimension(:,:,:) :: d_scalarDens_dsigma    ! partial derivative of the scalar density over sigma field (fm^-3/GeV)
  real,                  save, Allocatable, dimension(:,:,:) :: K_rmf, Diag_rmf, U_rmf    ! auxiliary arrays needed when RMF gradients are included

  real(singlePrecision), save, Allocatable, dimension(:,:,:,:), public :: mom4dens     ! four-momentum density field (not used in propagation)

  real(singlePrecision), save, Allocatable, dimension(:,:,:,:,:), public :: Tens      ! energy-momentum tensor (not used in propagation)


  !****************************************************************************
  !****g* densitymodule/nLargePoints
  ! SOURCE
  !
  integer, save, public :: nLargePoints=2
  ! PURPOSE
  ! Number of points which are considered to the left and right to smear
  ! density on
  !****************************************************************************

  integer, parameter,public :: nSmallPoints=2 !number of Gridpoints in small grid

  !normalized weights for smearing
  real, save,dimension(:,:),allocatable,public :: smearingWeights


  !****************************************************************************

  real, dimension(0:1), save :: pF_Box = -99.9
  type(histogram), dimension(0:1), save :: hP_Box
  type(dichte), save :: dichte_Box

  logical, save :: initFlag=.true.  ! Checks whether input was already read in


  !****************************************************************************
  !****f* densitymodule/DiracMass
  ! NAME
  ! real function DiracMass(Part)
  ! real function DiracMass(i1,i2,i3,barMass,id,charge,antiFlag)
  ! PURPOSE
  ! calculates the effective (Dirac) mass of a particle
  ! INPUTS
  ! * type(particle) :: Part
  !
  ! or:
  ! * integer :: i1,i2,i3 -- index in spatial grid
  ! * real :: barMass     -- bare mass  of the particle
  ! * logical :: antiFlag -- antiparticle
  !
  ! OUTPUT
  ! returns m^* = m + S
  ! NOTES
  ! returns m^* = m, if particle outside of grid or wrong particle id
  !****************************************************************************
  interface DiracMass
     module procedure DiracMass1, DiracMass2
  end interface DiracMass

  !****************************************************************************
  !****f* densitymodule/getGridIndex
  ! NAME
  ! logical function getGridIndex(pos,index)
  ! logical function getGridIndex(pos,index,add)
  ! logical function getGridIndex(pos,index,iSmall,small)
  ! PURPOSE
  ! Returns .false. if the point is outside the grid (then 'index' is invalid).
  ! 'index' is the index of the corresponding coordinate. It is even calculated
  ! if the vale '.false.' is returned. Then the validity of the results has
  ! to be checked elsewhere.
  !
  ! One may use the optional parameter 'add' to shift the test according
  ! being inside by a number of bins.
  !
  ! The third variant also calculates the index for a 'small' grid
  !
  ! INPUTS
  ! * real, dimension(1:3) :: pos
  ! * integer, optional :: add
  ! OUTPUT
  ! * integer, dimension(1:3) :: index
  ! * integer, dimension(1:3) :: iSmall
  ! * integer, optional :: small
  ! * return value signaling "inside the grid"
  ! NOTES
  ! * The first tests according the real grid size is necessary, since it
  !   may give integer overfloats, if one only checks the integer values.
  ! * Usually, we do not need the output of iSmall directly, but only small-
  !   The interface policy does not allow a distinction between 'add' and
  !   'small', so we add the redundant 'iSmall' to the function interface.
  !****************************************************************************
  interface getGridIndex
     module procedure getGridIndex1, getGridIndex2, getGridIndex3
  end interface getGridIndex

  !****************************************************************************
  !****s* densitymodule/getFieldRMF
  ! NAME
  ! subroutine getFieldRMF(pos, sigma,omega,rho)
  ! subroutine getFieldRMF(ix,iy,iz, sigma,omega,rho)
  !
  ! PURPOSE
  ! return the value of some fileds at the given position
  ! INPUTS
  ! * real, dimension(1:3) :: pos -- position to be considered
  ! * integer :: ix,iy,iz -- index of position
  ! OUTPUT
  ! * real, OPTIONAL :: sigma :: value of sigma field
  ! * real, dimension(0:3), OPTIONAL :: omega :: value of omega field
  ! * real, dimension(0:3), OPTIONAL :: rho :: value of rho field
  !
  ! NOTES
  ! * if index is given as coordinates, no check whether inside/outside is
  !   performed
  !****************************************************************************
  interface getFieldRMF
     module procedure getFieldRMF_i, getFieldRMF_p
  end interface getFieldRMF

  !****************************************************************************
  !****s* densitymodule/energyDeterminationRMF
  ! NAME
  ! subroutine energyDeterminationRMF(Part)
  ! subroutine energyDeterminationRMF(Parts)
  ! PURPOSE
  ! This subroutine determines the one-particle energy E^*,
  ! which is the zeroth component of the kinetic four-momentum
  ! in the frame, where the space components of the kinetic four-momentum
  ! are given.
  ! INPUTS
  ! * type(particle),intent(inOut) :: Part  -- Particle whose energy should
  !   be calculated.
  !
  ! or:
  ! * type((particle),dimension(:,:),intent(inOut) :: Parts
  !   -- Particles whose energy should be calculated.
  ! NOTES
  ! Should be used in RMF-mode. Please note, that not the full
  ! single-particle energy is computed here. The vector field contribution is
  ! missed in E^*.
  !****************************************************************************
  interface energyDeterminationRMF
     module procedure energyDeterminationRMF0, energyDeterminationRMF2
  end interface energyDeterminationRMF

  !****************************************************************************
  !****f* densitymodule/SelfEnergy_scalar
  ! NAME
  ! real function SelfEnergy_scalar(Part)
  ! real function SelfEnergy_scalar(i1,i2,i3,id,antiFlag)
  !
  ! PURPOSE
  ! calculates the scalar component of the RMF selfenergy of a particle
  !
  ! INPUTS
  ! * type(particle) :: Part
  !
  ! or:
  ! * integer :: i1,i2,i3 -- index in the grid
  ! * integer :: id -- id of the particle
  ! * logical :: antiFlag -- antiparticle or not
  !
  ! NOTES
  ! * for kaons selfenergy_scalar includes the Sigma_KN term only!
  ! * gs_kaon = Sigma_KN / f_pi^2
  !****************************************************************************
  interface SelfEnergy_scalar
     module procedure SelfEnergy_scalar1, SelfEnergy_scalar2
  end interface SelfEnergy_scalar


!!$  public :: hadronicPotEnergy
  public :: densityAt
  public :: getBaryonDensity
  public :: updateDensity
  public :: updateRMF, TensorRMF, TotalEnergyRMF
  public :: storeFieldsRMF
  public :: energyDeterminationRMF
  public :: Particle4MomentumRMF, true4MomentumRMF
  public :: SelfEnergy_scalar, SelfEnergy_vector, SelfEnergy_vector_old
  public :: DiracMass
  public :: boostToLRF
  public :: acceptGrid
  public :: getGridSpacing, getGridPoints, getGridSpacing0
  public :: fermiMomentum_sym, fermiMomentum_noIsospin
  public :: getDensitySwitch, setDensitySwitch, setDensityInput
  public :: cleanup
  public :: FermiMomAt
  public :: getGridIndex
  public :: writeDensityPlane
  public :: set_pF_Box
  public :: distMomBox
  public :: getFieldRMF

contains

!!$  !*************************************************************************
!!$  !****f* densitymodule/hadronicPotEnergy
!!$  ! NAME
!!$  ! real function hadronicPotEnergy(alpha,beta,tau)
!!$  ! PURPOSE
!!$  ! Evaluate the energy wich is stored in the hadronic potential for a
!!$  ! potential which is not momentum dependend.
!!$  ! INPUTS
!!$  ! * real alpha,beta,tau --  Potential parameters
!!$  ! RESULT
!!$  ! * Potential Energy in GeV
!!$  !*************************************************************************
!!$  function hadronicPotEnergy(alpha,beta,tau)
!!$    use constants, only : rhoNull
!!$
!!$    real :: hadronicPotEnergy
!!$    real, intent(in) :: alpha,beta,tau ! parameters of hadronic potential
!!$
!!$    ! Integrate the energy density
!!$    integer :: index_X,index_Y,index_Z
!!$    Do index_X=-gridPoints(1),gridPoints(1)
!!$       Do index_Y=-gridPoints(2),gridPoints(2)
!!$          Do index_Z=-gridPoints(3),gridPoints(3)
!!$             hadronicPotEnergy=hadronicPotEnergy &
!!$                  & +(alpha/2.*densField(index_X,index_Y,index_Z)%baryon(0)**2/rhoNull &
!!$                  & +beta/(tau+1.)*(densField(index_X,index_Y,index_Z)%baryon(0)**(tau+1.)) &
!!$                  & /(rhoNull**tau))*gridSpacing(1)*gridSpacing(2)*gridSpacing(3)
!!$          End do
!!$       End do
!!$    End do
!!$    hadronicPotEnergy=hadronicPotEnergy/1000. !in units of GeV
!!$  end function hadronicPotEnergy

  !****************************************************************************
  !****s* densitymodule/init
  ! NAME
  ! subroutine init
  !
  ! PURPOSE
  ! read the namelist, initalize (parts of) the module
  !****************************************************************************
  subroutine init
    use output, only: Write_InitStatus,Write_ReadingInput
    use RMF, only: getRMF_flag
    use nucleusDefinition
    use nucleus, only: getTarget
    use inputGeneral, only: eventType, numEnsembles
    use eventtypes, only: elementary, Box

    integer :: ios
    type(tNucleus), pointer :: targetNuc

    !**************************************************************************
    !****n* densitymodule/initDensity
    ! NAME
    ! NAMELIST /initDensity/
    ! PURPOSE
    ! Includes the input switches and variables:
    ! * densitySwitch
    ! * linearInterpolation
    ! * densityInput_proton
    ! * densityInput_neutron
    ! * gridSize
    ! * gridPoints
    ! * setnewsmearing
    ! * newsmearing
    ! * nLargePoints
    !**************************************************************************
    NAMELIST /initDensity/ densitySwitch, linearInterpolation, &
         densityInput_proton, densityInput_neutron, &
         gridSize, gridPoints, setnewsmearing, newsmearing, nLargePoints

    call Write_InitStatus('densitymodule',0)

    call Write_ReadingInput('initDensity',0)
    rewind(5)
    read(5,nml=initDensity,iostat=ios)
    call Write_ReadingInput('initDensity',0,ios)

    select case (eventType)
    case (elementary)
       densitySwitch = 0
       write(*,*) 'densitySwitch is set to 0 for elementary target'
    case (Box)
       if (densitySwitch .ne. 4) then
          densitySwitch = 0
          write(*,*) 'densitySwitch is set to 0 for Box'
       end if
    case default
       targetNuc => getTarget()
       if (targetNuc%mass==1) then
          densitySwitch = 0
          write(*,*) 'densitySwitch is set to 0 for nucleon target'
       end if
    end select

    write(*,'(A,I5,A)') ' Set densitySwitch to ',densitySwitch,'.'

    if (densitySwitch==3) then
       write(*,'(A)')        '   => We assume resting matter.'
       write(*,'(A,F5.3,A)') '      proton density:  ',densityInput_proton, '/fm^3.'
       write(*,'(A,F5.3,A)') '      neutron density: ',densityInput_neutron,'/fm^3.'
    end if
    write(*,'(A,L6,A)') ' Set linearInterpolation to  ',linearInterpolation,'.'
    if (setnewsmearing) then
       write(*,*) 'Smearingwidth set to ', newsmearing
    end if

    select case (densitySwitch)
    case (1,4)
       if (numEnsembles < 100) then
          write(*,*)
          write(*,*) 'Dynamic density needs smooth densities!'
          write(*,*) 'Increase numEnsembles in the jobcard to 100 or more.'
          call TraceBack()
       end if
    end select

    call Write_ReadingInput('initDensity',1)

    ! Set grid spacing and allocate vectors
    gridSpacing=gridSize/float(gridPoints)
    call alloc_fields

    initFlag=.false.

    call Write_InitStatus('densitymodule',1)
  end subroutine init

  !****************************************************************************
  !****************************************************************************
  subroutine alloc_Fields
    use RMF, only: getRMF_flag, grad_flag, Tens_flag

    write(*,'(A,2(F6.2,","),F6.2,A)') &
         ' The gridsize of the density grid is        = (', gridsize, ') fm'
    write(*,'(A,2(I6,","),I6,A)') &
         ' The number of gridpoints per dimension are = (', gridPoints, ') '
    write(*,'(A,2(F6.2,","),F6.2,A)') &
         ' The grid spacing is                        = (', gridSpacing, ') fm'

    allocate(densField(-gridPoints(1):gridPoints(1),&
         &             -gridPoints(2):gridPoints(2),&
         &             -gridPoints(3):gridPoints(3)))
    allocate(totalDens(-gridPoints(1):gridPoints(1),&
         &             -gridPoints(2):gridPoints(2),&
         &             -gridPoints(3):gridPoints(3)))

    if ( getRMF_flag() ) then
       allocate(sigmaField(-gridPoints(1):gridPoints(1),&
            &              -gridPoints(2):gridPoints(2),&
            &              -gridPoints(3):gridPoints(3)))
       allocate(omegaField(-gridPoints(1):gridPoints(1),&
            &              -gridPoints(2):gridPoints(2),&
            &              -gridPoints(3):gridPoints(3), 0:3))
       allocate(rhoField(-gridPoints(1):gridPoints(1),&
            &            -gridPoints(2):gridPoints(2),&
            &            -gridPoints(3):gridPoints(3), 0:3))
       allocate(scalarDens(-gridPoints(1):gridPoints(1),&
            &              -gridPoints(2):gridPoints(2),&
            &              -gridPoints(3):gridPoints(3)))
       allocate(d_scalarDens_dsigma(-gridPoints(1):gridPoints(1),&
            &                       -gridPoints(2):gridPoints(2),&
            &                       -gridPoints(3):gridPoints(3)))
       allocate(omegaField_old(-gridPoints(1):gridPoints(1),&
            &                  -gridPoints(2):gridPoints(2),&
            &                  -gridPoints(3):gridPoints(3), 0:3))
       allocate(rhoField_old(-gridPoints(1):gridPoints(1),&
            &                -gridPoints(2):gridPoints(2),&
            &                -gridPoints(3):gridPoints(3), 0:3))
       allocate(baryonCurrent_old(-gridPoints(1):gridPoints(1),&
            &                     -gridPoints(2):gridPoints(2),&
            &                     -gridPoints(3):gridPoints(3), 1:3))

       if (grad_flag) then
          allocate(K_rmf(-gridPoints(1):gridPoints(1),&
               &         -gridPoints(2):gridPoints(2),&
               &         -gridPoints(3):gridPoints(3)))
          allocate(Diag_rmf(-gridPoints(1):gridPoints(1),&
               &            -gridPoints(2):gridPoints(2),&
               &            -gridPoints(3):gridPoints(3)))
          allocate(U_rmf(-gridPoints(1):gridPoints(1),&
               &         -gridPoints(2):gridPoints(2),&
               &         -gridPoints(3):gridPoints(3)))
       end if

       if (Tens_flag) then
          allocate(Tens(-gridPoints(1):gridPoints(1),&
               &        -gridPoints(2):gridPoints(2),&
               &        -gridPoints(3):gridPoints(3), 0:3, 0:3))
          allocate(mom4dens(-gridPoints(1):gridPoints(1),&
               &            -gridPoints(2):gridPoints(2),&
               &            -gridPoints(3):gridPoints(3), 0:3))
          allocate(sigmaField_old(-gridPoints(1):gridPoints(1),&
               &                  -gridPoints(2):gridPoints(2),&
               &                  -gridPoints(3):gridPoints(3)))
       end if

    end if

  end subroutine alloc_Fields

  !****************************************************************************
  !****************************************************************************

  subroutine cleanUp
    call cleanup_fields
    if (allocated(SmearingWeights)) deallocate(SmearingWeights)
  end subroutine cleanUp

  !****************************************************************************

  subroutine cleanup_fields
    use RMF, only: getRMF_flag, grad_flag, Tens_flag

    if (allocated(densField)) deallocate(densField, totalDens)

    if ( getRMF_flag() ) then
       deallocate(sigmaField,omegaField,rhoField,scalarDens, &
            d_scalarDens_dsigma,omegaField_old, &
            & rhoField_old, baryonCurrent_old)
       if (grad_flag) deallocate(K_rmf,Diag_rmf,U_rmf)
       if (Tens_flag) deallocate(Tens,mom4dens,sigmaField_old)
    end if
  end subroutine cleanup_fields

  !****************************************************************************
  !****s* densitymodule/acceptGrid
  ! NAME
  ! subroutine acceptGrid(GridSpacing_x,GridSpacing_y,GridSpacing_z)
  !
  ! PURPOSE
  ! Accepts the grid spacings and recalculates
  ! the numbers of points in each direction. Needed for relativistic HIC,
  ! when initial nuclei are Lorentz-contracted. Called from the module
  ! initHeavyIon.
  !****************************************************************************
  subroutine acceptGrid(GridSpacing_x,GridSpacing_y,GridSpacing_z)
    use RMF, only: getRMF_flag
    use nucleusDefinition
    use nucleus, only: getTarget, getProjectile
    use inputGeneral, only: time_max

    type(tNucleus), pointer :: targetNuc, projectileNuc
    real, dimension(1:3) :: PosT, PosP, checkGrid

    real, intent(in) :: GridSpacing_x,GridSpacing_y,GridSpacing_z

    if (initFlag) call init

    if ( GridSpacing_x <= 0.) then
       write(*,*) 'Wrong grid spacing along x-axis:', GridSpacing_x
       call traceback()
    end if
    if ( GridSpacing_y <= 0.) then
       write(*,*) 'Wrong grid spacing along y-axis:', GridSpacing_y
       call traceback()
    end if
    if ( GridSpacing_z <= 0.) then
       write(*,*) 'Wrong grid spacing along z-axis:', GridSpacing_z
       call traceback()
    end if

    GridSpacing(1:3) = (/ GridSpacing_x, GridSpacing_y, GridSpacing_z /)

    ! Adjust gridSize for HIC;
    ! x- and y-directions: nuclei-radius + surface + 5fm
    ! z-direction in addition: (HIC)
    !       grid extension up to Pos_z = V_beam*time_max
    targetNuc => getTarget()
    projectileNuc => getProjectile()
    PosT(1:3) = targetNuc%pos(1:3) &
         + targetNuc%radius + targetNuc%surface + 5.
    PosP(1:3) = projectileNuc%pos(1:3) &
         + projectileNuc%radius + projectileNuc%surface + 5.

    PosT(3) = PosT(3) + targetNuc%vel(3)*time_max
    PosP(3) = PosP(3) + projectileNuc%vel(3)*time_max

    checkGrid(1:3) = aint( max(PosT(1:3),PosP(1:3)) )

    ! check that nuclei will stay inside the spatial grid
    if ( gridSize(1) < checkGrid(1) &
         .or. gridSize(2) < checkGrid(2) &
         .or. gridSize(3) < checkGrid(3) ) then

       write(*,*) ' nuclei partially outside of spatial grid!'
       write(*,'(A,2(F6.2,","),F6.2,A)') &
            &' checkGrid        = (', checkGrid, ') fm'
       write(*,*) ' readjust gridSize too'

       gridSize(1:3) = checkGrid(1:3)

    end if

    call cleanup_fields
    write(*,*) 'New grid parameters:'
    gridPoints = nint(gridSize/gridSpacing)
    call alloc_fields

  end subroutine acceptGrid

  !****************************************************************************
  !****f* densitymodule/getBaryonDensity
  ! NAME
  ! real function getBaryonDensity(pos)
  ! PURPOSE
  ! Evaluate the baryon density at a certain position.
  !****************************************************************************
  real function getBaryonDensity(pos)
    use idTable, only: nucleon
    use RMF, only: getRMF_flag, ModificationFactor
    real, dimension(1:3),intent(in) :: pos

    integer, dimension (1:3) :: i
    real :: factor

    if (getGridIndex(pos,i)) then

       if (getRMF_flag()) then
          ! Modification factor of the antinucleon vector coupling constants
          factor = ModificationFactor(nucleon,.true.,.true.)
       else
          ! Remember: w/o RMF totalDens is (baryon-antibaryon) density
          factor = 1.
       end if

       ! density of the baryons (excluding antibaryons)
       getBaryonDensity = (densField(i(1),i(2),i(3))%baryon(0) &
            + factor*totalDens(i(1),i(2),i(3))) / (1.+factor)

    else

       getBaryonDensity = 0.

    end if

  end function getBaryonDensity


  !****************************************************************************
  !****f* densitymodule/densityAt
  ! NAME
  ! type(dichte) function densityAt(r)
  ! PURPOSE
  ! Evaluate density at some space-point r.
  ! INPUTS
  ! * real, dimension(1:3) :: r --- position where density should be calculated
  !****************************************************************************
  type(dichte) function densityAt(r)
    use densityStatic, only: staticDensity
    use nucleus, only: getTarget
    use inputGeneral, only: continousBoundaries

    real, dimension(1:3), intent(in) :: r

    real,dimension(1:3) :: rNew
    integer :: ind(1:3)

    if (initFlag) call init

    select case (densitySwitch)
    case (0) ! Assume no density

       densityAt = dichte0

    case (1) ! dynamic density according to test-particle distribution

       if (getGridIndex(r,ind)) then   ! inside grid
          if (linearInterpolation) then
             densityAt = interpolate(r)
          else ! no interpolation
             densityAt = densField(ind(1),ind(2),ind(3))
          end if
       else                                  ! outside grid

          densityAt = dichte0

          if (continousBoundaries) then
             rNew=r
             where (r >  gridSize) rNew = r - 2*gridsize
             where (r < -gridSize) rNew = r + 2*gridsize

             if (getGridIndex(rNew,ind)) then  ! inside grid
                if (linearInterpolation) then
                   densityAt = interpolate(rNew)
                else                                 ! no interpolation
                   densityAt = densField(ind(1),ind(2),ind(3))
                end if
             end if
          end if
       end if

    case (2) ! static density of a resting target

       densityAt=staticDensity(r,getTarget())

    case (3) ! Use input density

       densityAt%baryon = (/densityInput_proton+densityInput_neutron,0.,0.,0./)
       densityAt%neutron= (/densityInput_neutron,0.,0.,0./)
       densityAt%proton = (/densityInput_proton,0.,0.,0./)

    case (4) ! Dynamic density in a Box

       densityAt = dichte_Box

    case default
       write(*,*) 'This density switch is not valid:',densitySwitch
       call TRACEBACK()

    end select



  contains

    !**************************************************************************
    ! For linear interpolation between the grid points.
    !**************************************************************************
    type(dichte) function interpolate(r)
      use inputGeneral, only: continousBoundaries

      real, dimension(1:3), intent(in) :: r  ! position where density should be calculated
      integer :: i,j,k,lowIndex(1:3)
      integer,dimension(0:7,1:3) :: grid     ! field with all corners of 3D-box where r is situated in
      real,dimension(1:3) :: gridPos         ! position of gridpoints
      real :: factor

      ! The point "r" where the density should be calculated
      ! is sitting in a 3D-box with lowest corner "LowIndex" on the
      ! density grid. We first construct this point. Then all other corners of
      ! the 3D-box are constructed.
      ! In the end a simple linear interpolation is used to make the density
      ! smooth inside the box.

      ! (1.) Construct lowest lying point (most negative point!)
      LowIndex=nint(r/gridSpacing-0.5)

      ! (2.) Define GridPoints on the corners of the 3D box in which the
      !      point "r" is situated
      do i=0,1
         do j=0,1
            do k=0,1
               grid(i+2*j+4*k,1:3)=LowIndex+(/i,j,k/)
            end do
         end do
      end do

      ! (3.) Do linear interpolation
      interpolate = dichte0

      do i=0,7
         gridPos(1:3)=grid(i,1:3)*gridSpacing(1:3) !position of grid point
         factor=1.  ! weight for linear interpolation for each grid point
         do j=1,3
            factor=factor*Abs(gridSpacing(j)-abs(r(j)-gridPos(j)))/gridspacing(j)
         end do
         if (continousBoundaries) then
            do j=1,3
               if (grid(i,j).lt.-gridPoints(j)) then
                  grid(i,j)=grid(i,j)+2*gridPoints(j)
               else if (grid(i,j).gt.gridPoints(j)) then
                  grid(i,j)=grid(i,j)-2*gridPoints(j)
               end if
            end do
         end if
         if ((abs(grid(i,1)).le.gridPoints(1)) &
              .and.(abs(grid(i,2)).le.gridPoints(2))&
              .and.(abs(grid(i,3)).le.gridPoints(3))) then
            interpolate=interpolate &
                 + densField(grid(i,1),grid(i,2),grid(i,3))*factor
         end if
      end do
    end function interpolate

  end function densityAt

  !****************************************************************************
  !****f* densitymodule/updateDensity
  ! NAME
  ! subroutine updateDensity(parts)
  ! PURPOSE
  ! Updates the vector densField which is used by densityAt and stores the
  ! density of the testparticles.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: parts
  !
  ! NOTE
  ! This routine does not gets faster when using OpenMP.
  !****************************************************************************
  subroutine updateDensity(Parts)

    use particledefinition
    use inputGeneral, only: continousboundaries
    use output, only: DoPr, Write_InitStatus
    use idtable, only: isBaryon

    type(particle),dimension(:,:),intent(in),target :: Parts

    !Flag that controls whether density weights are already calculated
    logical, save :: FirstTime=.true.

    integer :: bEdge ! binary edge code
    type(particle),pointer :: pPart
    integer :: iEns,iPart, nEns
    real :: SW
    integer :: i1,i2,i3
    integer, dimension(1:3) :: iPos
    integer, dimension(1:3) :: iSmall, indexC
    integer :: small, large
    logical :: flag

    if (initFlag) call init

    nEns = size(Parts(:,1))

    if (densitySwitch.eq.4) then
       call updateDensityBox
       return
    end if

    if (FirstTime) then
       call Write_InitStatus('updateDensity',0)
       call initDensityWeights
       call FillStaticDensity
       FirstTime=.false.
       call Write_InitStatus('updateDensity',1)
    end if

    if (densitySwitch.ne.1) return

    if (DoPr(2)) write(*,*) 'Updating density field'

    densField = dichte0
    totalDens = 0.

    do iEns=1,size(Parts,dim=1)
       do iPart=1,size(Parts,dim=2)

          pPart => Parts(iEns,iPart)

          if (pPart%id == 0) cycle
          if (pPart%Id <  0) exit

          ! position in large and small grid:
          flag = getGridIndex(pPart%pos, iPos, iSmall,small)
          if (.not.flag) then
             ! we are ignoring the error flag, since the validity
             ! if the returned values is checked below separately
             ! (and differently) (cf. continousBoundaries)

!!$             write(*,*) iEns,iPart
!!$             write(*,*) pPart%pos
!!$             write(*,*) iPos
!!$             write(*,*) iSmall
!!$             write(*,*) small
!!$             write(*,*) nint(pPart%pos/gridspacing)
!!$             call traceback("outside of grid")
          end if

          large=0

          ! Smearing particle over points in neighborhood:
          do i1=iPos(1)-nLargePoints,iPos(1)+nLargePoints
             do i2=iPos(2)-nLargePoints,iPos(2)+nLargePoints
                do i3=iPos(3)-nLargePoints,iPos(3)+nLargePoints
                   large=large+1
                   indexC = (/i1,i2,i3/)

                   if (continousBoundaries) then
                      ! Implement continous boundary conditions,
                      ! so every point outside the grid is attributed to a
                      ! point inside.

                      where (indexC >  gridPoints) &
                           indexC = indexC - 2*gridPoints
                      where (indexC < -gridPoints) &
                           indexC = indexC + 2*gridPoints

                   end if

                   SW = smearingWeights(small,large)

                   if ( (all(abs(indexC) <= gridPoints)) &   ! = inside grid
                        .and.(SW.gt.0.0) ) then

                      bEdge = 0

                      if (continousBoundaries) then
                         ! When having continous boundaries, then every
                         ! contribution to a point at the edge of the grid
                         ! must also be attributed to its opposing point
                         ! on the other side.
                         if (abs(indexC(1)) == gridPoints(1)) bEdge = bEdge+1
                         if (abs(indexC(2)) == gridPoints(2)) bEdge = bEdge+2
                         if (abs(indexC(3)) == gridPoints(3)) bEdge = bEdge+4
                      end if

                      call addTodensityfield( indexC )

                      select case (bEdge)
                      case (1)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                      case (2)
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                      case (3)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                         call addTodensityfield( indexC * (/ -1,-1, 1 /) )
                      case (4)
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                      case (5)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                         call addTodensityfield( indexC * (/ -1, 1,-1 /) )
                      case (6)
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                         call addTodensityfield( indexC * (/  1,-1,-1 /) )
                      case (7)
                         call addTodensityfield( indexC * (/ -1, 1, 1 /) )
                         call addTodensityfield( indexC * (/  1,-1, 1 /) )
                         call addTodensityfield( indexC * (/ -1,-1, 1 /) )
                         call addTodensityfield( indexC * (/  1, 1,-1 /) )
                         call addTodensityfield( indexC * (/ -1, 1,-1 /) )
                         call addTodensityfield( indexC * (/  1,-1,-1 /) )
                         call addTodensityfield( indexC * (/ -1,-1,-1 /) )
                      end select

                   end if

                end do ! loop I1
             end do ! loop I2
          end do ! loop I3

       end do
    end do

  contains

    !**************************************************************************
    subroutine addTodensityfield(ind)

      use idtable, only: nucleon, isBaryon
      use RMF, only: getRMF_flag, ModificationFactor

      integer, dimension(3), intent(in) :: ind
      real :: factor
      real,dimension(0:3) :: vec

      vec(0) = 1.
      vec(1:3) = pPart%vel(1:3)
      vec = vec*SW

      if (isBaryon(pPart%id)) then

         if (getRMF_flag()) then

            ! Modification factor for vector coupling constants (always positive !!!):
            factor = ModificationFactor(pPart%Id,pPart%anti,.true.)
            if ( pPart%anti ) factor = -factor

            totalDens(ind(1),ind(2),ind(3))=totalDens(ind(1),ind(2),ind(3))+SW

         else

            factor = 1.

            if ( .not.pPart%anti ) then
               totalDens(ind(1),ind(2),ind(3))=totalDens(ind(1),ind(2),ind(3))+SW
            else
               totalDens(ind(1),ind(2),ind(3))=totalDens(ind(1),ind(2),ind(3))-SW
            end if

         end if

         densField(ind(1),ind(2),ind(3))%baryon = densField(ind(1),ind(2),ind(3))%baryon &
              + vec*factor

         if ( pPart%ID.eq.nucleon ) then
            if ( pPart%charge.eq.0 ) then
               densField(ind(1),ind(2),ind(3))%neutron = densField(ind(1),ind(2),ind(3))%neutron &
                    + vec*factor
            else
               densField(ind(1),ind(2),ind(3))%proton = densField(ind(1),ind(2),ind(3))%proton &
                    & + vec*factor
            end if
         end if

      end if

      densField(ind(1),ind(2),ind(3))%charge = densField(ind(1),ind(2),ind(3))%charge &
           & + vec * pPart%charge

    end subroutine addTodensityfield


    !**************************************************************************
    !****f* updateDensity/initDensityWeights
    ! NAME
    ! subroutine initDensityWeights
    ! PURPOSE
    ! * Initializes weights which are used to evaluate the densities at the
    !   gridpoints.
    ! * Each particle has some coordinate off the gridpoints. Therefore it is
    !   needed to define weights which tell how much a particle at position r
    !   contributes to the ith grid point.
    !   The particle are smeared with a gaussian distribution.
    ! * The width is chosen such that it is equal to the maximum of the
    !   gridspacings.
    ! INPUTS
    ! * type(particle) Parts(:,:)
    ! USES
    ! * idTable, particleDefinition
    !**************************************************************************
    subroutine initDensityWeights
      use output

      !grid spacing in small grid:
      real, dimension(1:3) :: smallGridSpacing

      ! Gaussian width, used to smear the density:
      real, save, dimension(1:3) :: smearingWidth=0.5
      ! (cutoff for smearing)**2:
      real, save, dimension(1:3) :: smearingCutOff=5.

      real,dimension(1:3) :: rLarge        ! Point on large grid
      real,dimension(1:3) :: rSmall        ! Point on small grid
      integer :: iS1,iS2,iS3
      integer :: iL1,iL2,iL3
      integer :: small,large
      ! (distance of points on large and small grids)**2:
      real,dimension(1:3) :: rSquare
      real :: norm             ! Normalization
      integer :: i
      logical,save :: firstTime=.true.

      call Write_InitStatus("Smearing weights for density",0)

      if (.not.firstTime) &
           call TRACEBACK("WARNING: Allocating array again !!!",-1)

      allocate(smearingWeights(1:(2*nSmallPoints+1)**3,&
           &                   1:(2*nLargePoints+1)**3))

      ! Set smearing width in each direction equal to the grid spacing in
      ! this direction:

      if (setnewsmearing) then
         smearingWidth=gridSpacing*newsmearing
         smallGridSpacing=GridSpacing*newsmearing/(2.*float(nSmallPoints)+1.)
      else
         smearingWidth=gridSpacing
         smallGridSpacing=GridSpacing/(2.*float(nSmallPoints)+1.)
      end if


      ! Check that cut off is small than the maximal distance of a point on
      ! which I am smearing:
      do i = 1,3
         do
            if ( (float(nLargePoints)+0.5)*gridSpacing(i).lt.sqrt(smearingCutOff(i)) ) then
!!$            If(firstTime) then
!!$               write (*,'("     ************************************************************")')
!!$               write (*,'("     Error! Smearing cut off too large!!! Enlarge nLargePoints or decrease cut off!")')
!!$               write (*,'("     Setting smearing cut off to better value... ")')
!!$               write (*,'("     ************************************************************")')
!!$               firsttime=.false.
!!$            end if
               smearingCutOff(i)=smearingCutOff(i)*0.9
            else
               exit
            end if
         end do
      end do

      write(*,'(" Gridspacing in XYZ: ",3F8.5)') gridspacing(1:3)
      write(*,'(" Cut offs for smearing in XYZ: ",3F8.5," fm")') &
           sqrt(smearingCutOff)
      write(*,'(" Smearing with the function e^(-x^2/2/",F4.2,"-y^2/2/",F4.2,"-z^2/2/",F4.2,")")') &
           smearingWidth**2
      write(*,'(" Minimal distance of not-considered gridpoint: ",3F8.5)') &
           (float(nLargePoints)+0.5)*gridSpacing
      write(*,*) "points in real  grid: ",nLargePoints
      write(*,*) "points in small grid: ",nSmallPoints


      smearingWeights = 0.0

      !Loop over points on small grid
      do iS1=-nSmallPoints,nSmallPoints
         do iS2=-nSmallPoints,nSmallPoints
            do iS3=-nSmallPoints,nSmallPoints
               small=1+(nSmallPoints+iS3)&
                    & +(nSmallPoints+iS2)*(2*nSmallPoints+1) &
                    & +(nSmallPoints+iS1)*(2*nSmallPoints+1)**2
               rSmall=(/iS1,iS2,iS3/)*smallGridSpacing

               Norm=0.
               large=0

               ! Loop over points in the real grid:
               do iL1=-nLargePoints,nLargePoints
                  do iL2=-nLargePoints,nLargePoints
                     do iL3=-nLargePoints,nLargePoints
                        large=large+1
                        rLarge=(/iL1,iL2,iL3/)*GridSpacing

                        ! (Distance between points on small and real grid)**2:
                        rSquare(:) = (rLarge(:)-rSmall(:))**2

                        if (all(rSquare < SmearingCutoff)) then
                           smearingWeights(small,large)= &
                                exp(-sum(rSquare(:)/(2.*smearingWidth(:)**2)))
                        end if
                        Norm=Norm+smearingWeights(small,large)    !Normalization
                     end do
                  end do
               end do

               !Normalizing the weigths:
               smearingWeights(small,:) = smearingWeights(small,:)/Norm

            end do
         end do
      end do

      ! overall normalization factor:
      smearingWeights = smearingWeights &
           /(gridSpacing(1)*gridSpacing(2)*gridSpacing(3)) &
           /size(Parts(:,1))

      call Write_InitStatus("Smearing weights for density",1)
      firstTime = .false.

    end subroutine initDensityWeights

    !**************************************************************************
    !****f* updateDensity/FillStaticDensity
    ! NAME
    ! subroutine FillStaticDensity
    ! PURPOSE
    ! In the case of static density, this fills the fields densField and
    ! totalDens with values of the density parametrization. Needed in
    ! the case of RMF calculations and also for Coulomb potential.
    !**************************************************************************
    subroutine FillStaticDensity
      use densityStatic
      use nucleus, only: getTarget
      use output, only: Write_InitStatus

      real,dimension(1:3)    :: r  !position where density should be calculated

      integer :: I1,I2,I3

      if (densitySwitch.ne.2) return
      call Write_InitStatus("FillStaticDensity",0)
      do I1 = -gridPoints(1),gridPoints(1)
         r(1) = I1*gridSpacing(1)
         do I2 = -gridPoints(2),gridPoints(2)
            r(2) = I2*gridSpacing(2)
            do I3 = -gridPoints(3),gridPoints(3)
               r(3) = I3*gridSpacing(3)
               densField(I1,I2,I3) = staticDensity(r,getTarget())
               totalDens(I1,I2,I3) = densField(I1,I2,I3)%Baryon(0)
            end do
         end do
      end do
      call Write_InitStatus("FillStaticDensity",1)

    end subroutine FillStaticDensity

    !**************************************************************************
    !****f* updateDensity/updateDensityBox
    ! NAME
    ! subroutine updateDensityBox
    ! PURPOSE
    ! This updates the density in a box dynamically. It also calculates the
    ! momentum distribution
    !**************************************************************************
    subroutine updateDensityBox

      use constants, only: pi,hbarc
      use idtable, only: nucleon, isBaryon
      use hist

      integer :: nEns,iQ,iBin,iRun
      real :: fak,Vbox,mom
      real,dimension(0:3) :: vec

      if (FirstTime) then
         call Write_InitStatus('updateDensityBox',0)
         if (any(pF_Box < 0.)) &
              call TRACEBACK("pF_Box must be set")

         call CreateHist(hP_Box(0), "f_n(p)", 0.0, pF_Box(0)*5.0, pF_Box(0)/20)
         call CreateHist(hP_Box(1), "f_p(p)", 0.0, pF_Box(1)*5.0, pF_Box(1)/20)

         FirstTime=.false.
         call Write_InitStatus('updateDensityBox',1)
      end if

      call ClearHist(hP_Box(0))
      call ClearHist(hP_Box(1))

      Vbox = 8.*gridsize(1)*gridsize(2)*gridsize(3)

      nEns = size(Parts,dim=1)
      fak = 1.0/(nEns*Vbox)

      dichte_Box = dichte0

      vec(0) = 1.

      do iEns=1,nEns
         do iPart=1,size(Parts,dim=2)

            pPart => Parts(iEns,iPart)

            if (pPart%id == 0) cycle
            if (pPart%Id <  0) exit
            if (.not.isBaryon(pPart%id)) cycle

            if (pPart%anti) cycle

            vec(1:3) = pPart%vel(1:3)

            dichte_Box%baryon = dichte_Box%baryon + vec
            if (pPart%ID == nucleon) then
               mom = absMom(pPart)
               call addHist(hP_Box(pPart%charge), mom, 1./(mom**2)*fak*(pi**2*hbarc**3), 1.)
               select case (pPart%charge)
               case (0)
                  dichte_Box%neutron = dichte_Box%neutron + vec
               case (1)
                  dichte_Box%proton  = dichte_Box%proton  + vec
               end select
            end if

         end do
      end do
      dichte_Box = dichte_Box * fak

      call writeHist(hP_Box(0),file="hP_Box0.ORIG.dat")

      do iRun=1,9
         do iQ=0,1
            do iBin=1,19
               if ( ((hP_Box(iQ)%yVal(iBin,1)/hP_Box(iQ)%xBin)-1.0) &
                    *((hP_Box(iQ)%yVal(iBin+1,1)/hP_Box(iQ)%xBin)-1.0) < 0.0) then

                  hP_Box(iQ)%yVal(iBin,:) = (hP_Box(iQ)%yVal(iBin,:)+hP_Box(iQ)%yVal(iBin+1,:))/2
               end if
            end do
         end do

!         call writeHist(hP_Box(0),file="hP_Box0.SMOOTH"//achar(48+iRun)//".dat")
      end do

!      stop

    end subroutine updateDensityBox

  end subroutine updateDensity


  !****************************************************************************
  !****s* densitymodule/storeFieldsRMF
  ! NAME
  ! subroutine storeFieldsRMF
  ! PURPOSE
  ! Store the old values of the omega, rho and density fields.
  !****************************************************************************
  subroutine storeFieldsRMF
    use RMF, only: Tens_flag
    use output, only: DoPr

    integer :: k

    if (DoPr(2)) &
         write(*,*) 'Store the old values of the omega, rho and density fields'

    omegaField_old = omegaField
    rhoField_old   = rhoField
    do k=1,3
       baryonCurrent_old(:,:,:,k) = densField(:,:,:)%baryon(k)
    end do

    if (Tens_flag) sigmaField_old = sigmaField

  end subroutine storeFieldsRMF


  !****************************************************************************
  !****s* densitymodule/updateRMF
  ! NAME
  ! subroutine updateRMF(Parts,doInit)
  ! PURPOSE
  ! Updates the sigma, omega and rho fields needed when the propagation
  ! with the RMF is done.
  ! Updates also the baryon velocities, which are needed in the subsequent
  ! updating of the baryon 4-current by the subroutine updateDensity.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts
  ! * logical, optional :: doInit
  ! NOTES
  ! The scalarDens is computed as well.
  !****************************************************************************
  subroutine updateRMF(Parts,doInit_)

    use idTable
    use particleDefinition
    use RMF
    use ADI
    use constants, only: pi,hbarc,mN,mpi,f_pi
    use output, only: DoPr

    type(particle), dimension(:,:), intent(inOut) :: Parts
    logical, optional, intent(in) :: doInit_

    integer :: i,j,k,k_max
    integer :: niter
    integer :: I1,I2,I3

    real :: rho_lrf, pstar2, mstar, estar, pf_n, pf_p
    real :: funMax, fun, derfun, factor
    real :: eta_x, eta_y, tmp, lhs, error

    logical, parameter :: debug = .true.
    logical, save :: doInit=.true.

    integer, dimension(1:3) :: iPos, iSmall
    integer :: small, large

    if (initFlag) call init
    if (densitySwitch.eq.0) return

    if (present(doInit_)) then
       if (doInit_) doInit=.true.
    end if

    if (.not.doInit .and. densitySwitch.ne.1) then

       ! In the case of static density only
       ! update particle energies (E^*'s) and velocities:
       call energyDeterminationRMF( Parts )
       return

    end if

    if (DoPr(2)) write(*,*) 'Updating Relativistic Mean Fields'

    if (doInit) then

       ! Initialize sigmaField:

       do I1 = -gridPoints(1),gridPoints(1)
          do I2 = -gridPoints(2),gridPoints(2)
             do I3 = -gridPoints(3),gridPoints(3)
                rho_lrf = densField(I1,I2,I3)%baryon(0)**2 &
                     &- dot_product(densField(I1,I2,I3)%baryon(1:3),&
                     &densField(I1,I2,I3)%baryon(1:3))
                rho_lrf = sqrt(max(0.,rho_lrf))
                if (flagWalecka) then
                   sigmaField(I1,I2,I3) = -fshift(rho_lrf)/g_sigma
                else ! Feynman-Hellmann theorem, cf. T.D. Cohen et al., PRC 45, 1881 (1992)
                   sigmaField(I1,I2,I3) = f_pi*(1.-Sigma_N/(mpi*f_pi)**2*rho_lrf*hbarc**3)
                end if
             end do
          end do
       end do
       write(*,*) ' sigmaField initialized'

       do k = 0,3
          omegaField(:,:,:,k) = a_6/g_omega*densField(:,:,:)%baryon(k)
       end do
       write(*,*) ' omegaField initialized'

       if (g_rho /= 0.0) then
          do k = 0,3
             rhoField(:,:,:,k) = &
                  & a_7/g_rho*(densField(:,:,:)%proton(k)-densField(:,:,:)%neutron(k))
          end do
          write(*,*) ' rhoField initialized'
       else
          rhoField(:,:,:,:) = 0.0
       end if

       doInit = .false.

    end if


    if (debug.and.DoPR(2)) write(*,*) ' Iterations to find the sigma field:'


    niter = 0
    Loop_over_iterations : do

       niter = niter + 1

       scalarDens(:,:,:) = 0.
       d_scalarDens_dsigma(:,:,:) = 0.

       if (densitySwitch.eq.1) then

          do j=1,Size(Parts,dim=1)
             do i=1,Size(Parts,dim=2)

                if ( Parts(j,i)%id == 0 ) cycle
                if ( Parts(j,i)%id < 0 ) exit

                if ( Parts(j,i)%id >= pion) cycle ! Only baryons are accounted for presently

                ! Modification factor for the coupling constants:
                factor = ModificationFactor(Parts(j,i)%Id,Parts(j,i)%anti)

                pstar2 = dot_product(Parts(j,i)%mom(1:3),Parts(j,i)%mom(1:3))

                ! position in grid:
                if (.not.getGridIndex(Parts(j,i)%pos,iPos, iSmall,small)) cycle

                large=0

                ! Smearing particle over points in Neighborhood:
                do I1=iPos(1)-nLargePoints,iPos(1)+nLargePoints
                   do I2=iPos(2)-nLargePoints,iPos(2)+nLargePoints
                      do I3=iPos(3)-nLargePoints,iPos(3)+nLargePoints

                         large=large+1

                         if (         abs(I1).le.gridPoints(1) &
                              & .and. abs(I2).le.gridPoints(2) &
                              & .and. abs(I3).le.gridPoints(3) &
                              & .and. smearingWeights(small,large).gt.0. ) then

                            mstar = Parts(j,i)%mass + SelfEnergy_scalar(I1,I2,I3,Parts(j,i)%Id,Parts(j,i)%anti)
                            estar = sqrt( mstar**2 + pstar2 )

                            if (.not.flagWalecka) then
                               select case (Parts(j,i)%id)
                               case (S11_1535)
                                  factor = dmPD(sigmaField(I1,I2,I3),-1)/g_sigma
                               case default
                                  factor = dmPD(sigmaField(I1,I2,I3),1)/g_sigma
                               end select
                            end if

                            scalarDens(I1,I2,I3) = scalarDens(I1,I2,I3) &
                                     & + mstar/estar*smearingWeights(small,large) * factor
                            d_scalarDens_dsigma(I1,I2,I3) &
                                     & = d_scalarDens_dsigma(I1,I2,I3) &
                                     & + g_sigma*pstar2/estar**3*smearingWeights(small,large) * factor**2

                            if (.not.flagWalecka) then
                                d_scalarDens_dsigma(I1,I2,I3) = d_scalarDens_dsigma(I1,I2,I3) &
                                     & + mstar/estar*smearingWeights(small,large) * d2mPD(sigmaField(I1,I2,I3))/g_sigma
                            end if


                         end if

                      end do
                   end do
                end do

             end do
          end do

       else ! Static density profiles are used:

          do I1 = -gridPoints(1),gridPoints(1)
             do I2 = -gridPoints(2),gridPoints(2)
                do I3 = -gridPoints(3),gridPoints(3)

                   mstar = mN + SelfEnergy_scalar(I1,I2,I3,nucleon,.false.)
                   pf_n=(3.*pi**2*densField(I1,I2,I3)%neutron(0))**0.333333*hbarc
                   pf_p=(3.*pi**2*densField(I1,I2,I3)%proton(0))**0.333333*hbarc

                   if (pf_n.gt.1.e-06 .and. pf_p.gt.1.e-06) then

                      if (flagWalecka) then
                         factor = 1.
                      else
                         factor = dmPD(sigmaField(I1,I2,I3),1)/g_sigma
                      end if

                      scalarDens(I1,I2,I3) &
                           & = (  densField(I1,I2,I3)%neutron(0)*f(mstar/pf_n) &
                           &    + densField(I1,I2,I3)%proton(0)*f(mstar/pf_p) )*factor

                      d_scalarDens_dsigma(I1,I2,I3) &
                           & = (  densField(I1,I2,I3)%neutron(0)*g_sigma/pf_n*fprime(mstar/pf_n) &
                           &    + densField(I1,I2,I3)%proton(0)*g_sigma/pf_p*fprime(mstar/pf_p) )*factor**2

                      if (.not.flagWalecka) then
                         d_scalarDens_dsigma(I1,I2,I3) = d_scalarDens_dsigma(I1,I2,I3) &
                              & + (  densField(I1,I2,I3)%neutron(0)*f(mstar/pf_n) &
                              &       + densField(I1,I2,I3)%proton(0)*f(mstar/pf_p) )*d2mPD(sigmaField(I1,I2,I3))/g_sigma
                      end if

                   end if

                end do
             end do
          end do

       end if

       funMax = 0.

       if (grad_flag) then
          if (flagWalecka) then
             tmp = -(gridSpacing(3)*m_sigma/hbarc)**2/g_sigma
          else
             tmp = -(gridSpacing(3)*mu/hbarc)**2/g_sigma
          end if
          eta_x = (gridSpacing(3)/gridSpacing(1))**2
          eta_y = (gridSpacing(3)/gridSpacing(2))**2
       end if

       do I1 = -gridPoints(1),gridPoints(1)
          do I2 = -gridPoints(2),gridPoints(2)
             do I3 = -gridPoints(3),gridPoints(3)

                if (flagWalecka) then

                   fun = g_sigma*sigmaField(I1,I2,I3) &
                        & + a_1*scalarDens(I1,I2,I3) &
                        & + a_2*sigmaField(I1,I2,I3)**2 &
                        & + a_3*sigmaField(I1,I2,I3)**3   ! function we want to have = 0.

                   derfun = g_sigma*( 1. + a_4*sigmaField(I1,I2,I3) &
                        &             + a_5*sigmaField(I1,I2,I3)**2 ) &
                        & +a_1*d_scalarDens_dsigma(I1,I2,I3) ! d fun / d sigma

                else  ! Expression for the PDM

                   fun = -g_sigma*sigmaField(I1,I2,I3) &
                        & + a_1*scalarDens(I1,I2,I3) &
                        & + a_3*sigmaField(I1,I2,I3)**3 &
                        & - a_5*sigmaField(I1,I2,I3)**5 - a_0    ! function we want to have = 0.

                   derfun = -g_sigma + a_1*d_scalarDens_dsigma(I1,I2,I3) &
                        & + 3.*a_3*sigmaField(I1,I2,I3)**2 &
                        & - 5.*a_5*sigmaField(I1,I2,I3)**4    ! d fun / d sigma

                end if



                if (.not.grad_flag) then ! local fields:

                   if ( abs(fun) > funMax ) funMax = abs(fun)

                   if ( derfun /= 0. ) &
                        & sigmaField(I1,I2,I3) = sigmaField(I1,I2,I3) &
                        & - fun/derfun

                else ! gradients included:

                   K_rmf(I1,I2,I3) = tmp*( fun - derfun* sigmaField(I1,I2,I3) )

                   Diag_rmf(I1,I2,I3) = -tmp*derfun/3.

                   if (   abs(I1).lt.gridPoints(1) .and. &
                        & abs(I2).lt.gridPoints(2) .and. &
                        & abs(I3).lt.gridPoints(3) ) then

                      lhs = ( -eta_x * (  sigmaField(I1-1,I2,I3) &
                           &+ sigmaField(I1+1,I2,I3) ) &
                           & -eta_y * ( sigmaField(I1,I2-1,I3) &
                           &+ sigmaField(I1,I2+1,I3) ) &
                           & -sigmaField(I1,I2,I3-1) &
                           & -sigmaField(I1,I2,I3+1) &
                           & + 2.*(eta_x+eta_y+1.) * sigmaField(I1,I2,I3) ) / tmp

                      error = abs(fun - lhs)

                      if ( error .gt. funMax ) funMax = error

                   end if

                end if

             end do
          end do
       end do

       if (debug.and.DoPR(2)) &
            write(*,*) ' In updateRMF, niter, funMax: ', niter, funMax

       if (grad_flag) &
            & call ADI_solve_Douglas(sigmaField, K_rmf, Diag_rmf, eta_x, eta_y)

       if ( funMax <= 1.e-04 ) then
          exit Loop_over_iterations ! ==> success !
       else if ( niter == 50 ) then
          write(*,*) ' In updateRMF: bad convergence after 50 iterations:',&
               & funMax
          call TRACEBACK()
       end if

    end do Loop_over_iterations

    ! Update particle energies (E^*'s) and velocities:
    call energyDeterminationRMF( Parts )
    call updateDensity(Parts)

    if (lorentz_flag) then
       k_max=3    ! With space components of the omega- and rho- fields
    else
       k_max=0    ! W/o space components of the omega- and rho- fields
    end if

    !--------------------------------------------------------------------------
    if (debug.and.DoPR(2)) write(*,*) ' Iterations to find the omega field:'
    !--------------------------------------------------------------------------
    do k = 0,k_max

       if (.not.grad_flag) then ! local fields:

          omegaField(:,:,:,k) = a_6/g_omega*densField(:,:,:)%baryon(k)

       else ! gradients included:

          K_rmf(:,:,:) = gridSpacing(3)**2*g_omega*hbarc*densField(:,:,:)%baryon(k)

          if (k.eq.0) Diag_rmf(:,:,:) = (gridSpacing(3)*m_omega/hbarc)**2 / 3.

          U_rmf(:,:,:) = omegaField(:,:,:,k)

          call ADI_solve_Douglas(U_rmf, K_rmf, Diag_rmf, eta_x, eta_y)

          omegaField(:,:,:,k) = U_rmf(:,:,:)

       end if

    end do

    if (g_rho /= 0.0) then
       !-----------------------------------------------------------------------
       if (debug.and.DoPR(2)) write(*,*) ' Iterations to find the rho field:'
       !-----------------------------------------------------------------------

       do k = 0,k_max

          if (.not.grad_flag) then ! local fields:

             rhoField(:,:,:,k) = &
                  & a_7/g_rho*(densField(:,:,:)%proton(k)-densField(:,:,:)%neutron(k))

          else ! gradients included:

             K_rmf(:,:,:) = gridSpacing(3)**2*g_rho*hbarc*( densField(:,:,:)%proton(k)&
                  & -densField(:,:,:)%neutron(k) )

             if (k.eq.0) Diag_rmf(:,:,:) = (gridSpacing(3)*m_rho/hbarc)**2 / 3.

             U_rmf(:,:,:) = rhoField(:,:,:,k)

             call ADI_solve_Douglas(U_rmf, K_rmf, Diag_rmf, eta_x, eta_y)

             rhoField(:,:,:,k) = U_rmf(:,:,:)

          end if

       end do

    else
       do k = 0,k_max
          rhoField(:,:,:,k) = 0.0
       end do
    end if

  end subroutine updateRMF


  !****************************************************************************
  !****s* densitymodule/TensorRMF
  ! NAME
  ! subroutine TensorRMF(Parts)
  ! PURPOSE
  ! Calculates the energy-momentum tensor and the four-momentum density fields.
  ! INPUTS
  ! * type(particle), dimension(:,:) :: Parts      ! real particle vector
  ! NOTES
  ! Applicable in RMF mode (both Walecka and PDM).
  !****************************************************************************
  subroutine TensorRMF(Parts)

    use particleDefinition

    type(particle), dimension(:,:), intent(in) :: Parts

    real, dimension(0:3) :: momentum
    integer :: j,i,m,n
    integer, dimension(1:3) :: iPos
    integer, dimension(1:3) :: iSmall
    integer :: small, large
    integer :: i1,i2,i3

    Tens(:,:,:,:,:) = 0.
    mom4dens(:,:,:,:)=0.

    do j=1,Size(Parts,dim=1)
       do i=1,Size(Parts,dim=2)

          if ( Parts(j,i)%id == 0 ) cycle
          if ( Parts(j,i)%id < 0 ) exit

          ! position in large and small grid:
          if (.not.getGridIndex( Parts(j,i)%pos, iPos, iSmall,small)) &
               call traceback("outside of grid")


          large=0

          ! Smearing particle over points in Neighborhood:
          do I1=iPos(1)-nLargePoints,iPos(1)+nLargePoints
             do I2=iPos(2)-nLargePoints,iPos(2)+nLargePoints
                do I3=iPos(3)-nLargePoints,iPos(3)+nLargePoints

                   large=large+1

                   if (         abs(I1).le.gridPoints(1) &
                        & .and. abs(I2).le.gridPoints(2) &
                        & .and. abs(I3).le.gridPoints(3) &
                        & .and. smearingWeights(small,large).gt.0. ) then

                      do m=0,3
                         momentum(m) = Parts(j,i)%mom(m) &
                              & + SelfEnergy_vector(I1,I2,I3,m,Parts(j,i)%id,Parts(j,i)%charge,Parts(j,i)%anti)
                      end do

                      ! Particle contribution:
                      do n=0,3
                         do m=0,3
                            Tens(I1,I2,I3,n,m) = Tens(I1,I2,I3,n,m) &
                                 & + momentum(n)*Parts(j,i)%mom(m)/Parts(j,i)%mom(0)*smearingWeights(small,large)
                         end do
                      end do

                      momentum = true4MomentumRMF(Parts(j,i))
                      mom4dens(I1,I2,I3,0:3) = mom4dens(I1,I2,I3,0:3)&
                           & + momentum(0:3)*smearingWeights(small,large)

                   end if

                end do
             end do
          end do

       end do
    end do


    do I1 = -gridPoints(1),gridPoints(1)
       do I2 = -gridPoints(2),gridPoints(2)
          do I3 = -gridPoints(3),gridPoints(3)
             Tens(I1,I2,I3,:,:) = Tens(I1,I2,I3,:,:) + TensorRMF_fields(I1,I2,I3)
          end do
       end do
    end do


  end subroutine TensorRMF


  !****************************************************************************
  !****f* densitymodule/TensorRMF_fields
  ! NAME
  ! function TensorRMF_fields
  ! PURPOSE
  ! Calculates the fields contribution Tens_fields(0:3,0:3) to the
  ! energy-momentum tensor at the grid point specified by the indices I1,I2,I3.
  ! NOTES
  ! Applicable in RMF mode (both Walecka and PDM).
  !****************************************************************************
  function TensorRMF_fields(I1,I2,I3)

    use RMF
    use inputGeneral, only: delta_T
    use constants, only: hbarc

    integer, intent(in) :: I1,I2,I3
    real, dimension(0:3,0:3) :: TensorRMF_fields

    real, dimension(0:3,0:3) :: domega,drho  ! time-space derivatives of vector fields (GeV/fm),
    ! 1-st index -> Lorentz index of field component,
    ! 2-nd index -> Lorentz index of the derivative
    real,  dimension(0:3) :: dsigma  ! time-space derivatives of the scalar field (GeV/fm)
    real :: Vpot_omega,Vpot_rho,U,Pot
    integer :: m,n,iu

    TensorRMF_fields = 0.

    ! Local potential vector field contributions:
    Vpot_omega = omegaField(I1,I2,I3,0)**2
    Vpot_rho = rhoField(I1,I2,I3,0)**2
    if (lorentz_flag) then
       Vpot_omega = Vpot_omega - dot_product(omegaField(I1,I2,I3,1:3),omegaField(I1,I2,I3,1:3))
       Vpot_rho = Vpot_rho - dot_product(rhoField(I1,I2,I3,1:3),rhoField(I1,I2,I3,1:3))
    end if
    Vpot_omega = -0.5*m_omega**2*Vpot_omega/hbarc**3
    Vpot_rho = -0.5*m_rho**2*Vpot_rho/hbarc**3

    ! Local potential scalar field contributions:
    if (flagWalecka) then

       U = (  0.5*(m_sigma*sigmaField(I1,I2,I3))**2 &
            &  + g_2*sigmaField(I1,I2,I3)**3/3. &
            &  + g_3*sigmaField(I1,I2,I3)**4/4. ) / hbarc**3

    else  ! PDM

       U = ( -0.5*(mu*sigmaField(I1,I2,I3))**2 &
            &  + lambda*sigmaField(I1,I2,I3)**4/4. &
            &  - lambda6*sigmaField(I1,I2,I3)**6/6. &
            &  - epsilon*sigmaField(I1,I2,I3) ) / hbarc**3 - endens_vac

    end if

    ! Non-local contributions:
    if (grad_flag) then

       dsigma(0) = (sigmaField(I1,I2,I3)-sigmaField_old(I1,I2,I3))/delta_T
       domega(:,0) = (omegaField(I1,I2,I3,:)-omegaField_old(I1,I2,I3,:))/delta_T
       drho(:,0) = (rhoField(I1,I2,I3,:)-rhoField_old(I1,I2,I3,:))/delta_T

       if (I1.lt.gridPoints(1)) then
          iu = I1+1
       else
          iu = I1
       end if

       ! Minus because calculated are the  contravariant (upper Lorentz index) components:
       dsigma(1) = -(sigmaField(iu,I2,I3)-sigmaField(iu-1,I2,I3))/gridSpacing(1)
       domega(:,1) = -(omegaField(iu,I2,I3,:)-omegaField(iu-1,I2,I3,:))/gridSpacing(1)
       drho(:,1) = -(rhoField(iu,I2,I3,:)-rhoField(iu-1,I2,I3,:))/gridSpacing(1)

       if (I2.lt.gridPoints(2)) then
          iu = I2+1
       else
          iu = I2
       end if

       dsigma(2) = -(sigmaField(I1,iu,I3)-sigmaField(I1,iu-1,I3))/gridSpacing(2)
       domega(:,2) = -(omegaField(I1,iu,I3,:)-omegaField(I1,iu-1,I3,:))/gridSpacing(2)
       drho(:,2) = -(rhoField(I1,iu,I3,:)-rhoField(I1,iu-1,I3,:))/gridSpacing(2)

       if (I3.lt.gridPoints(3)) then
          iu = I3+1
       else
          iu = I3
       end if

       dsigma(3) = -(sigmaField(I1,I2,iu)-sigmaField(I1,I2,iu-1))/gridSpacing(3)
       domega(:,3) = -(omegaField(I1,I2,iu,:)-omegaField(I1,I2,iu-1,:))/gridSpacing(3)
       drho(:,3) = -(rhoField(I1,I2,iu,:)-rhoField(I1,I2,iu-1,:))/gridSpacing(3)

       U = U + 0.5*dot_product(dsigma(1:3),dsigma(1:3))/hbarc

       Vpot_omega = Vpot_omega - 0.5*(  dot_product(domega(0,1:3),domega(0,1:3)) &
            & -dot_product(domega(1,1:3),domega(1,1:3)) &
            & -dot_product(domega(2,1:3),domega(2,1:3)) &
            & -dot_product(domega(3,1:3),domega(3,1:3)) )/hbarc

       Vpot_rho = Vpot_rho - 0.5*(  dot_product(drho(0,1:3),drho(0,1:3)) &
            & -dot_product(drho(1,1:3),drho(1,1:3)) &
            & -dot_product(drho(2,1:3),drho(2,1:3)) &
            & -dot_product(drho(3,1:3),drho(3,1:3)) )/hbarc

       ! Kinetic field contribution:
       do n=0,3
          do m=1,3   ! m=0 component disappears for static Lagrangians
             tensorRMF_fields(n,m) = (  dsigma(n)*dsigma(m) &
                  &    - domega(0,n)*domega(0,m) + dot_product(domega(1:3,n),domega(1:3,m)) &
                  &    - drho(0,n)*drho(0,m) + dot_product(drho(1:3,n),drho(1:3,m))  )/hbarc
          end do
       end do

    end if

    Pot = U + Vpot_omega + Vpot_rho

    tensorRMF_fields(0,0) = tensorRMF_fields(0,0) + Pot
    tensorRMF_fields(1,1) = tensorRMF_fields(1,1) - Pot
    tensorRMF_fields(2,2) = tensorRMF_fields(2,2) - Pot
    tensorRMF_fields(3,3) = tensorRMF_fields(3,3) - Pot

  end function TensorRMF_fields


  !****************************************************************************
  !****s* densitymodule/TotalEnergyRMF
  ! NAME
  ! subroutine TotalEnergyRMF
  ! PURPOSE
  ! Calculates the total energy of the system (w/o Coulomb part).
  ! INPUTS
  ! * type(particle), dimension(:,:), optional :: Parts -- real particle vector
  ! OUTPUT
  ! *  real :: Etot   ! total energy (GeV)
  ! NOTES
  ! Applicable in RMF mode (both Walecka and PDM).
  !****************************************************************************
  subroutine TotalEnergyRMF(Etot,Parts)

    use particleDefinition

    real, intent(out) :: Etot   ! total energy (GeV)
    type(particle), dimension(:,:), optional, intent(in) :: Parts

    integer, parameter :: imode=2    ! 1 - particle contribution from energy-momentum tensor on the grid
                                     ! 2 - particle contribution by the direct summation

    integer :: j,i
    real, dimension(0:3) :: momentum
    real :: dV
    real, dimension(0:3,0:3) :: fields
    integer :: I1,I2,I3

    Etot = 0.

    dV = gridSpacing(1)*gridSpacing(2)*gridSpacing(3)

    select case (imode)

    case (1)

        do I1 = -gridPoints(1),gridPoints(1)
          do I2 = -gridPoints(2),gridPoints(2)
             do I3 = -gridPoints(3),gridPoints(3)
                Etot = Etot + Tens(I1,I2,I3,0,0)*dV
             end do
          end do
        end do

    case (2)

        if (.not.present(Parts)) then
           call TRACEBACK('Input particle vector in missed in TotalEnergyRMF')
        end if

        do j=1,Size(Parts,dim=1)
           do i=1,Size(Parts,dim=2)

              if ( Parts(j,i)%id == 0 ) cycle
              if ( Parts(j,i)%id < 0 ) exit

              momentum = Particle4MomentumRMF(Parts(j,i))
              Etot = Etot + momentum(0)

           end do
        end do

        Etot = Etot/Size(Parts,dim=1)

        do I1 = -gridPoints(1),gridPoints(1)
           do I2 = -gridPoints(2),gridPoints(2)
              do I3 = -gridPoints(3),gridPoints(3)
                 fields = TensorRMF_fields(I1,I2,I3)
                 Etot = Etot + fields(0,0)*dV
              end do
           end do
        end do

    end select

  end subroutine TotalEnergyRMF



  !****************************************************************************
  ! cf. interface energyDeterminationRMF
  !****************************************************************************
  subroutine energyDeterminationRMF0(Part)

    use particleDefinition
    use idTable, only: kaon, kaonBar, rho, omegaMeson, phi, isBaryon
    use RMF, only: ModificationFactor
    use dichteDefinition
    use minkowski, only: abs4
    use mesonPotentialMain, only: vecMes_massShift

    type(particle),intent(inOut) :: Part

    type(dichte):: density
    real    :: pstar2, mstar, factor, densityLRF, spot
    logical :: flagOk

    ! Check input:
    if ( Part%Id <= 0 ) then
       write(*,*) 'wrong input particle:', Part%Id
       call TRACEBACK()
    end if

    factor = ModificationFactor(Part%Id,Part%anti)

    pstar2 = dot_product(Part%mom(1:3),Part%mom(1:3))

    !for mesons mean-field optionally only for kaons and antikaons

    if ( (isBaryon(Part%Id) .or. Part%Id.eq.Kaon .or. Part%Id.eq.kaonBar) &
         & .and. factor.gt.0. ) then

       mstar = DiracMass(Part)
       Part%mom(0) = sqrt( mstar**2 + pstar2 )

    else

       select case (Part%id)
       case (rho,omegaMeson,phi)
          ! set scalar potential for vector mesons (mass shift!)
          density = densityAt(Part%pos)
          densityLRF=abs4(density%baryon,flagOk)
          if (flagOk) then
             spot = vecMes_massShift(Part%id,densityLRF)
          else
             spot=0.
          end if
       case default
          ! in all other cases: neglect potential
          spot = 0.
       end select

       Part%mom(0) = sqrt( (Part%mass+spot)**2 + pstar2 )

    end if

  end subroutine energyDeterminationRMF0
  !----------------------------------------------------------------------------
  subroutine energyDeterminationRMF2(Parts)
    use particleDefinition

    type(particle), dimension(:,:), intent(inout) :: Parts

    integer :: i,j

    do i=1,Size(Parts,dim=1)
       do j=1,Size(Parts,dim=2)
          if ( Parts(i,j)%id == 0 ) cycle
          if ( Parts(i,j)%id < 0 ) exit

          if (abs(Parts(i,j)%offshellPar).gt.1.e-06) cycle

          call energyDeterminationRMF0( Parts(i,j) )
          Parts(i,j)%vel(1:3) = Parts(i,j)%mom(1:3) / Parts(i,j)%mom(0)
       end do
    end do

  end subroutine energyDeterminationRMF2

  !****************************************************************************

  !****************************************************************************
  !****f* densitymodule/Particle4MomentumRMF
  ! NAME
  ! function Particle4MomentumRMF(Part) result(momentum)
  ! PURPOSE
  ! This subroutine determines the canonical four-momentum
  ! in computational frame (i.e. where mesonic mean fields are given).
  ! INPUTS
  ! * type(particle),intent(in) :: Part ---  Particle
  ! OUTPUT
  ! * real, dimension(0:3) :: momentum --- canonical 4-mom of the particle
  ! NOTES
  ! Should be used in RMF-mode. It is supposed, that the kinetic
  ! four-momentum of the particle is already determined by the subroutine
  ! energyDeterminationRMF earlier.
  !
  ! Electromagnetic part is not included.
  !****************************************************************************
  function Particle4MomentumRMF(Part) result(momentum)

    use particleDefinition
    use idTable, only: kaon,kaonBar,isBaryon
    use RMF, only: ModificationFactor

    real, dimension(0:3) :: momentum

    type(particle), intent(in) :: Part

    real :: fact
    integer, dimension(1:3) :: II

    ! Check input:
    if ( Part%Id <= 0 ) then
       write(*,*) 'wrong input particle', Part%Id
       call TRACEBACK()
    end if

    fact = ModificationFactor(Part%Id,Part%anti,.true.)

    momentum = Part%mom ! default return value

    if ( (isBaryon(Part%Id) .or. Part%id.eq.Kaon .or. Part%id.eq.kaonBar) &
         & .and. fact > 0. ) then

       if (getGridIndex(Part%pos, II)) then
          momentum(0:3) = momentum(0:3) &
               + SelfEnergy_vector(II(1),II(2),II(3),(/0,1,2,3/),&
               &  Part%id,Part%charge,Part%anti)
       end if

    end if

  end function Particle4MomentumRMF

  !****************************************************************************
  !****f* densitymodule/true4MomentumRMF
  ! NAME
  ! function true4MomentumRMF(Part,inside_grid_flag) result(momentum)
  ! PURPOSE
  ! This subroutine determines the "true" single-particle 4-momentum,
  ! i.e. the 4-momentum which is additive to produce the total
  ! 4-momentum of the system in the computational frame.
  ! INPUTS
  ! * type(particle) :: Part --- particle whose "true" 4-momentum should
  !   be calculated.
  ! OUTPUT
  ! * real, dimension(0:3) :: momentum ---  "true" 4-momentum of the particle
  ! * logical, optional :: inside_grid_flag --- .true. if the particle is
  !   inside grid, .false. otherwise
  ! NOTES
  ! * Should be used in RMF-mode. The input particle kinetic 4-momentum must
  !   be given in the computational frame (where the density field is defined).
  ! * Formula for the fieldenergy updated; gradient terms replaced the
  !   corresponding sources using the meson-field equations and partial
  !   integrations.
  !****************************************************************************
  function true4MomentumRMF(Part,inside_grid_flag) result(momentum)

    use constants, only: hbarc
    use particleDefinition
    use idTable, only: kaon,kaonBar,isMeson
    use RMF

    real, dimension(0:3) :: momentum

    type(particle),       intent(in) :: Part
    logical, optional,    intent(out) :: inside_grid_flag

    integer, dimension(1:3) :: I
    real    :: VectorFieldEnergy, ScalarFieldEnergy, fieldEnergy, bar_plus_antibar_density, factor

    real, dimension(0:3) :: BaryonCurrent, IsospinCurrent

    !--------------------------------------------------------------------------

    if (present(inside_grid_flag)) inside_grid_flag = .false.

    factor = ModificationFactor(Part%Id,Part%anti)

    if ( Part%ID <= 0 ) then
       write(*,*) 'wrong input particle Id', Part%ID
       call TRACEBACK()

    else if ( isMeson(Part%Id) .and. &
         &  ( Part%Id.ne.Kaon.and.Part%Id.ne.kaonBar .or. &
         &  ( Part%id.eq.Kaon.or.Part%Id.eq.kaonBar).and.factor.eq.0. ) ) then

       momentum(0:3) = Part%mom(0:3)
       return

    end if

    !--------------------------------------------------------------------------

    if (getGridIndex(Part%pos, I)) then

       if (present(inside_grid_flag)) inside_grid_flag = .true.

       ! Energy density of the fields:
       ! valid in RMF with and without gradient terms

       BaryonCurrent(:)  = densField(I(1),I(2),I(3))%baryon(:)
       IsospinCurrent(:) = densField(I(1),I(2),I(3))%proton(:)-densField(I(1),I(2),I(3))%neutron(:)


       if (lorentz_flag) then  ! With space components of the omega- & rho-fields:

          VectorFieldEnergy = &
               & 0.5*g_omega*( &
               &                 dot_product( omegaField(I(1),I(2),I(3),1:3),BaryonCurrent(1:3) )&
               &               - BaryonCurrent(0)*omegaField(I(1),I(2),I(3),0) ) &
               & + 0.5*g_rho*( &
               &                 dot_product(IsospinCurrent(1:3),rhoField(I(1),I(2),I(3),1:3))&
               &               - IsospinCurrent(0)*rhoField(I(1),I(2),I(3),0) )

       else ! W/o space components of the omega- & rho-fields:

          VectorFieldEnergy = - 0.5*g_omega*BaryonCurrent(0)*omegaField(I(1),I(2),I(3),0) &
                            & - 0.5*g_rho*IsospinCurrent(0)*rhoField(I(1),I(2),I(3),0)

       end if


       if (flagWalecka) then

          ScalarFieldEnergy = - 0.5*g_sigma*scalarDens(I(1),I(2),I(3))*sigmaField(I(1),I(2),I(3)) &
               & - ( g_2*sigmaField(I(1),I(2),I(3))**3/6. &
               &    +g_3*sigmaField(I(1),I(2),I(3))**4/4. ) / hbarc**3

       else ! PDM

          ScalarFieldEnergy = - 0.5*g_sigma*scalarDens(I(1),I(2),I(3))*sigmaField(I(1),I(2),I(3)) &
               &            - (  lambda*sigmaField(I(1),I(2),I(3))**4/4. - lambda6*sigmaField(I(1),I(2),I(3))**6/3. &
               &               + epsilon*sigmaField(I(1),I(2),I(3))/2. ) / hbarc**3 - endens_vac
       end if

       fieldEnergy = VectorFieldEnergy + ScalarFieldEnergy

       bar_plus_antibar_density = totalDens(I(1),I(2),I(3))

       if ( bar_plus_antibar_density > 1.e-06 ) then

          momentum(0) = Part%mom(0) + SelfEnergy_vector(I(1),I(2),I(3),0,Part%Id,Part%charge,Part%anti) &
               &      + fieldEnergy / bar_plus_antibar_density
       else

          momentum(0) = Part%mom(0)

       end if

       if (lorentz_flag) then
          momentum(1:3) = Part%mom(1:3) &
               + SelfEnergy_vector(I(1),I(2),I(3),(/1,2,3/),Part%Id,Part%charge,Part%anti)
       else
          momentum(1:3) = Part%mom(1:3)
       end if

    else

       momentum(0:3) = Part%mom(0:3)

    end if

  end function true4MomentumRMF


  !****************************************************************************
  ! cf. interface DiracMass
  !****************************************************************************
  real function DiracMass1(i1,i2,i3,barMass,id,charge,antiFlag)
    use IdTable, only: kaon,kaonBar
    use RMF, only: kaonpot_flag

    integer, intent(in) :: i1,i2,i3,id,charge
    logical, intent(in) :: antiFlag
    real,    intent(in) :: barMass

    real :: Masse, S
    real, dimension(0:3) :: V

    S = SelfEnergy_scalar(i1,i2,i3,id,antiFlag)

    if ( kaonpot_flag .and. (id==kaon .or. id==kaonBar) ) then
       V(0:3) = Selfenergy_vector(i1,i2,i3,(/0,1,2,3/),id,charge,antiFlag)
       Masse = sqrt(barMass**2-S+dot_product(V,V))
    else
       Masse = barMass + S
    end if

    DiracMass1 = Masse

  end function DiracMass1
  !----------------------------------------------------------------------------
  real function DiracMass2(Part)
    use particleDefinition, only: particle

    type(particle), intent(in) :: Part

    integer, dimension(1:3) :: II

    if (getGridIndex(Part%pos, II)) then
       DiracMass2 = DiracMass1(II(1),II(2),II(3),&
            Part%mass,Part%id,Part%charge,Part%anti)
    else
       DiracMass2 = Part%mass
    end if
  end function DiracMass2
  !****************************************************************************

  !****************************************************************************
  !****f* densitymodule/SelfEnergy_vector
  ! NAME
  ! real function SelfEnergy_vector(i1,i2,i3,k,id,charge,antiFlag)
  ! PURPOSE
  ! calculates the vector component of the RMF selfenergy of a particle
  !
  ! INPUTS
  ! * integer :: i1,i2,i3 -- index in the grid
  ! * integer :: k -- component of the vector (0..3)
  ! * integer :: id -- id of the particle
  ! * integer :: charge -- charge of the particle
  ! * logical :: antiFlag -- antiparticle or not
  !
  ! NOTES
  ! * gv_kaon = 3/(8 * f_pi*^2)
  !****************************************************************************
  elemental real function SelfEnergy_vector(i1,i2,i3,k,id,charge,antiFlag)
    use constants, only: hbarc
    use IdTable, only: nucleon,kaon,kaonBar
    use RMF, only: ModificationFactor, g_omega, g_rho, gv_kaon, kaonpot_flag

    integer, intent(in) :: i1,i2,i3,k,id,charge
    logical, intent(in) :: antiFlag

    real :: field, fac, isofac

    if ( kaonpot_flag .and.  (id == kaon .or. id==kaonBar) ) then
       fac = 1.
       if (id==kaonBar) fac = -1.
       field = densField(i1,i2,i3)%baryon(k) * gv_kaon * fac * hbarc**3 !GeV
    else
       fac = ModificationFactor(id,antiFlag,.true.)
       if (antiFlag) fac = -fac
       if (ID==nucleon) then
          if (charge==0) then
             isofac = -1.
          else
             isofac = 1.
          end if
       else
          isofac = 0.
       end if
       field = (omegaField(i1,i2,i3,k)*g_omega + rhoField(i1,i2,i3,k)*g_rho*isofac) * fac
    end if

    SelfEnergy_vector = field

  end function SelfEnergy_vector
  !****************************************************************************

  !****************************************************************************
  ! cf. interface SelfEnergy_scalar
  !****************************************************************************
  real function SelfEnergy_scalar1(i1,i2,i3,id,antiFlag)
    use constants, only: hbarc
    use IdTable, only: kaon,kaonBar,S11_1535
    use RMF, only: ModificationFactor, flagWalecka, g_sigma, gs_kaon, kaonpot_flag, mPD, m_nucleon, m_minus

    integer, intent(in) :: i1,i2,i3,id
    logical, intent(in) :: antiFlag

    real :: field, fac, sigma

    if (densitySwitch.eq.0) then
       SelfEnergy_scalar1 = 0.
       return
    end if

    if ( kaonpot_flag .and. (id == kaon .or. id==kaonBar) ) then
       field = scalarDens(i1,i2,i3) * gs_kaon * hbarc**3 !GeV**2 !!
    else ! baryons
       fac = ModificationFactor(id,antiFlag)
       sigma = sigmaField(i1,i2,i3)
       if (flagWalecka) then
          field = sigma*g_sigma*fac !GeV
       else
          select case (id)
          case (S11_1535)
             ! please note, that m_minus is maybe not hadron(S11_1535)%mass
             field = mPD(sigma,-1) - m_minus
          case default
             field = mPD(sigma, 1) - m_nucleon
          end select
       end if
    end if

    SelfEnergy_scalar1 = field

  end function SelfEnergy_scalar1
  !----------------------------------------------------------------------------
  real function SelfEnergy_scalar2(Part)
     use particleDefinition, only: particle

    type(particle), intent(in) :: Part

    integer, dimension(1:3) :: II

    if (getGridIndex(Part%pos, II)) then
       SelfEnergy_scalar2 = SelfEnergy_scalar1(II(1),II(2),II(3),&
            Part%id,Part%anti)
    else
       SelfEnergy_scalar2 = 0.
    end if
  end function SelfEnergy_scalar2
  !****************************************************************************

  !****************************************************************************
  !****f* densitymodule/SelfEnergy_vector_old
  ! NAME
  ! real function SelfEnergy_vector_old(i1,i2,i3,k,id,charge,antiFlag)
  ! PURPOSE
  ! calculates the vector component of the RMF selfenergy of a particle
  ! from the previous time step;
  ! needed for time derivatives in propagation_RMF
  !
  ! INPUTS
  ! * integer :: i1,i2,i3 -- index in the grid
  ! * integer :: k -- component of the vector (0..3)
  ! * integer :: id -- id of the particle
  ! * integer :: charge -- charge of the particle
  ! * logical :: antiFlag -- antiparticle or not
  !
  ! NOTES
  ! * gv_kaon = 3/(8 * f_pi*^2)
  !****************************************************************************
  real function SelfEnergy_vector_old(i1,i2,i3,k,id,charge,antiFlag)
    use constants, only: hbarc
    use IdTable, only: nucleon,kaon,kaonBar
    use RMF, only: ModificationFactor, g_omega, g_rho, gv_kaon, kaonpot_flag

    integer, intent(in) :: i1,i2,i3,k,id,charge
    logical, intent(in) :: antiFlag

    real :: field, fac, isofac

    if ( kaonpot_flag .and. (id == kaon .or. id==kaonBar) ) then
       fac = 1.
       if (id==kaonBar) fac = -1.
       field = baryonCurrent_old(i1,i2,i3,k) * gv_kaon * fac * hbarc**3 !GeV**2 !!!
    else
       fac = ModificationFactor(id,antiFlag,.true.)
       if (antiFlag) fac = -fac
       if (ID==nucleon) then
          if (charge==0) then
             isofac = -1.
          else
             isofac = 1.
          end if
       else
          isofac = 0.
       end if
       field = (omegaField_old(i1,i2,i3,k)*g_omega + rhoField_old(i1,i2,i3,k)*g_rho*isofac) * fac
    end if

    SelfEnergy_vector_old = field


  end function SelfEnergy_vector_old


  !****************************************************************************
  !****f* densitymodule/getDensitySwitch
  ! NAME
  ! function getDensitySwitch()
  ! PURPOSE
  ! If not yet done, reads input from jobcard and then returns densitySwitch.
  !****************************************************************************
  integer function getDensitySwitch()
    if (initFlag) call init
    getDensitySwitch=densitySwitch
  end function getDensitySwitch

  !****************************************************************************
  !****************************************************************************
  subroutine setDensitySwitch(d)
    integer,intent(in) :: d
    densitySwitch = d
  end subroutine setDensitySwitch

  !****************************************************************************
  !****************************************************************************
  subroutine setDensityInput(pro,neu)
    real, intent(in) :: pro, neu
    densityinput_proton = pro
    densityinput_neutron = neu
  end subroutine setDensityInput


  !****************************************************************************
  !****f* densitymodule/getGridSpacing0
  ! NAME
  ! function getGridSpacing0()
  ! PURPOSE
  ! returns gridSpacing
  ! RESULT
  ! real, dimension(1:3) :: gridSpacing -- in units of fermi
  !****************************************************************************
  function getGridSpacing0()
    real, dimension(1:3) :: getGridSpacing0
    if (initFlag) call init
    getGridSpacing0(1:3)=gridSpacing(1:3)
  end function getGridSpacing0


  !****************************************************************************
  !****f* densitymodule/getGridSpacing
  ! NAME
  ! function getGridSpacing()
  ! PURPOSE
  ! returns gridSpacing, but rescaled when linear interpolation is used.
  ! This is used for calculating derivatives.
  ! RESULT
  ! real, dimension(1:3) :: gridSpacing -- in units of fermi
  !****************************************************************************
  function getGridSpacing()
    real, dimension(1:3) :: getGridSpacing

    if (initFlag) call init
    if (linearInterpolation) then
       getGridSpacing(1:3)=gridSpacing(1:3)/4.
    else
       getGridSpacing(1:3)=gridSpacing(1:3)
    end if
  end function getGridSpacing


  !****************************************************************************
  !****f* densitymodule/getGridPoints
  ! NAME
  ! function getGridPoints()
  ! PURPOSE
  ! returns  GridPoints
  ! RESULT
  ! integer, dimension(1:3) :: gridPoints -- numbers
  !****************************************************************************
  function getGridPoints()
    integer, dimension(1:3) :: getGridPoints

    if (initFlag) call init
    getGridPoints(1:3)=gridPoints(1:3)
  end function getGridPoints

  !****************************************************************************
  ! cf. interface getGridIndex
  !****************************************************************************
  logical function getGridIndex1(pos,ind)
    real, dimension(1:3), intent(in) :: pos
    integer, dimension(1:3), intent(out) :: ind

    getGridIndex1 = .false.
    ind = nint(pos/gridSpacing)
    if (any(abs(pos).gt.gridSize+gridSpacing)) return ! --> failure
    if (any(abs(ind).gt.gridpoints)) return ! --> failure
    getGridIndex1 = .true.
  end function GetGridIndex1
  !----------------------------------------------------------------------------
  logical function getGridIndex2(pos,ind,add)
    real, dimension(1:3), intent(in) :: pos
    integer, dimension(1:3), intent(out) :: ind
    integer, intent(in) :: add

    getGridIndex2 = .false.
    ind = nint(pos/gridSpacing)
    if (any(abs(pos).gt.gridSize+gridSpacing)) return ! --> failure
    if (any(abs(ind).gt.gridpoints+add)) return ! --> failure
    getGridIndex2 = .true.
  end function getGridIndex2
  !----------------------------------------------------------------------------
  logical function getGridIndex3(pos,ind,iSmall,small)
    real, dimension(1:3), intent(in) :: pos
    integer, dimension(1:3), intent(out) :: ind
    integer, dimension(1:3), intent(out) :: iSmall
    integer, intent(out) :: small
    logical :: flag

    getGridIndex3 = .false.

    flag = getGridIndex1(pos,ind)

    ! position in small grid:
    iSmall=NINT((pos/gridSpacing-float(ind))*(2.*float(nSmallPoints)+0.999999))

    small=1+(nSmallPoints+iSmall(3))  &
         & +(nSmallPoints+iSmall(2))*(2*nSmallPoints+1)       &
         & +(nSmallPoints+iSmall(1))*(2*nSmallPoints+1)**2

    if (.not.flag) return ! --> failure
    if (any(abs(iSmall).gt.nSmallPoints)) return ! --> failure

    getGridIndex3 = .true.

  end function getGridIndex3
  !****************************************************************************

  !****************************************************************************
  !****s* densitymodule/boostToLRF
  ! NAME
  ! subroutine boostToLRF(Part, switch, density)
  ! PURPOSE
  ! Boosts particle between "Local Rest Frame (LRF)" and
  ! "calculation frame (CF)".
  !
  ! INPUTS
  ! * type(particle), intent(inout) :: Part  ! Particle to be boosted
  ! * integer, intent(in) :: switch          ! 1= boost from CF to LRF
  !                                          ! 2= boost from LRF to CF
  ! NOTES
  ! * The LRF is the frame in which the baryon current vanishes.
  ! * If the density is very small, then no boost takes place.
  ! * This boost may not work in RMF mode, since %baryon(0) may vanish,
  !   even if the vector components are non-zero. (In RMF mode this
  !   boost should not be necessary anyway.)
  !****************************************************************************
  subroutine boostToLRF(Part, switch, density_in)
    use ParticleDefinition
    use lorentzTrafo, only: lorentz

    type(particle), intent(inout) :: Part
    integer, intent(in) :: switch
    type(dichte), intent(in), optional :: density_in

    real, dimension(1:3) :: beta  ! beta of LRF in calculation frame
    type(dichte) :: density

    if (present(density_in)) then
       density=density_in
    else
       density=densityAt(Part%pos(1:3))
    end if

    if (abs(density%baryon(0))<1E-8) return  ! do nothing if density is too small

    beta(1:3) = density%baryon(1:3)/density%baryon(0)

    ! Boost particle's momentum
    select case (switch)
    case (1) ! The particle is boosted to LRF out of CF
       call lorentz( beta,Part%mom(0:3))
    case (2) ! The particle is boosted to calculation frame out of LRF.
       call lorentz(-beta,Part%mom(0:3))
    case default
       write(*,*) 'Wrong value for switch:', switch
       call TRACEBACK()
    end select

  end subroutine boostToLRF


  !****************************************************************************
  !****f* densitymodule/fermiMomentum_sym
  ! NAME
  ! real function fermiMomentum_sym(rho)
  ! PURPOSE
  ! Evaluate the fermi momentum for symmetric nuclear matter.
  ! Assumption: rho_p=rho_n !!!
  ! INPUTS
  ! * real :: rho ! density in fm^-3
  ! RESULT
  ! * Fermi momentum in GeV
  !****************************************************************************
  real function fermiMomentum_sym(rho)
    use constants, only: hbarc, pi

    real, intent(in) :: rho ! Density in fm^-3
    ! Take care of improper input by checking rho:
    if (rho.lt.-1E-20) then
       write(*,*) 'ERROR in fermiMomentum_sym: Rho less than zero!!!'
       stop 'density.f90: fermiMomentum_sym'
    else if (rho.lt.1E-20) then
       fermiMomentum_sym=0.
    else
       fermiMomentum_sym=(3./2.*pi**2*rho*hbarc**3)**(1./3.)
    end if
  end function fermiMomentum_sym

  !****************************************************************************
  !****f* densitymodule/fermiMomentum_noIsospin
  ! NAME
  ! real function fermiMomentum_noIsospin(rho)
  ! PURPOSE
  ! Evaluate the fermi momentum for a gas of fermions.
  ! INPUTS
  ! * real :: rho ! density in fm^-3
  !   ATTENTION: Do not use rho(nucleon) but rho(proton) or rho(neutron)
  ! RESULT
  ! * Fermi momentum in GeV
  !****************************************************************************
  real function fermiMomentum_noIsospin(rho)
    use constants, only: hbarc, pi

    real, intent(in) :: rho ! Density in fm^-3
    ! Take care of improper input by checking rho:
    if (rho.lt.-1E-20) then
       call TRACEBACK('ERROR in fermiMomentum_sym: Rho less than zero!!!')
    else if (rho.lt.1E-20) then
       fermiMomentum_noIsospin=0.
    else
       fermiMomentum_noIsospin=(3.*pi**2*rho*hbarc**3)**(1./3.)
    end if
  end function fermiMomentum_noIsospin

  !****************************************************************************
  !****f* densitymodule/FermiMomAt
  ! NAME
  ! real function FermiMomAt(pos,charge)
  !
  ! PURPOSE
  ! calculate the fermi momentum at some position. If no charge is given
  ! it uses the averaged proton/neutron densities, otherwise it uses
  ! the corresponding isospin channel
  !
  ! INPUTS
  ! * real, dimension(1:3) :: pos -- the space coordinates
  ! * integer, OPTIONAL :: charge -- the isospin channel
  !****************************************************************************
  real function FermiMomAt(pos,charge)
    use nucleusDefinition, only: tNucleus
    use minkowski, only: abs4
    use nucleus, only: getTarget
    use constants, only: pi,hbarc
    use densityStatic, only: staticDensity

    real,dimension(1:3), intent(in) :: pos
    integer, intent(in), optional :: charge

    real :: rho
    type(tNucleus),pointer :: nuc
    type(dichte) :: dens

    if (initFlag) call init
    select case (densitySwitch)
    case (2) !Static density of a resting target
       nuc => getTarget()
       dens = staticDensity(pos,nuc)
       if (nuc%ReAdjustForConstBinding) then
          dens%proton  = dens%proton/nuc%facP
          dens%neutron = dens%neutron/nuc%facN
          dens%baryon  = dens%proton + dens%neutron
       end if
    case default
       dens = densityAt(pos)
    end select

    if (present(charge)) then
       select case (charge)
       case (0)
          rho = abs4(dens%neutron)
       case (1)
          rho = abs4(dens%proton)
       end select
    else
       rho = abs4(dens%baryon)/2
    end if
    FermiMomAt = (3.*pi**2*rho)**(1./3.)*hbarc

  end function FermiMomAt


  !****************************************************************************
  !****s* densitymodule/writeDensityPlane
  ! NAME
  ! subroutine writeDensityPlane(fName,iPlane)
  !
  ! PURPOSE
  ! Write the density to a file as cuts according some planes
  !
  ! INPUTS
  ! * character*(*) :: fName  -- name of file to write to
  ! * integer :: iPlane -- selection of plane to write (see below)
  !
  ! possible values for iPlane:
  ! * 1: yz-plane (through origin) aka x = 0
  ! * 2: xz-plane (through origin) aka y = 0
  ! * 3: xy-plane (through origin) aka z = 0
  !
  !****************************************************************************
  subroutine writeDensityPlane(fName,iPlane)

    character*(*),intent(in) :: fName
    integer, intent(in) :: iPlane

    integer, parameter :: iFile = 645
    integer, dimension(1:3) :: indexOrigin
    integer :: i1, i2
    real, dimension(1:3) :: r

    if (.not.getGridIndex( (/ 0.,0.,0. /), indexOrigin, 0)) then
       call TraceBack('failure: Origin not in Grid!')
    end if

    open(iFile, file=fName, status='unknown')
    rewind(iFile)

    select case (iPlane)
    case (1)
       do i1 = -gridPoints(2),gridPoints(2)
          do i2 = -gridPoints(3),gridPoints(3)
             r = (/indexOrigin(1),i1,i2/) * gridSpacing
             write(iFile, '(3f12.4,1P,16e13.5)') &
                  r,densField(indexOrigin(1),i1,i2)
          end do
          write(iFile,*)
       end do

    case (2)
       do i1 = -gridPoints(1),gridPoints(1)
          do i2 = -gridPoints(3),gridPoints(3)
             r = (/i1,indexOrigin(2),i2/) * gridSpacing
             write(iFile, '(3f12.4,1P,16e13.5)') &
                  r,densField(i1,indexOrigin(2),i2)
          end do
          write(iFile,*)
       end do

    case (3)
       do i1 = -gridPoints(1),gridPoints(1)
          do i2 = -gridPoints(2),gridPoints(2)
             r = (/i1,i2,indexOrigin(3)/) * gridSpacing
             write(iFile, '(3f12.4,1P,16e13.5)') &
                  r,densField(i1,i2,indexOrigin(3))
          end do
          write(iFile,*)
       end do

    end select

    close(iFile)
  end subroutine writeDensityPlane

  !****************************************************************************
  !****s* densitymodule/set_pF_Box
  ! NAME
  ! subroutine set_pF_Box(pF)
  !
  ! PURPOSE
  ! set the internal values pF_Box.
  !
  ! These values are only used for adjusting the binning of the momentum
  ! distribution for densitySwitch=4
  !****************************************************************************
  subroutine set_pF_Box(pF)
    real, dimension(0:1), intent(in) :: pF
    pF_Box = pF
  end subroutine set_pF_Box

  !****************************************************************************
  !****f* densitymodule/distMomBox
  ! NAME
  ! real function distMomBox(mom,charge)
  !
  ! PURPOSE
  ! return the value of the distribution function for given momentum and
  ! charge of the nucleon
  !****************************************************************************
  real function distMomBox(mom,charge)

    use hist

    real, intent(in) :: mom
    integer, intent(in) :: charge

    if (initFlag) call init
    if (densitySwitch.ne.4) call TRACEBACK("wrong densitySwitch")

    distMomBox = GetHist(hP_Box(charge), mom, 1)
  end function distMomBox

  !****************************************************************************
  ! cf. interface getFieldRMF
  !****************************************************************************
  subroutine getFieldRMF_i(ix,iy,iz, sigma,omega,rho)

    integer, intent(in) :: ix,iy,iz
    real, optional, intent(out) :: sigma
    real, dimension(0:3), optional, intent(out) :: omega,rho

    if (present(sigma)) sigma = 0.
    if (present(omega)) omega(:) = 0.
    if (present(rho))   rho(:) = 0.

    if (present(sigma)) then
       if (allocated(sigmaField)) then
          sigma = sigmaField(ix,iy,iz)
       end if
    end if
    if (present(omega)) then
       if (allocated(omegaField)) then
          omega(:) = omegaField(ix,iy,iz,:)
       end if
    end if
    if (present(rho)) then
       if (allocated(rhoField)) then
          rho(:) = rhoField(ix,iy,iz,:)
       end if
    end if

  end subroutine getFieldRMF_i
  !----------------------------------------------------------------------------
  subroutine getFieldRMF_p(pos, sigma,omega,rho)
    use constants, only: f_pi

    real, dimension(1:3), intent(in) :: pos
    real, optional, intent(out) :: sigma
    real, dimension(0:3), optional, intent(out) :: omega,rho

    integer, dimension(1:3) :: ind

    if (present(sigma)) sigma = 0.
    if (present(omega)) omega(:) = 0.
    if (present(rho))   rho(:) = 0.

    if (getGridIndex1(pos,ind)) then
       if (present(sigma)) then
          if (allocated(sigmaField)) then
             sigma = sigmaField(ind(1),ind(2),ind(3))
          end if
       end if
       if (present(omega)) then
          if (allocated(omegaField)) then
             omega(:) = omegaField(ind(1),ind(2),ind(3),:)
          end if
       end if
       if (present(rho)) then
          if (allocated(rhoField)) then
             rho(:) = rhoField(ind(1),ind(2),ind(3),:)
          end if
       end if
    else
       if (present(sigma)) then
          if (allocated(sigmaField)) then
             sigma = f_pi
          end if
       end if
    end if

  end subroutine getFieldRMF_p
  !----------------------------------------------------------------------------

end module densitymodule
