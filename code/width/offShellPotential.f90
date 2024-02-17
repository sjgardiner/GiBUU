!******************************************************************************
!****m* /offShellPotential
! NAME
! module offShellPotential
! PURPOSE
! This module calculates the offshell potential.
!******************************************************************************
module offShellPotential

  use CallStack, only: TRACEBACK

  implicit none
  private

  !****************************************************************************
  !****g* offShellPotential/OffShell_debug
  ! PURPOSE
  ! Switch on or off whether the debug information shall be given.
  ! SOURCE
  logical, parameter :: OffShell_debug=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/useOffShellPotentialBaryons
  ! PURPOSE
  ! Switch on or off whether the offshellness should be used for baryons.
  ! SOURCE
  logical, save :: useOffShellPotentialBaryons=.false.
  !
  ! NOTES
  ! * must be set to "TRUE" if mediumSwitch_coll
  !   (see module BaryonWidthMedium) is .true.
  ! * if .true. then delta_T (see module inputGeneral) must be <=0.05
  !   AND delta_P (see module propagation)  must be <=0.002;
  !   AND delta_E (see module propagation)  must be <=0.002;
  !   slows down propagation by a factor of 10
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/useOffShellPotentialMesons
  ! PURPOSE
  ! Switch on or off whether the offshellness should be used for mesons.
  ! SOURCE
  logical, save :: useOffShellPotentialMesons=.false.
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/offshell_cutoff
  ! PURPOSE
  ! If abs(offshellParameter) is less than this value, then we treat
  ! it as zero -> getoffshellMass returns the pole mass!
  ! SOURCE
  real, parameter :: offshell_cutoff = 1E-5
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/max_offshellparameter
  ! PURPOSE
  ! The maximal value for the offshell parameter. Note: empirical value!
  ! This only applies to baryons. For mesons we have no restrictions
  ! on the offshell parameter.
  ! SOURCE
  real, save :: max_offshellparameter=5.
  !****************************************************************************

  !****************************************************************************
  !****g* offShellPotential/extrapolateBaryonWidth
  ! PURPOSE
  ! Whether to extrapolate the baryon width below minimal mass or not.
  ! SOURCE
  logical, save :: extrapolateBaryonWidth=.true.
  !****************************************************************************


  !****************************************************************************
  !****g* offShellPotential/relativistic
  ! PURPOSE
  ! * false: Use non-rel. off-shell parameter x=Delta m/Gamma,
  !   which obeys Stefan Leupold's non-rel. EOM.
  ! * true: Use rel. off-shell parameter x=Delta m^2/Gamma,
  !   which obeys Cassing's rel. EOM.
  ! SOURCE
  logical, save :: relativistic = .false.
  !****************************************************************************


  !****************************************************************************
  !****g* offShellPotential/SetOffShellEnergyFlag
  ! PURPOSE
  ! * false: the energy of off-shell particle is constant during time evolution
  !          (static nucleus)
  ! * true: the energy of off-shell particle varies during time evolution
  !          (dynamic case, e.g. heavy ion collision)
  ! SOURCE
  logical, save :: SetOffShellEnergyFlag = .false.
  !****************************************************************************


  logical, save :: initFlag    =.true.

  public :: get_useOffShellPotentialBaryons
  public :: get_useOffShellPotentialMesons
  public :: treatParticleOffShell
  public :: HamiltonFunc_offshell
  public :: get_offshell_debug
  public :: setOffShellParameter
  public :: setOffShellEnergy



  !****************************************************************************
  !****s* offShellPotential/setOffShellParameter
  ! NAME
  ! Interface setOffShellParameter
  ! PURPOSE
  ! This is an interface which can be called like a subroutine using
  ! "call setOffShellParameter(p,success)".
  !
  ! It calculates the offshell parameter for a single particle or a set of
  ! particles.
  ! If it fails for one of them, then it returns .false., otherwise .true. .
  ! The particles should be properly initialized, only the offshellparameter
  ! should be left-over to initialize.
  !
  ! INPUTS
  ! * type(particle), intent(inout), dimension (:) :: p
  ! or:
  ! * type(particle), intent(inout) :: p
  !
  ! OUTPUT
  ! * logical :: success -- .true. if the offshell parameter
  !   could be evaluated for all particles
  !****************************************************************************
  Interface setOffShellParameter
     Module Procedure setOffShellParameter_2dim, &
          setOffShellParameter_1dim, &
          setOffShellParameter_0dim
  end Interface



contains

  !****************************************************************************
  !****s* offShellPotential/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "OffShellPotential".
  !****************************************************************************
  subroutine readInput
    use output
    use baryonWidthMedium, only: get_MediumSwitch_coll
    use mesonWidthMedium, only: get_MediumSwitchMesons
    use inputGeneral
    use eventtypes


    integer :: ios

    !**************************************************************************
    !****n* offShellPotential/OffShellPotential
    ! NAME
    ! NAMELIST /OffShellPotential/
    ! PURPOSE
    ! Includes the switches:
    ! * useOffShellPotentialBaryons
    ! * useOffShellPotentialMesons
    ! * extrapolateBaryonWidth
    ! * max_offshellparameter
    ! * relativistic
    ! * SetOffShellEnergyFlag
    !**************************************************************************
    NAMELIST /offShellPotential/ &
         useOffShellPotentialBaryons, useOffShellPotentialMesons, &
         extrapolateBaryonWidth, max_offshellparameter, relativistic, &
         SetOffShellEnergyFlag

    call Write_ReadingInput('offShellPotential',0)
    rewind(5)
    read(5,nml=offShellPotential,IOSTAT=IOS)
    call Write_ReadingInput('offShellPotential',0,IOS)

    if (.not.get_MediumSwitch_coll() .and. useOffShellPotentialBaryons) then
       write(*,*) 'MediumSwitch_coll of baryons is switched off, therefore we now switch off useOffShellPotentialBaryons'
       useOffShellPotentialBaryons = .false.
    else if (get_MediumSwitch_coll() .and. .not.useOffShellPotentialBaryons) then
       call TRACEBACK('MediumSwitch_coll of baryons is switched on, therefore useOffShellPotentialBaryons must be set to TRUE!')
    end if

    if (useOffShellPotentialBaryons .and. .not.LRF_equals_CALC_frame) then
       call TRACEBACK('Baryonic off-shell potential works only in the "LRF=Calculation frame" assumption; The more general case is not yet implemented!!')
    end if

    write(*,*) 'use offShell potential for baryons:', useOffShellPotentialBaryons
    write(*,*) 'use offShell potential for mesons: ', useOffShellPotentialMesons
    write(*,*) 'extrapolate baryon width below minimal mass: ', extrapolateBaryonWidth
    write(*,*) 'max. offshell parameter for baryons:  ', max_offshellparameter
    write(*,*) 'use relativistic off-shell parameter: ', relativistic
    write(*,*) 'set off-shell energy: ', SetOffShellEnergyFlag

    call Write_ReadingInput('offShellPotential',1)

  end subroutine readInput


  !****************************************************************************
  !****f* offShellPotential/get_useOffShellPotentialBaryons()
  ! NAME
  ! logical function get_useOffShellPotentialBaryons()
  ! PURPOSE
  ! returns the value of useOffShellPotentialBaryons
  !****************************************************************************
  logical function get_useOffShellPotentialBaryons()
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    get_useOffShellPotentialBaryons=useOffShellPotentialBaryons
  end function get_useOffShellPotentialBaryons


  !****************************************************************************
  !****f* offShellPotential/get_useOffShellPotentialMesons()
  ! NAME
  ! logical function get_useOffShellPotentialMesons()
  ! PURPOSE
  ! returns the value of useOffShellPotentialMesons
  !****************************************************************************
  logical function get_useOffShellPotentialMesons()
    if (initFlag) then
       call readInput
       initFlag=.false.
    end if
    get_useOffShellPotentialMesons=useOffShellPotentialMesons
  end function get_useOffShellPotentialMesons


  !****************************************************************************
  !****f* offShellPotential/get_OffShell_debug()
  ! NAME
  ! logical function get_OffShell_debug()
  ! PURPOSE
  ! returns the value of OffShell_debug (which is a parameter)
  !****************************************************************************
  pure logical function get_OffShell_debug()
    get_OffShell_debug=OffShell_debug
  end function get_OffShell_debug


  !****************************************************************************
  !****f* offShellPotential/treatParticleOffShell(partID,partOffShellParameter)
  ! NAME
  ! logical function treatParticleOffShell(partID,partOffShellParameter)
  ! PURPOSE
  ! returns for a given particle, whether it is treated as
  ! offShell particle (e.g. in propagation).
  !****************************************************************************
  logical function treatParticleOffShell(partID,partOffShellParameter)
    use IdTable, only: isBaryon, isMeson
    integer, intent(in) :: partID
    real, intent(in) :: partoffShellParameter

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    treatParticleOffShell=.false.

    if (useOffShellPotentialBaryons) then
       if(isBaryon(partId).and.abs(partOffshellParameter).gt.offshell_cutoff) &
            treatParticleOffShell=.true.
    end if

    if (useOffShellPotentialMesons) then
       if(isMeson(partId).and.abs(partOffshellParameter).gt.offshell_cutoff) &
            treatParticleOffShell=.true.
    end if

  end function treatParticleOffShell


  !****************************************************************************
  !****f* offShellPotential/HamiltonFunc_offshell
  ! NAME
  ! real function HamiltonFunc_offshell(part,outOfBounds,massIsDetermined,full_offShell,flagOk)
  ! PURPOSE
  ! determines Hamilton function for the offshell potential prescription
  ! INPUTS
  ! * type(particle) :: part -- the particle
  ! * logical, OPTIONAL :: massIsDetermined -- ???
  ! * logical, OPTIONAL :: full_offShell    -- ???
  ! OUTPUT
  ! * function value
  ! * logical :: outOfBounds -- .true. if the width table is out of bounds
  ! * logical, OPTIONAL :: flagOk -- .true. if success
  !
  ! NOTES
  ! in the case of an internal failure (e.g. 'negative mass in massDet'),
  ! the return value is 99999.
  !****************************************************************************
  real function HamiltonFunc_offshell(part_in,outOfBounds,massIsDetermined_in,full_offShell_in,flagOk)
    use particleDefinition
    use potentialMain, only: potential_LRF,massDetermination
    use minkowski, only: abs4
    use output, only: writeParticle_debug, DoPR
    use IdTable, only: isMeson,rho,omegaMeson,phi
    use selfenergy_mesons, only: get_realPart
    use mediumDefinition
    use mediumModule, only: mediumAt
    use particleProperties, only: partName

    type(particle),intent(in) :: part_in
    logical, intent(out) :: outOfBounds
    logical, intent(out), optional :: flagOk
    logical, intent(in), optional :: massIsDetermined_in, full_offShell_in

    type(particle) :: part
    real :: mass,rp
    logical :: success
    logical :: massIsDetermined,full_offShell
    type(medium) :: med

    if (present(massIsDetermined_in)) then
       massIsDetermined=massIsDetermined_in
    else
       massIsDetermined=.false.
    end if
    if (present(full_offshell_in)) then
       full_offshell=full_offshell_in
    else
       full_offshell=.true.
    end if

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    part=part_in

    if (present(flagOk)) flagOk = .true.

    if (.not.treatParticleOffShell(part%ID,part%offshellPar)) then
       call offShellErrorMessage(part)
       call TraceBack('HAMILTONFUNC_offshell should not be called -> STOP')
    end if


    if (full_offshell) then
       !Offshell mass includes off shell potential!
       if (.not.massIsDetermined) then
          call massDetermination(part,success=success,verbose=.false.)
          if (.not. success) then
             hamiltonFunc_offshell=99999.
             if (present(flagOk)) flagOk = .false.
             if (DoPR(2)) write(*,'("WARNING in HamiltonFunc: ",A,A,i10)') &
                  'No success in massDet: ',partName(part),part%number
             ! call writeParticle_debug(part)
             return
          end if
       end if
       med = mediumAt(part%pos)
       med%useMedium = .true. ! to avoid threshold effects
       if (isMeson(part%ID)) then
          rp = get_realPart(part%ID, abs4(part%mom), med)
       else
          rp = 0.
       end if
       mass=getOffShellMass(part%ID,part%offshellPar,part%mom,part%mass,med,outOfBounds,success=success)
       if (.not. success) then
          hamiltonFunc_offshell=99999.
          if (present(flagOk)) flagOk = .false.
          if (DoPR(2)) write(*,'("WARNING in HamiltonFunc: ",A,A,i10)') &
                  'getOffShellMass: ',partName(part),part%number
          ! call writeParticle_debug(part)
          return
       end if
       select case (part%ID)
       case (rho,omegaMeson,phi)
          ! Optinal in-medium mass shift of vector mesons is contained in 'mass'.
          hamiltonFunc_offshell= sqrt(Dot_Product(part%mom(1:3),part%mom(1:3))+mass**2 + rp)
       case default
          hamiltonFunc_offshell= sqrt(Dot_Product(part%mom(1:3),part%mom(1:3))+mass**2 + rp) + potential_LRF(part)
       end select
    else
       outOfBounds=.false.
       ! no offshell potential
       hamiltonFunc_offshell= sqrt(Dot_Product(part%mom(1:3),part%mom(1:3))+part%mass**2) + potential_LRF(part)
    end if


  end function HamiltonFunc_offshell



  !****************************************************************************
  ! cf. Interface setOffShellParameter
  !****************************************************************************
  subroutine setOffShellParameter_0dim(p,success)
    use particleDefinition

    type(particle), intent(inout)  :: p
    logical, intent(out):: success

    p%offshellPar=getOffShellParameter(p%ID,p%mass,p%mom,p%pos,success)

  end subroutine setOffShellParameter_0dim

  !****************************************************************************
  ! cf. Interface setOffShellParameter
  !****************************************************************************
  subroutine setOffShellParameter_1dim(p,success)
    use particleDefinition

    type(particle), intent(inout), dimension (:) :: p
    logical, intent(out):: success
    logical :: flagOK
    integer :: i

    success=.true.
    do i= lbound(p,dim=1),ubound(p,dim=1)
       if (p(i)%ID.le.0) cycle
       call setOffShellParameter_0dim(p(i),flagOk)
       if (.not.flagOk) then
          success=.false.
       end if
    end do
  end subroutine setOffShellParameter_1dim

  !****************************************************************************
  ! cf. Interface setOffShellParameter
  !****************************************************************************
  subroutine setOffShellParameter_2dim(p,success)
    use particleDefinition

    type(particle), intent(inout), dimension (:,:) :: p
    logical, intent(out):: success
    logical :: flagOK
    integer :: iPart,iEns,nPart,nEns

    success=.true.

    nEns = size(p,dim=1)
    nPart = size(p,dim=2)
    do iEns = 1,nEns
       do iPart = 1,nPart
          if (p(iEns,iPart)%ID.le.0) cycle
          call setOffShellParameter_0dim(p(iEns,iPart),flagOk)
          if (.not.flagOk) then
             success=.false.
          end if
       end do
    end do
  end subroutine setOffShellParameter_2dim

  !****************************************************************************
  !****f* offShellPotential/getOffShellParameter
  ! NAME
  ! real function getOffShellParameter(partID,bareMass,momentum,position,
  ! success)
  !
  ! PURPOSE
  ! This function calculates the offshell parameter to be set into
  ! particle%offshellPar.
  ! When calling this routine, make sure, that momentum(0) is set correctly!
  !
  ! INPUTS
  ! * integer :: partID    -- ID of particle
  ! * real :: bareMass     -- bare mass= %mass
  ! * real, dimension(0:3) :: momentum  -- 4-momentum at production point
  ! * real, dimension(1:3) :: position  -- production point
  !
  ! OUTPUT
  ! * function value
  ! * logical :: success -- .true. if offshell parameter could be evaluated
  !****************************************************************************
  real function getOffShellParameter(partID,bareMass,momentum,position,success)
    use IDTable, only: rho, omegaMeson, phi, isBaryon
    use particleProperties, only: hadron
    use mediumDefinition
    use hist2D
    use mediumModule, only: mediumAt
    use minkowski, only: abs4,abs3
    use mesonWidthMedium, only: WidthMesonMedium

    integer, intent(in) :: partID
    real, intent(in) :: bareMass
    real, dimension(0:3), intent(in) :: momentum
    real, dimension(1:3), intent(in) :: position
    logical, intent(out):: success

    real :: width, poleMass
    type(medium) :: mediumAtPosition

    !for debugging
    integer,save :: numcalls=0
    type(histogram2D),save :: histo_nucl
    type(histogram2D),save :: histo_delta,histo_delta_rho
    type(histogram2D),save :: histo_rho
    logical, save :: initHistFirst=.true.
    logical :: outofBounds
    integer, save :: counter_too_large=0

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    numcalls=numcalls+1

    !set standard output
    getOffShellParameter=0.
    success=.false.

    if (partID.le.0) return

    mediumAtPosition = mediumAt(position)
    mediumAtPosition%useMedium = .true. ! to avoid threshold effects

    !BARYONS --------------------------------------------------
    if (useOffShellPotentialBaryons.and.isBaryon(partId)) then
       width=getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,outofBounds)
    !MESONS --------------------------------------------------
    else if (useOffShellPotentialMesons.and.(partID==rho .or. partID==omegaMeson .or. partID==phi)) then
       width = getMesonWidth(partID,baremass,momentum,mediumAtPosition)
    else
       success=.true.
       return
    end if

    poleMass=hadron(partID)%mass

    if (width.gt.0) then
       if (relativistic) then
         getOffShellParameter=(bareMass**2-poleMass**2)/(2.*abs4(momentum)*width)
       else
         getOffShellParameter=(bareMass-poleMass)/width
       end if
    else
       getOffShellParameter=0.
       success=.true.
       return
    end if


!    if (partID==103) then
!       write(742,*) bareMass,width,getOffShellParameter,abs4(momentum),abs3(momentum),mediumAtPosition%densityNeutron,mediumAtPosition%densityProton
!       flush(742)
!    end if

    !cutoff: we do not allow particles which are too far offshell
    if (isBaryon(partId) .and. abs(getoffshellparameter)>max_offshellparameter) then
       ! Declare failure: particle is too far off-shell
       if (counter_too_large.lt.100) then
          write(*,'(A,2G13.5)') 'WARNING: off-shell parameter too large!', &
               getoffshellparameter, max_offshellparameter
          write(*,'(12g13.5)') partID,bareMass,momentum,poleMass,width,&
               mediumAtPosition%density
          counter_too_large=counter_too_large+1
          if (counter_too_large==100) &
               write(*,*) 'This happended now 100 times. Stopping output!'
       end if
       getoffshellparameter=0.
       success=.false.
    else
       success=.true.
    end if

!    write(*,*)' Id, Mass, OSP : ', partID, bareMass, getOffShellParameter

    !for debugging:
    if (offshell_debug) then
       if (initHistFirst) then
          call CreateHist2D(histo_nucl, 'off shell parameter',&
               (/-20.,0.6/),(/20.,1.3/),(/0.05,0.01/))
          call CreateHist2D(histo_delta, 'off shell parameter',&
               (/-20.,0.9/),(/20.,1.6/),(/0.05,0.01/))
          call CreateHist2D(histo_delta_rho, 'off shell parameter',&
               (/-20.,0./),(/20.,0.2/),(/0.05,0.005/))
          call CreateHist2D(histo_rho, 'off shell parameter',&
               (/-20.,0./),(/20.,1.0/),(/0.05,0.01/))
          initHistFirst=.false.
       end if

       select case (partID)
       case (1)
          call AddHist2D(histo_nucl, (/getoffshellparameter,width/),1.)
       case (2)
          call AddHist2D(histo_delta, (/getoffshellparameter,width/),1.)
          call AddHist2D(histo_delta_rho, &
               (/getoffshellparameter,mediumAtPosition%density/),1.)
       case (103)
          call AddHist2D(histo_rho, (/getoffshellparameter,width/),1.)
       end select

       if (mod(numCalls,1000).eq.0) then
          call WriteHist2D_Gnuplot(histo_nucl,&
               file='offshell_params_nucl.dat',&
               mul=histo_nucl%xBin(1)*histo_nucl%xBin(2))
          call WriteHist2D_Gnuplot(histo_delta,&
               file='offshell_params_delta.dat',&
               mul=histo_delta%xBin(1)*histo_delta%xBin(2))
          call WriteHist2D_Gnuplot(histo_delta_rho,&
               file='offshell_params_delta_rho.dat',&
               mul=histo_delta_rho%xBin(1)*histo_delta_rho%xBin(2))
          call WriteHist2D_Gnuplot(histo_rho,&
               file='offshell_params_rho.dat',&
               mul=histo_rho%xBin(1)*histo_rho%xBin(2))
       end if
    end if


  end function getOffShellParameter


  !****************************************************************************
  !****f* offShellPotential/getOffShellMass
  ! NAME
  ! real function getOffShellMass(partID,offshellparameter,momentum,baremass,
  ! mediumAtPosition,outOfBounds,success)
  !
  ! FUNCTION
  ! This function calculates the offshellmass of a particle depending on its
  ! offshellparameter and the current kinematics
  ! -> full four-momentum needed!!!
  !
  ! INPUTS
  ! * integer :: partID                 -- ID of particle
  ! * real :: offshellparameter         --
  ! * real :: bareMass                  -- bare mass= %mass
  ! * real, dimension(0:3) :: momentum  -- 4-momentum
  ! * type(medium) :: mediumAtPosition  -- medium information
  !
  ! OUTPUT
  ! * real :: getOffShellMass
  ! * logical :: outOfBounds -- true if mass is not in grid
  ! * logical, OPTIONAL :: success -- indicates failure
  !****************************************************************************
  real function getOffShellMass(partID,offshellparameter,momentum,baremass,mediumAtPosition,outOfBounds,success)
    use particleDefinition
    use IdTable, only: isBaryon,isMeson,rho,omegaMeson,phi
    use particleProperties, only: hadron
    use baryonWidthMedium_tables, only: get_minMass
    use mediumDefinition
    use dichteDefinition
    use minkowski, only: abs4, SP
!    use mesonWidthMedium, only: WidthMesonMedium
    use mesonPotentialMain, only: vecMes_massShift

    integer, intent(in) :: partID
    real, intent(in) :: offshellparameter
    real, dimension(0:3),intent(in) :: momentum
    real, intent(in) :: baremass
    type(medium),intent(in) :: mediumAtPosition
    logical, optional, intent(out) :: success

    logical, intent(out) :: outOfBounds
    real :: width,spot

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    outOfBounds=.false.
    if (present(success)) success=.true.
    if (treatParticleOffShell(partID,OffShellParameter)) then

       ! Baryons
       if (isBaryon(partId)) then

          if (SP(momentum,momentum).le.0.) then
             call errorMessage_SPleZero()
             call traceback()
          end if
          width = getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,outOfBounds)
          if (relativistic) then
            getOffShellMass=sqrt(max(hadron(partID)%mass**2+offshellparameter*2.*abs4(momentum)*width,0.))
          else
            getOffShellMass=hadron(partID)%mass+offshellparameter*width
          end if
          if (.not.extrapolateBaryonWidth) then
             getOffShellMass=max(getOffShellMass,get_minMass(partID))
             if (outOfbounds) getOffShellMass=get_minMass(partID)
          end if

       ! Mesons
       else if (isMeson(partId)) then

          if (SP(momentum,momentum).le.0.) then
             if (present(success)) then
                success=.false.
             else
                call errorMessage_SPleZero()
             end if
             outOfBounds = .true.
             getOffShellMass = 0.
             return
          end if

          width = getMesonWidth(partID,baremass,momentum,mediumAtPosition)

          select case (partId)
          case (rho,omegaMeson,phi)
             ! set scalar potential for vector mesons (mass shift!)
             spot = vecMes_massShift(partId,mediumAtPosition%density)
          case default
             ! in all other cases: neglect potential
             spot = 0.
          end select

          if (relativistic) then
            getOffShellMass=sqrt(max((hadron(partID)%mass+spot)**2+offshellparameter*2.*abs4(momentum)*width,0.))
          else
            getOffShellMass=max(hadron(partID)%mass+spot+offshellparameter*width,0.)
          end if

       else
          write(*,*) 'Strange ID!! -> STOP', partID
          call TraceBack()
       end if

    else if (partID/=0) then
       write(*,*) 'CRITICAL ERROR: SHOULD NOT BE CALLED!! -> STOP',  &
            & partID,OffShellParameter
       call TraceBack()
    else
       write(*,*) 'WARNING: getOffShellMass: ID=0'
       getOffShellMass=0
    end if

  contains

    subroutine errorMessage_SPleZero()
      write(*,*) 'p^mu p_mu is less or equal zero in getOffShellMass:', SP(momentum,momentum)
      write(*,*) '  partID,bareMass:',partID,bareMass
!      write(*,*) 'momentum:'
      write(*,*) momentum
      write(*,*) '  offshellparameter:',offshellparameter
      write(*,*) '  density:',mediumAtPosition%density
      write(*,*)
    end subroutine errorMessage_SPleZero


  end function getOffShellMass


  !****************************************************************************
  !****s* offShellPotential/SetOffShellEnergy
  ! NAME
  ! subroutine SetOffShellEnergy(part,errCode)
  !
  ! PURPOSE
  ! Calculate the offshell energy and offshell mass of a particle depending on
  ! its offshellparameter and three-momentum.
  !
  ! INPUT
  ! * type(particle) :: part
  ! OUTPUT
  ! * type(particle) :: part -- with redefined part%mom(0) and part%mass
  ! * integer, OPTIONAL :: errCode -- 0 if okay, otherwise = 1,...5
  !****************************************************************************
  subroutine SetOffShellEnergy(part,errCode)
    use particleDefinition
    use IdTable, only: isBaryon,rho,omegaMeson,phi
    use particleProperties, only: hadron, partName
    use mediumDefinition
    use mediumModule, only: mediumAt
    use mesonWidthMedium, only: WidthMesonMedium
    use mesonPotentialMain, only: vecMes_massShift
    use constants, only: melec
    use output, only: WriteParticle, DoPr

    type(particle), intent(inout) :: part
    integer, intent(out), optional :: errCode

    real, parameter :: err=1.e-06, dE=1.e-04, mass2min=1.e-04

    real, dimension(0:3) :: momentum
    type(medium) :: mediumAtPosition
    real :: spot,pAbs2,massOld,Eold,PoleMass,PoleMass2,E,E0,Enew,derf,fE,f0,fnew,mass2
    integer :: iter
    logical :: flagTryDiv
    integer, parameter :: iterMax = 20  ! default: 10

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    if (present(errCode)) errCode = 0

    if(.not.SetOffShellEnergyFlag) return

    if(part%id.le.0) return

    if (.not.treatParticleOffShell(part%id,part%offshellPar)) return

    mediumAtPosition=mediumAt(part%pos)
    mediumAtPosition%useMedium = .true. ! to avoid threshold effects

    select case (part%id)
    case (rho,omegaMeson,phi)
       ! set scalar potential for vector mesons (mass shift!)
       spot = vecMes_massShift(part%id,mediumAtPosition%density)
    case default
       ! in all other cases: neglect potential
       spot = 0.
    end select

    momentum=part%mom

    pAbs2=momentum(1)**2+momentum(2)**2+momentum(3)**2

    massOld=part%mass
    Eold=sqrt((massOld+spot)**2+pAbs2)

    PoleMass=hadron(part%id)%mass + spot
    PoleMass2=PoleMass**2

    if(momentum(0).lt.0.) then
       write(*,*) ' Particle with negative energy: ', &
            part%Id, part%number, momentum
       call TraceBack()
    end if

    E=momentum(0)

    fE=f(E)

    if(abs(fE).lt.err) then
       if(mass2.gt.mass2min) then
          part%mass = sqrt(mass2) - spot
       else
          if (DoPR(2)) write(*,'(A,A,i10,G12.3)') &
               "WARNING in SetOffShellEnergy (1): too small (inv.mass)^2: ",&
               partName(part), part%number, mass2
          part%mass=massOld
          part%mom(0)=Eold
          if (present(errCode)) errCode = 1
       end if
       return
    end if

    flagTryDiv=.false.

    iter=0
    do

       iter=iter+1

       derf=(f(E+dE)-fE)/dE

       if(derf.ne.0.) then
          E0=E
          f0=fE
          E = E - fE/derf
       else  ! just shift from the extremum
          E = E + dE
       end if

       E=abs(E)

       fE=f(E)

!       write(*,*)' iter-secant, E, fE : ', iter, E, fE

       if(abs(fE).lt.err) then
          if(mass2.gt.mass2min) then
             part%mass = sqrt(mass2) - spot
             part%mom(0)=E
          else
             if (DoPR(2)) write(*,'(A,A,i10,G12.3)') &
                  "WARNING in SetOffShellEnergy (2): too small (inv.mass)^2: ",&
                  partName(part), part%number, mass2
             part%mass=massOld
             part%mom(0)=Eold
             if (present(errCode)) errCode = 2
          end if
          return
       else if(iter.gt.3 .and. f0*fE.lt.0.) then
          flagTryDiv=.true.
          exit
       else if(iter.eq.iterMax) then
          if (DoPR(2)) write(*,'(A,A,i10,G12.3)') &
               "WARNING in SetOffShellEnergy (3): poor convergence: ",&
               partName(part), part%number, abs(fE)
          part%mass=massOld
          part%mom(0)=Eold
          if (present(errCode)) errCode = 3
          return
       end if

    end do

    if(flagTryDiv) then
       ! Use the method of division by two
       iter=0
       do
          iter=iter+1

          Enew = 0.5*(E+E0)
          fnew=f(Enew)
          if(fnew*fE.lt.0.) then
             E0=Enew
          else
             E=Enew
             fE=fnew
          end if
!          write(*,*)' iter-division, E, fE : ', iter, Enew, fnew
          if(abs(fnew).lt.err) then
              if(mass2.gt.mass2min) then
                 part%mass = sqrt(mass2) - spot
                 part%mom(0)=Enew
              else
                 if (DoPR(2)) write(*,'(A,A,i10,G12.3)') &
                      "WARNING in SetOffShellEnergy (4): too small (inv.mass)^2: ",&
                      partName(part), part%number, mass2
                 part%mass=massOld
                 part%mom(0)=Eold
                 if (present(errCode)) errCode = 4
              end if
              return
           else if(iter.eq.100) then
              if (DoPR(2)) write(*,'(A,A,i10,G12.3)') &
                   "WARNING in SetOffShellEnergy (5): poor convergence: ",&
                   partName(part), part%number, abs(fnew)
              part%mass=massOld
              part%mom(0)=Eold
              if (present(errCode)) errCode = 5
              return
          end if
       end do
    end if

  contains

    real function f(p0)

      real, intent(in) :: p0    ! particle energy (GeV)

      real :: mass,width
      logical :: outofBounds

      momentum(0)=p0

      mass2=momentum(0)**2-pAbs2

!      if (mass2.le.0.) then
!         write(*,*) 'space-like particle in offshellPotential/SetOffShellEnergy'
!         call WriteParticle(6,99,1,part)
!         stop
!      end if

      mass=sqrt(max(4.*melec**2,mass2))

      if (isBaryon(part%id)) then
          width = getBaryonWidth(part%id,mass,momentum,mediumAtPosition,outOfBounds)
      else
          width = getMesonWidth(part%id,mass,momentum,mediumAtPosition)
      end if

      if (relativistic) then
         f = mass2 - PoleMass2 - part%offshellPar*2.*mass*width
      else
         f = mass - PoleMass - part%offshellPar*width
      end if

!      write(*,*)' id, mass, density, width : ', part%id, mass, mediumAtposition%density, width
!      write(*,*)' momentum: ', momentum

    end function f

  end subroutine SetOffShellEnergy


  !****************************************************************************
  !****if* offShellPotential/getBaryonWidth
  ! NAME
  ! real function getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,
  ! outofBounds)
  !
  ! PURPOSE
  ! returns the baryonWidth. In case, extrapolateBaryonWidth is set to .true.
  ! the baryon width is extrapolated to masses below minimalmass.
  !
  ! IMPORTANT: *** Baryon Width is calculated in the Lab Frame! ***
  !
  ! INPUTS
  ! * integer :: partID   -- ID of particle
  ! * real :: bareMass     -- bare mass= %mass
  ! * real, dimension(0:3) :: momentum  -- 4-momentum
  ! * type(medium) :: mediumAtPosition   -- medium information
  !
  ! OUTPUT
  ! * function value
  !****************************************************************************
  real function getBaryonWidth(partID,bareMass,momentum,mediumAtPosition,outofBounds)
    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use baryonwidthmedium, only: WidthBaryonMedium
    use mediumDefinition
    use minkowski, only: abs4
    use twoBodyTools, only: pCM
    use constants, only: mN
    use lorentzTrafo, only: lorentz

    integer, intent(in) :: partID
    real, intent(in) :: bareMass
    real, dimension(0:3), intent(in) :: momentum
    type(medium),intent(in) :: mediumAtPosition
    logical,intent(out) :: outofBounds

    real, dimension(0:3) :: momLRF
    real :: width
    real :: x1,x2,y1,y2,m,b,s_cut,s_real
    integer, parameter :: Method=3
    real :: p_squared

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    if (.not.isBaryon(partId)) then
       write(*,*) 'partID = ',partID
       call TRACEBACK('particle is no baryon!')
    end if

    ! Determine momentum in LRF
    momLRF=momentum
    call lorentz(mediumAtPosition%betaLRF,momLRF)


    if (extrapolateBaryonWidth.and.baremass.lt.(hadron(partid)%minmass+0.015)) then
       x1=hadron(partid)%minmass+0.015
       y1=WidthBaryonMedium(partID,x1,momLRF,mediumATposition,outofBounds)
       select case (Method)
       case (1)
          ! tina's method
          x2=hadron(partid)%minmass+0.014
          y2=WidthBaryonMedium(partID,x2,momLRF,mediumATposition,outofBounds)
          m=(y1-y2)/(x1-x2)
          b=y1-m*x1
          width=m*baremass+b
          if (offshell_debug) write(*,*) 'fakewidth=',Width
          width=max(0.001,width)  !minimalwidth as set in baryonWidthMedium_tables
          outofBounds=.false.
       case (2)
          ! Constant width
          width=y1
          getBaryonWidth=width
          outofBounds=.true.
       case (3)
          ! Assume that width is due to NN absorption below threshold
          p_squared=Dot_Product(momentum(1:3),momentum(1:3))
          s_cut =sqrt((mN+sqrt(x1**2      +p_squared) )**2-p_Squared )
          s_real=sqrt((mN+sqrt(baremass**2+p_squared) )**2-p_squared )
          width=y1*pcm(s_real,mN,mN)/pcm(s_cut,mN,mN)
          outofBounds=.false.
       end select
    else
       width=WidthBaryonMedium(partID,baremass,momLRF,mediumATposition,outofBounds)
       if (offshell_debug) write(*,*) 'width=',Width,partID,baremass,momentum,mediumATposition,outofBounds
    end if

    ! Gamma_lab = Gamma_restframe / gamma
    getBaryonWidth = width * abs4(momentum)/momentum(0)

  end function getBaryonWidth



  !****************************************************************************
  !****f* offShellPotential/getMesonWidth
  ! NAME
  ! real function getMesonWidth
  !
  ! PURPOSE
  ! Returns the meson width in the meson rest frame.
  !
  ! INPUTS
  ! * integer :: partID   -- ID of particle
  ! * real :: bareMass     -- bare mass= %mass
  ! * real, dimension(0:3) :: momentum  -- 4-momentum
  ! * type(medium) :: mediumAtPosition   -- medium information
  !
  ! OUTPUT
  ! * function value
  !****************************************************************************
  real function getMesonWidth(partID,bareMass,momentum,mediumAtPosition)
    use mediumDefinition
    use mesonWidthMedium, only: WidthMesonMedium
    use lorentzTrafo, only: lorentz

    integer, intent(in) :: partID
    real, intent(in) :: bareMass
    real, dimension(0:3), intent(in) :: momentum
    type(medium),intent(in) :: mediumAtPosition

    real, dimension(0:3) :: momLRF

    ! Determine momentum in LRF
    momLRF=momentum
    call lorentz(mediumAtPosition%betaLRF,momLRF)

    getMesonWidth = WidthMesonMedium(partID,bareMass,momLRF,mediumAtposition)

  end function getMesonWidth



  !****************************************************************************
  !****s* offShellPotential/offShellErrorMessage
  ! NAME
  ! subroutine offShellErrorMessage(partID,bareMass,momentum,position)
  ! PURPOSE
  ! Routine prints particle information.
  ! INPUTS
  ! * type(particle) :: part
  !****************************************************************************
  subroutine offShellErrorMessage(part)
    use particleDefinition
    use minkowski, only: SP

    type(particle), intent(in) :: part
    character(20) :: form
    form='(A,4G12.5)'

    write(*,'(A)')'offShellErrorMessage: printing particle properties...:'
    write(*,form)'Particle ID:                ', part%ID
    write(*,form)'Particle mass:              ', part%mass
    write(*,form)'Particle offshellparameter: ', part%offshellPar
    write(*,form)'Particle number:            ', part%number
    write(*,form)'Particle firstevent:        ', part%firstevent
    write(*,form)'Particle history:           ', part%history
    write(*,form)'Particle momentum:          ', part%mom
    write(*,form)'Particle perturbative?:     ', part%pert
    write(*,form)'Particle charge:            ', part%Charge
    write(*,form)'Particle position:          ', part%pos
    write(*,form)'SP(momentum,momentum):      ', SP(part%mom,part%mom)
    write(*,form)'SQ(abs(SP(momentum,momentum):',sqrt(abs(SP(part%mom,part%mom)))
    write(*,*)

  end subroutine offShellErrorMessage



end module offShellPotential
