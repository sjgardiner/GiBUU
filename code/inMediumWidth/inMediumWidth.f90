!******************************************************************************
!****m* /inMediumWidth
! NAME
! module inMediumWidth
! PURPOSE
! * Implements the routines for the medium width of the baryons
! * Generates tabulated input files
! * see also code/inMediumWidth/testRun/tabulateImagPart.f90 which
!   creates an executable using this routine to tabulate the width
!******************************************************************************
module inMediumWidth

  implicit none

  private

  !****************************************************************************
  !****g* inMediumWidth/debugFlag
  ! SOURCE
  !
  logical,save  :: debugFlag=.false.
  !
  ! PURPOSE
  ! * Switch for debug information
  !****************************************************************************

  !****************************************************************************
  !****g* inMediumWidth/writeLocal
  ! SOURCE
  !
  logical,save  :: writeLocal=.true.
  !
  ! PURPOSE
  ! * If .true. then output files are written to ./"filename", else to buuinput.
  !****************************************************************************

  !****************************************************************************
  !****g* inMediumWidth/maxRes
  ! SOURCE
  !
  integer,save  :: maxRes=1000
  !
  ! PURPOSE
  ! Read the data table up to a maximum resonance ID. ONLY FOR TESTING!!!
  !****************************************************************************

  !****************************************************************************
  !****g* inMediumWidth/minRes
  ! SOURCE
  !
  integer,save  :: minRes=-1000
  !
  ! PURPOSE
  ! Read the data table starting at this minimal resonance ID.
  ! ONLY FOR TESTING!!!
  !****************************************************************************


  ! flags for mesons
  integer, save :: maxMes = 1000      ! max. meson ID to tabulate
  integer, save :: minMes = -1000     ! min. meson ID to tabulate
  real,    save :: max_absP_Mes = 3.0 ! max. momentum for tabulation
  integer, save :: num_MonteCarlo_Points_mesons = 250 ! number of monte-carlo points
  integer, save :: cmin = 0, cmax = 0 ! bounds of charge loop

  logical, save :: readInput_flag=.true.

  public :: tabulate_inMediumWidth_baryons
  public :: tabulate_inMediumWidth_mesons
  public :: evaluateCollisionBroadening_mesons

contains


  !****************************************************************************
  !****s* inMediumWidth/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "inMediumWidth".
  !****************************************************************************
  subroutine readInput

    use output

    integer :: ios ! checks file behavior

    NAMELIST /inMediumWidth/ debugFlag,maxRes,minRes,writeLocal, &
         maxMes,minMes,max_absP_Mes,num_MonteCarlo_Points_mesons,cmin,cmax

    call Write_ReadingInput('inMediumWidth',0)
    rewind(5)
    read(5,nml=inMediumWidth,IOSTAT=IOS)
    call Write_ReadingInput('inMediumWidth',0,IOS)

    write(*,*) 'Debugging information        ?   ', debugFlag
    write(*,*) 'maxRes                       ?   ', maxRes
    write(*,*) 'minRes                       ?   ', minRes
    write(*,*) 'write locally                ?   ', writeLocal

    write(*,*) 'maxMes                       ?   ', maxMes
    write(*,*) 'minMes                       ?   ', minMes
    write(*,*) 'max_absP_Mes                 ?   ', max_absP_Mes
    write(*,*) 'num_MonteCarlo_Points_mesons ?   ', num_MonteCarlo_Points_mesons
    write(*,*) 'cmin                         ?   ', cmin
    write(*,*) 'cmax                         ?   ', cmax

    call Write_ReadingInput('inMediumWidth',1)
  end subroutine readInput


  !****************************************************************************
  !****s* inMediumWidth/tabulate_inMediumWidth_baryons
  ! NAME
  ! subroutine tabulate_inMediumWidth_baryons()
  !
  ! PURPOSE
  ! * Tabulates the in-medium width of baryons according to
  !   Gamma=Gamma_free*Pauli_Blocking+sigma*rho*v *Pauli_Blocking
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming
  !   particles in the vacuum
  ! * Simplification: An average over the charge of the baryon is performed.
  !
  ! INPUTS
  ! * NONE
  !
  ! OUTPUT
  ! * Files "inMediumWidth/InMediumWidth.particleID.1.dat.bz2" and
  !   "inMediumWidth/InMediumWidth.particleID.2.dat.bz2"
  !****************************************************************************
  subroutine tabulate_inMediumWidth_baryons()
    use idTable, only: nres,nucleon
    use mediumDefinition
    use particleProperties, only: hadron
    use inputGeneral, only: path_To_Input
    use output, only: intToChar,line
    use baryonWidthMedium_tables, only: numSteps_rhoN,numSteps_rhoP, &
         numSteps_absP,numsteps_mass,max_absP,max_rhoP,max_rhoN, &
         num_MonteCarlo_Points, get_min_charge_loop, get_max_charge_loop
    use bzip

    real, save ::  delta_rhoN,delta_rhoP,delta_absP
    real, save, dimension(1:nres+1) :: delta_mass,min_mass,max_mass
    real, save, dimension(:,:,:,:,:) , Allocatable:: widthTable
    ! (:,:,:,:,:,1) = collisional broadening
    ! (:,:,:,:,:,2) = pauli blocked vacuum width

    integer :: ID,iAbsP,iMass, iRhoN,iRhoP,dummy,i
    integer, parameter :: showColl=1,showPauli=2
    character(20) :: format
    character(1000) :: fileName
    real :: rhoN, rhoP,absP, mass
    type(bzFile) :: f
    character(len=100) :: buf

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.

       write(*,*) 'tabulation parameters:'
       write(*,*) 'numSteps_rhoN= ',numSteps_rhoN
       write(*,*) 'numSteps_rhoP= ',numSteps_rhoP
       write(*,*) 'numSteps_absP= ',numSteps_absP
       write(*,*) 'numSteps_mass= ',numSteps_mass
       write(*,*) 'max_absP= ',max_absP
       write(*,*) 'max_rhoP= ',max_rhoP
       write(*,*) 'max_rhoN= ',max_rhoN
       write(*,*) 'num_MonteCarlo_Points= ',num_MonteCarlo_Points
       ! The following construction is necessary to avoid a recursive I/O operation:
       dummy=get_min_charge_loop()
       write(*,*) 'min_charge_loop= ',dummy
       dummy=get_max_charge_loop()
       write(*,*) 'max_charge_loop= ',dummy
    end if

    allocate(widthTable(0:numSteps_absP,0:numSteps_mass,0:numSteps_rhoN,0:numSteps_rhoP,1:2))
    delta_rhoN=max_rhoN/float(numSteps_rhoN)
    delta_rhoP=max_rhoP/float(numSteps_rhoP)
    delta_absP=max_absP/float(numSteps_absP)

    format='(' // intToChar(numSteps_rhoP+1) // 'E12.4)'

    ! Set the bounds for the parametrization of the mass
    do ID=1,nres+1
       !         min_mass(ID)=max(minimalMass(ID)+0.02,hadron(ID)%mass-max(0.1,hadron(ID)%width*3.))
       !         max_mass(ID)=hadron(ID)%mass+max(0.1,hadron(ID)%width*3.)
       if (ID.eq.nucleon) then
          min_mass(ID)= 0.700 ! 0.80
       else
          min_mass(ID)=hadron(ID)%minmass+0.01
       end if
       max_mass(ID)= 4. ! 7.
       delta_mass(ID)=(max_mass(ID)-min_mass(ID))/float(numSteps_mass)
    end do

    write(*,*) 'Tabulating IN-MEDIUM WIDTH table for BARYONS'
    do ID=max(minRes,1),min(maxRes,nres+1)
       do iAbsP=0,numSteps_absP
          absP=float(iAbsP)*delta_absP
          write(*,*) (ID-1)*(numSteps_absP+1)+iAbsP,'/',(nres+1)*(numSteps_absP+1)
          write(*,'(A,F15.4,2(A,I8))') 'absP=', absP,' steps:',iAbsP,'/',numsteps_absP
          do iMass=0,numSteps_mass
             mass=min_mass(ID)+float(iMass)*delta_mass(ID)
             write(*,'(A,F15.4,2(A,I8))') 'mass=', mass,' steps:',iMass,'/',numsteps_mass
             do iRhoN=0,numSteps_rhoN
                rhoN=float(iRhoN)*delta_rhoN
                do iRhoP=0,numSteps_rhoP
                   if (mass.lt.hadron(ID)%minmass) then
                      widthTable(iAbsP,iMass,iRhoN,iRhoP,showPauli)= 0.
                      widthTable(iAbsP,iMass,iRhoN,iRhoP,showColl)=  0.
                   else
                      rhoP=float(iRhoP)*delta_rhoP
                      widthTable(iAbsP,iMass,iRhoN,iRhoP,showPauli)=&
                           & get_pauliBlockedDecayWidth(ID,absP,mass,rhoN,rhoP)
                      widthTable(iAbsP,iMass,iRhoN,iRhoP,showColl)=&
                           & evaluateCollisionBroadening_baryons(ID,absP,mass,rhoN,rhoP)
                   end if
                end do
             end do
          end do
       end do

       ! Write the table to files:
       do i=1,2
          if (writelocal) then
             fileName='./InMediumWidth.'//trim(intToChar(ID))&
                  & //'_'//trim(intToChar(i))//'.dat.bz2'
          else
             fileName=trim(path_to_Input)//'/inMediumWidth/InMediumWidth.'//trim(intToChar(ID))&
                  & //'_'//trim(intToChar(i))//'.dat.bz2'
          end if
          f = bzOpenW(trim(fileName))
          write(buf,'(A,7I6)')'#', numSteps_absP, numSteps_mass , numSteps_rhoN, numSteps_rhoP,num_MonteCarlo_Points,&
               & get_min_charge_loop(),get_max_charge_loop()
          call bzWriteLine(f,buf)
          write(buf,'(A,5E15.4)')'#', max_absP, max_rhoP , max_rhoN,min_mass(ID),max_mass(ID)
          call bzWriteLine(f,buf)
          do iAbsP=0,numSteps_absP
             do iMass=0,numSteps_mass
                do iRhoN=0,numSteps_rhoN
                   write(buf,format) widthTable(iAbsP,iMass,iRhoN,0:numSteps_rhoP,i)
                   call bzWriteLine(f,buf)
                end do
             end do
          end do
          call bzCloseW(f)
       end do
       write(*,*) '... finished ID=',ID
    end do
    write(*,*) '... finished tabulating IN-MEDIUM WIDTH table for BARYONS'
    write(*,*) line
    write(*,*)

  end subroutine tabulate_inMediumWidth_baryons



  !****************************************************************************
  !****s* inMediumWidth/tabulate_inMediumWidth_mesons
  ! NAME
  ! subroutine tabulate_inMediumWidth_mesons()
  !
  ! PURPOSE
  ! * Tabulates the in-medium width of mesons according to
  !   Gamma=Gamma_free+sigma*rho*v
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming
  !   particles in the vacuum
  ! * Simplification: An average over the charge of the meson is performed.
  !
  ! INPUTS
  ! * NONE
  !
  ! OUTPUT
  ! * Files "inMediumWidth/InMediumWidth.particleID.dat.bz2"
  !****************************************************************************
  subroutine tabulate_inMediumWidth_mesons()
    use idTable, only: nmes,pion
    use inputGeneral, only: path_To_Input
    use output, only: intToChar,line
    use mesonWidthMedium_tables, only: numSteps_rhoN,numSteps_rhoP,numSteps_absP,numsteps_mass,max_rhoP,max_rhoN
    use bzip

    real, parameter :: min_mass=0.01
    real, parameter :: max_mass=3.0
    real ::  delta_rhoN,delta_rhoP,delta_absP,delta_mass
    real, dimension(:,:,:,:), allocatable :: widthTable
    integer :: iAbsP,iMass,iRhoN,iRhoP,ID
    character(20) :: format
    character(1000) :: fileName
    real :: rhoN, rhoP, absP, mass
    type(bzFile) :: f
    character(len=100) :: buf

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.

       write(*,*) 'tabulation parameters:'
       write(*,*) 'numSteps_rhoN= ',numSteps_rhoN
       write(*,*) 'numSteps_rhoP= ',numSteps_rhoP
       write(*,*) 'numSteps_absP= ',numSteps_absP
       write(*,*) 'numSteps_mass= ',numSteps_mass
       write(*,*) 'max_rhoP= ',max_rhoP
       write(*,*) 'max_rhoN= ',max_rhoN
    end if

    allocate(widthTable(0:numSteps_absP,0:numSteps_mass,0:numSteps_rhoN,0:numSteps_rhoP))
    delta_rhoN=max_rhoN/float(numSteps_rhoN)
    delta_rhoP=max_rhoP/float(numSteps_rhoP)
    delta_absP=max_absP_Mes/float(numSteps_absP)
    delta_mass=(max_mass-min_mass)/float(numSteps_mass)

    format='(1P,'//intToChar(numSteps_rhoP+1)//'E11.4)'

    write(*,*) 'Tabulating IN-MEDIUM WIDTH table for MESONS'
    do ID=max(pion,minMes),min(pion+nmes-1,maxMes)
       write(*,*) 'particleID=',ID
       do iAbsP=0,numSteps_absP
          absP=float(iAbsP)*delta_absP

          do iMass=0,numSteps_mass
             mass=min_mass+float(iMass)*delta_mass

             write(*,'(A,F15.4,2(A,I6),A,A,F15.4,2(A,I6),A)') &
                  'absP=',absP,' (',iAbsP,'/',numsteps_absP,')    ', &
                  'mass=',mass,' (',iMass,'/',numsteps_mass,') '

             do iRhoN=0,numSteps_rhoN
                rhoN=float(iRhoN)*delta_rhoN
                do iRhoP=0,numSteps_rhoP
                   rhoP=float(iRhoP)*delta_rhoP
                   widthTable(iAbsP,iMass,iRhoN,iRhoP)=&
                            & evaluateCollisionBroadening_mesons(ID,absP,mass,rhoN,rhoP)
                end do
             end do
             write(*,*) '     w: ',widthTable(iAbsP,iMass,1,1),'...', &
                        widthTable(iAbsP,iMass,numSteps_rhoN,numSteps_rhoP)

          end do
       end do

       ! Write the table to file:
       if (writelocal) then
          fileName='./InMediumWidth.'//trim(intToChar(ID))//'.dat.bz2'
       else
          fileName=trim(path_to_Input)//'/inMediumWidth/InMediumWidth.'//trim(intToChar(ID))//'.dat.bz2'
       end if
       f = bzOpenW(trim(fileName))
       write(buf,'(A,7I6)')'#', numSteps_absP,numSteps_mass,numSteps_rhoN,numSteps_rhoP,num_MonteCarlo_Points_mesons,cmin,cmax
       call bzWriteLine(f,buf)
       write(buf,'(A,1P,5E15.4)')'#', max_absP_Mes, max_rhoP , max_rhoN,min_mass,max_mass
       call bzWriteLine(f,buf)
       do iAbsP=0,numSteps_absP
          do iMass=0,numSteps_mass
             do iRhoN=0,numSteps_rhoN
                write(buf,format) widthTable(iAbsP,iMass,iRhoN,0:numSteps_rhoP)
                call bzWriteLine(f,buf)
             end do
          end do
       end do
       call bzCloseW(f)

       write(*,*) '... finished particleID=',ID
    end do
    write(*,*) '... finished tabulating IN-MEDIUM WIDTH table for MESONS'
    write(*,*) line
    write(*,*)

    deallocate(widthTable)

  end subroutine tabulate_inMediumWidth_mesons



  !****************************************************************************
  !****f* inMediumWidth/get_pauliBlockedDecayWidth
  ! NAME
  ! function get_pauliBlockedDecayWidth(particleID,momentumLRF,mass,rhoN,rhoP) result(gammaDecay)
  !
  ! PURPOSE
  ! * Returns the in-medium width of baryons according to
  !   Gamma=Gamma_free*Pauli_Blocking
  ! * An average over the Fermi sea is performed
  ! * Simplification: An average over the charge of the baryon is performed.
  !
  ! INPUTS
  ! * integer, intent(in) :: particleID  -- ID of baryon
  ! * real, intent(in)    :: absP        -- absolute Momentum
  ! * real, intent(in)    :: mass        -- Mass of baryon
  ! * real, intent (in)   :: rhoN,rhoP   -- proton and neutron density in fm^-3
  !
  ! OUTPUT
  ! * real :: gammaDecay
  !****************************************************************************
  function get_pauliBlockedDecayWidth(particleID,momLRF,mass,rhoN,rhoP) result(gammaDecay)

    use random
    use master_1Body, only: decayParticle
    use particleProperties, only: validcharge
    use idTable, only: nucleon
    use particleDefinition
    use baryonWidthMedium_tables, only: num_MonteCarlo_Points,get_min_charge_loop,get_max_charge_loop

    integer, intent(in)              :: particleID
    real, intent(in)                 :: momLRF
    real, intent(in)                 :: mass

    real :: gammaDecay

    type(particle)                   :: part
    type(particle),dimension(1:10)   :: finalState


    logical:: pauliBlocking

    logical :: collisionFlag, finalFlag
    real :: rhoN,rhoP ! neutron and proton density
    real :: gamma, fermi_n, fermi_p

    integer :: i,j

    call setToDefault(part)
    part%ID=particleID
    part%mass=mass
    part%mom(1:2)= 0
    part%mom(3)=momLRF
    part%mom(0) = freeEnergy(part)
    part%vel=part%mom(1:3)/part%mom(0)
    part%pos=999.   ! outside nucleus

    fermi_p=fermiMom(rhoP)
    fermi_n=fermiMom(rhoN)
    finalFlag=.true.

    gammaDecay=0.
    monteCarloLoop :  do i=1,num_MonteCarlo_Points

       !***********************************************************************
       ! PRELIMINARY : AVERAGING OVER CHARGE
       do
          part%charge=NINT(float(get_min_charge_loop())+float(get_max_charge_loop()-get_min_charge_loop())*rn())
          if (validCharge(part%ID,part%charge)) exit
       end do
       !***********************************************************************

       finalState%ID=0

       call decayParticle(part,finalState,collisionFlag,finalFlag,0.,gamma)
       if (.not.collisionFlag) cycle monteCarloLoop

       ! Check Pauli-Blocking

       pauliBlocking=.false.
       pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)
          if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%anti) then
             if (finalState(j)%charge.eq.0) then
                if (AbsMom(finalState(j)).lt.fermi_n) then
                   pauliBlocking=.true.
                   exit pauliLoop
                end if
             else if (finalState(j)%charge.eq.1) then
                if (AbsMom(finalState(j)).lt.fermi_p) then
                   pauliBlocking=.true.
                   exit pauliLoop
                end if
             else
                write(*,*) 'error in inMediumWidth.charge.', finalstate(j)
                stop 'error in inMediumWidth.charge.'
             end if
          end if
       end do pauliLoop

       if (pauliBlocking) cycle monteCarloLoop

       gammaDecay=gammaDecay+gamma
    end do monteCarloLoop
    gammaDecay=gammaDecay/float(num_MonteCarlo_Points)
  end function get_pauliBlockedDecayWidth




  !****************************************************************************
  !****f* inMediumWidth/evaluateCollisionBroadening_baryons
  ! NAME
  ! function evaluateCollisionBroadening_baryons(particleID,momentumLRF,mass,rhoN,rhoP) RESULT(gcoll)
  !
  ! PURPOSE
  ! * Returns the in-medium width of baryons according to
  !   Gamma=sigma*rho*v *Pauli_Blocking
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming
  !   particles in the vacuum
  ! * Simplification: An average over the charge of the baryon is performed.
  !
  ! INPUTS
  ! * integer, intent(in) :: particleID  -- ID of baryon
  ! * real, intent(in)    :: absP        -- absolute Momentum
  ! * real, intent(in)    :: mass        -- Mass of baryon
  ! * real, intent (in)   :: rhoN,rhoP   -- proton and neutron density in GeV^3
  ! * integer, optional, intent (in) :: monte_in
  !   -- if this input is given, then the number of Monte Carlo points is chosen
  !   according to this input; otherwise it is equal to "num_MonteCarlo_Points".
  !
  ! OUTPUT
  ! * real :: gcoll
  ! * real, optional, intent(out) :: gcoll_elastic_out
  !   -- elastic collisional width
  !****************************************************************************
  function evaluateCollisionBroadening_baryons(particleID,momLRF,mass,rhoN,rhoP,monte_in,gcoll_elastic_out) RESULT(gcoll)
    use particleDefinition
    use random
    use constants, only: GeVSquared_times_mb, mN
    use particleProperties, only: validCharge
    use master_2Body, only: generateFinalState
    use IDTable, only: nucleon
    use baryonWidthMedium_tables, only: get_min_charge_loop,get_max_charge_loop
    use twoBodyTools, only: pCM

    integer, intent(in)              :: particleID
    real, intent(in)                 :: momLRF
    real, intent(in)                 :: mass
    real, intent(in)                 :: rhoN,rhoP
    integer, optional, intent(in)    :: monte_in ! Input for Monte carlo points
    real, optional, intent(out)      :: gcoll_elastic_out
    real :: gColl

    type(particle)                   :: part,nuc

    integer :: i,j,charge,num_MonteCarlo_Points_local

    real    :: stringFactor
    integer :: numEnsembles
    type(particle), dimension(1:2)  :: pair
    type(particle), dimension(1:100) :: finalState
    real    :: time
    logical :: collisionFlag
    logical :: HiEnergyFlag  ! .true. if fritiof was used
    integer :: HiEnergyType  ! 0:LowEnergy, 1:Fritiof, 2:Pythia
    real    :: sigmaTot

    real :: fermiMomentum, pnuc,cost,srts,vrel, fermi_p, fermi_n,rho

    real :: gammaLorentz

    logical :: pauliBlocking
    logical :: usePauli
    logical,parameter :: vac_check=.false. ! only for debugging!!!

    !real ::   imsig2,imsig3,imsigq,absP
    real :: gcoll_elastic
    logical :: pauliIsUsedforXsection!,elastic
    integer :: pauliblocked

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    if (present(monte_in)) then
       num_MonteCarlo_Points_local=monte_in
    else
       num_MonteCarlo_Points_local= 250 !for speed up, was "num_MonteCarlo_Points"
    end if

    stringFactor=1.
    numEnsembles=1
    time=999.

    call setToDefault(part)
    part%ID=particleID

    if (particleID.eq.nucleon) then
       ! We don't know the off-shell cross sections. Therefore we
       ! assume that the nucleon width is independent of mass:
       part%mass=mN
    else
       part%mass=mass
    end if
    part%mom(1:2)= 0
    part%mom(3)=momLRF
    part%mom(0) = freeEnergy(part)
    part%vel=part%mom(1:3)/part%mom(0)
    part%pos=999.   ! outside nucleus


    call setToDefault(nuc)
    nuc%ID=nucleon
    nuc%mass=mN
    nuc%pos=part%pos

    gColl=0.
    gColl_elastic=0.

    fermi_p=fermiMom(rhoP)
    fermi_n=fermiMom(rhoN)



    nucleonChargeLoop: do charge=0,1
       if (vac_check.and.charge.eq.0) cycle
       if (charge.eq.1) then
          fermiMomentum=fermi_p
          rho=rhoP
       else if (charge.eq.0) then
          fermiMomentum=fermi_n
          rho=rhoN
       end if
       if (rho.lt.1E-10) cycle

       nuc%charge=charge


       usePauli=.true.
       pauliblocked=0

       !check if srts>4 (much above smooth transition region) and might run into Pythia
       !if yes, reduce number of Monte-Carlo-points to 1 and neglect Pauli blocking
       !(only cross section needed)
       !take minimal srts: both momenta parallel

       nuc%mom(1) = 0.
       nuc%mom(2) = 0.
       nuc%mom(3) = fermiMomentum
       nuc%mom(0) = freeEnergy(nuc)
       srts=sqrtS(nuc,part)
       if (debugFlag) write(*,*) 'min srts', srts
       if (srts.gt.4.) then
          num_MonteCarlo_Points_local=1
          usePauli=.false.
          if (debugFlag) write(*,*) 'speedup for srts', srts
       end if


       monteCarloLoop :  do i=1,num_MonteCarlo_Points_local
          ! Setting up the incoming nucleon
          !  * momentum of incoming nucleon:

          !********************************************************************
          ! PRELIMINARY : AVERAGING OVER CHARGE
          do
             if (vac_check) then
                part%charge=1 ! Set resonance charge to 1 for debugging
             else
                part%charge=NINT(float(get_min_charge_loop())+float(get_max_charge_loop()-get_min_charge_loop())*rn())
             end if
             if (validCharge(part%ID,part%charge)) exit
          end do
          !********************************************************************

          pnuc=(rn())**(1./3.)*fermiMomentum
          cost=(rn()-0.5)*2.
          nuc%mom(1) = pnuc*sqrt(max(1.-cost**2,0.))
          nuc%mom(2) = 0.
          nuc%mom(3) = pnuc*cost

          if (vac_check) nuc%mom(1:3) = 0.!set nuc momentum to 0 for debugging

          nuc%mom(0) = freeEnergy(nuc)
          nuc%vel    = nuc%mom(1:3)/nuc%mom(0)

          pair(1)=part
          pair(2)=nuc

          call setToDefault(finalState)
          finalState%ID=0

          call generateFinalState(pair, finalState, stringFactor, numEnsembles,&
               time, collisionFlag, HiEnergyFlag, HiEnergyType, &
               sigTot_out=sigmaTot, pauliIncluded_out=pauliIsUsedforXsection)
          if (.not. collisionFlag) cycle monteCarloLoop

          ! Check Pauli-Blocking
          pauliBlocking=.false.
          if (usePauli.and.(.not.pauliIsUsedforXsection)) then
             pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)

                if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%anti) then
                   if (finalState(j)%charge.eq.0) then
                      if (AbsMom(finalState(j)).lt.fermi_n) then
                         pauliBlocking=.true.
                         exit pauliLoop
                      end if
                   else if (finalState(j)%charge.eq.1) then
                      if (AbsMom(finalState(j)).lt.fermi_p) then
                         pauliBlocking=.true.
                         exit pauliLoop
                      end if
                   else
                      write(*,*) 'error in inMediumWidth.charge.',finalState(j)
                      stop 'error in inMediumWidth.charge.'
                   end if
                end if
             end do pauliLoop
             if (pauliBlocking) then
                pauliblocked=pauliblocked+1
                if (debugFlag) write(*,*) 'Pauli blocked'
                cycle monteCarloLoop
             end if
          end if
          ! Evaluate relative velocity:

          ! vrel=sqrt(Dot_Product(nuc%vel-part%vel,nuc%vel-part%vel)) !!! THIS IS A BUG
          srts=sqrtS(pair)
          vrel= pcm(srts,pair(1)%mass,pair(2)%mass)*srts/(pair(1)%mom(0)*pair(2)%mom(0))


          if (debugFlag) then
             write(*,*) 'In:',pair%ID
             write(*,*) 'OUT:',finalState(1:3)%ID
             write(*,*) 'sigma:',sigmaTot
             write(*,*) 'vrel:',vrel
             write(*,*) 'rho:',rho
             write(*,*) 'flags:',collisionFlag,HiEnergyFlag,HiEnergyType
             write(*,*) 'Gamma:',vrel*rho*sigmaTot* GeVSquared_times_mb
             write(*,'(4G18.3)') absmom(nuc), absMom(part),momLRF
             write(*,'(4G18.3)') part%mass, nuc%mass
             write(*,'(4G18.3)') finalstate(1)%mass, finalstate(2)%mass
             write(*,'(4G18.3)') AbsMom(finalState(1)),AbsMom(finalState(2)),fermi_n,fermi_p
          end if

          ! Evaluate width in GeV:
          gColl=gColl+ vrel*rho*sigmaTot* GeVSquared_times_mb

          if (particleID.ne.nucleon) then
             if ((finalState(1)%ID.eq.nucleon.or.finalState(2)%ID.eq.nucleon)&
                  & .and.(finalState(1)%ID.eq.particleID.or.finalState(2)%ID.eq.particleID) &
                  & .and.(finalState(3)%ID.le.0) ) then
                ! Elastic event
                gColl_elastic=gColl_elastic+ vrel*rho*sigmaTot* GeVSquared_times_mb
             end if
          else
             if (finalState(1)%ID.eq.nucleon.and.finalState(2)%ID.eq.nucleon.and.finalState(3)%ID.le.0) then
                ! Elastic event
                gColl_elastic=gColl_elastic+ vrel*rho*sigmaTot* GeVSquared_times_mb
             end if
          end if

       end do monteCarloLoop

       if (debugFlag) write(*,*) 'Pauli blocked events', pauliblocked

    end do nucleonChargeLoop
    gColl=gColl/float(num_MonteCarlo_Points_local)
    gColl_elastic=gColl_elastic/float(num_MonteCarlo_Points_local)
    ! Transform the widht into the particles rest-frame
    gammaLorentz = 1./sqrt( 1. - dot_product(part%vel(1:3),part%vel(1:3)) )
    if (vac_check) gammaLorentz=1.
    gColl=gColl*gammaLorentz
    if (present(gcoll_elastic_out)) gcoll_elastic_out=gcoll_elastic*gammaLorentz

  end function evaluateCollisionBroadening_baryons



  !****************************************************************************
  !****f* inMediumWidth/evaluateCollisionBroadening_mesons
  ! NAME
  ! function evaluateCollisionBroadening_mesons(particleID,momentumLRF,mass,rhoN,rhoP) RESULT(gcoll)
  !
  ! PURPOSE
  ! * Returns the in-medium width of mesons according to Gamma=sigma*rho*v
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming
  !   particles in the vacuum
  ! * Simplification: An average over the charge of the meson is performed.
  !
  ! INPUTS
  ! * integer, intent(in)  :: particleID  -- ID of particle
  ! * real, intent(in)     :: absP        -- absolute Momentum
  ! * real, intent(in)     :: mass        -- Mass of meson
  ! * real, intent(in)     :: rhoN,rhoP   -- proton and neutron density in GeV^3
  !
  ! OUTPUT
  ! * real :: gcoll
  !****************************************************************************
  function evaluateCollisionBroadening_mesons(particleID,momLRF,mass,rhoN,rhoP) RESULT(gcoll)
    use particleDefinition
    use random
    use constants, only: GeVSquared_times_mb, mN
    use particleProperties, only: validCharge
    use master_2Body, only: generateFinalState
    use IDTable, only: nucleon
    use twoBodyTools, only: pCM

    integer, intent(in)              :: particleID
    real, intent(in)                 :: momLRF
    real, intent(in)                 :: mass
    real, intent(in)                 :: rhoN,rhoP
    real :: gColl

    type(particle)                   :: part,nuc

    integer :: i,j,charge,num_MonteCarlo_Points_local

    real    :: stringFactor=1.
    integer :: numEnsembles=1
    type(particle), dimension(1:2)  :: pair
    type(particle), dimension(1:100) :: finalState
    real    :: time=999.
    logical :: collisionFlag
    logical :: HiEnergyFlag  ! .true. if fritiof was used
    integer :: HiEnergyType  ! 0:LowEnergy, 1:Fritiof, 2:Pythia
    real    :: sigmaTot

    real :: fermiMomentum, pnuc,cost,srts,vrel, fermi_p, fermi_n,rho

    real :: gammaLorentz

    logical :: pauliBlocking
    logical :: usePauli

    logical :: pauliIsUsedforXsection

    if (readInput_flag) then
       call readinput
       readInput_flag=.false.
    end if

    num_MonteCarlo_Points_local = num_MonteCarlo_Points_mesons

    call setToDefault(part)
    part%ID=particleID

    part%mass=mass
    part%mom(1:2)= 0
    part%mom(3)=momLRF
    part%mom(0) = freeEnergy(part)
    part%vel=part%mom(1:3)/part%mom(0)
    part%pos=999.   ! outside nucleus

    call setToDefault(nuc)
    nuc%ID=nucleon
    nuc%mass=mN
    nuc%pos=part%pos

    gColl=0.

    fermi_p=fermiMom(rhoP)
    fermi_n=fermiMom(rhoN)

    nucleonChargeLoop: do charge=0,1
       if (charge.eq.1) then
          fermiMomentum=fermi_p
          rho=rhoP
       else if (charge.eq.0) then
          fermiMomentum=fermi_n
          rho=rhoN
       end if
       if (rho.lt.1E-10) cycle

       nuc%charge=charge

       usePauli=.true.

       !check if srts>4 (much above smooth transition region) and might run into Pythia
       !if yes, reduce number of Monte-Carlo-points to 1 and neglect Pauli blocking
       !(only cross section needed)
       !take minimal srts: both momenta parallel

       nuc%mom(1) = 0.
       nuc%mom(2) = 0.
       nuc%mom(3) = fermiMomentum
       nuc%mom(0) = freeEnergy(nuc)
       srts=sqrtS(nuc,part)
       if (debugFlag) write(*,*) 'min srts', srts
       if (srts>4.) then
          num_MonteCarlo_Points_local=1
          usePauli=.false.
          if (debugFlag) write(*,*) 'speedup for srts', srts
       end if


       monteCarloLoop :  do i=1,num_MonteCarlo_Points_local

          !********************************************************************
          ! PRELIMINARY : AVERAGING OVER CHARGE
          do
             part%charge=NINT(float(cmin)+float(cmax-cmin)*rn())
             if (validCharge(part%ID,part%charge)) exit
          end do
          !********************************************************************

          ! Setting up the incoming nucleon
          !  * momentum of incoming nucleon:
          pnuc=(rn())**(1./3.)*fermiMomentum
          cost=(rn()-0.5)*2.
          nuc%mom(1) = pnuc*sqrt(max(1.-cost**2,0.))
          nuc%mom(2) = 0.
          nuc%mom(3) = pnuc*cost
          nuc%mom(0) = freeEnergy(nuc)
          nuc%vel    = nuc%mom(1:3)/nuc%mom(0)

          pair(1)=part
          pair(2)=nuc
          call setToDefault(finalState)
          finalState%ID=0
          call generateFinalState(pair, finalState, stringFactor, numEnsembles,&
               time, collisionFlag, HiEnergyFlag, HiEnergyType, &
               sigTot_out=sigmaTot, pauliIncluded_out=pauliIsUsedforXsection)
          if (.not.collisionFlag) cycle monteCarloLoop

          ! Check Pauli-Blocking
          pauliBlocking=.false.
          if (usePauli.and.(.not.pauliIsUsedforXsection)) then
             pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)
                if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%anti) then
                   select case(finalState(j)%charge)
                   case (0)
                      if (AbsMom(finalState(j))<fermi_n) then
                         pauliBlocking=.true.
                         exit pauliLoop
                      end if
                   case (1)
                      if (AbsMom(finalState(j))<fermi_p) then
                         pauliBlocking=.true.
                         exit pauliLoop
                      end if
                   end select
                end if
             end do pauliLoop
             if (pauliBlocking) then
                cycle monteCarloLoop
             end if
          end if
          ! Evaluate relative velocity:
          ! vrel=sqrt(Dot_Product(nuc%vel-part%vel,nuc%vel-part%vel)) !!! THIS IS A BUG
          srts=sqrtS(pair)
          vrel= pcm(srts,pair(1)%mass,pair(2)%mass)*srts/(pair(1)%mom(0)*pair(2)%mom(0))

          if (debugFlag) then
             write(*,*) 'In:',pair%ID
             write(*,*) 'OUT:',finalState(1:3)%ID
             write(*,*) 'sigma:',sigmaTot
             write(*,*) 'vrel:',vrel
             write(*,*) 'rho:',rho
             write(*,*) 'flags:',collisionFlag,HiEnergyFlag,HiEnergyType
             write(*,*) 'Gamma:',vrel*rho*sigmaTot* GeVSquared_times_mb
             write(*,'(4G18.3)') absmom(nuc), absMom(part),momLRF
             write(*,'(4G18.3)') part%mass, nuc%mass
             write(*,'(4G18.3)') finalstate(1)%mass, finalstate(2)%mass
             write(*,'(4G18.3)') AbsMom(finalState(1)),AbsMom(finalState(2)),fermi_n,fermi_p
          end if

          ! Evaluate width in GeV:
          gColl=gColl+ vrel*rho*sigmaTot* GeVSquared_times_mb

       end do monteCarloLoop

    end do nucleonChargeLoop
    gammaLorentz = 1./sqrt( 1. - dot_product(part%vel(1:3),part%vel(1:3)) )
    gColl=gColl*gammaLorentz/float(num_MonteCarlo_Points_local)

  end function evaluateCollisionBroadening_mesons



  !****************************************************************************
  !****if* inMediumWidth/fermiMom
  ! NAME
  ! real function fermiMom(rho)
  !
  ! PURPOSE
  ! * Returns the fermi momentum.
  !
  ! INPUTS
  ! * real, intent(in)   :: rho     -- proton or neutron density in GeV^-3
  !
  ! OUTPUT
  ! * Fermi momentum in GeV
  !****************************************************************************
  real function fermiMom(rho)
    use constants, only: pi

    real, intent(in) :: rho

    fermiMom=(3.*pi**2*rho)**(1./3.)
  end function fermiMom

end module inMediumWidth
