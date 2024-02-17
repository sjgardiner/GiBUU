!******************************************************************************
!****m* /inMediumWidth_rho
! NAME
! module inMediumWidth_rho
! PURPOSE
! * Implements the routines for the medium width of the rho meson
!   at finite baryon density and temperature
! * Generates tabulated input files
!******************************************************************************
module inMediumWidth_rho

  use constants

  implicit none

  private

  integer, parameter :: n_mass = 100   ! number of mass points
  integer, parameter :: n_p= 50         ! number of momentum points
  integer, parameter :: n_dens = 60    ! number of density points
  integer, parameter :: n_temp = 20    ! number of temperature points

  integer, save :: n_dens_min = 0,   n_dens_max = n_dens   ! for jobcard
  integer, save :: n_temp_min = 0, n_temp_max = n_temp     ! for jobcard
  integer, save :: n_p_min = 0, n_p_max = n_p              ! for jobcard

  real,    parameter :: min_mass=2.*mElec  ! min. mass (GeV)
  real,    parameter :: max_mass=3.0   ! max. mass (GeV)
  real,    parameter :: p_max = 3.0    ! max. momentum for tabulation (GeV/c)
  real,    parameter :: dens_max = 0.96 ! max. baryon density for tabulation (fm^-3)
  real,    parameter :: temp_max = 0.2 ! max. temperature for tabulation (GeV)

  integer, parameter :: num_mc = 1000 ! number of monte-carlo points

  logical, save :: usePauli=.true.   ! if .false. -- no Pauli blocking

  integer, save :: imode=1      ! 1 -- final state generation
                                ! 2 -- full resonance contribution in the medium (w/o Pauli blocking)
                                ! 3 -- general expression for cross section (w/o Pauli blocking)

  logical, save  :: writeLocal=.true.  ! if .false., the output bz2-file
                                      ! will be written in the buuinput/inMediumWidth directory

  logical, save  :: readInp=.true.


  public :: tabulate_GammaColl_rho, GammaColl

contains


  !****s* inMediumWidth_rho/readInput
  ! NAME
  ! subroutine readInput
  ! PURPOSE
  ! Reads input in jobcard out of namelist "inMediumWidth_rho".
  !****************************************************************************
  subroutine readInput

    use output

    implicit none
    integer :: ios ! checks file behavior

    NAMELIST /inMediumWidth_rho/ n_dens_min, n_dens_max, n_temp_min, n_temp_max,&
                               & n_p_min, n_p_max, usePauli, imode

    call Write_ReadingInput('inMediumWidth_rho',0)
    rewind(5)
    read(5,nml=inMediumWidth_rho,IOSTAT=IOS)
    call Write_ReadingInput('inMediumWidth_rho',0,IOS)

    write(*,*)' n_p_min, n_p_max : ', n_p_min, n_p_max
    write(*,*)' n_dens_min, n_dens_max : ', n_dens_min, n_dens_max
    write(*,*)' n_temp_min, n_temp_max : ', n_temp_min, n_temp_max
    write(*,*)' usePauli : ', usePauli
    write(*,*)' imode : ', imode

    call Write_ReadingInput('inMediumWidth_rho',1)
  end subroutine readInput


  !****************************************************************************
  !****s* inMediumWidth_rho/tabulate_GammaColl_rho
  ! NAME
  ! subroutine tabulate_GammaColl_rho
  !
  ! PURPOSE
  ! * Tabulates the in-medium collisional width of rho meson according to Gamma_coll=sigma*rho*v
  ! * An average over the Fermi distribution at finite temperature is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  !
  ! INPUTS
  ! * NONE
  !
  ! OUTPUT
  ! * File "GammaColl.103.dat.bz2"
  !****************************************************************************
  subroutine tabulate_GammaColl_rho

    use mesonWidthVacuum, only: vacuumWidth
    use IDTable
    use output
    use inputGeneral, only: path_To_Input
    use bzip
!    use inMediumWidth, only: evaluateCollisionBroadening_mesons

    real :: delta_mass,delta_p,delta_dens,delta_temp,mass,p,dens,temp
    integer :: partId,index_mass,index_p,index_dens,index_temp
    real, dimension(:,:,:,:), allocatable :: widthTable
    character(25) :: format
    character(1000) :: fileName
    type(bzFile) :: f
    character(len=200) :: buf

    allocate(widthTable(0:n_mass,0:n_p,0:n_dens,0:n_temp))

    delta_mass=(max_mass-min_mass)/float(n_mass)
    delta_p=p_max/float(n_p)
    delta_dens=dens_max/float(n_dens)
    delta_temp=temp_max/float(n_temp)

    format='(' // intToChar(n_temp_max-n_temp_min+1) // '(1x,E13.7))'

    partId=rho

    if(readInp) then
       call readinput
       readInp=.false.
    end if

    write(*,*) 'Tabulating collisional width for rho-meson'

    open(1,file='GammaColl.'//inttochar(partID)//'.tst',status='unknown')

    do index_temp=n_temp_min,n_temp_max
        temp=float(index_temp)*delta_temp
        write(*,'(A,F15.4,2(A,I8))') 'temp=', temp,' steps:',index_temp,'/',n_temp
        do index_dens=n_dens_min,n_dens_max
            dens=float(index_dens)*delta_dens
            write(*,'(A,F15.4,2(A,I8))') 'dens=', dens,' steps:',index_dens,'/',n_dens
            write(1,'(2(A,F15.4))')'# temp=', temp,' dens=', dens
            write(1,'(A)')'# p:             mass:          Gamma_coll:' !    Gamma_vac:'
            do index_p=n_p_min,n_p_max
                p=float(index_p)*delta_p
                write(*,'(A,F15.4,2(A,I8))') 'p=', p,' steps:',index_p,'/',n_p
                do index_mass=0,n_mass
                    mass=min_mass+float(index_mass)*delta_mass
                    write(*,'(A,F15.4,2(A,I8))') 'mass=', mass,' steps:',index_mass,'/',n_mass
                    widthTable(index_mass,index_p,index_dens,index_temp)= &
                         & GammaColl(partId,mass,p,dens,temp)
!                    & evaluateCollisionBroadening_mesons(partID,p,mass,dens/2.*0.197**3,dens/2.*0.197**3)
                    write(1,'(4(2x,e13.7))') p,mass,widthTable(index_mass,index_p,index_dens,index_temp)!,vacuumWidth(mass,partID)
                end do
            end do
        end do
    end do

    close(1)

    write(*,*) '... finished tabulating collisional width for rho-meson'

    ! Write the table to file:
    if (writelocal) then
        fileName='./GammaColl.'//trim(intToChar(partID))//'.dat.bz2'
    else
        fileName=trim(path_to_Input)//'/inMediumWidth/GammaColl.'//trim(intToChar(partID))//'.dat.bz2'
    end if
    f = bzOpenW(trim(fileName))
    write(buf,'(A,4I6)')'#', n_temp, n_dens, n_p, n_mass
    call bzWriteLine(f,buf)
    write(buf,'(A,5(1x,E13.7))')'#', temp_max, dens_max, p_max, min_mass, max_mass
    call bzWriteLine(f,buf)
    write(buf,'(A,6I6)')'#', n_p_min, n_p_max, n_dens_min, n_dens_max, n_temp_min, n_temp_max
    call bzWriteLine(f,buf)
    do index_dens=n_dens_min,n_dens_max
        do index_p=n_p_min,n_p_max
            do index_mass=0,n_mass
                write(buf,format) widthTable(index_mass,index_p,index_dens,n_temp_min:n_temp_max)
                call bzWriteLine(f,buf)
            end do
        end do
    end do
    call bzCloseW(f)

    deallocate(widthTable)

  end subroutine tabulate_GammaColl_rho




  !******************************************************************************
  !****f* inMediumWidth_rho/GammaColl
  ! NAME
  ! real function GammaColl(partId,mass,p,dens,temp)
  !
  ! PURPOSE
  ! * Returns the collisional width of a meson in its rest frame according to GammaColl_rho=(sigma*rho*v) * gamma_Lor
  ! * An average over the Fermi sea is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  !
  ! INPUTS
  ! * integer, intent(in)            :: partId        ! id of meson
  ! * real, intent(in)               :: mass          ! off-shell mass of meson (GeV)
  ! * real, intent(in)               :: p             ! absolute Momentum (GeV/c)
  ! * real, intent(in)               :: dens          ! baryon (proton+neutron) density (fm^-3)
  ! * real, intent(in)               :: temp          ! temperature (GeV)
  !
  ! OUTPUT
  ! * real :: GammaColl (GeV)
  !****************************************************************************
  real function GammaColl(partId,mass,p,dens,temp)

    use particleDefinition
    use twoBodyTools, only: pCM
    use IDTable
    use random, only: rn, rnOmega
    use constants, only: hbarc, mN
    use master_2Body, only: generateFinalState, XsectionMaster
    use thermoDyn, only: muAt
    use distributions, only: Fermi

    use mediumDefinition
    use dichtedefinition
    use resonanceCrossSections, only: barMes2resonance
    use densitymodule, only: densityAt
    use mediumModule, only: mediumAt
    use preEventDefinition

    integer, intent(in)            :: partId        ! id of meson
    real, intent(in)               :: mass          ! off-shell mass of meson (GeV)
    real, intent(in)               :: p             ! absolute Momentum (GeV/c)
    real, intent(in)               :: dens          ! baryon (proton+neutron) density (fm^-3)
    real, intent(in)               :: temp          ! temperature (GeV)

    type(particle)                   :: part,nuc
    type(particle), dimension(1:2)  :: pair
    type(particle), dimension(1:100) :: finalState
    real, save    :: stringFactor=1.
    integer, save :: numEnsembles=1
    real, save    :: time=999.
    logical :: collisionFlag,pauliIsUsedforXsection
    logical :: HiEnergyFlag  ! .true. if fritiof was used
    integer :: HiEnergyType  ! 0:LowEnergy, 1:Fritiof, 2:Pythia
    real    :: mu,E_max,p_max,Fmax,gcoll,p_nuc,E_nuc,sigmaTot,Pacc,srts,vrel,gammaLorentz
    integer :: charge,i,isamp,imode_local,j

    type(dichte) :: density
    type(medium) :: mediumAtColl
    real, dimension(0:3) :: momLRF
    real, dimension(Delta:nbar) :: massRes       !  Resonance masses
    real, dimension(Delta:nbar) :: sigmaRes      ! rho N -> R cross section
    real, dimension(Delta:nbar) :: gcollRes      ! partial resonance contributions to the collisional width of the meson

    type(preEvent),dimension(1:4) :: chosenEvent
    real, dimension(0:7) :: sigs

    GammaColl = 0.

    if (dens.lt.1.E-06) return

    if(readInp) then
       call readinput
       readInp=.false.
    end if

    mu=muAt(dens*hbarc**3,temp)               ! chemical potential (GeV)
!    write(*,*)' dens, temp, mu : ', dens, temp, mu
    E_max = mu + 10.*temp                     ! arbitrary energy cutoff
    p_max = sqrt(E_max**2 - mN**2)            ! corresponding momentum
    Fmax=Fermi(mn,mu,temp)                    ! maximum occupation number
!    write(*,*)' Fmax : ', Fmax

    call setToDefault(part)
    part%ID=partId
    part%charge=0       ! assume rho^0

    part%mass=mass
    part%mom(1:2)= 0.
    part%mom(3)=p
    part%mom(0)=freeEnergy(part)
    part%vel=part%mom(1:3)/part%mom(0)
    part%pos=999.   ! outside nucleus

    call setToDefault(nuc)
    nuc%ID=nucleon
    nuc%mass=mN
    nuc%pos=part%pos

    gcoll=0.
    gcollRes=0.
    nucleonChargeLoop: do charge=0,1

        nuc%charge=charge

        monteCarloLoop :  do i=1,num_mc

           ! Sample nucleon momentum:
            isamp=0
            do
              isamp=isamp+1
               p_nuc=p_max*rn()**0.333333
               E_nuc=sqrt(p_nuc**2+mn**2)
               if(rn()*Fmax .lt. Fermi(E_nuc,mu,temp)) exit  ! check Fermi distribution
!               if(isamp.eq.100) then
!                  write(*,*)' can not sample nucleon momentum'
!                  write(*,*)' mu, temp :', mu,temp
!                  stop
!               end if
            end do
            if(isamp.gt.10) write(*,*)' isamp = ', isamp
            nuc%mom(1:3) = p_nuc * rnOmega()

            nuc%mom(0) = freeEnergy(nuc)
            nuc%vel    = nuc%mom(1:3)/nuc%mom(0)

            pair(1)=part
            pair(2)=nuc

            srts=sqrtS(pair)

            density=densityAt(pair(1)%pos)
            mediumAtColl=mediumAt(density,pair(1)%pos)
            momLRF=pair(1)%mom+pair(2)%mom
!            write(*,*)' mediumAtColl%density, mediumAtColl%useMedium :', mediumAtColl%density, mediumAtColl%useMedium
!            write(*,*)' mediumAtColl%betaLRF : ', mediumAtColl%betaLRF

            imode_local=imode
            if(imode.eq.1 .and. srts.gt.2.0) imode_local=3

            select case(imode_local)

            case(1)

                  call setToDefault(finalState)
                  call generateFinalState (pair, finalState, stringFactor, numEnsembles, time, collisionFlag,&
                    &HiEnergyFlag, HiEnergyType, sigTot_out=sigmaTot, pauliIncluded_out=pauliIsUsedforXsection)
                  if (.not.collisionFlag) cycle monteCarloLoop

            case(2)

                  sigmaRes=0.
                  sigmaRes = barMes2resonance (pair(1)%id,pair(2)%id,pair(1)%charge,pair(2)%charge,.true.,mediumAtColl, &
                         & momLRF,massRes,pair(1)%Mass,pair(2)%Mass,pair(1)%pos,.false.,srts)
                  sigmaTot=sum(sigmaRes)

            case(3)

                  call XsectionMaster(srts, pair, mediumAtColl, momLRF, &
                                      & chosenEvent, sigs, HiEnergyFlag,&
                                      & pauliIncluded_out=pauliIsUsedforXsection)
                  sigmaTot=sigs(0)

            end select

            if (usePauli.and.(.not.pauliIsUsedforXsection).and.imode_local.eq.1) then ! Check Pauli-Blocking
               Pacc=1.
               pauliLoop : do j=lbound(finalState,dim=1),ubound(finalState,dim=1)
                   if (finalState(j)%ID.eq.nucleon.and..not.finalstate(j)%anti) then
                      Pacc=Pacc*(1.-Fermi(finalState(j)%mom(0),mu,temp))
                   end if
                end do pauliLoop
                if(rn().gt.Pacc) cycle monteCarloLoop
            end if

            ! Evaluate relative velocity:
            vrel= pcm(srts,pair(1)%mass,pair(2)%mass)*srts/pair(1)%mom(0)/pair(2)%mom(0)

!            write(*,*)' vrel, sigmaTot: ', vrel, sigmaTot

            gcoll=gcoll+vrel*sigmaTot

            if(imode.eq.2) gcollRes=gcollRes+vrel*sigmaRes

        end do monteCarloLoop
    end do nucleonChargeLoop

    gcoll=gcoll/float(num_mc)*dens/2. * 0.1*hbarc   ! GeV

    gammaLorentz = 1./sqrt( 1. - dot_product(part%vel(1:3),part%vel(1:3)) )

    GammaColl=gcoll*gammaLorentz

    if(imode.eq.2) then
       gcollRes= gcollRes/float(num_mc)*dens/2. * 0.1*hbarc * gammaLorentz
       write(58,'(32(1x,e13.6))') p, mass, gcollRes(2:31)
    end if

  end function GammaColl

end module inMediumWidth_rho
