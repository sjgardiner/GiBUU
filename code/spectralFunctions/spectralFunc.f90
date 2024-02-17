!******************************************************************************
!****m* /spectralFunc
! NAME
! module spectralFunc
!
! PURPOSE
! Includes the spectral functions of the baryons.
!******************************************************************************
module spectralFunc

  use CallStack, only: TRACEBACK

  implicit none
  private

  public :: specFunc
  public :: propagator_nenner

  !****************************************************************************
  !****g* spectralFunc/inMed_collTerm
  ! SOURCE
  !
  logical, save :: inMed_collTerm=.false.
  !
  ! PURPOSE
  ! Use in-medium spectral functions according to sigma*rho*v,
  ! where sigma is taken according to the collision term.
  ! The correct normalisation is included here.
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFunc/inMed_flag
  ! SOURCE
  !
  logical, save :: inMed_flag=.false.
  !
  ! PURPOSE
  ! Use in-medium width according to settings in baryonWidthMedium.
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFunc/relativistic
  ! SOURCE
  !
  logical, save :: relativistic=.true.
  !
  ! PURPOSE
  ! Use either relativistic or non-relativistic spectral functions.
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFunc/which_nuclwidth
  ! SOURCE
  !
  integer, save :: which_nuclwidth=1
  !
  ! PURPOSE
  ! This flag decides what is used for the nucleon width.
  !
  ! Possible values:
  ! * 1 : use constant width given in const_nuclwidth
  ! * 2 : use width increasing linear with density;
  !   Gamma=const*rho/rho0 with const given in nuclwidth_dens
  ! * 3 : use toy model (constant NN cross section)
  ! * 4 : use realistic width (cf. diploma thesis of D. Kalok)
  ! * 5 : use realistic width: width based on our collision term
  !
  ! NOTES
  ! The correct normalisation has not been included here!!
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFunc/nuclwidth
  ! SOURCE
  !
  real, save :: nuclwidth=0.001
  !
  ! PURPOSE
  ! if which_nuclwidth=1, nuclwidth gives the width used in the Breit-Wigner
  ! for the nucleon
  !****************************************************************************


  !****************************************************************************
  !****g* spectralFunc/nuclwidth_dens
  ! SOURCE
  !
  real, save :: nuclwidth_dens=0.006
  !
  ! PURPOSE
  ! if which_nuclwidth=2, nuclwidth_dens gives the width used in density
  ! dependent width
  !
  ! 6 MeV are motivated in F. Froemel dissertation
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFunc/nuclwidth_sig
  ! SOURCE
  !
  real, save :: nuclwidth_sig=5.5
  !
  ! PURPOSE
  ! if which_nuclwidth=3, nuclwidth_sig gives the NN cross section in fm^2
  !****************************************************************************

  !****************************************************************************
  !****g* spectralFunc/widthMass
  ! SOURCE
  !
  integer, save :: widthMass = 1
  !
  ! PURPOSE
  ! select which mass is used to calculate the width:
  ! * 1: bare mass
  ! * 2: invariant mass
  ! * 3: invariant mass + mass shift of nucleon (for PDM)
  !   (= 1 for RMF and Skyrme)
  !****************************************************************************


  !****************************************************************************
  !****g* spectralFunc/debugFlag
  ! SOURCE
  logical, parameter :: debugFlag=.false.
  ! PURPOSE
  ! Switch for debug information
  !****************************************************************************

  logical, save :: initflag=.true.
!   integer, save :: scapotFailure_counter = 0  ! number of times errors will be displayed



contains

  subroutine readinput
    use output, only: Write_ReadingInput
    use baryonWidthMedium, only: get_MediumSwitch_coll, get_MediumSwitch

    integer :: IOS

    !**************************************************************************
    !****n* spectralFunc/spectralFunction
    ! NAME
    ! NAMELIST /spectralFunction/
    ! PURPOSE
    ! Includes the switches:
    ! * which_nuclwidth
    ! * nuclwidth
    ! * nuclwidth_dens
    ! * nuclwidth_sig
    ! * relativistic
    ! * widthMass
    !**************************************************************************
    NAMELIST /spectralFunction/ which_nuclwidth, nuclwidth, &
         nuclwidth_dens, nuclwidth_sig, relativistic, &
         widthMass

    inMed_collTerm=get_MediumSwitch_coll()
    !in-medium width according to our collision term is used, is determined by
    !the flag mediumSwitch_coll in module baryonWidthMedium
    inMed_flag=get_MediumSwitch()
    !if in-medium width is calculated or not is determined by
    !the flag mediumSwitch in module baryonWidthMedium


    call Write_ReadingInput('spectralFunction',0)
    rewind(5)
    read(5,nml=spectralFunction,IOSTAT=IOS)
    call Write_ReadingInput('spectralFunction',0,IOS)
    if (debugflag) write(*,'(A)') ' Debugging mode'
    if (inMed_flag) then
       if (inMed_collTerm) then
          write(*,'(A)') ' use in-med width based on our collision term!!!!!'
       else
          write(*,'(A)') ' use in-med width for the resonances according to settings in baryonWidthMedium.'
          write(*,'(A,I4)') ' for the nucleon: which_nuclwidth=', which_nuclwidth
          select case(which_nuclwidth)
          case (1)
             write(*,'(A,F12.5)') ' nuclWidth=',nuclWidth
          case (2)
             write(*,'(A,F12.5)') ' nuclwidth_dens=',nuclwidth_dens
          case (3)
             write(*,'(A,F12.5)') ' nuclwidth_sig=',nuclwidth_sig
             if (nuclwidth_sig.le.0.) then
                call TRACEBACK('nuclwidth_sig is smaller than 0')
             end if
          end select

          write(*,'(A,L8)') ' Relativistic spectral functions?',relativistic
       end if
    else
       write(*,'(A,F10.4,A)') ' use free width for the resonances and ',nuclwidth,' GeV for the nucleon'
       write(*,'(A,L8)') ' Relativistic spectral functions?',relativistic
    end if

    write(*,'(A,I4)') ' mass used for the width: ',widthMass
    if (widthMass<1 .or. widthMass>3) call TRACEBACK("wrong value")

    call Write_ReadingInput('spectralFunction',1)
  end subroutine readinput


  !****************************************************************************
  !****f* spectralFunc/specFunc
  ! NAME
  ! real function specFunc(ID,charge,mom,pos,bareMass,width)
  !
  ! PURPOSE
  ! Returns the spectral function of a baryon as a function of
  ! absolute momentum, charge and position
  !
  ! INPUTS
  ! * integer                 :: id       -- Id of baryon
  ! * integer                 :: charge   -- Charge of baryon
  ! * real, dimension(0:3)    :: mom      -- momentum of baryon in LRF
  ! * real, dimension(1:3)    :: pos      -- position of baryon
  !
  ! OUTPUT
  ! * real, OPTIONAL :: bareMass -- bare mass of the particle
  ! * real, OPTIONAL :: width -- width used for the spectral function
  !
  ! NOTES
  ! The input parameter mom has to be the in-medium momentum, since
  ! the Breit-Wigner uses the pole mass shifted by the scalar potential.
  !****************************************************************************
  real function specFunc(ID,charge,mom,pos, bareMass_out,width_out)
    use constants, only: pi
    use IdTable, only: isBaryon, nucleon, S11_1535
    use particleProperties, only: hadron
    use minkowski, only: abs4
    use potentialMain, only: scapot
    use RMF, only: getRMF_flag, mDiracNucleon_Approx, mDirac1535_Approx
    use constants, only: mN
    use dichteDefinition
    use densityModule, only: densityAt
    use baryonWidthMedium_tables, only : get_inMediumWidth

    integer, intent(in)         :: ID
    integer, intent(in)         :: Charge
    real, dimension(0:3), intent(in) :: mom
    real, dimension(1:3), intent(in) :: pos
    real, optional, intent(out) :: bareMass_out
    real, optional, intent(out) :: width_out
	!real, optional, intent(out) :: mEff

    real :: width
    real :: invMass, baremass, mEff, mCorr
    real :: dummy,scalarPotential
    logical :: flagOK
    type(dichte) :: density
    real :: dens

    logical, parameter :: addCollWidth = .true.

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    ! Default return values:
    specFunc=0.
    if (present(baremass_out)) baremass_out=0.
    if (present(width_out)) width_out=0.

    ! Check input
    if (.not.isBaryon(id)) then
       write(*,*) 'ID=',id
       call TRACEBACK('only designed for baryons')
    end if

    ! Invariant mass of the particle:
    invMass=abs4(mom,flagOK)
    if (.not.flagOK) return ! ==> failure

    if (inMed_collTerm) then

       if (getRMF_flag()) then
          call TRACEBACK("RMF mode not yet implemented")
       end if

       dummy= scapot(ID,charge,mom,pos,baremass,flagOK)  !dummy call needed to get baremass
       if (.not.flagOK) return ! ==> failure
       specFunc=specFunc_Medium(ID,mom(0),sqrt(max(1.E-10,Dot_Product(mom(1:3),mom(1:3)))),pos,charge)

    else

       if (invMass.lt.hadron(ID)%minmass) return ! ==> failure

       if (.not.getRMF_flag()) then

          mEff = hadron(ID)%mass+scapot(ID,charge,mom,pos,baremass,flagOK)
          if (.not.flagOK) return ! ==> failure
          mCorr = baremass

       else

          density = densityAt(pos)
          dens = density%baryon(0)

          select case(ID)
          case (nucleon)
             mEff = mDiracNucleon_Approx(dens)
             scalarPotential = mEff - mN
          case (S11_1535)
             mEff = mDirac1535_Approx(dens)
             scalarPotential = mEff - hadron(ID)%mass
          case default
             mEff = mDiracNucleon_Approx(dens) - mN + hadron(ID)%mass
             scalarPotential = mEff - hadron(ID)%mass
          end select

          bareMass = invMass - scalarPotential
          mCorr = invmass - (mDiracNucleon_Approx(dens) - mN)

       end if

       select case (widthMass)
       case (1)
          dummy = baremass ! old default
       case (2)
          dummy = invmass
       case (3)
          dummy = mCorr ! for PDM (= baremass for RMF,Skyrme)
       end select
       width=getWidth(ID,charge,mom,pos,dummy)

       ! for temporary use only:
!!$       if (addCollWidth) then
!!$          width = width + get_inMediumWidth(ID,mom,&
!!$               invmass,density%neutron(0),density%proton(0),1)
!!$
!!$       end if


       if (present(width_out)) width_out=width

       if (relativistic) then
          specfunc=width/pi&
               *invMass/((invMass**2-mEff**2)**2 +invMass**2*width**2)
       else
          specfunc=1./(4.*pi)/mEff*width/  &
               & ((invMass-mEff)**2+width**2/4.)
       end if
    end if

!!$    if (.not.success) then
!!$       if (scapotFailure_counter.lt.100) then
!!$         scapotFailure_counter = scapotFailure_counter + 1
!!$         write(*,*) 'SpecFunc: No succes in scapot',ID,charge,p,pos,baremass
!!$         write(*,*) 'Setting specFunc=0.!'
!!$         if (scapotFailure_counter.eq.100) write(*,'(A)') &
!!$          & 'WARNING : Error displayed 100 times.. Will suppress error output from now on'
!!$       endif
!!$       specfunc=0.
!!$    end if
    if (present(baremass_out)) baremass_out=baremass
  end function specFunc


  real function specFunc_Medium(ID,E,p,pos,charge)
    use mediumDefinition
    use mediumModule, only: mediumAt
    use particleProperties, only: hadron
    !use idTable, only : nucleon
    use constants, only: pi
    use selfenergy_baryons, only: get_RealPart,selfenergy_Imag

    real, dimension(1:3),intent(in) :: pos
    integer, intent(in) :: ID
    integer, intent(in) :: charge
    real, intent(in) :: p,E
    real :: mass
    real :: imPart
    real :: realPart!,gamma,imagPart,width
    type(medium) :: med

    med=mediumAt(pos)

    !Check Pauli blocking
    !   if(id.eq.nucleon) then
    !     Select case(charge)
    !    case(0)
    !      if(p.lt.fermiMom(rho_P*0.197**3)) then
    !        specfunc_Medium=0.
    !       return
    !   end if
    !case(1)
    !  if(p.lt.fermiMom(rho_N*0.197**3)) then
    !    specfunc_Medium=0.
    !   return
    !end if
    !end select
    !end if

    ImPart = selfenergy_Imag(ID,p,E,med,pos)
    if (abs(imPart).lt.1E-10) then
       specfunc_Medium=0.
       !write(*,*) 'SpecFunc_Medium=0 ! ID=',ID, '   p=',p, '    E=',E, '    pos=', pos
    else
       mass=sqrt(E**2-p**2)
       !if(mass.lt.minimalMass(ID)) then
       !   specfunc_Medium=0.
       !   return
       !end if
       realPart=get_RealPart(ID,p,mass,med,pos)
       specfunc_Medium=-1./pi*imPart/((E**2-p**2-hadron(Id)%mass**2-realPart)**2 + imPart**2)
    end if

  end function specFunc_Medium



!!$  real function fermiMom(rho)
!!$    use constants, only : pi
!!$    real :: rho
!!$    fermiMom=(3.*pi**2*rho)**(1./3.)
!!$  end function fermiMom


  !****************************************************************************
  !****f* spectralFunc/getWidth
  ! NAME
  ! real function getWidth(ID,charge,mom,pos,bareMass)
  !
  ! PURPOSE
  ! * Returns the width of a baryon as a function of absolute momentum in [GeV],
  !   charge and position
  !
  ! INPUTS
  ! * integer                 :: id       -- Id of baryon
  ! * integer                 :: charge   -- Charge of baryon
  ! * real, dimension(0:3)    :: mom      -- momentum of baryon in LRF
  ! * real, dimension(1:3)    :: pos      -- position of baryon
  ! * real                    :: bareMass
  !
  !****************************************************************************
  real function getWidth(ID_out,charge_out,mom,pos,bareMass)
    use particleDefinition
    use particleProperties, only: hadron
    use IDtable, only: nucleon
    use constants, only: rhoNull, hbarc, pi, mN
    use baryonwidthmedium, only: WidthBaryonMedium
    use mediumDefinition
    use mediumModule, only: mediumAt
    use minkowski, only: abs3,abs4
    use pn_medium_width, only: proton_width_medium
    use baryonWidthMedium_tables, only: get_inMediumWidth
    use potentialMain, only: scapot

    integer, intent(in) :: ID_out
    integer, intent(in) :: charge_out
    real, dimension(0:3), intent(in) :: mom
    real, dimension(1:3), intent(in) :: pos
    real, intent(in) :: bareMass
    real :: invmass

    type(particle) :: finalstate
    type(medium) :: mediumAtPosition

    real, dimension(0:3) :: p_fermi
    real :: omega_fermi

    !define outgoing particle
    call setToDefault(finalstate)
    finalstate%ID=ID_out
    finalstate%charge=charge_out
    finalstate%mom=mom
    finalstate%pos=pos
    finalState%pert=.true.

    mediumAtPosition=mediumAt(finalState%pos)

    invMass=abs4(finalstate%mom)

    ! if baryon is strongly off-shell, like in u-channel, then width=0,
    if (invMass .le. hadron(ID_out)%minmass) then
       getWidth=0
       return
    end if

    if (ID_out.eq.nucleon) then
       select case (which_nuclwidth)

       case (1)
          getWidth=nuclwidth

       case (2) !ONLY FOR TESTING PURPOSES
          getWidth=nuclwidth_dens*mediumAtPosition%density/rhoNull

       case (3) !ONLY FOR TESTING PURPOSES
          getWidth=mediumAtPosition%density*abs3(mom)/mN*hbarc*nuclwidth_sig

       case (4)
          if (charge_out.eq.1) then
             p_fermi(0:3)=(/hadron(ID_out)%mass,&
                  hbarc*(3*pi**2*mediumAtPosition%densityProton)**(1./3.),0.,0./)
          else
             p_fermi(0:3)=(/hadron(ID_out)%mass,&
                  hbarc*(3*pi**2*mediumAtPosition%densityNeutron)**(1./3.),0.,0./)
          end if
          omega_fermi=sqrt((hadron(ID_out)%mass+scapot(ID_out,charge_out,p_fermi,pos))**2 +p_fermi(1)**2)
          getWidth=proton_width_medium(mom,&
               & mediumAtPosition%densityProton,&
               & mediumAtPosition%densityNeutron,omega_fermi)

       case (5)
          getWidth=get_inMediumWidth(ID_out,mom,bareMass,&
               & mediumAtPosition%densityNeutron, &
               & mediumAtPosition%densityProton,3)

       case default
          write(*,*) 'which_nuclwidth=',which_nuclwidth
          call TRACEBACK('not implemented')
       end select

    else
       !write(*,'(A,I4,A,g11.5,A,4g11.5)') 'In getWidth: ID=',finalstate%ID, '   baremass=',bareMass, '   momentum=',finalstate%mom
       getWidth=WidthBaryonMedium(finalstate%ID,bareMass,finalstate%mom,mediumATposition)
    end if

  end function getWidth

  !****************************************************************************
  !****f* spectralFunc/propagator_nenner
  ! NAME
  ! real function propagator_nenner(ID,charge,mom,pos,bareMass)
  !
  ! PURPOSE
  ! Returns the denominator of the baryon propagator as a function of
  ! absolute momentum, charge and position
  !
  ! useful when one has to sum up different diagrams and take into account
  ! interference terms
  !
  ! INPUTS
  ! * integer                 :: id       -- Id of baryon
  ! * integer                 :: charge   -- Charge of baryon
  ! * real, dimension(0:3)    :: mom      -- momentum of baryon in LRF
  ! * real, dimension(1:3)    :: pos      -- position of baryon
  ! OUTPUT
  ! * real, OPTIONAL :: bareMass -- bare mass of the particle
  !****************************************************************************
  complex function propagator_nenner(ID,charge,mom,pos,bareMass_out)
    use constants, only: ii
    use IdTable, only: isBaryon
    use particleProperties, only: hadron
    use minkowski, only: sp
    use potentialMain, only: scapot

    integer, intent(in)         :: ID
    integer, intent(in)         :: Charge
    real, dimension(0:3)        :: mom
    real, dimension(1:3)        :: pos
    real, optional, intent(out) :: bareMass_out
    real :: invMass, invMass2, baremass, mass, width
    real :: scalarPotential
    logical :: success

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    ! Check input
    if (.not.isBaryon(id)) then
       write(*,*) 'ID=',id
       call TRACEBACK('only designed for baryons')
    end if

    baremass=hadron(ID)%mass
    mass=baremass

    ! Invariant mass of the particle:
    invMass2=SP(mom,mom)
    if (invMass2.gt.0) then
       invMass=sqrt(invMass2)
       if (invMass.le.hadron(ID)%minmass) then
          invMass=baremass
          scalarPotential=0
          width=0
       else
          scalarPotential=scapot(ID,charge,mom,pos,baremass,success) ! baremass=OUT
          if (.not.success) scalarPotential=0
          width=getWidth(ID,charge,mom,pos,baremass)
       end if
    else
       ! this happens in the u-channel
       invMass=baremass
       scalarPotential=0
       width=0  !hadron(ID)%width
    end if

    propagator_nenner=SP(mom,mom)-(mass+scalarPotential)**2 + ii*invMass*width

    if (present(baremass_out)) baremass_out=baremass
  end function propagator_nenner




!   complex function propagator_nenner_Medium(ID,E,p,pos,charge)
!     use mediumDefinition
!     use mediumModule, only: mediumAt
!     use particleProperties, only: hadron
!     use minkowski, only : ii
!     use selfenergy_baryons, only : get_RealPart,selfenergy_Imag
!
!     real, dimension(1:3),intent(in) :: pos
!     integer,intent(in) :: ID
!     integer , intent(in) :: charge
!     real,intent(in) :: p,E
!
!     real :: mass, imPart, realPart
!     type(medium) :: med
!
!     med=mediumAt(pos)
!
!     ImPart  =    selfenergy_Imag(ID,p,E,med)
!     iF(abs(imPart).lt.1E-10) then
!        propagator_nenner_Medium=1.e12
!     else
!        mass=sqrt(E**2-p**2)
!        realPart=get_RealPart(ID,p,mass,med,pos)
!        propagator_nenner_Medium=1./(E**2-p**2-hadron(Id)%mass**2-realPart + ii*imPart)
!     end if
!
!   end function propagator_nenner_Medium



end module spectralFunc
