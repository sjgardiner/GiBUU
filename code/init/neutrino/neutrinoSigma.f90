!******************************************************************************
!****m* /neutrinoSigma
! NAME
! module neutrinoSigma
!
! PURPOSE
! This module bundles routines which are actually called by initNeutrino
! and represent some intermediate layer above neutrinoXsection
!******************************************************************************
module neutrinoSigma

  use particleDefinition
  use eN_eventDefinition
  use callstack, only: traceback

  implicit none

  private

  public :: Xsec_integratedSigma
  public :: Xsec_dSigmadCosThetadElepton
  public :: Xsec_dSigmadQ2dElepton
  public :: Xsec_dSigmadCosTheta
  public :: Xsec_dSigmadElepton
  public :: Xsec_SigmaMC
  public :: Xsec_SigmaMC_Q2
  public :: Xsec_SigmaMC_W
  public :: SetXsecMC
  public :: get_sigma_namelist

  !****************************************************************************
  !****g* neutrinoSigma/enu
  ! SOURCE
  real, save :: enu=-10.
  ! PURPOSE
  ! neutrino energy, read in by namelist
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/delta_enu
  ! SOURCE
  real, save :: delta_enu=-10.
  ! PURPOSE
  ! value by which the neutrino energy is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/Q2
  ! SOURCE
  real, save :: Q2=-10.
  ! PURPOSE
  ! momentum transfer squared
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/delta_Q2
  ! SOURCE
  real, save :: delta_Q2=-10.
  ! PURPOSE
  ! value by which Q2 is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/W
  ! SOURCE
  real, save :: W=-10.
  ! PURPOSE
  ! invariant mass defined as (p+q)^2
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/delta_W
  ! SOURCE
  real, save :: delta_W=-10.
  ! PURPOSE
  ! value by which W is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/costheta
  ! SOURCE
  real, save :: costheta=-10.
  ! PURPOSE
  ! cosine of the angle between the neutrino (z-direction) and the
  ! outgoing lepton
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/delta_costheta
  ! SOURCE
  real, save :: delta_costheta=-10.
  ! PURPOSE
  ! value by which costheta is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/elepton
  ! SOURCE
  real, save :: elepton=-10.
  ! PURPOSE
  ! energy of the outgoing lepton
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/delta_elepton
  ! SOURCE
  real, save :: delta_elepton=-10.
  ! PURPOSE
  ! value by which elepton is increased
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/eps
  ! SOURCE
  real, save :: eps=-10.
  ! PURPOSE
  ! polarisation in the case of EM interaction
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/MC_xmax
  ! SOURCE
  real,save :: MC_xmax = 2.0
  ! PURPOSE
  !
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/integralPrecision
  ! SOURCE
  integer, save ::  integralPrecision=3
  ! PURPOSE
  ! precision for the Gauss integration
  ! (reduce it for nuXsectionMode.eq.0 (sigma) to e.g. 2)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoSigma/integralPrecisionQE
  ! SOURCE
  integer, save ::  integralPrecisionQE=500
  ! PURPOSE
  ! precision for the Gauss integration over the QE peak
  ! (reduce it for nuXsectionMode.eq.0 (sigma) to e.g. 300)
  !****************************************************************************


  real,save :: MC_x,MC_y
  real,save :: xyJacob  ! Jacobian from dx dy -> dcos(theta) dE'

  public :: MC_x,MC_y

  logical, save :: initflag=.true.

contains

  subroutine readinput(nuXsectionmode, isExp)

    use output, only: Write_ReadingInput, Write_InitStatus
    use neutrino_IDTable
    use neutrinoParms, only: readInput_neutrino
    use neutrinoXsection, only: readinput_XS => readinput

    integer :: IOS
    integer, intent(in) :: nuXsectionmode
    logical, intent(in), optional :: isExp

    !**************************************************************************
    !****n* neutrinoSigma/nl_integratedSigma
    ! NAME
    ! NAMELIST /nl_integratedSigma/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=integratedsigma and
    ! includes:
    ! * enu
    ! * delta_enu
    ! * integralPrecision
    ! * integralPrecisionQE
    !**************************************************************************
    NAMELIST /nl_integratedSigma/ enu, delta_enu, &
         integralPrecision, integralPrecisionQE

    !**************************************************************************
    !****n* neutrinoSigma/nl_dSigmadCosThetadElepton
    ! NAME
    ! NAMELIST /nl_dSigmadCosThetadElepton/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadCosThetadElepton
    ! and includes:
    ! * enu
    ! * costheta
    ! * elepton
    ! * delta_elepton
    !**************************************************************************
    NAMELIST /nl_dSigmadCosThetadElepton/ enu, costheta, elepton, &
         delta_elepton

    !**************************************************************************
    !****n* neutrinoSigma/nl_dSigmadQ2dElepton
    ! NAME
    ! NAMELIST /nl_dSigmadQ2dElepton/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadQ2dElepton
    ! and includes:
    ! * enu
    ! * Q2
    ! * elepton
    ! * delta_elepton
    !**************************************************************************
    NAMELIST /nl_dSigmadQ2dElepton/ enu, Q2, elepton, delta_elepton

    !**************************************************************************
    !****n* neutrinoSigma/nl_dSigmadQ2
    ! NAME
    ! NAMELIST /nl_dSigmadQ2/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadQ2 and includes:
    ! * enu
    ! * Q2
    ! * delta_Q2
    !**************************************************************************
    NAMELIST /nl_dSigmadQ2/ enu, Q2, delta_Q2

    !**************************************************************************
    !****n* neutrinoSigma/nl_dSigmadW
    ! NAME
    ! NAMELIST /nl_dSigmadW/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadW and includes:
    ! * enu
    ! * W
    ! * delta_W
    !**************************************************************************
    NAMELIST /nl_dSigmadW/ enu, W, delta_W

    !**************************************************************************
    !****n* neutrinoSigma/nl_dSigmadcostheta
    ! NAME
    ! NAMELIST /nl_dSigmadcostheta/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadcostheta and
    ! includes:
    ! * enu
    ! * costheta
    ! * delta_costheta
    ! * integralPrecision
    ! * integralPrecisionQE
    !**************************************************************************
    NAMELIST /nl_dSigmadCosTheta/ enu, costheta, delta_costheta, &
         integralPrecision, integralPrecisionQE

    !**************************************************************************
    !****n* neutrinoSigma/nl_dSigmadElepton
    ! NAME
    ! NAMELIST /nl_dSigmadElepton/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmadElepton and
    ! includes:
    ! * enu
    ! * elepton
    ! * delta_elepton
    ! * integralPrecision
    ! * integralPrecisionQE
    !**************************************************************************
    NAMELIST /nl_dSigmadElepton/ enu, elepton, delta_elepton, &
         integralPrecision, integralPrecisionQE

    !**************************************************************************
    !****n* neutrinoSigma/nl_SigmaMC
    ! NAME
    ! NAMELIST /nl_SigmaMC/
    ! PURPOSE
    ! This Namelist is read in when nuXsectionmode=dSigmaMC and
    ! includes:
    ! * enu
    ! * MC_xmax
    !**************************************************************************
    NAMELIST /nl_SigmaMC/ enu,MC_xmax

    if (.not.initFlag) return

    call Write_InitStatus('neutrinodSigma',0)

    call readinput_XS()

    write(*,*)
    write(*,*) 'nuXsectionMode =',sXsectionMode(nuXsectionMode)
    write(*,*)

    select case (nuXsectionMode)

    case (integratedSigma)

       call Write_ReadingInput('nl_integratedSigma',0)
       rewind(5)
       read(5,nml=nl_integratedSigma,IOSTAT=IOS)
       call Write_ReadingInput('nl_integratedSigma',0,IOS)
       if (enu.lt.0..or.delta_enu.lt.0.) then
          write(*,*) 'input error in neutrinoXsection, enu= ', enu
          call traceback()
       end if
       write(*,*) 'integral precicions: ',integralPrecision,integralPrecisionQE
       call Write_ReadingInput('nl_integratedSigma',1)

    case (dSigmadCosThetadElepton)

       call Write_ReadingInput('nl_dSigmadCosThetadElepton',0)
       rewind(5)
       read(5,nml=nl_dSigmadCosThetadElepton,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadCosThetadElepton',0,IOS)
       if (enu.lt.0..or.abs(costheta).gt.1..or.elepton.lt.0..or.delta_elepton.lt.0.) then
          call traceback()
       end if
       call Write_ReadingInput('nl_dSigmadCosThetadElepton',1)

    case (dSigmadQ2dElepton)

       call Write_ReadingInput('nl_dSigmadQ2dElepton',0)
       rewind(5)
       read(5,nml=nl_dSigmadQ2dElepton,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadQ2dElepton',0,IOS)
       if (enu.lt.0..or.Q2.lt.0..or.elepton.lt.0..or.delta_elepton.lt.0.) then
          call traceback()
       end if
       call Write_ReadingInput('nl_dSigmadQ2dElepton',1)

    case (dSigmadCosTheta)

       call Write_ReadingInput('nl_dSigmadCosTheta',0)
       rewind(5)
       read(5,nml=nl_dSigmadCosTheta,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadCosTheta',0,IOS)
       if (enu.lt.0..or.abs(costheta).gt.1..or.delta_costheta.lt.0.) then
          call traceback()
       end if
       write(*,*) 'integral precicions: ',integralPrecision,integralPrecisionQE
       call Write_ReadingInput('nl_dSigmadCosTheta',1)

    case (dSigmadElepton)

       call Write_ReadingInput('nl_dSigmadElepton',0)
       rewind(5)
       read(5,nml=nl_dSigmadElepton,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadElepton',0,IOS)
       if (enu.lt.0..or.elepton.lt.0..or.delta_elepton.lt.0.) then
   !      write(*,*) 'enu=',enu,'elepton=',elepton,'delta_elepton=',delta_elepton
          call traceback()
       end if
       write(*,*) 'integral precicions: ',integralPrecision,integralPrecisionQE
       call Write_ReadingInput('nl_dSigmadElepton',1)

    case (dSigmaMC)

       if (present(isExp)) then
          if (isExp) then
             enu = 99.9 ! dummy
          end if
       end if

       call Write_ReadingInput('nl_SigmaMC',0)
       rewind(5)
       read(5,nml=nl_SigmaMC,IOSTAT=IOS)
       call Write_ReadingInput('nl_SigmaMC',0,IOS)

       write(*,'(" Enu      =",f12.3)') enu
       write(*,'(" xmax     =",f12.3)') MC_xmax

       if (enu.lt.0.) then
          call traceback('enu < 0')
       end if

       call Write_ReadingInput('nl_SigmaMC',1)

    case (dSigmaMC_dW)

       call Write_ReadingInput('nl_dSigmadW',0)
       rewind(5)
       read(5,nml=nl_dSigmadW,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadW',0,IOS)
       if (enu.lt.0. .or. W.lt.0. .or. delta_W.lt.0.) then
          write(*,*) 'input error in nl_dSigmadW: enu or W or delta_W <0'
          call traceback()
       end if
       call Write_ReadingInput('nl_dSigmadW',1)

    case (dSigmaMC_dQ2)

       call Write_ReadingInput('nl_dSigmadQ2',0)
       rewind(5)
       read(5,nml=nl_dSigmadQ2,IOSTAT=IOS)
       call Write_ReadingInput('nl_dSigmadQ2',0,IOS)
       if (enu.lt.0. .or. Q2.lt.0. .or. delta_Q2.lt.0.) then
          write(*,*) 'input error in nl_dSigmadQ2: enu or Q2 or delta_Q2 <0'
          call traceback()
       end if
       call Write_ReadingInput('nl_dSigmadQ2',1)

    case default
       write(*,*) 'error in case nuXsectionMode', nuXsectionMode
       call traceback()
    end select

    call Write_InitStatus('neutrinodSigma',1)

    initFlag=.false.

  end subroutine readinput

  !****************************************************************************
  !****s* neutrinoSigma/Xsec_SigmaMC
  ! NAME
  ! subroutine Xsec_SigmaMC(eNev,IP,raiseFlag,raiseVal,OutPart,sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates the integrated cross section depending
  ! on the input variables. It applies a MC integration technique.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_SigmaMC(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out

    if (initFlag) then
       call readInput(dSigmaMC, present(flux_enu))

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if

    end if

    raiseVal=eNev%lepton_in%mom(0)

    ! set default output
    call setToDefault(OutPart)
    sig=0.


    if (xyJacob .eq. 0.0) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)
    ! Monte-Carlo integration over x gives * MC_xmax; for the integration over y the  the factor is 1
    ! xyJacob is Jacobian from x,y to cos(theta) dE'
    sig = sig * xyJacob * MC_xmax

  end subroutine Xsec_SigmaMC

  !****************************************************************************
  !****s* neutrinoSigma/SetXsecMC
  ! NAME
  ! subroutine SetXsecMC()
  ! PURPOSE
  ! set the values of the integration variables in the MC integration
  ! called from subroutine initNeutrino in file initNeutrino.f90
  !****************************************************************************
  subroutine SetXsecMC(eNev, flux_enu, XsectionMode)
    use random, only: rn
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3c
    use neutrino_IDTable
    use minkowski, only: SP, abs4Sq

    type(electronNucleon_event),  intent(inout) :: eNev
    real, intent(in) :: flux_enu
    integer, intent(in) :: XsectionMode

    real :: Ein, Eprime
    real :: PK, PP
    logical :: flagOK

    if (initFlag) then
       select case (MOD(XsectionMode,10))
       case (dsigmaMC)
          call readInput(dSigmaMC, (XsectionMode>10))
       case (dsigmaMC_dW)
          call readInput(dSigmaMC_dW)
       case (dsigmaMC_dQ2)
          call readInput(dSigmaMC_dQ2)
       end select

       if (flux_enu.gt.0.0) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if

    Ein = enu
    if (flux_enu.gt.0.0) Ein = flux_enu

    call eNev_init_nuStep2(eNev,Ein) ! we ignore threshold checks
    PK = SP(eNev%lepton_in%mom,eNev%nucleon%mom)
    PP = abs4Sq(eNev%nucleon%mom)

    select case (MOD(XsectionMode,10))
    case (dsigmaMC)
      MC_x = rn()*MC_xmax
      MC_y = rn()
    case (dsigmaMC_dQ2)
      MC_y=rn()
      MC_x=Q2/(2.*MC_y*PK)
    case (dsigmaMC_dW)
      MC_y=rn()
      MC_x=1.-(W**2-PP)/(2.*MC_y*PK)
      !write(*,*) 'MC_y=',MC_y, '    MC_x=',MC_x
    end select

    call eNev_init_nuStep3c(eNev,MC_x,MC_y,flagOK)
    if (.not.flagOK) then
       xyJacob = 0.0
       !write(*,*) 'after eNev_init_nuStep3c flakOK is false. xyJacob = 0.0'
    else
       Eprime=eNev%lepton_out%mom(0)
       ! Jacobian xyJacob from dsigma/dE1dcostheta to dsigma/dxdy:
       xyJacob = (MC_y*PK**2)/(Ein *sqrt(Eprime**2-eNev%lepton_out%mass**2)*eNev%nucleon%mom(0))
       !write(*,*) 'xyJacob=',xyJacob
    end if

    ! MC integration is transferred to the subroutine  Xsec_SigmaMC

  end subroutine SetXsecMC

  !****************************************************************************
  !****s* neutrinoSigma/Xsec_SigmaMC_Q2
  ! NAME
  ! subroutine Xsec_SigmaMC_Q2(eNev,IP,raiseFlag,raiseVal,OutPart,sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates the integrated cross section depending
  ! on the input variables. It applies a MC integration technique.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_SigmaMC_Q2(eNev,IP,raiseFlag,raiseVal,OutPart,sig,flux_enu)

    use neutrino_IDTable
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out


    if (initFlag) then
       call readInput(dSigmaMC_dQ2)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if

    if (raiseFlag) then
       Q2=Q2+delta_Q2
       write(*,'(3(A,g12.5))') 'Enu=', enu, &
            & '      Q2 is raised by ...', delta_Q2, '  to  Q2=', Q2
    end if
    raiseVal=Q2


    ! set default output
    call setToDefault(OutPart)
    sig=0.


    if (xyJacob .eq. 0.0) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)

    !factor /Q2 *MC_x is   from dsi/x to dsi/dQ2
    sig = sig * xyJacob /Q2 *MC_x

  end subroutine Xsec_SigmaMC_Q2



  !****************************************************************************
  !****s* neutrinoSigma/Xsec_SigmaMC_W
  ! NAME
  ! subroutine Xsec_SigmaMC_W(eNev,IP,raiseFlag,raiseVal,OutPart,sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates the integrated cross section depending
  ! on the input variables. It applies a MC integration technique.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_SigmaMC_W(eNev,IP,raiseFlag,raiseVal,OutPart,sig,flux_enu)

    use neutrino_IDTable
    use minkowski, only: SP
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real    :: PK

    if (initFlag) then
       call readInput(dSigmaMC_dW)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if


    if (raiseFlag) then
       W=W+delta_W
       write(*,'(3(A,g12.5))') 'Enu=', enu, &
            & '      W is raised by ...', delta_W, '  to  W=', W
    end if
    raiseVal=W

    ! set default output
    call setToDefault(OutPart)
    sig=0.


    if (xyJacob .eq. 0.0) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)
    !write(*,*) 'sig=',sig
    !factor W /(PK *MC_y)  is   from dsi/x to dsi/dW
    PK = SP(eNev%lepton_in%mom,eNev%nucleon%mom)
    !write(*,*) 'PK=',PK, '   W=',W, '   MC_y=',MC_y, '    xyJacob=',xyJacob
    sig = sig * xyJacob *W /(PK *MC_y)

  end subroutine Xsec_SigmaMC_W




  !****************************************************************************
  !****s* neutrinoSigma/Xsec_integratedSigma
  ! NAME
  ! subroutine Xsec_integratedSigma(eNev,IP,raiseFlag,raiseVal,OutPart,sig,
  ! flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates the integrated cross section depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_integratedSigma(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use idtable, only: nucleon
    use random
    use gauss_integration
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a
    !!! use lepton_kinematics_free, only: minmaxE1_costheta
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge


    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min

    integer :: intprec,intprec1,l,j,n2,n1
    real :: eleptonint,costhetaint
    real :: sigmaximum,sigrd
    real, dimension(:),allocatable :: y,x,yy,xx
    logical :: flagOK


    if (initFlag) then
       call readInput(integratedSigma)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
    end if


    ! set default output
    call setToDefault(OutPart)
    sig=0.

    if (raiseFlag.and..not.present(flux_enu)) then
       enu=enu+delta_enu
       write(*,'(a,F12.4,a,F12.4)') 'Enu is raised by ...', &
            & delta_enu, ' to ', enu
    end if

    if (present(flux_enu)) enu=flux_enu
    raiseVal=enu

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return


    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    !!! call getCosthetaLimits(costheta_min,costheta_max,enu)
    call getCosthetaLimits(costheta_min,costheta_max)
    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)

    ! set integral precision
    intprec=integralPrecision
    intprec1=integralPrecision
    if (IP.eq.nucleon) intprec=integralPrecisionQE  ! due to small Breit-Wigner

    allocate (x(20*intprec1))
    allocate (y(20*intprec1))
    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    sigmaximum=0.

    call sg20r(costheta_min,costheta_max,intprec1,x,n1)
    call sg20r(elepton_min, elepton_max, intprec, xx, n2)

    do j=1,n1
       costhetaint=x(j)

       ! in this version the elepton_min depends on the costheta
       ! thus integration is performed only over kinematically allowed region
       !!! call minmaxE1_costheta(IP,enu,eNev%nucleon%mass,eNev%lepton_out%mass,costhetaint,elepton_min,elepton_max,success=flagOK)
       !!! if (.not.flagOK) cycle
       !!! if (debugflag) write(*,'(3(A,f12.4))') ' costheta=', costhetaint, '      E1min=', elepton_min, '    E1max=', elepton_max
       !!! call sg20r(elepton_min, elepton_max, intprec, xx, n2)

       yy = 0.0
       do l=1,n2
          eleptonint=xx(l)
          call eNev_init_nuStep3a(eNev,eleptonint,costhetaint,flagOK)
          if (.not.flagOK) cycle
          call XsecdCosthetadElepton(eNev,IP,OutPart,yy(l))
          if (yy(l).gt.sigmaximum) sigmaximum=yy(l)
       end do
       call rg20r(elepton_min,elepton_max,intprec,yy,y(j))
    end do
    call rg20r(costheta_min,costheta_max,intprec1,y,sig)

    ! set kinematics !!!! TODO: avoid infinite loop !!!!
    sigrd=0.
    if (sigmaximum.gt.0.) then
       do
          eleptonint=elepton_min+rn()*(elepton_max-elepton_min)
          costhetaint=costheta_min+rn()*(costheta_max-costheta_min)

          call eNev_init_nuStep3a(eNev,eleptonint,costhetaint,flagOK)
          if (.not.flagOK) cycle
          call XsecdCosthetadElepton(eNev,IP,OutPart, sigrd)
          if (sigmaximum*rn().le.sigrd) exit
       end do
    end if

    deallocate(x,y,xx,yy)

  end subroutine Xsec_integratedSigma


  !****************************************************************************
  !****s* neutrinoSigma/Xsec_dSigmadCosThetadElepton
  ! NAME
  ! subroutine Xsec_dSigmadCosThetadElepton(eNev,IP,raiseFlag,raiseVal,
  ! OutPart,sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadCosThetadElepton depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadCosThetadElepton(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min
    logical :: flagOK

    if (initFlag) then
       call readInput(dSigmadCosThetadElepton)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'costheta= ', costheta
       write(*,'(a,F12.4)') 'elepton=  ', elepton
    end if


    !set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if
    raiseVal=Elepton

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    !check threshold for costheta and elepton
    call getCosthetaLimits(costheta_min,costheta_max)
    if (.not.checkLimits(costheta_min,costheta_max,costheta)) return

    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    if (.not.checkLimits(elepton_min,elepton_max,elepton)) return

    call eNev_init_nuStep3a(eNev,elepton,costheta,flagOK)
    if (.not.flagOK) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)

  end subroutine Xsec_dSigmadCosThetadElepton


  !****************************************************************************
  !****s* neutrinoSigma/Xsec_dSigmadQ2dElepton
  ! NAME
  ! subroutine Xsec_dSigmadQ2dElepton(eNev,IP,raiseFlag,raiseVal,OutPart,
  ! sig,flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadQ2dElepton depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadQ2dElepton(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3b
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: elepton_max, elepton_min
    real :: Q2_min, Q2_max
    logical :: flagOK

    if (initFlag) then
       call readInput(dSigmadQ2dElepton)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'Q2=      ', Q2
       write(*,'(a,F12.4)') 'elepton= ', elepton
    end if


    !set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if
    raiseVal=Elepton

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    !check threshold for Q2 and elepton
    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    if (.not.checkLimits(elepton_min,elepton_max,elepton)) return

    call getQ2Limits(eNev,elepton,Q2_min,Q2_max)
    if (.not.checkLimits(Q2_min,Q2_max,Q2)) return

    call eNev_init_nuStep3b(eNev,elepton,Q2,flagOK)
    if (.not.flagOK) return

    call XsecdCosthetadElepton(eNev,IP,OutPart, sig)

    ! correct cross section for Jacobian:
    sig = sig/(2.*eNev%lepton_in%mom(0)*sqrt(elepton**2-eNev%lepton_out%mass**2))

  end subroutine Xsec_dSigmadQ2dElepton



  !****************************************************************************
  !****s* neutrinoSigma/Xsec_dSigmadcostheta
  ! NAME
  ! subroutine Xsec_dSigmadcostheta(eNev,IP,raiseFlag,raiseVal,OutPart,sig,
  ! flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadcostheta depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadcostheta(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use idtable, only: nucleon
    use random
    use gauss_integration
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min

    integer :: intprec,l,n2
    real :: eleptonint
    real :: sigmaximum,sigrd
    real, dimension(:),allocatable :: yy,xx
    logical :: flagOK

    if (initFlag) then
       call readInput(dSigmadcostheta)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'costheta= ', costheta
    end if

    ! set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       costheta=costheta+delta_costheta
       write(*,*) 'costheta is raised by ...', delta_costheta, ' to ', costheta
    end if
    raiseVal=costheta

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return


    call getCosthetaLimits(costheta_min,costheta_max)
    if (.not.checkLimits(costheta_min,costheta_max,costheta)) return

    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)

    !set integral precision
    intprec=integralPrecision
    if (IP.eq.nucleon) intprec=integralPrecisionQE ! due to small Breit-Wigner

    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    call sg20r(elepton_min,elepton_max,intprec,xx,n2)
    sigmaximum=0.0
    yy = 0.0

    do l=1,n2
       eleptonint=xx(l)

       call eNev_init_nuStep3a(eNev,eleptonint,costheta,flagOK)
       if (.not.flagOK) cycle

       call XsecdCosthetadElepton(eNev,IP,OutPart, yy(l))
       if (yy(l).gt.sigmaximum) sigmaximum=yy(l)
    end do

    call rg20r(elepton_min,elepton_max,intprec,yy,sig)

    ! set kinematics
    if (sigmaximum.gt.0.) then
       do
          eleptonint=elepton_min+rn()*(elepton_max-elepton_min)

          call eNev_init_nuStep3a(eNev,eleptonint,costheta,flagOK)
          if (.not.flagOK) cycle

          call XsecdCosthetadElepton(eNev,IP,OutPart, sigrd)

          if (sigmaximum*rn().le.sigrd) exit
       end do
    end if

    deallocate(xx,yy)

  end subroutine Xsec_dSigmadcostheta


  !****************************************************************************
  !****s* neutrinoSigma/Xsec_dSigmadElepton
  ! NAME
  ! subroutine Xsec_dSigmadElepton(eNev,IP,raiseFlag,raiseVal,OutPart,sig,
  ! flux_enu)
  !
  ! PURPOSE
  ! This subroutine calculates dSigmadElepton depending
  ! on the input variables:
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eNev -- the Lepton-Nucleon Event
  ! * integer             :: IP           -- ID of outgoing hadron/process
  ! * logical             :: raiseFlag    -- shall energy etc. be increased?
  ! * real, optional      :: flux_enu     -- if present, use this value as
  !   neutrino energy (used for experiments)
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * real                         :: raiseVal -- the actual value of the
  !   running variable (i.e. elepton, enu, costheta, ...)
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !****************************************************************************
  subroutine Xsec_dSigmadElepton(eNev, IP, raiseFlag, raiseVal, &
       & OutPart,sig,flux_enu)

    use neutrino_IDTable
    use idtable, only: nucleon
    use random
    use gauss_integration
    use eN_event, only: eNev_init_nuStep2,eNev_init_nuStep3a
    use neutrinoXsection, only: XsecdCosthetadElepton, SetHadronCharge

    type(electronNucleon_event),  intent(inout) :: eNev
    integer,                      intent(in)    :: IP
    logical,                      intent(in)    :: raiseFlag
    real,                         intent(out)   :: raiseVal
    type(particle), dimension(:), intent(out)   :: OutPart
    real,                         intent(out)   :: sig
    real, optional,               intent(in)    :: flux_enu

    integer :: charge_out, pion_charge_out
    real :: costheta_max, costheta_min
    real :: elepton_max, elepton_min

    integer :: intprec,l,n2
    real :: costhetaint
    real :: sigmaximum,sigrd
    real,dimension(:),allocatable :: yy,xx
    logical :: flagOK


    if (initFlag) then
       call readInput(dSigmadElepton)

       if (present(flux_enu)) then
          write(*,*) 'enu from experimental flux', flux_enu
       else
          write(*,'(a,F12.4)') 'enu= ', enu
       end if
       write(*,'(a,F12.4)') 'elepton= ', elepton
    end if


    ! set default output
    call setToDefault(OutPart)
    sig=0.

    if (present(flux_enu)) enu=flux_enu

    if (raiseFlag) then
       elepton=elepton+delta_elepton
       write(*,*) 'Elepton is raised by ...', delta_elepton, ' to ', elepton
    end if
    raiseVal=Elepton

    ! set neutrino 4-vector: incoming neutrino in z-direction
    ! plus threshold check for enu

    call eNev_init_nuStep2(eNev,enu,IP,flagOK)
    if (.not.flagOK) return

    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return

    call getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    if (.not.checkLimits(elepton_min,elepton_max,elepton)) return

    call getCosthetaLimits(costheta_min,costheta_max,enu)

    ! set integral precision
    intprec=integralPrecision
    if (IP.eq.nucleon) intprec=integralPrecisionQE  ! due to small Breit-Wigner

    allocate (xx(20*intprec))
    allocate (yy(20*intprec))

    call sg20r(costheta_min,costheta_max,intprec,xx,n2)
    sigmaximum=0.0
    yy = 0.0

    do l=1,n2
       costhetaint=xx(l)
       call eNev_init_nuStep3a(eNev,elepton,costhetaint,flagOK)
       if (.not.flagOK) cycle

       call XsecdCosthetadElepton(eNev,IP,OutPart, yy(l))
       if (yy(l).gt.sigmaximum) sigmaximum=yy(l)
    end do

    call rg20r(costheta_min,costheta_max,intprec,yy,sig)

    ! set kinematics
    if (sigmaximum.gt.0.) then
       do
          costhetaint=costheta_min+rn()*(costheta_max-costheta_min)

          call eNev_init_nuStep3a(eNev,elepton,costhetaint,flagOK)
          if (.not.flagOK) cycle

          call XsecdCosthetadElepton(eNev,IP,OutPart, sigrd)

          if (sigmaximum*rn().le.sigrd) exit
       end do
    end if

    deallocate(xx,yy)

  end subroutine Xsec_dSigmadElepton

  !****************************************************************************
  !****s* neutrinoSigma/getCosthetaLimits
  ! NAME
  ! subroutine getCosthetaLimits(costheta_min,costheta_max)
  !
  ! PURPOSE
  ! This subroutine returns the limits in costheta.
  !
  ! OUTPUT
  ! * real             :: costheta_min
  ! * real             :: costheta_max
  !****************************************************************************
  subroutine getCosthetaLimits(costheta_min,costheta_max,Enu)

    real, intent(out) :: costheta_min,costheta_max
    real, intent(in), optional :: Enu
    costheta_min=-1.
    costheta_max=1.
    if (present(Enu)) then
       ! practical cut for high neutrino energies -
       ! the formula is pure "educated guess"
       ! the reason for cut is physical --- above Q^2>4 everything dies,
       ! so for high neutrino energies only forward scattering is possible
       if (Enu.ge.2) costheta_min=1.-4./Enu/Enu
    end if

  end subroutine getCosthetaLimits


  !****************************************************************************
  !****s* neutrinoSigma/getEleptonLimits
  ! NAME
  ! subroutine getEleptonLimits(IP,eNev,elepton_min,elepton_max)
  !
  ! PURPOSE
  ! This subroutine returns the limits in elepton depending on the neutrino
  ! energy.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * integer                      :: IP   -- ID of outgoing hadron/process
  !
  ! OUTPUT
  ! * real             :: elepton_min
  ! * real             :: elepton_max
  !****************************************************************************
  subroutine getEleptonLimits(IP,eNev,elepton_min,elepton_max)
    use idtable, only: nucleon,delta
    use constants, only: mPi

    integer, intent(in) :: IP
    real, intent(out) :: elepton_min,elepton_max
    type(electronNucleon_event), intent(in) :: eNev

    real :: enu, ml_out

    enu = eNev%lepton_in%mom(0)
    ml_out = eNev%lepton_out%mass

    elepton_min=ml_out
    if (IP.eq.nucleon) elepton_max=enu
    if (IP.ge.delta) elepton_max=enu-mPi
    if (IP.eq.35) elepton_max=enu ! IP=35 = 2p2h process
  end subroutine getEleptonLimits


  !****************************************************************************
  !****s* neutrinoSigma/getQ2Limits
  ! NAME
  ! subroutine getQ2Limits(eNev,elepton,Q2_min,Q2_max)
  !
  ! PURPOSE
  ! This subroutine returns the limits in Q2 depending on the neutrino energy
  ! and the energy of the outgoing lepton.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * real                         :: elepton -- energy of outgoing lepton
  !
  ! OUTPUT
  ! * real             :: Q2_min
  ! * real             :: Q2_max
  !****************************************************************************
  subroutine getQ2Limits(eNev,elepton,Q2_min,Q2_max)

    type(electronNucleon_event), intent(in) :: eNev
    real, intent(in) :: elepton
    real, intent(out) :: Q2_min,Q2_max

    real :: enu,ml_out

    enu = eNev%lepton_in%mom(0)
    ml_out = eNev%lepton_out%mass

    !! Attention! you should have checked before, that elepton>ml_out
    !! is guaranteed. Therefore we could skip the "sqrt(max(" stuff !!

    Q2_max=-ml_out**2+2.*enu*(elepton+sqrt(max((elepton**2-ml_out**2),0.)))
    Q2_min=-ml_out**2+2.*enu*(elepton-sqrt(max((elepton**2-ml_out**2),0.)))
  end subroutine getQ2Limits

  !****************************************************************************
  !****f* neutrinoSigma/checkLimits
  ! NAME
  ! logical function checkLimits(X_min,X_max,X)
  !
  ! PURPOSE
  ! This function checks whether X is out of bounds.
  ! If not, checkLimits=.true., if yes checkQ2Limits=.false.
  !
  ! INPUTS
  ! * real  :: X_min,X_max,X
  !****************************************************************************
  pure logical function checkLimits(X_min,X_max,X)

    real, intent(in) :: X_min,X_max,X

    checkLimits = .false.
    if (X.lt.X_min) return
    if (X.gt.X_max) return
    checkLimits = .true.

  end function checkLimits

  !****************************************************************************
  !****s* neutrinoSigma/get_sigma_namelist
  ! NAME
  ! subroutine get_sigma_namelist(XsectionMode, Genu, Gelepton, Gcostheta)
  !
  ! PURPOSE
  ! This function returns global parameters. It first checks, whether it has
  ! to read the jobcard.
  ! Please note, that all arguments are optional and the name has to be given
  ! explicitely.
  !
  ! INPUTS
  ! * integer, optional :: XsectionMode
  !
  ! OUTPUT
  ! * real,optional :: Genu
  ! * real,optional :: Gcostheta,
  ! * real,optional :: Gelepton
  !****************************************************************************
  subroutine get_sigma_namelist(XsectionMode, Genu, Gelepton, Gcostheta)

    integer, optional, intent(in)  :: XsectionMode
    real,    optional, intent(out) :: Genu, Gelepton, Gcostheta

    if (present(XsectionMode)) then
       if (initflag) then
          call readInput(MOD(XsectionMode,10),(XsectionMode>10))
       end if
    end if

    if (present(Genu)) Genu=enu
    !    if (present(Gdelta_enu)) Gdelta_enu= delta_enu
    !    if (present(GQ2)) GQ2 =Q2
    !    if (present(Gdelta_Q2)) Gdelta_Q2 =delta_Q2
    if (present(Gcostheta)) Gcostheta =costheta
    !    if (present(Gdelta_costheta)) Gdelta_costheta =delta_costheta
    if (present(Gelepton)) Gelepton =elepton
    !    if (present(Gdelta_elepton)) Gdelta_elepton =delta_elepton



  end subroutine get_sigma_namelist

end module neutrinoSigma
