!******************************************************************************
!****m* /neutrinoXsection
! NAME
! module neutrinoXsection
!
! PURPOSE
! This module calculates various neutrino cross sections, depending on the
! choice of the user:
! * integrated or differential or double differential cross sections
!   (according to the value of nuXsectionmode, the corresponding
!   subroutines are called in the neutrino init file)
! * EM, CC or NC (according to the value of process_ID given as input to the
!   corresponding subroutines)
! * muon, electron or tau flavor (according to the value of flavor_ID given
!   as input to the corresponding subroutines)
!******************************************************************************
module neutrinoXsection

  use particleDefinition
  use eN_eventDefinition
  use callstack, only: traceback
  use constants, only: hbarc

  implicit none
  private


  !****************************************************************************
  !****g* neutrinoXsection/debugflag
  ! SOURCE
  logical, parameter :: debugflag=.false.
  ! PURPOSE
  ! to switch on/off debug information
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoXsection/nuclear_phasespace
  ! SOURCE
  logical, parameter :: nuclear_phasespace=.true.
  ! PURPOSE
  ! to change between different phasespace factors.
  ! (change only for debugging purposes):
  ! * .false. = 1/SP(k_in,p_in)  i.e. phasespace factor for each nucleon
  ! * .true.  = 1/k_in(0)*p_in(0) i.e. global phasespace factor for the nucleus
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/singlePiModel
  ! SOURCE
  integer, save :: singlePiModel=1
  ! PURPOSE
  ! to change between different models for the pion nucleon cross section:
  ! * 0 = pi N according to Nieves et al (hep-ph/0701149)
  ! * 1 = MAID-like model
  ! * 2 = Bosted-Christy
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoXsection/invariantMassCut
  ! SOURCE
  real, save :: invariantMassCut=100.
  ! PURPOSE
  ! cut events with invariant Mass above this value (in GeV);
  ! cut pion production from Delta and DIS on Wrec = Sqrt[M^2 + 2*M*nu - Q^2]
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/invariantMassCut_BG
  ! SOURCE
  real, save :: invariantMassCut_BG=100.
  ! PURPOSE
  ! cut MAID-like background events with invariantMass_BG above this value
  ! (in GeV);
  ! cut 1pi BG on Wrec = Sqrt[M^2 + 2*M*nu - Q^2]
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DIScutW
  ! SOURCE
  real, save :: DIScutW = 3.0
  ! PURPOSE
  ! W-cut for sigmoid onset of DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DIScutwidth
  ! SOURCE
  real, save :: DIScutwidth = 0.2
  ! PURPOSE
  ! width for sigmoid onset of DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/REScutW
  ! SOURCE
  real, save :: REScutW = 2.0
  ! PURPOSE
  ! W-cut for end of 1pi,2pi BGs
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISformfakEM
  ! SOURCE
  integer,save :: DISformfakEM = 2
  ! PURPOSE
  ! Introduce an additional form factor for the DIS cross section, when
  ! processed via a photon:
  ! * 0: no form factor
  ! * 1: Q^2/(mcutDIS^2+Q^2)
  ! * 2: Q^4/(mcutDIS^2+Q^2)^2
  !
  ! In case of electron induced events, we need choose 2 in order to be
  ! compatible with Pythia's electron machinery.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISformfakNCCC
  ! SOURCE
  integer,save :: DISformfakNCCC = 1
  ! PURPOSE
  ! Introduce an additional form factor for the DIS cross section, when
  ! processed via W or Z boson:
  ! * 0: no form factor
  ! * 1: Q^2/(mcutDIS^2+Q^2)
  ! * 2: Q^4/(mcutDIS^2+Q^2)^2
  !
  ! In case of electron induced events, we need choose 2 in order to be
  ! compatible with Pythia's electron machinery.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/mcutDIS
  ! SOURCE
  real, save :: mcutDIS = 0.6
  ! PURPOSE
  ! parameter to control Q^2 dependence of DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISrespectHad
  ! SOURCE
  logical, save :: DISrespectHad = .true.
  ! PURPOSE
  ! Flag to indicate, whether hadronization failures should be respected and
  ! affect the overall DIS cross section
  !
  ! Pythia is run to generate the DIS cross section. But not every of the
  ! generated events may lead to a correct hadronic final state.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/DISdoMSTP23
  ! SOURCE
  logical, save :: DISdoMSTP23 = .true.
  ! PURPOSE
  ! Flag to indicate, whether in Pythia for neutrino-DIS the value MSTP(23)=1
  ! should be used or not
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/new2piBG
  ! SOURCE
  logical, save :: new2piBG = .true.
  ! PURPOSE
  ! Flag to turn on the new treatment of 2pi BG for electrons and neutrinos
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoXsection/indBG
  ! SOURCE
  integer, save :: indBG = 3
  ! PURPOSE
  ! Index to choose Bloom-Gilman like BG parametrization
  ! 1 : original Bloom Gilman
  ! 2 : Niculescu fit
  ! 3 : nonresonant BG fit from Christy-Bosted
  !****************************************************************************

  logical, save :: initflag=.true.

  real, parameter :: unit_factor=hbarc**2 * 10.**12 ! factor for sigma unit conversion to 10^(-38) cm^2, unit : GeV^2*fm^2 = GeV^2 * 10^(-26)cm^2

  public :: XsecdCosthetadElepton
  public :: SetHadronCharge
!  public :: get_xsection_namelist
  public :: readinput


contains

  !****************************************************************************
  !****s* neutrinoXsection/readinput
  ! NAME
  ! subroutine readinput
  ! PURPOSE
  ! This subroutine reads input out of jobcard.
  !****************************************************************************
  subroutine readinput
    use output, only: Write_ReadingInput, Write_InitStatus
    use neutrino_IDTable
    use neutrinoParms, only: readInput_neutrino, new_eN

    integer :: IOS

    ! dummies for reading:
    integer :: integralPrecision=-99, integralPrecisionQE=-99

    !**************************************************************************
    !****n* neutrinoXsection/nl_neutrinoxsection
    ! NAME
    ! NAMELIST /nl_neutrinoxsection/
    ! PURPOSE
    ! This Namelist includes:
    ! * singlePiModel
    ! * invariantMassCut
    ! * invariantMassCut_BG
    ! * DIScutW
    ! * DIScutwidth
    ! * REScutW
    ! * DISformfakEM
    ! * DISformfakNCCC
    ! * mcutDIS
    ! * DISrespectHad
    ! * DISdoMSTP23
    ! * new2piBG
    ! * indBG
    !**************************************************************************
    NAMELIST /nl_neutrinoxsection/ integralPrecision,&
         & integralPrecisionQE,singlePiModel,&
         & invariantMasscut,invariantMasscut_BG, &
         & DIScutW, DIScutwidth, REScutW, &
         & DISformfakEM, DISformfakNCCC,mcutDIS, &
         & DISrespectHad, DISdoMSTP23,new2piBG,indBG

    if (.not.initFlag) return

    call Write_InitStatus('neutrinoXsection',0)

    call readInput_neutrino()

    ! reset some default values before reading the jobcard:
    if (new_en) then
       singlePiModel = 2
    end if

    call Write_ReadingInput('nl_neutrinoxsection',0)
    rewind(5)
    read(5,nml=nl_neutrinoxsection,IOSTAT=IOS)
    call Write_ReadingInput('nl_neutrinoxsection',0,IOS)

    if (integralPrecision > 0. .or. integralPrecisionQE>0) then
       write(*,*) 'integralPrecision or integralPrecisionQE without function!'
       write(*,*) 'use replacements in module neutrinoSigma.'
       call Traceback()
    end if

    select case (singlePiModel)
    case (0)
       write(*,*) 'pi N cross section (Delta + interfering background) according to Nieves et al'
    case (1)
       write(*,*) 'pi N background (resonances are subtracted) is taken MAID-like'
    case (2)
       write(*,*) 'pi N background from Bosted-Christy parametrization'
    case default
       write(*,*) 'no valid input for singlePiModel -> stop', singlePiModel
       call traceback()
    end select

    write(*,'(a,F12.4)') ' cut         events with invariant masses above', &
         & invariantMassCut
    write(*,'(a,F12.4)') ' cut MAID BG events with invariant masses above', &
         & invariantMassCut_BG

    write(*,'(A,3F8.3)') ' DIS W-cut = ',DIScutW,DIScutwidth
    write(*,'(A,3F8.3)') ' RES W-cut = ',REScutW
    select case (DISformfakEM)
    case (0)
       write(*,'(A)') ' DIS form factor (EM)   : 0,  = 1'
    case (1)
       write(*,'(A)') ' DIS form factor (EM)   : 1,  = Q^2/(mcutDIS^2+Q^2)'
    case (2)
       write(*,'(A)') ' DIS form factor (EM)   : 2,  = Q^4/(mcutDIS^2+Q^2)^2'
    case default
       write(*,*) 'value ',DISformfakEM,' not valid for DISformfakEM!'
       call traceback()
    end select

    select case (DISformfakNCCC)
    case (0)
       write(*,'(A)') ' DIS form factor (CC,NC): 0,  = 1'
    case (1)
       write(*,'(A)') ' DIS form factor (CC,NC): 1,  = Q^2/(mcutDIS^2+Q^2)'
    case (2)
       write(*,'(A)') ' DIS form factor (CC,NC): 2,  = Q^4/(mcutDIS^2+Q^2)^2'
    case default
       write(*,*) 'value ',DISformfakNCCC,' not valid for DISformfakNCCC!'
       call traceback()
    end select

    write(*,*) 'mcutDIS = ', mcutDIS

    write(*,*) 'DIS respect Hadronization: ',DISrespectHad
    write(*,*) 'DIS use MSTP(23)=1:        ',DISdoMSTP23

    write(*,*) 'new2piBG= ',new2piBG
    write(*,*) 'indBG   = ',indBG

    call Write_ReadingInput('nl_neutrinoxsection',1)

    call Write_InitStatus('neutrinoXsection',1)

    initFlag=.false.

  end subroutine readinput

  !****************************************************************************
  !****s* neutrinoXsection/XsecdCosthetadElepton
  ! NAME
  ! subroutine XsecdCosthetadElepton(eNev,IP,OutPart,sig)
  !
  ! PURPOSE
  ! This subroutine is the basic subroutine, which does all the job of
  ! calculating the cross section dSigma/dCost dElepton and generating
  ! a corresponding final state particle vector.
  !
  ! This routine is called by all routines called "Xsec_...".
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * integer                      :: IP   -- ID of outgoing hadron/process
  !
  ! OUTPUT
  ! * type(electronNucleon_event)  :: eNev     -- the Lepton-Nucleon Event
  ! * type(particle), dimension(:) :: OutPart  -- FinalState particles
  ! * real                         :: sig      -- calculated cross section
  !
  ! NOTES
  ! For the outgoing particles OutPart, the following entries have to
  ! be set afterwards:
  ! * OutPart%firstEvent
  ! * OutPart%event(1:2)
  ! * OutPart%pert
  ! * OutPart%vel
  ! * OutPart%offshellPar
  ! * OutPart%perweight
  ! All the other values are set in this routine.
  !
  ! Returned cross section is in 10^-38 cm^2/GeV, except for EM, where the
  ! units are nb/GeV
  !****************************************************************************
  subroutine XsecdCosthetadElepton(eNev,IP,OutPart,sig)

    use constants, only: pi, mN, mPi, twopi
    use minkowski, only: SP, op_ang, abs4
    use NeutrinoMatrixElement
    use spectralFunc, only:specfunc
    use ParticleProperties, only: hadron
    use leptonicID
    use idTable, only: nucleon,pion
    use neutrino_IDTable
    use singlePionProductionMAIDlike
    use singlePionProductionNHVlike, only: Nieves1piN_elepton_ct
    use Coll_nuN
    use lepton2p2h, only: lepton2p2h_DoQE,lepton2p2h_DoDelta
    use neutrino2piBack, only: DoNu2piBack,BGstruct,BGDualXS,BC2piBG
    use distributions, only: sigmoid
    use eN_event, only: nuclearFluxFactor_correction
    use eventGenerator_eN_lowEnergy, only: init_2Pi
    use ParamEP_BC, only: sigma_lNBC
    use lorentzTRafo
    use electronPionProduction_kine, only: getKinematics_eN
    use degRad_conversion, only: degrees
    use random, only: rn, rnCos
    use neutrinoParms, only: new_eN, new_eNres,normRes, normpiBG, normBC, &
        & Wtrans,twopiBG_on, ME_version

    integer,             intent(in)  :: IP
    real,                intent(out) :: sig
    type(electronNucleon_event), intent(inout) :: eNev
    type(particle), dimension(:), intent(out) :: OutPart ! FinalState particles

    ! for internal purposes:

    integer :: process_ID
    real, dimension(0:3) :: k_in,k_out,p_in,p_out, pion_momentum_out
    real                 :: mass_out, ml_out
    integer              :: charge_out, pion_charge_out
    real, dimension(1:3) :: position
    integer              :: charge_in

    real :: kin_factor,invMassSQ,Wrec,Q2,nu,Wfree,W,QE_corrfactor,qvec2
    real :: phi_pi,theta_pi
    logical :: tworoots
    real, dimension (1:8) :: BCIP = (/1,2,4,7,16,10,3,31/) ! resonances in Bosted-Christy analysis

    logical :: success,Q2cutoff = .FALSE.

    real :: plep
    real :: sigBC, sigDIS,res_sigma,nr_sigma,SigSIS

    if (initFlag) call Traceback('not yet initialized')

    !set default output
    call setToDefault(OutPart)
    sig= 0.

    ! The charges are maybe already calculated earlier, but we do it here again
    if (.not.SetHadronCharge(eNev,IP,charge_out,pion_charge_out)) return


    ! set some abbreviations:
    process_ID = eNev%idProcess   ! 1: electron, 2: CC neutrino, 3: NC neutrino

    k_in  = eNev%lepton_in%mom
    k_out = eNev%lepton_out%mom
    nu = k_in(0) - k_out(0)
    Q2 = eNev%Q2
    p_in  = eNev%nucleon%mom
    p_out = p_in+k_in-k_out               ! momentum of outgoing particle

    position = eNev%nucleon%pos
    charge_in = eNev%nucleon%charge
    ml_out = eNev%lepton_out%mass
    plep = sqrt(max((eNev%lepton_out%mom(0)**2-ml_out**2),0.))

 !  Now invariant mass cut on reconstructed Wrec = sqrt(mN**2 + 2*mN*nu - Q2)
 !   = invariant mass in incoming channel for free nucleon  at rest
 !  used for all production mechanisms
    W = eNev%W
    Wrec = eNev%W_rec
    Wfree = eNev%W_free
    if (Wrec .gt. invariantMasscut) return

    if (Wrec < 0.7) return

    select case (IP)     ! IP = reaction type: QE, Delta, N*1,N*2, ......
                         ! reaction types defined in initNeutrino, line ~ 500

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! QE and RES production
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (:31)

       ! cutting the resonances off at high W
       if (Wrec .gt. REScutW + 0.1) return

       invMassSQ=SP(p_out,p_out)
       if (invMassSQ.le.0.) return !reaction not possible

       if (-SP(p_out-p_in,p_out-p_in).lt.0.) then
          if (debugflag) write(*,*) 'neutrinoXsec, -SP(pf-pi,pf-pi).lt.0.', &
               & -SP(p_out-p_in,p_out-p_in), &
               & -SP(k_in-k_out,k_in-k_out), IP,charge_in,k_in,k_out
          return !reaction not possible
       end if

       kin_factor=plep/(pi*16.)/(SP(k_in,p_in))*specfunc(IP,charge_out,p_out,position,mass_out)

       if (nuclear_phasespace) kin_factor=kin_factor*nuclearFluxFactor_correction(p_in,k_in)

       if (mass_out.eq.0.0) return ! failure in specfunc

       !avoid the production of res. below threshold:
       if (mass_out.le.hadron(IP)%minMass) then
          if (debugflag) write(*,*) 'less than threshold -> sig=0'
          return
       end if

       ! first only QE:
       if (IP == 1 ) then   !QE done in the old (Leitner) way, for all versions

          sig=kin_factor*nuMaEl(process_ID,IP,charge_in,k_in,k_out,p_in,p_out,mass_out,position)
          !multiply cross section by unit_factor to obtain results in 10^(-38) cm^2
          sig=unit_factor*sig

          if (new_eN .and. ME_version==6) then
             ! multiply with qvec^2 QE quenching factor for electrons
             qvec2 = Q2 + nu**2
             QE_corrfactor = 1.0/(1. + 0.2445*qvec2)
             sig = sig * QE_corrfactor
          end if


       ! now also higher resonances:
       else if (IP > 1) then
          if (.not. new_eN ) then

             sig=kin_factor*nuMaEl(process_ID,IP,charge_in,k_in,k_out,p_in,p_out,mass_out,position)
             !multiply cross section by unit_factor to obtain results in 10^(-38) cm^2
             sig=unit_factor*sig

          else if (new_eNres) then  ! new treatment of resonances only for neutrinos

             call sigma_lNBC(process_ID,IP,Q2,charge_in,k_in,k_out,p_in,p_out,res_sigma,nr_sigma,position)
             ! lNBC gives lepton cross section, based on Bosted-Cristy parametrization
             ! cross section parametrization for eN collisions in Lab frame
             sig = unit_factor * res_sigma
             mass_out = abs4(p_out)

          else    ! use Leitner resonances

             if (.not. ANY(BCIP==IP)) return
             sig=kin_factor*nuMaEl(process_ID,IP,charge_in,k_in,k_out,p_in,p_out,mass_out,position)
             !multiply cross section by unit_factor to obtain results in 10^(-38) cm^2
             sig=unit_factor*sig

          end if
       end if

       if (IP .ge. 3 .and. abs(process_ID) > 1) then
          sig = sig * normRES
       end if

       ! smooth switching off the resonances,
       sig = sig * Sigmoid(W,REScutW,-0.05)

       ! Setting the outgoing resonance properties:

       OutPart%pos(1)=position(1)
       OutPart%pos(2)=position(2)
       OutPart%pos(3)=position(3)

       OutPart%formTime = -999

       OutPart(1)%mass = mass_out
       OutPart(1)%ID   = IP
       OutPart(1)%Charge = Charge_Out
       OutPart(1)%anti=.false.
       OutPart(1)%mom = p_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 2p2h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (chQE2p2h)

       if (Wrec > 1.5) return

       call lepton2p2h_DoQE(eNev,outPart,sig)
       ! returned XS is dsigma/dE'dOmega in mb/GeV/sr/A

       sig = sig * twopi * 1.e11 ! cross section dsigma/dE'dcost in 10^-38 cm^2

    case (chDelta2p2h)

       call lepton2p2h_DoDelta(eNev,outPart,sig)
       ! returned XS is dsigma/dE'dOmega in mb/GeV/A
       sig = sig * twopi * 1.e11 ! cross section dsigma/dE'dcost in 10^-38 cm^2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! single pi N background
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Single and Double pi N backgrounds are obtained only for the W region between 1.2 and REScutW GeV.
! Their cross section sum is given by the background in the Bosted-Christy parametrization

    case (chOnePionN,chOnePionP)

       ! cutting the background off at high and low W
       if (W < 0.8 .or. W .gt. REScutW + 0.1) return

       ! invariant mass cut for pion BG contribution
       if (Wrec .gt. invariantMasscut_BG) return

       ! remarks: the routines to be called here should return:
       ! - dsigma/(dElepton dcostheta_lepton) in units of GeV**4 for CC,EM,NC
       ! - p_out, mass_out
       ! - pion_momentum_out

       select case (singlePiModel)

       case (0) !HNV model

          !if (debugflag) write(*,*) 'singlePiModel nieves started'
          call Nieves1piN_elepton_ct( process_ID,k_in,k_out,.true., &
               & p_in,position,charge_in, &
               & p_out,charge_out,pion_momentum_out,pion_charge_out,sig)


       case (1) !MAID like

          !check on outgoing particles
          invMassSQ=SP(p_out,p_out)
          if (invMassSQ.le.0.) return !reaction not possible

          sig=MAIDlike_singlePi(eNev,charge_out,pion_charge_out, &
               & p_out, pion_momentum_out,nuclear_phasespace)
          !multiply cross section by factor to obtain results in 10^(-38) cm^2
          sig=unit_factor*sig


       case (2) !Bosted-Christy background

          call sigma_lNBC(process_ID,IP,Q2,charge_in,k_in,k_out,p_in,p_out,res_sigma,nr_sigma,position)
          ! lNBC gives Bosted-Cristy cross section parametrization for eN collisions in Lab frame
          
      !  if(nr_sigma > 200.0) write(*,*) '598 Process_ID=',Process_ID,'IP=',IP,'W=',Wrec,'Q2=',Q2,'nu=',nu,'nr_sigma=',nr_sigma        

          ! Now set kinematics of final pi-N state:
          ! Assume here no pion potential.
          ! If pion potential is desired, then use coding from electronPionProduction_medium_eN
          theta_pi = degrees(rnCos())
          phi_pi = degrees(rn()*2*pi)
          call getkinematics_eN(eNev,pion_Charge_out,charge_out,phi_pi,theta_pi,&
               & pion_momentum_out,p_out,twoRoots,success,pionNucleonSystem=2)

          if (.not.success) then
             return
          end if

            select case (process_ID)

            case (1,-1)
               sig = nr_sigma * 0.5  ! factor 0.5 because BC background contains both final states in 1pi decay
            case (2)
               if (charge_in == 1) sig = nr_sigma
               if (charge_in == 0) sig = nr_sigma * 0.5

            case (-2)
               if (charge_in == 1) sig = nr_sigma * 0.5
               if (charge_in == 0) sig = nr_sigma

            case (3,-3)
               if ( IP == chOnePionN ) sig = sig * 1.15
               if ( IP == chOnePionP ) sig = sig * 0.85
               sig = 0.5 * sig       ! reduction for NC because less dofs of outgoing lepton
            end select

          sig = unit_factor * sig      ! convert to proper units

       case default

          write(*,*) 'wrong choice for singlePiModel -> stop',singlePiModel
          call traceback()

       end select

       mass_out = mN

       ! smooth switching off the 1pi background to allow for 2pi contribution
       ! 0.13 chosen such that 0.13/sqrt(W-1.25) = 1
       if (new_en .and. W > 1.267) then
          sig = sig * 1./3. *(1 + 0.13/Sqrt(W - 1.25)*2.)
       end if

       ! smooth switching off the pi background for transition to PYTHIA:
       sig = sig * Sigmoid(W,REScutW,-0.05)

       sig = sig * NormpiBG
       if ( IP == chOnePionN ) sig = sig * 1.15
       if ( IP == chOnePionP ) sig = sig * 0.85

       ! Setting the outgoing particles:

       OutPart%pos(1)=position(1)
       OutPart%pos(2)=position(2)
       OutPart%pos(3)=position(3)

       OutPart%formTime = -999

       OutPart(1:2)%ID = (/nucleon,pion/)
       OutPart(1:2)%Mass=(/mass_out,mPi/)
       OutPart(1:2)%Charge=(/charge_out,pion_charge_out/)
       OutPart(1:2)%anti=.false.

       OutPart(1)%mom = p_out
       OutPart(2)%mom = pion_momentum_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! two pion background
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (chTwoPion)

       ! cutting the background off at high W
       if (W < 1.26 .or.  W > REScutW + 0.1) return
       if (Q2 > 10.) return

       invMassSQ=SP(p_out,p_out)
       if (invMassSQ.le.0.) return ! reaction not possible

       if (new_eN) then

          call BC2piBG(eNev,outPart,sig)
          sig = unit_factor * sig
          if (abs(process_iD) == 3) sig = 0.5 * sig    ! reduction for NC because less dofs of outgoing lepton
          sig = sig * NormpiBG

       else

          select case (abs(process_ID))
          case (1)                 ! Electrons
             if (new2piBG) then
                call DoNu2piBack(eNev,outPart,sig,new2piBG,indBG)
                ! from DoNu2piBack X-section d2sigma/dE' dOmega in mb.
                ! now converted to d2sigma/dE' dcost in 10^(-38) cm^2/GeV
                sig = sig * twopi * 1.e11
                ! The following check rejects events where in-medium correction
                ! did not find a solution
                if (abs(sig) < 1.e-16) return

             else
                call init_2Pi(eNev,OutPart,sig, .true.)
                sig = sig * twopi * 1e11 ! cross section dsigma/dE'dcost in 10^(-38) cm^2/GeV
             end if

          case (2:3)                 ! Neutrinos
             call DoNu2piBack(eNev,outPart,sig,new2piBG,indBG)
             ! from DoNu2piBack X-section d2sigma/dE' dOmega in mub/GeV.
             ! now converted to d2sigma/dE' dcost in 10^(-38) cm^2/GeV

             ! The following check rejects events where the in-medium correction
             ! did not find a solution
             if (abs(sig) < 1.e-16) return

             sig = unit_factor * sig * twopi
             sig = sig * NormpiBG
          end select

       end if

       OutPart%pos(1)=position(1)
       OutPart%pos(2)=position(2)
       OutPart%pos(3)=position(3)

       if( new_eN .and. W > 1.267) then
         sig = sig * (1. - 0.13/Sqrt(W - 1.25)) * 2./3.
       end if


       ! smooth switching off the 2pi background for transition to PYTHIA:
       sig = sig * Sigmoid(W,REScutW,-0.05)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (chDIS)

       ! Here the transition from RES to DIS is made such that first up to
       ! DIScutW the total cross sectio is given by the background contribution
       ! of the Bosted-Christy fit while the final state particles are determined
       ! by Pythia. From DIScutW onwards also the cross section is obtained from PYTHIA

       if (W .lt. REScutW - 0.1) return

       !check on outgoing particles
       invMassSQ=SP(p_out,p_out)
       if (invMassSQ.le.0.) return !reaction not possible

       OutPart%pos(1)=position(1)
       OutPart%pos(2)=position(2)
       OutPart%pos(3)=position(3)
       sigDIS = 0.0
       sigBC = 0.0
       sig = 0.0

       ! first get final state particles from PYTHIA

       call DoColl_nuN_Py(eNev,OutPart,success,sigDIS,DISrespectHad,DISdoMSTP23)
       if (.not.success) then
          sigDIS=0.0
          sig = 0.0
        ! write(*,*) 'W =', W, 'PYTHIA not success'
          return
       end if
       sigDIS = sigDIS * 1e11 ! cross section in 10^-38 cm^2/GeV (from mb=10^{-27} to 10^{-38})

       ! If W is smaller than DIScutW then get the total cross section from the BC parametrization

       if (W < DIScutW + DIScutwidth)  then ! overwrite the X-section with the BC BG X-section
          call sigma_lNBC(process_ID,34,Q2,charge_in,k_in,k_out,p_in,p_out,res_sigma,sigSIS,position)
          sigSIS = unit_factor * sigSIS
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (Q2cutoff) then
          select case (eNev%idProcess)
          case (-1,1) !=== EM
             select case (DISformfakEM)
             case (1)
                sigDIS = sigDIS * eNev%Q2/(eNev%Q2+mcutDIS**2)
             case (2)
                sigDIS = sigDIS * (eNev%Q2/(eNev%Q2+mcutDIS**2))**2
             end select
          case default !=== CC, NC
             select case (DISformfakNCCC)
             case (1)
                sigDIS = sigDIS * eNev%Q2/(eNev%Q2+mcutDIS**2)
             case (2)
                sigDIS = sigDIS * (eNev%Q2/(eNev%Q2+mcutDIS**2))**2
             end select
          end select
       end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       ! correction for the nuclear phase space
       if (nuclear_phasespace) sigDIS=sigDIS*nuclearFluxFactor_correction(p_in,k_in)

       if (new_eN) then
          sig = sigSIS * Sigmoid(W,DIScutW,-DIScutwidth) + sigDIS * Sigmoid(W,DIScutW,DIScutwidth)     !smooth transition from SIS to DIS

       else

          sig = sigDIS
          if (Wrec < Wtrans + 0.4) then
             !for W < Wtrans GeV parametrization for BG from Christy-Bosted is used
             call BGDualXS(eNev,sigBC,indBG)
             ! cross section comes in millibarns
             sigBC = sigBC * twopi *1.e11       ! X-section in 10^(-38) cm^2/GeV
             ! correction for the nuclear phase space
             if (nuclear_phasespace) sigBC=sigBC*nuclearFluxFactor_correction(p_in,k_in)
             if (abs(process_Id) > 1) sigBC = sigBC * normBC ! tune for neutrinos
             sig = SigBC * Sigmoid(Wrec,Wtrans,-0.1) + sigDIS* Sigmoid(Wrec,Wtrans,0.1)
             sig = sig * Sigmoid(Wrec,DIScutW,DIScutwidth)
          end if
       end if

       ! smooth switching on the transition to PYTHIA:
       sig = sig * Sigmoid(W,REScutW,0.05)

    case default
       write(*,*) 'strange IP, STOP:',IP
       call traceback()

    end select


    if (process_ID.eq.EM .or. process_ID.eq.antiEM) sig=1e-5*sig    !cross section in nanobarn for el-m reactions


  end subroutine XsecdCosthetadElepton



  !****************************************************************************
  !****f* neutrinoXsection/SetHadronCharge
  ! NAME
  ! logical function SetHadronCharge(eNev,IP,Q_R,Q_pi)
  !
  ! PURPOSE
  ! This function sets the charge of the outgoing hadrons depending on the
  ! reaction process. If successful, SetHadronCharge=.true., if not (when
  ! the reaction is not possible) SetHadronCharge=.false.
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eNev -- the Lepton-Nucleon Event
  ! * integer                      :: IP   -- ID of reaction
  !
  ! OUTPUT
  ! * integer                      :: Q_R  -- charge of outgoing baryon
  ! * integer                      :: Q_pi -- charge of outgoing pion
  ! * function value -- indicates possible failures
  !****************************************************************************
  logical function SetHadronCharge(eNev,IP,Q_R,Q_pi)
    use neutrino_IDTable
    use leptonicID
    use ParticleProperties, only: hadron

    type(electronNucleon_event), intent(in) :: eNev
    integer, intent(in) :: IP
    integer, intent(out) :: Q_R, Q_pi

    !set default output
    Q_R=0
    Q_pi=0

    SetHadronCharge=.true.

    select case (IP)

    case (1:31)
       !=======================================================================
       !===== QE and RES production
       !=======================================================================
       Q_pi = 0

       select case (eNev%idProcess)
       case default
          ! == EM, NC, antiEM, antiNC
          Q_R = eNev%nucleon%charge
       case (2)
          ! == CC
          Q_R = eNev%nucleon%charge + 1
          if ((Q_R.eq.2).and.(hadron(IP)%isoSpinTimes2.ne.3)) &
               & SetHadronCharge=.false. ! reaction not possible
       case (-2)
          ! == antiCC
          Q_R = eNev%nucleon%charge - 1
          if ((Q_R.eq.-1).and.(hadron(IP)%isoSpinTimes2.ne.3)) &
               & SetHadronCharge=.false. ! reaction not possible

       end select

    case (chOnePionN)
       !=======================================================================
       !===== pi + n production
       !=======================================================================

       !CC:    nu + n -> l- + pi+ + n  ;  nu~ + n -> l+  + pi- + n
       !CC:    nu + p ->  ---          ;  nu~ + p -> l+  + pi0 + n
       !NC/EM: nu + n -> nu + pi0 + n  ;  nu~ + n -> nu~ + pi0 + n
       !NC/EM: nu + p -> nu + pi+ + n  ;  nu~ + p -> nu~ + pi+ + n

       Q_R = 0
       select case (eNev%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          Q_pi = eNev%nucleon%charge
       case (2)    ! == CC
          Q_pi = eNev%nucleon%charge+1
          if (Q_pi.eq.2) SetHadronCharge=.false. ! reaction not possible
       case (-2)   ! == antiCC
          Q_pi = eNev%nucleon%charge-1
       end select

    case (chOnePionP)
       !=======================================================================
       !===== pi + p production
       !=======================================================================

       !CC:    nu + n -> l- + pi0 + p  ;  nu~ + n ->  ---
       !CC:    nu + p -> l- + pi+ + p  ;  nu~ + p -> l+  + pi- + p
       !NC/EM: nu + n -> nu + pi- + p  ;  nu~ + n -> nu~ + pi- + p
       !NC/EM: nu + p -> nu + pi0 + p  ;  nu~ + p -> nu~ + pi0 + p

       Q_R = 1
       select case (eNev%idProcess)
       case default ! == EM, NC, antiEM, antiNC
          Q_pi = eNev%nucleon%charge-1
       case (2)    ! == CC
          Q_pi = eNev%nucleon%charge
       case (-2)   ! == antiCC
          Q_pi = eNev%nucleon%charge-2
          if (Q_pi.eq.-2) SetHadronCharge=.false. ! reaction not possible
       end select

    case (chDIS,chQE2p2h,chDelta2p2h)
       !=======================================================================
       !===== DIS, 2p2h
       !=======================================================================
       !... everything as the defaults

    case (chTwoPion)
       !=======================================================================
       !===== 2 pion backround
       !=======================================================================
       !... everything as the defaults

    case default
       !=======================================================================
       !===== DEFAULT
       !=======================================================================
       write(*,*) 'wrong IP in SetHadronCharge:',IP,' STOP!'
       call traceback()

    end select

  end function SetHadronCharge

  !****************************************************************************
  !****s* neutrinoXsection/get_xsection_namelist
  ! NAME
  ! subroutine get_xsection_namelist(XsectionMode, ...)
  !
  ! PURPOSE
  ! This subroutine returns variables that are set in the according xsection
  ! namelist.
  !
  ! INPUTS:
  ! * integer, optional :: XsectionMode
  !
  ! OUTPUT
  ! * logical,optional :: Gdebugflag
  ! * logical,optional :: Gnuclear_phasespace
  ! * integer,optional :: GsinglePiModel
  ! * integer,optional :: GintegralPrecision
  ! * integer,optional :: GintegralPrecisionQE
  ! * real,optional :: Genu
  ! * real,optional :: Gdelta_enu
  ! * real,optional :: GQs
  ! * real,optional :: Gdelta_Qs
  ! * real,optional :: Gcostheta,
  ! * real,optional :: Gdelta_costheta
  ! * real,optional :: Gelepton
  ! * real,optional :: Gdelta_elepton
  ! * real,optional :: GinvariantMasscut
  !
  !****************************************************************************
  subroutine get_xsection_namelist(XsectionMode,Gdebugflag,Gnuclear_phasespace,&
       GsinglePiModel,GintegralPrecision,GintegralPrecisionQE, &
       Genu,Gdelta_enu,GQs,Gdelta_Qs,Gcostheta,Gdelta_costheta, &
       Gelepton,Gdelta_elepton,GinvariantMasscut)

    integer,optional,intent(in) ::  XsectionMode
    logical,optional,intent(out) :: Gdebugflag,Gnuclear_phasespace
    integer,optional,intent(out) :: GsinglePiModel,GintegralPrecision,GintegralPrecisionQE
    real,optional,intent(out) :: Genu,Gdelta_enu,GQs,Gdelta_Qs,Gcostheta,Gdelta_costheta,Gelepton,Gdelta_elepton,GinvariantMasscut

    if (present(XsectionMode)) then
       if (initflag) then
!!!!          call readInput(MOD(XsectionMode,10),(XsectionMode>10))
          call Traceback("not yet implemented anymore.")
          initflag = .false.
       end if
    end if

    if (present(Gdebugflag)) Gdebugflag=debugflag
    if (present(Gnuclear_phasespace)) Gnuclear_phasespace=nuclear_phasespace
    if (present(GsinglePiModel)) GsinglePiModel=singlePiModel
    if (present(GinvariantMasscut)) GinvariantMasscut =invariantMasscut

!    if (present(GintegralPrecision)) GintegralPrecision =integralPrecision
!    if (present(GintegralPrecisionQE)) GintegralPrecisionQE =integralPrecisionQE
!    if (present(Genu)) Genu=enu
!    if (present(Gdelta_enu)) Gdelta_enu= delta_enu
!    if (present(GQs)) GQs =Qs
!    if (present(Gdelta_Qs)) Gdelta_Qs =delta_Qs
!    if (present(Gcostheta)) Gcostheta =costheta
!    if (present(Gdelta_costheta)) Gdelta_costheta =delta_costheta
!    if (present(Gelepton)) Gelepton =elepton
!    if (present(Gdelta_elepton)) Gdelta_elepton =delta_elepton

    if ( present(GintegralPrecision).or. &
         present(GintegralPrecisionQE).or. &
         present(Genu).or. &
         present(Gdelta_enu).or. &
         present(GQs).or. &
         present(Gdelta_Qs).or. &
         present(Gcostheta).or. &
         present(Gdelta_costheta).or. &
         present(Gelepton).or. &
         present(Gdelta_elepton)) then
       call Traceback("not yet implemented anymore.")
    end if

  end subroutine get_xsection_namelist


end module neutrinoXsection
