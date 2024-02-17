!******************************************************************************
!****m* /resProd_lepton
! NAME
! module resProd_lepton
! PURPOSE
! * Evaluates cross sections for gamma N -> R
! * For details see the notes about this in the work of Oliver Buss
!******************************************************************************
module resProd_lepton
  implicit none
  private

  public :: dSigmadOmega_fdE_f_resProd_eN
  public :: dSdO_fdE_fdO_k_med_res_eN
  public :: sigma_resProd
  public :: dSdOmega_k_med_res
!  public :: sigma_pipi_res
  public :: sigma_pipi_res_vac
  public :: sigma_barMes_res_vac

  !****************************************************************************
  !****g* resProd_lepton/debug
  ! SOURCE
  logical, parameter :: debug=.false.
  ! PURPOSE
  !
  !****************************************************************************

  logical, save :: initFlag=.true.

contains

  !****************************************************************************
  !****f* resProd_lepton/dSigmadOmega_fdE_f_resProd_eN
  ! NAME
  ! real function dSigmadOmega_fdE_f_resProd_eN(eN,resID,pout,bareMass,processID)
  !
  ! PURPOSE
  ! * Evaluates cross sections for e N -> R
  ! * For details see the notes about this in the work of Oliver Buss
  !
  ! INPUTS
  ! * type(electronNucleon_event)  :: eN   -- underlying electron nucleon event
  ! * integer                      :: resID    -- resonance ID
  ! * integer, optional            :: processID
  !
  ! OUTPUT
  ! * function value -- dsigma/dOmega/dE in units of mb/(GeV sr)
  ! * real                         :: bareMass    -- bare resonance mass
  ! * real, dimension(0:3)         :: pout        -- resonance momentum
  !
  !****************************************************************************
  function dSigmadOmega_fdE_f_resProd_eN(eN,resID,pout,bareMass,processID) result(xSection)
    use leptonicID, only: EM, CC, antiCC
    use spectralFunc, only: specFunc
    use constants, only: GeVSquared_times_mb, pi
    use minkowski, only: contract,SP
    use hadronTensor_ResProd
    use leptonTensor
    use degRad_conversion
    use particleDefinition
    use eN_eventDefinition, only: electronNucleon_event

    type(electronNucleon_event) , intent(in)  :: eN
    integer                     , intent(in)  :: resID
    integer, optional           , intent(in)  :: processID
    real                        , intent(out) :: bareMass
    real, dimension(0:3)        , intent(out) :: pout

    ! local variables
    integer                    :: process_ID
    real, dimension(0:3)       :: pin, lin,q,lout
    complex,dimension(0:3,0:3) :: hadronTensor,leptonTens
    complex                    :: matrixElement_Squared
    real                       :: kinematics, Xsection
    real                       :: spec

    if (initFlag) then
       call readInput
       initFlag=.false.
    end if

    ! ********** Kinematics ******************

    ! Incoming nucleon
    pin=eN%nucleon%mom

    ! Initial lepton: Assume lepton in z-direction
    lin=eN%lepton_in%mom

    !Outgoing lepton
    lout=eN%lepton_out%mom

    q=lin-lout
    pout=pin+q

    ! ********** Cross section ******************

    process_ID=EM

    if (present(processID)) then
       process_ID=processID
       if (processID.eq.CC) process_ID=99
       if(processID.eq.antiCC) process_ID=100
    !  Needed in function hadronTensor_R in hadron Tensor_ResProd
    end if

    leptonTens=l_munu_matrix(lin,lout)
    spec=specFunc(resID,en%nucleon%Charge,pout,en%nucleon%pos,baremass)
    if (hadronTensor_R(pin,pout,resID,en%nucleon%Charge,process_ID,hadronTensor,baremass) ) then
       matrixElement_Squared=Contract( leptonTens, hadronTensor)
    else
       matrixElement_Squared=0.
       xsection=0.
       return
    end if

    if (abs(AIMAG(matrixElement_Squared)).gt.0.0000001) then
       write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', matrixElement_Squared
    end if


    kinematics=1./SP(lin,pin)*sqrt(Dot_Product(lout(1:3),lout(1:3)))/(32.*pi**2)

    Xsection=kinematics  *REAL(matrixElement_Squared)*spec


    ! ********** UNIT CONVERSION ******************

    ! Convert from 1/(sr GeV^3) to mb/(GeV sr)
    ! 1/(GeV^2)=1/(1000 MeV)^2=1/(1000/(197 fm)^2) =0.197**2 *10 mb
    Xsection=Xsection/GeVSquared_times_mb

    ! Convert from mb/(GeV sr) to mb/(MeV sr)
    !  Xsection=Xsection/1000.

  end function dSigmadOmega_fdE_f_resProd_eN



  !****************************************************************************
  !****f* resProd_lepton/sigma_pipi_res_vac
  ! NAME
  ! function sigma_pipi_res_vac(targetN,q) result(sigma)
  !
  ! PURPOSE
  ! Evaluates the resonance contribution to double-pion production in
  ! gamma N -> R -> N 2Pi scattering.
  ! The return value is sigma in [mb]. Converts target nucleon first to
  ! vacuum nucleon!!!
  !
  ! dsigma/dOmega_electron/dE_electron/dOmega_pion
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetN   -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! OUTPUT
  ! * sigma(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production
  ! * sigma(1) -> nucleon piMinus piPlus
  ! * sigma(2) -> nucleon piPlus piNull or nucleon piMinus piNull
  ! * sigma(3) -> nucleon piNull piNull
  ! * sigma(0) Total Xsection into nucleon+2 Pions
  !****************************************************************************
  function sigma_pipi_res_vac(targetN,q) result(sigma_pipi)
    use particleDefinition
    use constants, only: mN

    type(particle), intent(in)        :: targetN
    real, dimension (0:3), intent(in) :: q
    real,dimension(0:3)               :: sigma_pipi

    type(particle) :: targetN_vac     ! Target nucleon, converted to vacuum

    targetN_vac=targetN

    targetN_vac%mass=mN
    targetN_vac%pos=1000.
    targetN_vac%mom(0)=freeEnergy(targetN)

    sigma_pipi=sigma_pipi_res(targetN_vac,q)

  end function Sigma_pipi_res_vac

  !****************************************************************************
  !****f* resProd_lepton/sigma_pipi_res
  ! NAME
  ! function sigma_pipi_res(targetN,q) result(sigma)
  !
  ! PURPOSE
  ! Evaluates the resonance contribution to double-pion production in
  ! gamma N -> R -> N 2Pi scattering.
  ! The return value is sigma in [mb].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetN   -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! OUTPUT
  ! * sigma(0:3) Cross sections for gamma+nucleon->nucleon+2Pi production
  ! * sigma(1) -> nucleon piMinus piPlus
  ! * sigma(2) -> nucleon piPlus piNull or nucleon piMinus piNull
  ! * sigma(3) -> nucleon piNull piNull
  ! * sigma(0) Total Xsection into nucleon+2 Pions
  !****************************************************************************
  function sigma_pipi_res(targetN,q) result(sigma_pipi)
    use particleDefinition
    use idTable, only: nres, nucleon,pion,P11_1440,rho,sigmaMeson,Delta
    use baryonWidth, only: partialWidthBaryon, FullWidthBaryon
    use clebschGordan, only: CG
    use particleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in)       :: targetN
    real, dimension(0:3), intent(in) :: q
    real, dimension(0:3)             :: sigma_pipi

    ! dsigma/dOmega_electron/dE_electron:
    real, dimension(nucleon+1:nucleon+nres,0:3) :: sigma


    real, dimension(1:3,1:4) :: clebsch
    integer :: qPi1, qPi2    ! Charges of outgoing pions
    integer :: resID
    !real :: branchingRatio
    integer :: i
    real :: sigma_tot
    integer :: res_I2, res_Iz2, nuc_Iz2
    real :: mass, baremass, fullWidth
    real, dimension(1:4) :: cc
    real, dimension(1:4) :: gamma_Out

    if (sqrt(Dot_Product(targetN%pos,targetN%pos)).lt.50) then
       write(*,*) 'In medium not yet implemented in Sigma_pipi_res. STOP!'
       write(*,*) 'position=',targetN%pos
       call TRACEBACK()
       stop
    end if

    sigma(:,:) = 0.

    resIDLoop: do resID=nucleon+1,nucleon+nres
!       pf_res=targetN%mom+q

       sigma_tot=sigma_resProd(targetN,resId,q,baremass)
       mass=bareMass

       fullWidth=FullWidthBaryon(resID,mass)
       if (fullWidth.lt.1E-10) then
          sigma(resID,:)=0.
          cycle resIdLoop
       end if

       ! Get gamma N -> pi pi N by summing just over the ...
       ! gamma N -> pion Delta, gamma N -> N rho, gamma N -> N sigma, gamma N -> pion P11_1440
       ! ... channels. All resonances but the P_11(1440) decay fully into N pi.
       ! We correct for the P_11(1440) later.
       gamma_Out(1)=partialwidthBaryon(resID,mass,.false.,pion,Delta)
       gamma_Out(2)=partialwidthBaryon(resID,mass,.false.,rho,nucleon)
       gamma_Out(3)=partialwidthBaryon(resID,mass,.false.,sigmaMeson,nucleon)
       ! Here we make the simplifying assumption that the decay ratio of P11_1440 to NPi
       ! is constant in mass:
       gamma_Out(4)=partialwidthBaryon(resID,mass,.false.,pion,P11_1440)  &
                    & *hadron(P11_1440)%decays(1)

       gamma_Out(:) = gamma_Out(:)/fullWidth

       res_I2 = hadron(resId)%isoSpinTimes2
       res_Iz2 = 2*targetN%charge-1

       clebsch(:,:) = 0.
       loop1: do qPi1=-1,1
          loop2: do qPi2=-1,1
             if (abs(qPi1+qPi2).gt.1) cycle loop2

             cc(:) = 0.
             nuc_Iz2 = res_Iz2 - 2*(qPi1+qPi2)

             ! R ->  pion Delta -> pion pion N
             if (abs(nuc_Iz2).le.1) then
                cc(1) = CG(2,3,res_I2, 2*qPi1,res_Iz2-2*qPi1)**2 &
                     *  CG(2,1,3, 2*qPi2,nuc_Iz2)**2
             end if

             ! R ->  rho N -> pion pion N
             if (abs(nuc_Iz2).le.1) then
                cc(2) = CG(2,1,res_I2, 2*(qPi1+qPi2),nuc_Iz2)**2 &
                     *  CG(2,2,2, 2*qPi2,2*qPi1)**2
             end if

             ! R ->  sigma N -> pion pion N
             if (qPi1+qPi2.eq.0) then
                cc(3) = CG(2,2,0, 2*qPi2,2*qPi1)**2
             end if

             ! R ->  pion P_11(1440) -> pion pion N
             if (abs(nuc_Iz2).le.1 .and. abs(res_Iz2-2*qPi1).le.1) then
                cc(4) = CG(2,1,res_I2, 2*qPi1,res_Iz2-2*qPi1)**2 &
                     *  CG(2,1,1, 2*qPi2,nuc_Iz2)**2
             end if

             if (((qPi1.eq.1).and.(qPi2.eq.-1)) &
                  .or.((qPi1.eq.-1).and.(qPi2.eq.1))) then
                clebsch(1,:) = clebsch(1,:)+cc(:)
             else if (((qPi1.eq.1).and.(qPi2.eq.0)) &
                  .or.((qPi1.eq.0).and.(qPi2.eq.1)) &
                  .or.((qPi1.eq.-1).and.(qPi2.eq.0))  &
                  .or.((qPi1.eq.0).and.(qPi2.eq.-1))) then
                clebsch(2,:) = clebsch(2,:)+cc(:)
             else if ((qPi1.eq.0).and.(qPi2.eq.0)) then
                clebsch(3,:) = clebsch(3,:)+cc(:)
             end if
          end do loop2
       end do loop1

       !write(*,'(5E15.4)') sigma_tot, gamma_Out
       do i=1,4
          sigma(resID,1:3)=sigma(resID,1:3)+sigma_tot*gamma_out(i)*clebsch(1:3,i)
       end do
       sigma(resID,0)=sum(sigma(resID,1:3))

    end do resIDLoop

    sigma_pipi(:) = sum(sigma(:,:),dim=1)

  end function Sigma_pipi_res


  !****************************************************************************
  !****f* resProd_lepton/sigma_barMes_res_vac
  ! NAME
  ! function sigma_barMes_res_vac(targetN,q,IDbar,IDmes) result (sigma)
  !
  ! PURPOSE
  ! Evaluates the resonance contribution of gamma N -> R -> B m^0 scattering
  ! (where X may be a nucleon or Delta, while m^0 is a neutral meson).
  ! The return value is sigma in [mb]. Converts target nucleon first to
  ! vacuum nucleon!!!
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetN -- Target nucleon
  ! * real, dimension (0:3)  :: q       -- Virtual photon 4-momentum
  ! * integer     :: IDbar     -- ID of produced baryon N (nucleon or Delta)
  ! * integer     :: IDmes     -- array containing the IDs of produced mesons m
  ! OUTPUT
  ! * sigma -- Cross sections for gamma N -> R -> B m^0 production, for all mesons you asked for
  !****************************************************************************
  function sigma_barMes_res_vac(targetN,q,IDbar,IDmes) result (sigma)
    use particleDefinition
    use constants, only: mN

    type(particle), intent(in)        :: targetN
    real, dimension (0:3), intent(in) :: q
    integer, intent(in)               :: IDbar, IDmes(:)
    real :: sigma(lbound(IDMes,dim=1):ubound(IDMes,dim=1))

    type(particle) :: targetN_vac  ! Target nucleon, converted to vacuum

    targetN_vac=targetN

    targetN_vac%mass=mN
    targetN_vac%pos=1000.
    targetN_vac%mom(0)=freeEnergy(targetN)

    sigma = sigma_barMes_res(targetN_vac,q,IDbar,IDmes)

  end function Sigma_barMes_res_vac



  !****************************************************************************
  !****f* resProd_lepton/sigma_barMes_res
  ! NAME
  ! function sigma_barMes_res(targetN,q,IDbar,IDmes) result(sigma)
  !
  ! PURPOSE
  ! Evaluates the resonance contribution of gamma N -> R -> B m^0 scattering
  ! (where B may be a nucleon or Delta, while m^0 is a neutral meson).
  ! The return value is sigma in [mb]. The cross section is calculated
  ! separately for all mesons which are passed in IDmes.
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)        :: targetN -- Target nucleon
  ! * real, dimension(0:3)  :: q       -- Virtual photon 4-momentum
  ! * integer :: IDbar     -- ID of produced baryon B (nucleon or Delta)
  ! * integer :: IDmes(:)  -- array containing the IDs of produced mesons m
  ! OUTPUT
  ! * sigma(:) -- Cross sections for gamma N -> R -> B m^0 production,
  !   for all mesons you asked for (cf. IDmes)
  !****************************************************************************
  function sigma_barMes_res(targetN,q,IDbar,IDmes) result(sigma_VM)
    use particleDefinition
    use idTable, only: nres,nucleon
    use baryonWidth, only: partialWidthBaryon, FullWidthBaryon
    use particleProperties, only: hadron
    use clebschGordan, only: CG
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in)       :: targetN
    real, dimension(0:3), intent(in) :: q
    integer, intent(in)              :: IDbar, IDmes(:)
    real, dimension(lbound(IDmes,dim=1):ubound(IDmes,dim=1)) :: sigma_VM

    real, dimension(nucleon+1:nucleon+nres,lbound(IDmes,dim=1):ubound(IDmes,dim=1)) :: sigma
    integer :: resID,i
    real :: sigma_tot, mass, baremass, fullWidth
    integer :: res_I2, res_Iz2, bar_I2


    if (sqrt(Dot_Product(targetN%pos,targetN%pos)).lt.50) then
       write(*,*) 'In medium not yet implemented in Sigma_VM_res. STOP!'
       write(*,*) 'position=',targetN%pos
       call TRACEBACK()
    end if

    sigma = 0.

    bar_I2 = hadron(IDbar)%isospinTimes2

    ! loop over intermediate resonances
    do resID=nucleon+1,nucleon+nres

      sigma_tot=sigma_resProd(targetN,resId,q,baremass)
      mass=bareMass

      fullWidth=FullWidthBaryon(resID,mass)
      if (fullWidth.lt.1E-10) cycle

      res_I2 = hadron(resId)%isoSpinTimes2
      res_Iz2= 2*targetN%charge-1

      ! loop over meson final states
      do i = lbound(IDmes,dim=1),ubound(IDmes,dim=1)
         ! Branching ratios into V N
         sigma(resID,i) = &
              partialwidthBaryon(resID,mass,.false.,IDmes(i),IDbar) &
              * CG(hadron(IDmes(i))%isospinTimes2,bar_I2,res_I2, 0,res_Iz2)**2
      end do
      sigma(resID,:) = sigma(resID,:) * sigma_tot/fullWidth

    end do

    ! sum over all resonances
    sigma_VM(:) = sum(sigma,dim=1)

  end function Sigma_barMes_res

  !****************************************************************************
  !****f* resProd_lepton/dSdOmega_k_med_res
  ! NAME
  ! function dSdOmega_k_med_res(targetN,q,k,pf) result(sigma_dOmega)
  !
  ! PURPOSE
  ! Evaluates the resonance contribution to pion production in
  ! gamma R->eNPi scattering. The return value
  ! is dsigma/dOmega(pion) in [mub/Sr].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  ! * Isotropic decay of the resonance in its rest-frame.
  ! * In the vacuum.
  !
  ! INPUTS
  ! * type(particle)        :: targetN -- Target nucleon
  ! * real, dimension(0:3)  :: q       -- Virtual photon 4-momentum
  ! * real, dimension(0:3)  :: k       -- pion 4-momentum
  ! * real, dimension(0:3)  :: pf      -- Outgoing nucleon 4-momentum
  ! OUTPUT
  ! * logical               :: success -- flag
  ! * real, dimension(-1:1) :: sigma_dOmega -- dsigma/dOmega_pion; Index: qPion
  !****************************************************************************
  function dSdOmega_k_med_res(targetN,q,k,pf,success) result(sigma_dOmega)
    use constants, only: pi
    use particleDefinition
    use mediumDefinition, only: vacuum
    use idTable, only: nres, nucleon,pion
    use baryonWidthMedium, only: partialWidthBaryonMedium, WidthBaryonMedium
    use clebschGordan, only: CG
    use particleProperties, only: hadron
    use CALLSTACK, only: TRACEBACK

    type(particle), intent(in)            :: targetN
    real, dimension(0:3), intent(in)      :: q
    real, dimension(-1:1,0:3), intent(in) :: pf
    real, dimension(-1:1,0:3), intent(in) :: k
    real, dimension(-1:1) :: sigma_dOmega  ! dsigma/dOmega_electron/dE_electron/dOmega_pion

    real, dimension(-1:1,nucleon+1:nucleon+nres) :: sigma  ! dsigma/dOmega_electron/dE_electron

    real, dimension (0:3) :: pf_res!,q_res
    integer :: qPi    ! Charge of outgoing pion
    integer :: resID
    real :: branchingRatio,mass,fullWidth
    real :: sigma_tot, baremass
    logical :: success

    if (sqrt(Dot_Product(targetN%pos,targetN%pos))<=50.) then
       write(*,*) 'In medium not yet implemented in dSdOmega_k_med_res. STOP!'
       write(*,*) 'position=',targetN%pos
       call TRACEBACK()
    end if

    sigma = 0.
    success=.true.

    resIDLoop: do resID=nucleon+1,nucleon+nres
       pf_res=targetN%mom+q

       sigma_tot=sigma_resProd(targetN,resId,q,baremass)
       mass=bareMass

       fullWidth=WidthBaryonMedium(resID,mass,pf_res,vacuum)
       if (fullWidth.lt.1E-10) then
          sigma(:,resID)=0.
          cycle resIdLoop
       end if

       branchingRatio=partialWidthBaryonMedium(resID,mass,.false.,pion,nucleon,pf_res,vacuum)&
            &        /fullWidth

       sigma_tot=sigma_tot*branchingRatio/4./pi
       do qPi=targetN%charge-1,targetN%charge
          sigma(qPi,resID) = sigma_tot &
               * dOmegaCM_dOmega(pf_res,k(qPi,:),pf(qPi,:),success) &
               * CG(2,1,hadron(resId)%isoSpinTimes2, 2*qPi,2*targetN%charge-1-2*qPi)**2
          if (.not.success) then
             write(*,*)
             write(*,*) "k=",k(qPi,:)
             write(*,*) "pf=", pf(qPi,:)
             write(*,*) "pf of res=", pf_res(:)
             write(*,*) "resID, pion charge ", resID,qPi
             sigma_dOmega=0.
             return
          end if
       end do
    end do resIDLoop

    sigma_dOmega(:) = sum(sigma(:,:),dim=2)

  end function dSdOmega_k_med_res




  !****************************************************************************
  !****f* resProd_lepton/dSdO_fdE_fdO_k_med_res_EN
  ! NAME
  ! function dSdO_fdE_fdO_k_med_res_EN(eN,k,pf,processID) result(sigma_dOmega)
  !
  ! PURPOSE
  ! Evaluates the resonance contribution to pion production in
  ! eN->eR->eNPi scattering. The return value
  ! is dsigma/dOmega(electron)/dE(electron)/dOmega(pion) in [mb/GeV/Sr**2].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  ! * Isotropic decay of the resonance in its rest-frame.
  ! * In the vacuum.
  !
  ! INPUTS
  ! * type(electronNucleon_event) :: eN  -- electron-nucleon scattering event
  ! * real, dimension(0:3)        :: pf  -- Outgoing nucleon 4-momentum
  ! * real, dimension(0:3)        :: k   -- Outgoing pion 4-momentum
  ! * integer, optional           :: processID -- See module leptonicID for
  !   usage
  ! * integer, optional           :: pionNucleonSystem --
  !   If this parameter is set to 1, then we evaluate dOmega_pion in the
  !   calculation frame. If it's 2 then it is evaluated
  !   in the cm frame of the outgoing pion  and nucleon.
  ! OUTPUT
  ! * real,dimension(-1:1):: sigma_dOmega --
  !   dsigma/dOmega_electron/dE_electron/dOmega_pion; Index: qPion
  !
  ! NOTES
  ! * Enhances dSdO_fdE_fdO_k_med_res by allowing arbitrary electron
  !   momentum directions
  !****************************************************************************
  function dSdO_fdE_fdO_k_med_res_EN(eN,k,pf,processID_IN,pionNucleonSystem)  &
       result(sigma_dOmega)
    use particleDefinition
    use particleProperties, only: hadron
    use constants, only: pi
    use mediumDefinition, only: vacuum
    use idTable, only: nres, nucleon,pion
    use baryonWidthMedium, only: partialWidthBaryonMedium, WidthBaryonMedium
    use clebschGordan, only: CG
    use leptonicID, only: EM,CC,antiCC
    use eN_eventDefinition, only: electronNucleon_event,write_electronNucleon_event
    use CALLSTACK, only: TRACEBACK

    ! Input
    type(electronNucleon_event), intent(in)       :: eN
    real, dimension(0:3),        intent(in)       :: pf
    real, dimension(0:3),        intent(in)       :: k
    integer, optional,           intent(in)       :: processID_IN
    integer, optional,           intent(in)       :: pionNucleonSystem
    ! Result:
    real,dimension(-1:1)                          :: sigma_dOmega
    ! dsigma/dOmega_electron/dE_electron/dOmega_pion

    real, dimension(-1:1,nucleon+1:nucleon+nres)  :: sigma
    ! dsigma/dOmega_electron/dE_electron

    real, dimension(0:3)  :: pf_res
    integer               :: qPi    ! Charge of outgoing pion
    integer               :: resID
    real                  :: branchingRatio,mass,fullWidth
    integer               :: processID
    real                  :: sigma_tot
    real                  :: baremass
    logical               :: piN_inCM_frame


    if (present(pionNucleonSystem)) then
       piN_inCM_frame=(pionNucleonSystem.eq.2)
    end if

    if (sqrt(Dot_Product(eN%nucleon%pos,eN%nucleon%pos))<=50.) then
       write(*,*) 'In medium not yet implemented in dSdO_fdE_fdO_k_med_res_EN. STOP!'
       write(*,*) 'position=',eN%nucleon%pos
       call write_electronNucleon_event(eN,.FALSE.,.FALSE.)
       call TRACEBACK()
    end if

    sigma(:,:)=0.

    processID=EM
    if (present(processID_IN)) processID=processID_IN

    resIDLoop: do resID=nucleon+1,nucleon+nres
       sigma_tot=dSigmadOmega_fdE_f_resProd_eN(eN,resID,pf_res,baremass,processID)

       mass=bareMass

       fullWidth=WidthBaryonMedium(resID,mass,pf_res,vacuum)
       if (fullWidth.lt.1E-10) then
          sigma(:,resID)=0.
          cycle resIdLoop
       end if

       branchingRatio=partialWidthBaryonMedium(resID,mass,.false.,pion,nucleon,pf_res,vacuum)&
            &        /fullWidth

       if (piN_inCM_frame) then
          sigma_tot=sigma_tot*branchingRatio/(4.*pi)
       else
          sigma_tot=sigma_tot*branchingRatio*dOmegaCM_dOmega(pf_res,k,pf)/(4.*pi)
       end if

       ! NC not yet implemented for pi bg

       select case (processID)
       case (EM)
          do qPi=eN%nucleon%charge-1,eN%nucleon%charge
             sigma(qPi,resID)=sigma_tot &
                  * CG(2,1,hadron(resId)%isoSpinTimes2, 2*qPi,2*eN%nucleon%charge-1-2*qPi)**2
          end do

       case (CC)
          do qPi=eN%nucleon%charge-1,eN%nucleon%charge
             sigma(qPi,resID)=sigma_tot &
                  * CG(2,1,hadron(resId)%isoSpinTimes2, 2*qPi,2*eN%nucleon%charge+1-2*qPi)**2
          end do

       case (antiCC)
          do qPi=eN%nucleon%charge-2,eN%nucleon%charge-1
             sigma(qPi,resID)=sigma_tot &
                  * CG(2,1,hadron(resId)%isoSpinTimes2, 2*qPi,2*eN%nucleon%charge-3-2*qPi)**2
          end do

       end select

    end do resIDLoop

    sigma_dOmega(:) = sum(sigma(:,:),dim=2)


  end function dSdO_fdE_fdO_k_med_res_eN


  !****************************************************************************
  !****f* resProd_lepton/dOmegaCM_dOmega
  ! NAME
  ! real function dOmegaCM_dOmega()
  ! PURPOSE
  ! Evaluates the Jacobian for dOmega_CM(pion)/dOmega_lab(pion)
  !****************************************************************************
  real function dOmegaCM_dOmega(pf_res,k,pf,success)
    use minkowski, only: abs4
    use twoBodyTools, only: pcm
    use constants, only: mN, mPi

    real :: abs_k_vec,abs_pf_vec!,abs_q_vec
    real :: cos_theta_k
    real :: kcm
    real, dimension (0:3),intent(in) :: k,pf,pf_res
    logical, optional :: success

    if (present(success)) success=.true.
    kcm= pCM(abs4(pf_res), mN, mPi)
    if (kcm.lt.1E-10) then
       write(*,*) 'WARNING: Trying to produce resting pion in CM=', k
       write(*,*) 'kcm=', kcm
       dOmegaCM_dOmega=0.
       if (present(success)) success=.false.
       return
    end if

    abs_k_vec=sqrt(Dot_Product(k(1:3),k(1:3)))
    !abs_q_vec=sqrt(Dot_Product(q(1:3),q(1:3)))
    abs_pf_vec=sqrt(Dot_Product(pf(1:3),pf(1:3)))

    ! Make case study for the evaluation of cos(theta)
    if (abs_k_vec.lt.1E-10) then
       write(*,*) 'WARNING: Trying to produce resting pion=', k
       write(*,'(A,4E15.3)') 'k=', k
       !write(*,'(A,4E15.3)') 'energy_li, energy_lf=', energy_li, energy_lf
       dOmegaCM_dOmega=0.
       if (present(success)) success=.false.
       return
    else if (abs_pf_vec.lt.1E-10) then
       dOmegaCM_dOmega=abs4(pf_res)*abs_K_vec**2/kCM/(abs_K_vec*pf(0))

    else
       cos_theta_k=Dot_Product(k(1:3),pf(1:3))/abs_k_vec/abs_pf_vec
       dOmegaCM_dOmega=abs4(pf_res)*abs_K_vec**2/kCM/   &
                       & (abs_K_vec*pf(0)-abs_pf_vec*k(0)*cos_theta_k)
    end if

  end function dOmegaCM_dOmega




  !****************************************************************************
  !****f* resProd_lepton/sigma_resProd
  ! NAME
  ! function sigma_resProd(targetN,resID,q,baremass) result(xSection)
  !
  ! PURPOSE
  !
  ! Evaluates the cross section for gamma N -> R scattering. The return value
  ! is sigma in [mb].
  !
  ! Assumptions:
  ! * No interferences among resonances.
  !
  ! INPUTS
  ! * type(particle)         :: targetN   -- Target nucleon
  ! * real, dimension (0:3)  :: q         -- Virtual photon 4-momentum
  ! * integer                :: resId     -- ID of resonance
  ! OUTPUT
  ! *  real                  :: baremass  --
  !    bare mass of resonance (mass without scalar potential)
  ! *  real                  :: xSection  -- sigma
  !****************************************************************************
  function sigma_resProd(targetN,resID,q,baremass_res) result(xSection)
    use leptonicID, only: EM
    use spectralFunc, only: specFunc
    use constants, only: GeVSquared_times_mb, pi
    use minkowski, only: SP, abs4
    use hadronTensor_ResProd
    use degRad_conversion
    use particleDefinition
    use constants, only: electronChargeSQ
    use RMF
    use densitymodule
    use dichteDefinition

    type(particle),       intent(in)     :: targetN
    integer,              intent(in)     :: resID
    real, dimension(0:3), intent(in)     :: q
    real,                 intent(out)    :: baremass_res

    real                                 :: xSection
    real, dimension(0:3)                 :: pin, pout ! kinetic
    real, dimension(0:3)                 :: hin, hout ! canonical
    complex, dimension(0:3,0:3)          :: HT ! hadron tensor
    complex                              :: matrixElement_Squared
    real                                 :: kinematics
    logical, parameter                   :: debug_this=.false.
    real                                 :: spec
    real                                 :: fakMedium

    ! the default is F,F,F, the other option is T,T,T:
    logical, parameter :: useCanonical = .false.
    logical, parameter :: useMediumModif = .false.
    logical, parameter :: useV = .false.


    if (debug_this) then
       write(*,'(A,4E15.4)')      'q=',q
       write(*,'(A,I8)')          'resID=',resID
       write(*,'(A,4E15.4)')      'p=',targetN%mom
    end if

    baremass_res= 0.
    xSection    = 0.
    fakMedium   = 1.0

    pin = targetN%mom
    pout = pin + q

    hin = pin
    hin(0) = FreeEnergy(targetN)
    ! other way to subtract the potential:
    ! (please note, that this does not guarantee that abs4(hin)=0.938 !
    ! The difference may be up to -50 MeV)
    if (useV.and.getRMF_flag()) hin = Particle4MomentumRMF(targetN)

    hout = hin + q

    if (useMediumModif) then
       fakMedium = (pin(0)-pin(3))*abs4(hin)*abs4(hout)*pout(0) &
            /((hin(0)-hin(3))*abs4(pin)*abs4(pout)*hout(0))
    end if

!!$    call PrintSomeStuff(targetN,resID,q)


    ! ********** Cross section ******************

    ! The spectral function has to be calculated with the kinetic momentum,
    ! since the pole mass is shifted by the scalar density
    Spec=specFunc(resID,targetN%Charge,pout,targetN%pos,baremass_res)

    if (useCanonical) then
       pin = hin
       pout = hout
    end if

    if (.not.hadronTensor_R(pin,pout,resID,targetN%Charge,EM,HT,baremass_res) ) return ! ==> failure

    matrixElement_Squared=1./2.*electronChargeSQ &
         & * (-HT(0,0)+HT(1,1)+HT(2,2)+HT(3,3))

    if (debug_this) then
       write(*,'(A,4E15.3)') 'M=',&
            Real(HT(0,0)),Real(HT(1,1)), &
            Real(HT(2,2)),Real(HT(3,3))
       write(*,*) real(matrixElement_Squared), &
            1./2.*(-real(HT(0,0))+real(HT(1,1))+real(HT(2,2))+real(HT(3,3)))
    end if

    if (abs(AIMAG(matrixElement_Squared))>1E-7) &
         write(*,*) 'ERROR: (Matrix element)^2  is imaginary:', &
         matrixElement_Squared

    kinematics=1./(4.*abs(SP(q,pin)))

    Xsection=kinematics *REAL(matrixElement_Squared)*2.*pi * Spec
    if (useMediumModif) Xsection = Xsection * fakMedium

    ! ********** UNIT CONVERSION ******************
    ! Convert from 1/GeV**2 to mb
    Xsection=Xsection/GeVSquared_times_mb

    if (Xsection<-1E-9) then
       write(*,*) 'Xsection less than zero'
       write(*,*) targetN%mom
       write(*,*) targetN%charge
       write(*,*) resID
       write(*,*) q
       stop
    end if

    if (debug_this) then
       write(*,'(A,E15.4)') 'Sigma=',xsection
       write(77,'(3E15.4)') q(0),xsection,baremass_res
       write(*,*)
    end if

  end function sigma_resProd

  subroutine PrintSomeStuff(targetN,resID,q)
    use particleDefinition
    use RMF
    use densitymodule
    use dichteDefinition
    use spectralFunc, only: specFunc
    use leptonicID, only: EM
    use hadronTensor_ResProd
    use degRad_conversion
    use particleDefinition
    use constants, only: electronChargeSQ
    use minkowski, only: SP, abs4

    type(particle),       intent(in)     :: targetN
    integer,              intent(in)     :: resID
    real, dimension(0:3), intent(in)     :: q


    real, dimension(0:3,1:3) :: p, r
    real, dimension(1:3) :: ME, kin
    type(dichte) :: density
    real :: mDirac1, mDirac2, spec, baremass_res, sigma
    integer :: i, resQ
    complex, dimension(0:3,0:3) :: HT ! hadron tensor


    resQ = targetN%Charge ! as abbrev

    p(:,1) = targetN%mom ! kinetic

    p(:,3) = p(:,1)
    p(0,3) = FreeEnergy(targetN) ! free (vacuum)

    if (getRMF_flag()) then
       p(:,2) = Particle4MomentumRMF(targetN) ! canonical
    else
       p(:,2) = p(:,3)
    end if

    do i=1,3
       r(:,i) = p(:,i) + q
       kin(i) = 1./(4.*abs(SP(q(:),p(:,i))))
    end do



    density = densityAt(targetN%pos)
    mDirac1 = mDiracNucleon_Approx(density%baryon(0))
    mDirac2 = mDirac1535_Approx(density%baryon(0))
    spec = specFunc(resID,resQ,r(:,1),targetN%pos,baremass_res)
    call getFieldRMF(targetN%pos, sigma=sigma)

    write(829,*) abspos(targetN),abs4(p(:,1)),abs4(r(:,1)),mDirac1,p(0,1),r(0,1)

    ME(:) = 0.
    do i=1,3
       if (hadronTensor_R(p(:,i),r(:,i),resID,resQ,EM, HT,baremass_res) ) &
            ME(i) = REAL(-HT(0,0)+HT(1,1)+HT(2,2)+HT(3,3))
    end do
    ME(:) = ME(:) * 1./2.*electronChargeSQ

    write(827,*) abspos(targetN), resQ, spec, baremass_res, & ! 1,2,3,4
         (abs4(p(:,i)),i=1,3), & ! 5,6,7
         (abs4(r(:,i)),i=1,3), & ! 8,9,10
         kin(1:3), & ! 11,12,13
         ME(1:3), & ! 14,15,16
         mDirac1, mDirac2, & ! 17,18
         density%baryon(0), & ! 19
         sigma ! 20


  end subroutine PrintSomeStuff


  subroutine readInput

!!$    use output
!!$
!!$    integer :: ios
!!$    !*************************************************************************
!!$    !****n* resProd_lepton/resonanceProd_electron
!!$    ! NAME
!!$    ! NAMELIST /resonanceProd_electron/
!!$    ! PURPOSE
!!$    ! Namelist for moduleresProd_lepton includes:
!!$    ! * debug
!!$    !*************************************************************************
!!$    NAMELIST /resonanceProd_electron/ debug
!!$
!!$    call Write_ReadingInput('resonanceProd_electron',0)
!!$    rewind(5)
!!$    read(5,nml=resonanceProd_electron,IOSTAT=ios)
!!$    call Write_ReadingInput("resonanceProd_electron",0,ios)
!!$
!!$    write(*,*) 'debug?', debug
!!$
!!$    call Write_ReadingInput('resonanceProd_electron',1)

  end subroutine readInput


end module resProd_lepton
