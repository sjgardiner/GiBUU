!******************************************************************************
!****m* /eventprops
! NAME
! module eventprops
!
! PURPOSE
! ?
!******************************************************************************
module eventprops

  use CallStack, only: Traceback

  implicit none
  private

  ! module variables as storage:
  real, dimension(1:3), save, public :: pos_
  real, dimension(0:3), public :: momIn_,momOut_
  integer, save, public :: charge_
  logical, save, public :: IsSet_ = .false.
  real, save :: rho_


  public :: set_evprops
  public :: get_evprops


  public :: Delta_inmedcorr
  public :: nucleon_inmedcorr





contains

  !****************************************************************************
  !****s* eventprops/set_evprops
  ! NAME
  ! subroutine set_evprops(charge,in_momentum,out_momentum,position)
  !****************************************************************************
  subroutine set_evprops(charge,in_momentum,out_momentum,position)
    ! sets e-N event properties

    use densityModule, only: densityAt
    use dichteDefinition

    real, dimension(1:3), intent(in) :: position
    real, dimension(0:3), intent(in) :: in_momentum, out_momentum   ! momenta for nucleon
    integer, intent(in) :: charge

    type(dichte) :: D

    charge_ = charge
    momIn_ = in_momentum
    momOut_ = out_momentum
    pos_ = position

    ! We use the density at the point of the delta to evaluate gamma.
    !dens_atDelta=densityAt(pos_delta)
    !rho=dens_atDelta%proton(0)+dens_atDelta%neutron(0)

    D = densityAt(position)
    rho_ = D%baryon(0)

    IsSet_ = .true.

  end subroutine set_evprops

  subroutine get_evprops(charge,in_momentum,out_momentum,position,density)

    integer, intent(out) :: charge
    real, dimension(1:3), intent(out) :: position
    real, dimension(0:3), intent(out) :: in_momentum,out_momentum
    real, intent(out) :: density

    if (.not. IsSet_) call Traceback('event properties not set!')

    charge = charge_
    in_momentum = momIn_
    out_momentum = momOut_
    position = pos_
    density = rho_

  end subroutine get_evprops

  !****************************************************************************
  !****s* eventprops/Delta_inmedcorr
  ! NAME
  ! subroutine Delta_inmedcorr(invmass,charge,mom,pos,rho,inmedwidth_corr,inmedmass_corr)
  !
  ! PURPOSE
  ! Provide the in-medium mass-shift and collisional width for the Delta resonance
  !
  !
  ! INPUTS
  ! *
  ! * integer                 :: charge       -- charge of Delta
  ! * real                    :: invmass      -- invariant mass of Delta
  ! * real                    :: mom          -- momentum of Delta
  ! * real                    :: pos          -- position of Delta (needed for density-dependence)
  ! * real                    :: rho          -- baryon density at position of Delta
  !
  ! OUTPUT
  ! * real                    :: inmedwidth_corr  -- inmedium change of Delta width (from Salcedo-Oset)
  ! * real                    :: inmedmass_corr   -- inmedium change of Delta mass
  !
  !****************************************************************************
  subroutine Delta_inmedcorr(invmass,charge,mom,pos,rho,inmedwidth_corr,inmedmass_corr)

    use spectralfunc
    use potentialMain, only: scapot
    use deltaWidth
    use baryonWidthMedium, only: get_mediumSwitch_Delta

    implicit none

    real, intent(in) :: rho,invmass
    integer, intent(in) :: charge
    real, intent(out) :: inmedwidth_corr, inmedmass_corr
    real :: imsig1,imsig3,imsigq,baremass
    real, intent(in),dimension(0:3)    :: mom  ! momentum of baryon in LRF
    real, intent(in), dimension(1:3)    :: pos  ! position of baryon
    logical :: flagOK,mediumswitch_Delta

    mediumswitch_Delta = get_mediumSwitch_Delta()
    !write(*,*) 'Delta_inmedcorr mediumswitch_Delta=',mediumswitch_Delta
    if(mediumSwitch_Delta) then
       call deloset(invmass,rho,imsig1,imsig3,imsigq)
       inmedwidth_corr = 2.*(imsig1 + imsig3 + imsigq)
    end if

    inmedmass_corr = scapot(2,charge,mom,pos,baremass,flagOK)
    if (.not.flagOK) inmedmass_corr = 0.

  end subroutine Delta_inmedcorr

  subroutine Nucleon_inmedcorr(charge,mom,pos,nucleonmass_corr)

    use potentialMain, only: scapot

    implicit none
    integer, intent(in) :: charge
    real, intent(in), dimension(0:3) ::  mom
    real, intent(in), dimension(1:3) ::  pos
    real, intent(out) :: nucleonmass_corr
    real :: baremass
    logical :: flagOK

    nucleonmass_corr = scapot(1,charge,mom,pos,baremass,flagOK)
    if(.not.flagOK) nucleonmass_corr = 0.

  end subroutine Nucleon_inmedcorr


end module eventprops
