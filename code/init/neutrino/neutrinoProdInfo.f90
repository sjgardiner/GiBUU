!******************************************************************************
!****m* /neutrinoProdInfo
! NAME
! module neutrinoProdInfo
!
! PURPOSE
! This module stores information about the initial neutrino event.
! (structure taken from EventInfo_HiLep)
!
!******************************************************************************
module neutrinoProdInfo
  implicit none

  private

  !****************************************************************************
  !****t* neutrinoProdInfo/tneutrinoProdInfo
  ! SOURCE
  !
  type tneutrinoProdInfo
     real,dimension(0:3) :: mom_lepIn  = 0.
     real,dimension(0:3) :: mom_lepOut = 0.
     real,dimension(0:3) :: mom_bos    = 0.
     real,dimension(0:3) :: mom_nuc    = 0.
     logical             :: flagOK     = .false.
     integer             :: prod_id    = 0
     real                :: perweight  = 0.
     integer             :: chrg_nuc   = 0
  end type tneutrinoProdInfo
  !
  ! PURPOSE
  ! This holds all information we want to store connected to
  ! neutrino induced reactions.
  !****************************************************************************

  type(tneutrinoProdInfo),save,dimension(:),allocatable :: nuProdInfo

  public :: neutrinoProdInfo_Init
  public :: neutrinoProdInfo_Store
  public :: neutrinoProdInfo_Get
  public :: neutrinoProdInfo_Dump

contains

  !****************************************************************************
  !****s* neutrinoProdInfo/neutrinoProdInfo_Init
  ! NAME
  ! subroutine neutrinoProdInfo_Init(NumInitialEvents)
  !
  ! PURPOSE
  ! allocate memory and reset the corresponding arrays.
  !
  ! INPUTS
  ! * integer :: NumInitialEvents -- number of possible initial events
  !
  ! OUTPUT
  ! ---
  !
  ! NOTES
  ! The current implementation is very 'expensive'
  !****************************************************************************
  subroutine neutrinoProdInfo_Init(NumInitialEvents)

    integer,intent(in) :: NumInitialEvents

    if (allocated(nuProdInfo)) then
       deallocate(nuProdInfo)
    end if
    allocate(nuProdInfo(1:NumInitialEvents))

  end subroutine NeutrinoProdInfo_Init


  !****************************************************************************
  !****s* neutrinoProdInfo/neutrinoProdInfo_Store
  ! NAME
  ! subroutine neutrinoProdInfo_Store(i,prod_id,perweight,Mom_LepIn,
  ! Mom_LepOut,Mom_Bos,Mom_Nuc,Chrg_Nuc)
  !
  ! PURPOSE
  ! Store the event info connected with number "i":
  !
  !
  ! INPUTS
  ! * integer :: i -- actual number of event
  ! * integer :: prod_id (1=N, 2=Delta, 3=P_11(1440)  32=pi-neutron-background 34=DIS  and so on)
  ! * real    :: perweight
  ! * real,dimension(0:3)    :: Mom_LepIn
  ! * real,dimension(0:3)    :: Mom_LepOut
  ! * real,dimension(0:3)    :: Mom_Bos
  ! * real,dimension(0:3)    :: Mom_Bos
  ! * integer :: Chrg
  ! OUTPUT
  ! ---
  !****************************************************************************
  subroutine neutrinoProdInfo_Store(i,prod_id,perweight,Mom_LepIn,Mom_LepOut,Mom_Bos,Mom_Nuc,Chrg_Nuc)

    integer,intent(in)             :: i
    integer,intent(in)             :: prod_id
    real,intent(in)                :: perweight
    real,dimension(0:3),intent(in) :: Mom_LepIn
    real,dimension(0:3),intent(in) :: Mom_LepOut
    real,dimension(0:3),intent(in) :: Mom_Bos
    real,dimension(0:3),intent(in) :: Mom_Nuc
    integer,intent(in) :: Chrg_Nuc

    nuProdInfo(i)%flagOK=.true.
    nuProdInfo(i)%prod_id=prod_id
    nuProdInfo(i)%perweight=perweight
    nuProdInfo(i)%mom_LepIn=mom_LepIn
    nuProdInfo(i)%mom_LepOut=mom_LepOut
    nuProdInfo(i)%mom_Bos=mom_Bos
    nuProdInfo(i)%mom_nuc=Mom_Nuc
    nuProdInfo(i)%chrg_nuc=Chrg_Nuc

  end subroutine NeutrinoProdInfo_Store


  !****************************************************************************
  !****f* neutrinoProdInfo/neutrinoProdInfo_Get
  ! NAME
  ! logical function neutrinoProdInfo_Get(i,prod_id,perweight,Mom_LepIn,
  ! Mom_LepOut,Mom_Bos,Mom_Nuc,Chrg_Nuc
  !
  ! PURPOSE
  ! Get the event info stored connected with number "i".
  !
  ! INPUTS
  ! * integer :: i -- actual number of event
  !
  ! OUTPUT
  ! * integer :: prod_id
  ! * real    :: perweight
  ! * real,dimension(0:3)    :: Mom_LepIn
  ! * real,dimension(0:3)    :: Mom_LepOut
  ! * real,dimension(0:3)    :: Mom_Bos
  ! * real,dimension(0:3), OPTIONAL    :: Mom_Bos
  ! * integer, OPTIONAL :: Chrg_Nuc
  !****************************************************************************
  logical function neutrinoProdInfo_Get(i,prod_id,perweight,Mom_LepIn,Mom_LepOut,Mom_Bos,Mom_Nuc, Chrg_Nuc)

    integer,intent(in)              :: i
    integer,intent(out)             :: prod_id
    real,intent(out)                :: perweight
    real,dimension(0:3),intent(out) :: Mom_LepIn
    real,dimension(0:3),intent(out) :: Mom_LepOut
    real,dimension(0:3),intent(out) :: Mom_Bos
    real,dimension(0:3),intent(out),OPTIONAL :: Mom_Nuc
    integer,intent(out),OPTIONAL             :: Chrg_Nuc


    neutrinoProdInfo_Get = .FALSE.

    if (.not.allocated(nuProdInfo)) return

    prod_id=nuProdInfo(i)%prod_id
    perweight=nuProdInfo(i)%perweight
    mom_LepIn=nuProdInfo(i)%mom_LepIn
    mom_LepOut=nuProdInfo(i)%mom_LepOut
    mom_Bos=nuProdInfo(i)%mom_Bos
    if (present(Mom_Nuc)) mom_Nuc=nuProdInfo(i)%mom_Nuc
    if (present(Chrg_Nuc)) Chrg_Nuc = nuProdInfo(i)%chrg_nuc

    neutrinoProdInfo_Get = nuProdInfo(i)%flagOK

  end function NeutrinoProdInfo_Get

  !****************************************************************************
  !****s* neutrinoProdInfo/neutrinoProdInfo_Dump
  ! NAME
  ! subroutine neutrinoProdInfo_Dump(i)
  !
  ! PURPOSE
  ! Dump the event info stored connected with number "i" to stdout.
  !
  ! INPUTS
  ! * integer :: i -- actual number of event
  !****************************************************************************
  subroutine neutrinoProdInfo_Dump(i)

    integer,intent(in)               :: i

    integer              :: prod_id
    real                 :: perweight
    real,dimension(0:3)  :: mom_LepIn, mom_LepOut,mom_Bos, mom_Nuc
    integer :: Chrg_Nuc

    if (neutrinoProdInfo_Get(i,prod_id,perweight,mom_LepIn,mom_LepOut,mom_Bos, mom_Nuc, Chrg_Nuc)) then
       write(*,'("** neutrinoProdInfo(",i7,")    prod_id=",i7,"   perweight=",1P,e13.4,0P)') i,prod_id,perweight

    else
       write(*,'("** ",A,i7," not found!")') 'neutrinoProdInfo:',i
    end if

  end subroutine neutrinoProdInfo_Dump

end module NeutrinoProdInfo
