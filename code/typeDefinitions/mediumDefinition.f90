!******************************************************************************
!****m* /mediumDefinition
! NAME
! module mediumDefinition
!
! PURPOSE
! Includes a type which defines the nuclear medium
!******************************************************************************

module mediumDefinition

  implicit none

  !****************************************************************************
  !****t* mediumDefinition/medium
  ! SOURCE
  !
  Type medium
     real :: density=0.        ! sqrt(j_mu j^mu) of baryons in units 1/fm**3
     real :: densityProton=0.  ! sqrt(j_mu j^mu) of protons in units 1/fm**3
     real :: densityNeutron=0. ! sqrt(j_mu j^mu) of neutrons in units 1/fm**3
     real :: temperature=0.    ! Temperature of medium
     real, dimension(1:3) :: betaLRF = (/0.,0.,0./) ! boost from calc frame to LRF
     logical  :: useMedium=.false.     ! whether to use medium informations
  End Type medium
  !****************************************************************************


  type(medium), parameter :: vacuum = medium(0.,0.,0.,0.,(/0.,0.,0./),.false.)


contains

  !****************************************************************************
  !****s* mediumDefinition/writeMedium
  ! NAME
  ! subroutine writeMedium(med)
  !
  ! PURPOSE
  ! Prints the content of an instance of type medium to the screen.
  !
  ! INPUTS
  ! * type(medium) :: med
  ! OUTPUT
  ! * NONE
  !****************************************************************************
  subroutine writeMedium(med)
    type(medium) :: med
    write(*,'(A,1P,3E15.3,3E15.4,E15.4,L8)') 'Medium:', &
         med%density, med%densityProton, med%densityNeutron, &
         med%betaLRF, &
         med%temperature, med%useMedium
  end subroutine writeMedium

end module mediumDefinition
