program testPotNucleus

  ! this program prints a table of the potential of the nucleons and the deltas
  ! as function of the radius

  ! The code prints values for two values of the momentum (0, pF), while these
  ! are just chosen as examples. The User may choose other values ad libitum.

  use inputGeneral
  use version, only: printVersion
  use particleDefinition
  use particleProperties, only: initParticleProperties, hadron
  use nucleusDefinition
  use nucleus, only: getTarget
  use initNucleus_in_PS, only: initNucPhaseSpace
  use densityModule, only: updateDensity, FermiMomAt
  use coulomb, only: updateCoulomb
  use yukawa, only: updateYukawa
  use propagation, only: updateVelocity
  use energyCalc, only: updateEnergies

  use CallStack, only: traceBack
  use potentialModule, only: potential_LRF


  implicit none

  type(tNucleus), pointer :: targetNuc
  type(particle), allocatable :: realParticles(:,:)
  integer :: lengthReal = 0

  call printVersion
  call readInputGeneral
  call initParticleProperties

  targetNuc => getTarget()
  lengthReal = targetNuc%mass
  allocate(realparticles(1:numEnsembles,1:lengthReal))

  call initNucPhaseSpace(realparticles,targetNuc)
  call updateDensity(realParticles)
  call updateCoulomb
  call updateYukawa(.true.)
  call updateVelocity(realParticles)
  call updateEnergies(realParticles)


  call doPlot


contains

  subroutine doPlot

    real, parameter :: maxR = 15.0, dR = 0.0001
    integer, parameter :: nR = maxR/dR
    integer :: iR, iQ, iP
    real, dimension (1:3) :: pos=0.
    real, dimension (0:3) :: mom=0.
    real, parameter :: nP = 2
    real, dimension (1:nP) :: aP

    type(particle) :: teilchen

    real, dimension(100) :: aVal
    integer :: iVal
    real :: pF


    do iR=0,nR
       write(*,*) 'Doing r=',iR*dR

       pos = (/ 0.0, 0.0, iR*dR /)

       pF = FermiMomAt(pos) ! averaged value

       aP = (/ 0.0, pF /)

       teilchen%pos = pos

       aVal = 0.0
       iVal = 0

       !===== NUKLEON =====
       teilchen%ID = 1
       teilchen%mass = hadron(1)%mass

       !...charges:
       do iQ=0,1
          teilchen%charge = iQ

          !...momenta:
          do iP=1,nP
             teilchen%mom(1:3) = (/ aP(iP), 0.0, 0.0 /)
             teilchen%mom(0)   = FreeEnergy(teilchen)

             iVal = iVal+1
             aVal(iVal) = potential_LRF(teilchen)

          end do
       end do


       !===== DELTA =====
       teilchen%ID = 2
       teilchen%mass = hadron(2)%mass

       !...charges:
       do iQ=-1,2
          teilchen%charge = iQ

          !...momenta:
          do iP=1,nP
             teilchen%mom(1:3) = (/ aP(iP), 0.0, 0.0 /)
             teilchen%mom(0)   = FreeEnergy(teilchen)

             iVal = iVal+1
             aVal(iVal) = potential_LRF(teilchen)

          end do
       end do


       !===== OUTPUT =====

       write(111,*) iR*dR,aVal(1:iVal)


    end do

  end subroutine doPlot

end program testPotNucleus
