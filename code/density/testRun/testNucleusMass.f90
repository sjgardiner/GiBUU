program testNucleusMass

  use inputGeneral
  use version, only: PrintVersion
  use particleDefinition
  use particleProperties, only: initParticleProperties
  use nucleusDefinition
  use nucleus, only : getProjectile,getTarget
  use dichteDefinition
  use densityStatic
  use output, only: Write_ReadingInput
  use CallStack, only: traceBack
  use initNucleus_in_PS, only: initNucPhaseSpace
  use densityModule, only: updateDensity, updateRMF, storeFields
  use coulomb, only: updateCoulomb
  use yukawa, only: updateYukawa
  use propagation, only: updateVelocity
  use energyCalc, only: updateEnergies
  use RMF, only: getRMF_flag


  implicit none

  type(tNucleus),pointer :: targetNuc
  integer :: ios

  type(particle),Allocatable :: realParticles(:,:)
  integer :: lengthReal = 0

  integer :: iEns, iPart
  real, dimension(0:3) :: Mom0
  real, dimension(1:2) :: aveM, aveP, aveE, avePF
  real :: h
  integer :: iAnneal, nAnneal

  call PrintVersion

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


!  nAnneal = 10
!  nAnneal = 500/lengthReal
!  nAnneal = 2000/lengthReal
!  do iAnneal=1,nAnneal
!     call Anneal()
!  end do


  aveM = 0.0
  aveP = 0.0
  aveE = 0.0
  avePF = 0.0
  do iEns=1,numEnsembles
     Mom0 = 0.0
     do iPart=1,lengthReal
        Mom0 = Mom0 + realparticles(iEns,iPart)%mom
        h = absMom(realparticles(iEns,iPart))
        AvePF = AvePF + (/ h,h**2 /)
     end do

     aveE = aveE + (/ Mom0(0), Mom0(0)**2 /)
     h = sum(Mom0(1:3)**2)
     aveP = aveP + (/ sqrt(h), h /)

     h = Mom0(0)**2-h
     aveM = aveM + (/ sqrt(h), h /)

     AvePF = AvePF/lengthReal

     write(501,'(i5,1P,5e13.4)') iEns, sqrt(Mom0(0)**2-sum(Mom0(1:3)**2)), Mom0
     write(503,'(i5,1P,5e13.4)') iEns, AvePF
  end do

  aveP = aveP/numEnsembles
  aveM = aveM/numEnsembles
  aveE = aveE/numEnsembles

  write(502,'(i5,1P,6e13.4)') lengthReal, aveM, aveP, aveE

contains

  subroutine Anneal()

    use random, only: rnOmega

    real, dimension(0:3) :: MomNew, Mom1
    real :: hMom0, hMom1

    do iEns=1,numEnsembles
       Mom0 = 0.0
       do iPart=1,lengthReal
          Mom0 = Mom0 + realparticles(iEns,iPart)%mom
       end do
       hMom0 = sum(Mom0(1:3)**2)

       do iPart=1,lengthReal
          h = absMom(realparticles(iEns,iPart))
          MomNew(1:3) = h*rnOmega()
          MomNew(0) = realparticles(iEns,iPart)%mom(0)
          Mom1 = Mom0 - realparticles(iEns,iPart)%mom + MomNew
          hMom1 = sum(Mom1(1:3)**2)

          if (hMom1 < hMom0) then
             realparticles(iEns,iPart)%mom = MomNew
             Mom0 = Mom1
             hMom0 = hMom1
          end if

       end do

    end do
  end subroutine Anneal


end program testNucleusMass
