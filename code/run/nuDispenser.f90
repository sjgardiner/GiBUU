module nuDispenser

  use nucleusDefinition, only: tNucleus
  use particleDefinition, only: particle, setToDefault, setNumbersToDefault

  implicit none
  private

  public :: initNuDispenser

  type(particle),Allocatable :: realParticles(:,:)
  type(particle),Allocatable :: pertParticles(:,:)
  type(tNucleus),pointer :: targetNuc => NULL()

contains

  subroutine initNuDispenser(jobcard)

    use baryonPotentialModule, only: HandPotentialToDensityStatic
    use checks, only: ChecksSetDefaulSwitches
    use coulomb, only: updateCoulomb
    use densityModule, only: updateDensity
    use deuterium_PL, only: deuteriumPL_assign
    use energyCalc, only: updateEnergies
    use initNucleus_in_PS, only: initNucPhaseSpace
    use inputGeneral, only: readinputGeneral, eventType, numEnsembles
    use insertion, only: GarbageCollection
    use nucleus, only: getTarget
    use particleProperties, only: initParticleProperties
    use PILCollected, only: PILCollected_ZERO
    use propagation, only: updateVelocity
    use random, only: setRandom
    use version, only: printVersion
    use yukawa, only: updateYukawa


    character(*), intent(in) :: jobcard

    integer :: lengthReal = 0  ! max number of real particles per ensemble
    integer :: lengthPert = 0  ! max number of pert particles per ensemble

    open(5,file=jobcard,status="old")

    call printVersion
    call readInputGeneral
    call initParticleProperties
    call ChecksSetDefaulSwitches(EventType)
    call setRandom()

    targetNuc => getTarget()    ! set up target resting at 0
    lengthReal = targetNuc%mass ! Real particles fixed by nucleons in target

    if (targetNuc%ReAdjustForConstBinding) then
       write(*,*) 'we now have to readjust the density'
       call handPotentialToDensityStatic(targetNuc)
    end if

    lengthPert =  max(100,10*targetNuc%mass)

    allocate(realparticles(1:numEnsembles,1:lengthReal))
    allocate(pertparticles(1:numEnsembles,1:lengthPert))

    call setNumbersToDefault

    call GarbageCollection(pertParticles,.true.)
    call GarbageCollection(realParticles)

    call PILCollected_ZERO ! reset all particle info lists

    !=== setUpTarget ===

    call initNucPhaseSpace(realparticles,targetNuc)
    call updateDensity(realParticles)
    call updateCoulomb
    call updateYukawa(.true.)
    call updateVelocity(realParticles)
    call updateEnergies(realParticles)

    call deuteriumPL_assign(realParticles)

  end subroutine initNuDispenser


end module nuDispenser
