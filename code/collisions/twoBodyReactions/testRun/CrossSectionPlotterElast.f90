!***************************************************************************
!****m* /CrossSectionPlotter2
! NAME
! program CrossSectionPlotter2
!
! PURPOSE
! This is a standalone main routine. It is used to plot the cross sections
! of specific channels, e.g. p p -> p p omega  or
!                            p p -> p Delta+
!
! INPUTS
! The Namelist "Plotter2" in the Jobcard and additional (usual namelists)
!
! NOTES
! particle 2 is the moving one, i.e. everything is in the rest frame of
! particle 1.
!
! The cross sections plotted here will either come from the low-energy
! resonance model or the high-energy string model or a mixture of both,
! according to the high-energy switching parameters in the namelist
! 'master_2Body', such as 'HiEnergySchwelleBarBar' etc.
!
! The 'elastic' contribution is defined to have the same particles in the
! outgoing channel as in the incomming. This includes intermediate resonance
! states, which are usually called 'inelastic'.
!
!***************************************************************************
program CrossSectionPlotterElast

  use version, only: PrintVersion
  use inputGeneral, only: readInputGeneral
  use particleProperties, only: initParticleProperties, hadron, PartName, validCharge
  use master_2Body, only: generateFinalState, HiEnergyContrib
  use particleDefinition
  use idTable, only: nucleon, pion, isBaryon, isMeson
  use preEventDefinition
  use mediumDefinition
  use mediumModule, only: mediumAt
  use energyCalc, only: energyDetermination
  use lorentzTrafo, only: lorentzCalcBeta
  use collisionTerm, only: forceDecays
  use twoBodyTools, only: isElastic

  implicit none

  type(particle), dimension(1:2) :: pair               ! incoming pair of particles
  type(particle), dimension(1,1:20) :: fState          ! produced final state
  type(particle), dimension(1,1) :: realVectorDummy
  logical :: HiFlag,collisionFlag,anti1,anti2
  integer :: i,j,k,q1,q2,id1,id2,id_out,q_out = 0
  real, dimension(0:3) :: momentum_calc
  real, dimension(1:3) :: betaToCM,betaToCF
  type(medium) :: tMedium
  real :: srts, sigTot, rHiEnergy
  real :: sigmaTotal, sigmaElast
  character(len=15) :: name(2)
  integer :: charge(1:3)

  integer, save :: Nevents  = 1000     ! number of events
  real, save    :: srts_min = 0.0      ! minimum sqrt(s)
  real, save    :: srts_max = 4.0      ! maximum sqrt(s)
  real, save    :: dp       = 0.1      ! momentum interval

  integer :: hiEnergyType, fail, totfail, pc1, pc2, pc3
  integer :: N_tot, N_nuc, N_pi, N_prongs, N_out_tot, N_nuc_in
  integer :: N_out(-1:2)                ! index = charge

  call PrintVersion
  call readInputGeneral
  call initParticleProperties

  call readinput

  betaToCF = 0.

  ! make ID_OUT particle stable
  hadron(id_out)%stability = 0

  pair%Id=(/id1,id2/)
  pair%charge=(/q1,q2/)
  pair%anti=(/anti1,anti2/)
  pair%pert=.false.
  pair(1)%event=1
  pair(2)%event=2

  pair(1)%mass=hadron(id1)%mass
  pair(2)%mass=hadron(id2)%mass

  pair(1)%pos = (/ 0., 0., 0.1 /)
  pair(2)%pos = (/ 0., 0., 0.0 /)

  pair(1)%mom(1:3)=0.

  !pair(1)%mom(0)=sqrt(pair(1)%mass**2+Dot_Product(pair(1)%mom(1:3),pair(1)%mom(1:3)))
  call energyDetermination(pair(1),(/0.,0.,0./))
  call energyDetermination(pair(2),(/0.,0.,0./))

  pair(1)%vel=pair(1)%mom(1:3)/pair(1)%mom(0)
  pair(2)%vel=pair(2)%mom(1:3)/pair(2)%mom(0)

  if (.not. validCharge(pair(1))) then
     write(*,'(A,I3,A,I3)') "Error: particle 1 with ID ", id1, " has invalid charge", q1
     stop
  end if
  if (.not. validCharge(pair(2))) then
     write(*,'(A,I3,A,I3)') "Error: particle 2 with ID ", id2, " has invalid charge", q2
     stop
  end if

  tMedium = mediumAt(pair(2)%pos)

!  if (isMeson(id_out)) srts_min = max(srts_min, hadron(pair(1)%ID)%mass+hadron(pair(1)%ID)%mass+hadron(id_out)%minmass)

  write(*,*) '***********************'
  write(*,*) 'positions:'
  write(*,*) pair(1)%pos
  write(*,*) pair(2)%pos
  write(*,*) 'momenta:'
  write(*,*) pair(1)%mom
  write(*,*) pair(2)%mom
  write(*,*) 'sqrt(s)=',sqrts(pair(1),pair(2))
  write(*,*) 'Total momentum=',pair(1)%mom+pair(2)%mom
  write(*,*) 'sqrts_min = ', srts_min
  write(*,*) '***********************'

  betaToCM = (/0.,0.,0./)


  open(23,file='xs.dat')
  write(23,'("#  >> ",A," + ",A," <<")') trim(name(2)),trim(name(1))
  write(23,'("# 1: pLab, 2: sqrt(s), 3: sigma_tot, 4: sigma_el, 5,6: (lo energy), 7,8: (hi energy), 9: rHi, 10: Ekin_lab")')
  write(23,'("#")')



  do i=1,200000
     pair(1)%mom(1:3) = (/ 0., 0., float(i)*dp /)
!     pair(1)%mom(3)= pair(1)%mom(3) * 10.**(1./50.)
!     pair(1)%mom(3)= pair(1)%mom(3) * 10.**(1./100.)
     pair(1)%mom(3)= pair(1)%mom(3) * 10.**(1./500.)
     pair(1)%mom(0)   = sqrt(pair(1)%mass**2+Dot_Product(pair(1)%mom(1:3),pair(1)%mom(1:3)))

     call energyDetermination(pair(1),betaToCF)
     pair(1)%vel=pair(1)%mom(1:3)/pair(1)%mom(0)

     momentum_calc = pair(1)%mom(0:3)+pair(2)%mom(0:3)

!!$     betaToCM = lorentzCalcBeta (momentum_calc)
!!$     call energyDetermination (pair(1), betaToCM)
!!$     pair(1)%vel=pair(1)%mom(1:3)/pair(1)%mom(0)

     srts=sqrts(pair(1),pair(2))

     if (srts < srts_min) cycle
     if (srts > srts_max) exit

     write(*,*) '--- srts = ',srts

     rHiEnergy = HiEnergyContrib(srts,pair%ID,pair%anti)

     sigmaTotal = 0.
     sigmaElast = 0.
     totfail = 0


     evtloop: do j=1,Nevents  ! event loop to produce final states

        fail = 0

        ! (1) generate final state
        do
           call setToDefault(fState)
           fState%ID = -1
           call generateFinalState(pair, fState(1,:),1.,1,0.,collisionFlag,HiFlag,HiEnergyType,sigTot_out=sigTot)
           if (collisionFlag) then
              exit
           else
              fail = fail + 1
              if (fail == 100) then
                 exit evtloop
              end if
           end if
        end do

        totfail = totfail + fail

        ! (2) force decays
        call ForceDecays(fState, realVectorDummy, 0.)

        ! (3) analysis
        N_tot = 0
        do k=lbound(fState,2),ubound(fState,2)

           if (fState(1,k)%ID < 0) exit
           if (fState(1,k)%ID == 0) cycle

           N_tot = N_tot + 1
        end do

        sigmaTotal = sigmaTotal + sigTot/Nevents

        if (N_tot == 2) then
           if (isElastic(pair,fState(1,1:2))) &
                sigmaElast = sigmaElast + sigTot/Nevents
        end if


     end do evtloop

     if (totfail > 0) write(*,*) '   fail = ',(1.0*totfail)/Nevents

     write(23,'(10G17.9)') absmom(pair(1)), srts, &
          sigmaTotal, sigmaElast, &
          0.0, 0.0, 0.0, 0.0, &
          rHiEnergy, kineticEnergy(pair(1))
     flush(23)

  end do

  close(23)

contains

  subroutine readInput
    use output, only: Write_ReadingInput

    NAMELIST /Plotter2/ id1,id2,q1,q2,anti1,anti2,id_out,q_out,Nevents,srts_max,srts_min,dp

    call Write_ReadingInput('Plotter2',0)
    rewind(5)
    read(5,nml=Plotter2)
    write(*,*) '  Id of first particle      : ', id1
    write(*,*) '  Id of second particle     : ', id2
    write(*,*) '  Charge of first particle  : ', q1
    write(*,*) '  Charge of second particle : ', q2
    write(*,*) '  antiparticle 1            : ', anti1
    write(*,*) '  antiparticle 2            : ', anti2
!!$    write(*,*) '  Id of outgoing particle   : ', id_out
!!$    write(*,*) '  Charge of outgoing part.  : ', q_out
    write(*,*)
    write(*,*) '  Number of events          : ', Nevents
    write(*,*) '  minimum sqrt(s)           : ', srts_min
    write(*,*) '  maximum sqrt(s)           : ', srts_max
    write(*,*) '  momentum interval dp      : ', dp

    name(1) = PartName(id1,q1,anti1)
    name(2) = PartName(id2,q2,anti2)

    write(*,*)
    write(*,*) '  >> ',trim(name(2)),' --> ',trim(name(1)),' <<'
    write(*,*)

    call Write_ReadingInput('Plotter2',1)


  end subroutine readInput


end program CrossSectionPlotterElast
