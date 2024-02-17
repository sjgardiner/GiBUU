!******************************************************************************
!****m* /initExternal
! NAME
! module initExternal
! PURPOSE
! Initializes a hadronic system according to an external data file.
!******************************************************************************
module initExternal

  implicit none
  private

  !****************************************************************************
  !****g* initExternal/inputFile
  ! SOURCE
  !
  character*1000, save :: inputFile='./source.inp'
  ! PURPOSE
  ! the absolute name of the input file with hadrons to be propagated.
  !
  ! possible values:
  ! * if not set, default is './source.inp'
  ! * if given, but does not contain '/':
  !   default is './[inputFile]'
  ! * otherwise: filename is absolute, including path
  !
  ! NOTE
  ! if you want to use the file 'XXX.inp' in the actual directory,
  ! give it as './XXX.inp'
  !****************************************************************************

  !****************************************************************************
  !****g* initExternal/DoPerturbative
  ! SOURCE
  logical, save :: DoPerturbative = .false.
  ! PURPOSE
  ! if true, the particles will be inserted into the perturbative particle
  ! vector, the real particles have to be initialized via some nucleus
  ! definition
  !****************************************************************************

  !****************************************************************************
  !****g* initExternal/NumberingScheme
  ! SOURCE
  integer, save :: NumberingScheme = 1
  ! PURPOSE
  ! The way, how particles%event will be numbered:
  ! * 1: event = iPart, i.e. the particle number in the ensemble
  !   (historical, but does not work for fullensemble)
  ! * 2: event = pert_numbering() or real_numbering()
  !   (good both for perturbative and real mode)
  !****************************************************************************

  !****************************************************************************
  !****g* initExternal/posSRC
  ! SOURCE
  logical, save :: posSRC = .false.
  ! PURPOSE
  ! If true, the position vectors of the proton and neutron from SRC will
  ! be sampled by Monte-Carlo. Relevant when the target nucleus was initialized
  ! before calling initializeExternal and if there are only proton and neutron
  ! in the external source.
  !****************************************************************************


  !****************************************************************************
  !****g* initExternal/flagPH
  ! SOURCE
  logical, save :: flagPH = .false.
  ! PURPOSE
  ! If true, a particle-hole excitation will be added to the target nucleus.
  ! The momentum of the particle is obtained by adding transverse momentum transfer along x-axis
  ! and from "-" momentum conservation.
  !****************************************************************************

  !****************************************************************************
  !****g* initExternal/pt
  ! SOURCE
  real, save :: pt = 0.
  ! PURPOSE
  ! Transverse momentum transfer to the struck nucleon (GeV/c).
  ! Relevant when flagPH=.true.
  !****************************************************************************

  public :: initializeExternal, ExternalIsPerturbative

  logical, save :: initFlag=.true.

contains

  !****************************************************************************
  !****f* initExternal/ExternalIsPerturbative
  ! NAME
  ! logical function ExternalIsPerturbative()
  ! PURPOSE
  ! Returns the value of DoPerturbative
  !****************************************************************************
  logical function ExternalIsPerturbative()
    ExternalIsPerturbative = DoPerturbative
  end function ExternalIsPerturbative

  !****************************************************************************
  !****s* initExternal/initInput
  ! NAME
  ! subroutine initInput
  ! PURPOSE
  ! Reads input out of jobcard. Namelist 'externalSystem'.
  !****************************************************************************
  subroutine initInput
    use output, only: Write_ReadingInput

    !**************************************************************************
    !****n* initExternal/externalSystem
    ! NAME
    ! NAMELIST externalSystem
    ! PURPOSE
    ! Includes the switches:
    ! * inputFile
    ! * DoPerturbative
    ! * NumberingScheme
    ! * posSRC
    ! * flagPH
    ! * pt
    !**************************************************************************
    NAMELIST /externalSystem/ inputFile, DoPerturbative, NumberingScheme, posSRC, flagPH, pt

    call Write_ReadingInput('externalSystem',0)
    rewind(5)
    read(5,nml=externalSystem)
    if (len_trim(inputFile)>0) then
       if (index(inputFile,"/")>0) then
          inputFile = trim(inputFile)
       else
          inputFile = './'//trim(inputFile)
       end if
    else
       inputFile = './source.inp'
    end if



    write(*,*) ' inputFile:  ', trim(inputFile)
    write(*,*) ' DoPerturbative : ', DoPerturbative
    write(*,*) ' NumberingScheme: ', NumberingScheme
    write(*,*) ' posSRC : ', posSRC
    write(*,*) ' flagPH : ', flagPH
    if (flagPH) then
       write(*,*) ' pt, GeV/c : ', pt
    end if
    call Write_ReadingInput('externalSystem',1)

    if(posSRC .and. flagPH) then
       write(*,*)' posSRC and flagPH can not be .true. simultaneoulsly !'
       stop
    end if

    if(DoPerturbative .and. flagPH) then
       write(*,*)' DoPerturbative and flagPH both .true. not implemented'
       stop
    end if

    initFlag = .false.

  end subroutine initInput

  !****************************************************************************
  !****s* initExternal/initializeExternal
  ! NAME
  ! subroutine initializeExternal(PartsReal,PartsPert)
  ! PURPOSE
  ! Read the particles from the file
  !****************************************************************************
  subroutine initializeExternal(PartsReal,PartsPert,targetNuc)
    use particleDefinition
    use insertion, only: GarbageCollection
    use output, only: Write_ReadingInput, Write_InitStatus
    use checks, only: ChecksSwitchRealRun
    use initSRC, only: DoInit_pnSRC, writeHistoSRC, writeProbSRC
    use nucleusDefinition
    use random, only: rn
    use residue, only: InitResidue, ResidueAddPH
    use collisionNumbering, only: real_numbering,pert_Numbering
    use constants, only : mn
    use RMF, only: getRMF_flag, g_rho
    use baryonPotentialMain, only: getsymmetryPotFlag_baryon
    use densitymodule, only: FermiMomAt

    type(particle), dimension(:,:), intent(inOut), target :: PartsReal
    type(particle), dimension(:,:), intent(inOut), target :: PartsPert
    type(tNucleus), pointer :: targetNuc

    type(particle), dimension(:,:), pointer :: pParts
    type(particle) :: StruckNuc

    real, dimension(3) :: r1,r2
    integer :: i,j,index,k,iens,IOS,count
    real :: alpha,pF
    logical :: flagSuccess

    write(*,*)
    call Write_InitStatus('External',0)
    if (initFlag) call initInput

    if(.not.flagPH) then

        call Write_ReadingInput(trim(inputFile),0)

        open(1,file=trim(inputFile),status='old')
        if(posSRC) open(2,file=trim(inputFile)//'_pos',status='unknown')

    end if

    call ChecksSwitchRealRun(.not.DoPerturbative)
    if (DoPerturbative) then
       pParts => PartsPert
    else
       pParts => PartsReal
    end if

    if(flagPH) then

        ! Initialize target residue determination (for analysis only):
        call InitResidue(size(pParts,dim=1),1,targetNuc%mass,targetNuc%charge)

        ensemble_loop_PH : do i=1,size(pParts,dim=1)

                 flagSuccess=.false.

                 attempt_loop : do j=1,100

                      index = int(rn()*targetNuc%mass)+1

                      StruckNuc=pParts(i,index)

                      alpha=(freeEnergy(StruckNuc)-StruckNuc%mom(3))/mn

                      StruckNuc%mom(1)=StruckNuc%mom(1)+pt

                      StruckNuc%mom(3)=(mn**2+StruckNuc%mom(1)**2+StruckNuc%mom(2)**2-(alpha*mn)**2)/(2.*alpha*mn)

                      if (getRMF_flag()) then
                         if (g_rho/=0.) then
                            pF=FermiMomAt(StruckNuc%pos,StruckNuc%charge)
                         else
                            pF=FermiMomAt(StruckNuc%pos)
                         end if
                      else
                         if (getsymmetryPotFlag_baryon()) then
                            pF=FermiMomAt(StruckNuc%pos,StruckNuc%charge)
                         else
                            pF=FermiMomAt(StruckNuc%pos)
                         end if
                      end if

                      if(absMom(StruckNuc) .gt. pF) then
                         ! Add a hole in the target nucleus (analysis only):
                         call ResidueAddPH(1000*i+1,pParts(i,index))
                         ! Create a particle excitation:
                         pParts(i,index)=StruckNuc
                         pParts(i,index)%event=real_numbering()
                         flagSuccess=.true.
                         exit attempt_loop
                      end if

                 end do attempt_loop

                 if(.not.flagSuccess) then
                    write(*,*)' Particle-hole configuration was not created after ',j,' attempts'
                    stop
                 end if

                 pParts(i,:)%firstevent=1000*i+1   ! needed for real particles to create holes in collisions

        end do ensemble_loop_PH

        return

    end if

    count = 0

    ensemble_loop : do i=1,size(pParts,dim=1)

       index=0

       if(posSRC) call DoInit_pnSRC(r1,r2)

       backspace(1)

       particle_loop : do

          index=index+1

          if(pParts(i,index)%id.gt.0) cycle   ! Find empty space

          if (index.gt.size(pParts,dim=2)) then
            write(*,*) 'Particle vector too small. Stop in initializeExternal.'
            stop
          end if

          call setToDefault(pParts(i,index)) ! set particle to its default values
          call setNumber(pParts(i,index)) ! give each particle a unique number

          read(1,*,IOSTAT=IOS) pParts(i,index)%id,pParts(i,index)%charge,&
                            &pParts(i,index)%mass,&
                            &(pParts(i,index)%pos(k), k=1,3),&
                            &(pParts(i,index)%mom(k), k=1,3),iens

          if (IOS.lt.0) then ! E.o.f. is reached
            call setToDefault(pParts(i,index))
            exit ensemble_loop
          end if

          if (iens.ne.i) then
            call setToDefault(pParts(i,index))
            exit particle_loop
          end if

          if(posSRC) then
             select case (pParts(i,index)%charge)
             case(1)
                pParts(i,index)%pos=r1
             case(0)
                pParts(i,index)%pos=r2
             case default
                write(*,*)' UNEXPECTED CHARGE OF PARTICLE IN SRC : ', pParts(i,index)%charge
             end select
             ! Modified input with positions:
             write(2,FMT=55) pParts(i,index)%id,pParts(i,index)%charge,&
                            &pParts(i,index)%mass,&
                            &(pParts(i,index)%pos(k), k=1,3),&
                            &(pParts(i,index)%mom(k), k=1,3),iens
55           format(i4,1x,i2,1x,f5.3,3(1x,f8.3),3(1x,f8.3),1x,i5)
          end if

          if (pParts(i,index)%id.lt.0) then
             pParts(i,index)%id=abs(pParts(i,index)%id)
             pParts(i,index)%anti=.true.
          end if

          select case (NumberingScheme)
          case (1)
             pParts(i,index)%event=index ! This has to be changed in the full-ensemble mode
          case (2)
             if(DoPerturbative) then
                pParts(i,index)%event=pert_numbering()
             else
                pParts(i,index)%event=real_numbering()
             end if
          case (3)
             pParts(i,index)%event=count
          end select

          pParts(i,index)%pert = DoPerturbative

          count = count + 1

       end do particle_loop

    end do ensemble_loop

    close(1)
    if(posSRC) close(2)

    call Write_ReadingInput(trim(inputFile),1)
    write(*,*) 'Particles read: ', count

    call GarbageCollection(pParts)

    call Write_InitStatus('External',1)

    if(posSRC) then
       call writeHistoSRC
       call writeProbSRC
    end if

  end subroutine initializeExternal


end module initExternal
