!******************************************************************************
!****m* /BoxAnalysis
! NAME
! module BoxAnalysis
! PURPOSE
! Do some analysis for the box init.
!
! NOTES
! * At the moment, only BoxAnalysis_FinalSum_Mult_SetX.dat,
!   E_SetX_finalSum.dat, p_SetX_finalSum.dat are able to accumulate statistics
!   by doing multiple runs per energy.
!   All other output is just from the last run.
!******************************************************************************
module BoxAnalysis

  use hist
  use histMP
  use histMC
  use Averager
  use TmunuDefinition, only: tTmunuNmu, fillTmunu, headTmunu

  implicit none

  private

  public :: DoBoxAnalysisTime


  type(histogram), save :: hMassRho
  type(histogramMC), save :: hMassDelta
  type(histogramMC), save :: hMCMomPion

  integer, parameter :: nSet = 5
  !****************************************************************************
  !****g* BoxAnalysis/useSet
  ! SOURCE
  logical, dimension(nSet), save :: useSet = (/ .false., .false., .false., .true., .true. /)
  ! PURPOSE
  ! Array to indicate, which particle set will be used for output
  !****************************************************************************
  type(histogramMP), dimension(nSet), save :: &
       hMP_ESet, hMP_pSet, hMP_ESet_final, hMP_pSet_final


  type(tTmunuNmu), dimension(:), allocatable, save :: arrTmunuNmu
  type(tTmunuNmu), dimension(2), save :: arrTmunuNmu_hadr

  type(tAverager), dimension(10,2) :: AveNpi
  type(tAverager), dimension(10) :: AveQ, AveNch, AveR, AveF, AveFpi

  logical, save :: initFlag=.true.
  logical, parameter :: do_Tmunu_pirho=.false.
  logical, parameter :: do_Tmunu_sub=.false.

  real, parameter :: widthZ = 0.1d0


  integer, save :: nCallsFinal = 0

  ! List of files handles:
  ! * 101?: BoxAnalysis_Mult_Set"//achar(48+iSet)//".dat, iSet=1,nSet
  ! * 1021: BoxAnalysis_Tmunu.dat
  ! * 1022: BoxAnalysis_Tmunu_Sub.dat
  ! * 1023: BoxAnalysis_Tmunu.pion.dat
  ! * 1024: BoxAnalysis_Tmunu.rho.dat
  ! * ????: BoxAnalysis_Tmunu."//intToChar4(i)//".dat, i=1,nEns
  ! * 1025: BoxAnalysis_Mult.dat
  ! * 1026: BoxAnalysis_MultDiff.dat
  ! * ????: p_Set'//achar(48+iSet)//'_'//intTochar4(timestep)//'.dat'
  ! * ????: E_Set'//achar(48+iSet)//'_'//intTochar4(timestep)//'.dat
  ! * ????: massRho_'//intTochar4(timestep)//'.dat
  ! * ????: massDelta_'//intTochar4(timestep)//'.dat
  ! * ????: BoxAnalysis_Final_Mult_Set"//achar(48+iSet)//".dat
  ! * ????: p_Set'//achar(48+iSet)//'_final.dat
  ! * ????: E_Set'//achar(48+iSet)//'_final.dat
  ! * ????: BoxAnalysis_FinalSum_Mult_Set"//achar(48+iSet)//".dat
  ! * ????: p_Set'//achar(48+iSet)//'_finalSum.dat
  ! * ????: E_Set'//achar(48+iSet)//'_finalSum.dat
  ! * 1031: BoxAnalysis_Tmunu.dat.bin
  ! * 1032: BoxAnalysis_Tmunu_Sub.dat.bin
  ! * 1033: BoxAnalysis_Tmunu.pion.dat.bin
  ! * 1034: BoxAnalysis_Tmunu.rho.dat.bin


  !****************************************************************************
  !****g* BoxAnalysis/do_Tmunu
  ! SOURCE
  logical,save :: do_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output. default: Only one file for all ensemble!
  ! you may change this with the flag perEnsemble_Tmunu
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/perEnsemble_Tmunu
  ! SOURCE
  logical,save :: perEnsemble_Tmunu=.false.
  ! PURPOSE
  ! Switch for Tmunu output. One file per ensemble!
  !
  ! NOTES
  ! this may slow down the execution dramatically, since huge output to the
  ! hard drive is induced.
  ! You may observe this, if e.g the cpu load drops permanently to 30%.
  ! Thus: switch it on, only if you want it!
  ! NOTES
  ! This is mutually exclusive with do_TmunuZ
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_TmunuZ
  ! SOURCE
  logical,save :: do_TmunuZ = .false.
  ! PURPOSE
  ! switch for Tmunu for every z-coordinate
  ! NOTES
  ! This is mutually exclusive with perEnsemble_Tmunu
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/selectTmunuFormat
  ! SOURCE
  integer,save :: selectTmunuFormat = 2
  ! PURPOSE
  ! select output format of Tmunu (binary encoded):
  ! * 1: ASCII
  ! * 2: Binary
  ! * 3: ASCII + Binary
  !****************************************************************************


  !****************************************************************************
  !****g* BoxAnalysis/do_P
  ! SOURCE
  logical,save :: do_P=.false.
  ! PURPOSE
  ! Switch for dN/p^2 dp output
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_velrel
  ! SOURCE
  logical,save :: do_velrel=.false.
  ! PURPOSE
  ! Switch for calculating velrel
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_posrel
  ! SOURCE
  logical,save :: do_posrel = .false.
  ! PURPOSE
  ! Switch for calculating relative distance
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_Cumulants
  ! SOURCE
  logical,save :: do_Cumulants=.false.
  ! PURPOSE
  ! Switch for calculating cumulants
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_Mult
  ! SOURCE
  logical,save :: do_Mult = .false.
  ! PURPOSE
  ! Switch for counting multiplicities
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/factorSubBoxVolume
  ! SOURCE
  real,save :: factorSubBoxVolume = 0.5
  ! PURPOSE
  ! Volume of the sub box relative to the full box
  !****************************************************************************

  !****************************************************************************
  !****g* BoxAnalysis/do_DumpPartVec
  ! SOURCE
  logical,save :: do_DumpPartVec = .false.
  ! PURPOSE
  ! Switch for writing the whole particle vector to file at the end of the run
  !****************************************************************************


  ! like in initBox, but here also starting values:
  real,dimension(101:122), save  :: densMes0 =-1, densMes =0
  real,dimension(1:2,1:61), save :: densBar0 =-1, densBar =0 ! Baryon, Anti-


contains

  !****************************************************************************
  !****s* BoxAnalysis/readInput
  ! NAME
  ! subroutine readInput
  !
  ! PURPOSE
  ! Reads input in jobcard out of namelist "BoxAnalysis"
  !****************************************************************************
  subroutine readInput

    use output, only: Write_ReadingInput
    use CallStack, only: TRACEBACK

    !**************************************************************************
    !****n* BoxAnalysis/BoxAnalysis
    ! NAME
    ! NAMELIST BoxAnalysis
    ! PURPOSE
    ! Includes the switches:
    ! * do_Tmunu
    ! * do_TmunuZ
    ! * perEnsemble_Tmunu
    ! * selectTmunuFormat
    ! * do_P
    ! * do_velrel
    ! * do_posrel
    ! * do_Cumulants
    ! * useSet
    ! * factorSubBoxVolume
    ! * do_DumpPartVec
    !**************************************************************************
    NAMELIST /BoxAnalysis/ &
         do_Tmunu, do_TmunuZ, perEnsemble_Tmunu, selectTmunuFormat, &
         do_P, do_velrel, do_posrel, &
         do_Cumulants, do_Mult, useSet, &
         factorSubBoxVolume, &
         do_DumpPartVec
    integer :: ios

    call Write_ReadingInput('BoxAnalysis',0)
    rewind(5)
    read(5,nml=BoxAnalysis,IOSTAT=ios)
    call Write_ReadingInput('BoxAnalysis',0,ios)

    write(*,*) '  do Tmunu:    ',do_Tmunu,'   perEnsemble: ',perEnsemble_Tmunu
    write(*,*) '  do TmunuZ:   ',do_TmunuZ
    write(*,*) '  Tmunu format:',selectTmunuFormat

    if (selectTmunuFormat==0) then
       call Traceback("selectTmunuFormat must not == 0")
    end if

    if (do_TmunuZ.and.perEnsemble_Tmunu) then
       call Traceback("do_TmunuZ and perEnsemble_Tmunu are mutually exclusive.")
    end if

    if (perEnsemble_Tmunu.and.selectTmunuFormat==2) then
       call Traceback("perEnsemble_Tmunu purely binary output not possible.")
    end if

    write(*,*) '  do P:        ',do_P
    write(*,*) '  do velrel:   ',do_velrel
    write(*,*) '  do posrel:   ',do_posrel
    write(*,*) '  do cumlants: ',do_cumulants
    write(*,*) '  do mult:     ',do_Mult
    write(*,*) '  use MP set : ',useSet
    write(*,*) '  factor SubBox volume =',factorSubBoxVolume
    write(*,*)
    write(*,*) '  DumpPartVec: ',do_DumpPartVec


    call Write_ReadingInput('BoxAnalysis',1)
  end subroutine readInput

  !****************************************************************************
  !****s* BoxAnalysis/DoBoxAnalysisTime
  ! NAME
  ! subroutine DoBoxAnalysisTime(realPart,timestep)
  ! PURPOSE
  !****************************************************************************
  subroutine DoBoxAnalysisTime(realPart,timestep)

    use CallStack, only: TRACEBACK
    use densityModule, only: gridsize
    use history, only: history_getParents
    use output, only: Write_InitStatus, intToChar4 !, WriteParticleVector
    use particleDefinition
    use collisionReporter, only: CR_write

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep

    real, save :: boxVol ! the volume of the box in fm^3

    integer, save :: nEns,nPart
    real, save :: mulfak

    integer :: i,j, iCh, iSet
    real :: mom0,mom, mass
    type(particle), POINTER :: pPart
    integer :: parents(1:3)

    type(tTmunuNmu) :: TmunuNmu, TmunuNmuSub

    real, dimension(10) :: countQ, countNch
    real, dimension(10,2) :: countNpi

    if (initFlag) then
       call Write_InitStatus('BoxAnalysis',0)
       initFlag=.false.

       call readInput

       if ((do_Tmunu) .and. (perEnsemble_Tmunu) .and. (nEns > 9999)) then
          call TRACEBACK("BoxAnalysis: Tmunu not prepared for more than 9999 ensembles.")
       end if

       nEns  = size(realPart,dim=1)
       nPart = size(realPart,dim=2)
       boxVol = 8.*gridsize(1)*gridsize(2)*gridsize(3)
       mulfak = 1.0/(nEns*boxVol)

       !----- mass distribution: -----
       call CreateHist(hMassRho, "mass(rho)", 0.0, 2.5, 0.01)
       call CreateHistMC(hMassDelta, "mass(Delta)", 1.0, 2.5, 0.005,4)
       hMassDelta%yDesc(1:4) = (/ "Delta- ", "Delta0 ", "Delta+ ", "Delta++" /)
       call CreateHistMC(hMCMomPion, "momentum(pion)", 0.0, 2.5, 0.02, 6)
       hMCMomPion%yDesc(1:6) = (/ "original  ",  &
            "rho       ", "sigma     ", "other dec ", &
            "pi pi     ", "other coll" /)

       !----- multiplicities: -----

       do iSet=1,nSet
          if (.not.useSet(iSet)) cycle
          call CreateHistMP(hMP_ESet(iSet), "dN/pE dE",  0.0, 3.0, 0.01, iSet)
          call CreateHistMP(hMP_pSet(iSet), "dN/p^2 dp", 0.0, 3.0, 0.01, iSet)
          call CreateHistMP(hMP_ESet_final(iSet), "dN/pE dE",  0.0, 3.0, 0.01, iSet)
          call CreateHistMP(hMP_pSet_final(iSet), "dN/p^2 dp", 0.0, 3.0, 0.01, iSet)
          open(1010+iSet,file="BoxAnalysis_Mult_Set"//achar(48+iSet)//".dat", status="unknown")
          call WriteHistMP_Names(iSet,1010+iSet)
          flush(1010+iSet)
       end do

       !----- hydro tensors: -----
       if (do_Tmunu) then

          if (iand(selectTmunuFormat,1)==1) then
             open(1021,file="BoxAnalysis_Tmunu.dat", status="unknown")
             write(1021,'(A)') headTmunu
             flush(1021)
          end if
          if (iand(selectTmunuFormat,2)==2) then
             open(1031,file="BoxAnalysis_Tmunu.dat.bin", status="unknown",form="unformatted")
             rewind(1031)
             flush(1031)
          end if

          if (do_Tmunu_sub) then
             if (iand(selectTmunuFormat,1)==1) then
                open(1022,file="BoxAnalysis_Tmunu_Sub.dat", status="unknown")
                write(1022,'(A)') headTmunu
                flush(1022)
             end if
             if (iand(selectTmunuFormat,2)==2) then
                open(1032,file="BoxAnalysis_Tmunu_Sub.dat.bin", status="unknown",form="unformatted")
                rewind(1032)
                flush(1032)
             end if
          end if

          if (do_Tmunu_pirho) then
             if (iand(selectTmunuFormat,1)==1) then
                open(1023,file="BoxAnalysis_Tmunu.pion.dat", status="unknown")
                write(1023,'(A)') headTmunu
                flush(1023)
             end if
             if (iand(selectTmunuFormat,2)==2) then
                open(1033,file="BoxAnalysis_Tmunu.pion.dat.bin", status="unknown",form="unformatted")
                rewind(1033)
                flush(1033)
             end if

             if (iand(selectTmunuFormat,1)==1) then
                open(1024,file="BoxAnalysis_Tmunu.rho.dat", status="unknown")
                write(1024,'(A)') headTmunu
                flush(1024)
             end if
             if (iand(selectTmunuFormat,2)==2) then
                open(1034,file="BoxAnalysis_Tmunu.rho.dat.bin", status="unknown",form="unformatted")
                rewind(1034)
                flush(1034)
             end if
          end if

          if (perEnsemble_Tmunu) then
             allocate( arrTmunuNmu(nEns) )
             do i=1,nEns
                open(123,file="BoxAnalysis_Tmunu."//intToChar4(i)//".dat", status="unknown")
                write(123,'(A)') headTmunu
                close(123)
             end do
          end if
       end if

       if (do_TmunuZ) then
          i = ceiling(2*gridsize(3)/widthZ)
          allocate( arrTmunuNmu(1:i) )
       end if

       !----- multiplicities: -----
       if (do_Mult) then
          open(1025,file="BoxAnalysis_Mult.dat", status="unknown")
          write(1025,'(A)') "# (to be implemented)"
          flush(1025)

          open(1026,file="BoxAnalysis_MultDiff.dat", status="unknown")
          write(1026,'(A)') "# (to be implemented)"
          flush(1026)
       end if

       call Write_InitStatus('BoxAnalysis',1)
    end if

    do iSet=1,nSet
       if (.not.useSet(iSet)) cycle
       call ClearHistMP(hMP_pSet(iSet))
       call ClearHistMP(hMP_ESet(iSet))
    end do

    call ClearHistMC(hMCMomPion)
    call ClearHist(hMassRho)
    call ClearHistMC(hMassDelta)

    if (do_Tmunu) then
       if (perEnsemble_Tmunu) arrTmunuNmu = TmunuNmu ! set all values to zero
       arrTmunuNmu_hadr = TmunuNmu ! set all values to zero
    end if

    if (do_TmunuZ) arrTmunuNmu = TmunuNmu ! set all values to zero

    if (do_cumulants) then
       countQ = 0.0
       countNch = 0.0
       countNpi = 0.0
    end if

    !
    ! accumulate data:
    !

    do i=1,nEns
       do j=1,nPart
          pPart => realPart(i,j)
          if (pPart%Id <  0) exit
          if (pPart%Id <= 0) cycle

          mom = absMom(pPart)
          mom0 = pPart%mom(0)

          select case (pPart%ID)
          case (101)
             parents = history_getParents(pPart%history)
             if (parents(2) == 0) then
                select case (parents(1))
                case (0)
                   iCh = 1
                case (103)
                   iCh = 2
                case (104)
                   iCh = 3
                case default
                   iCh = 4
                end select
             else
                if (parents(1)==101 .and. parents(2)==101) then
                   iCh = 5
                else
                   iCh = 6
                end if
             end if
             call AddHistMC(hMCMomPion, mom, iCh, 1.0/(mom**2))
          case (103)
             mass = sqrtS(pPart)
             call AddHist(hMassRho, mass, 1.0)
          case (2)
             if (.not.pPart%anti) then
                mass = sqrtS(pPart)
                call AddHistMC(hMassDelta, mass, pPart%charge+2, 1.0)
             end if
          end select

          do iSet=1,nSet
             if (.not.useSet(iSet)) cycle
             call AddHistMP(hMP_pSet(iSet), pPart, mom, 1.0/(mom**2), 1.0)
             call AddHistMP(hMP_ESet(iSet), pPart, mom0, 1.0/(mom0*mom), 1.0)
          end do

          ! fill Tmunu and Jmu:
          if (do_Tmunu) then
             call fillTmunu(TmunuNmu, pPart)
             if (do_Tmunu_sub) then
                if (isInSubBox(pPart)) then
                   call fillTmunu(TmunuNmuSub, pPart)
                end if
             end if

             if (do_Tmunu_pirho) then
                select case (pPart%ID)
                case (101)
                   call fillTmunu(arrTmunuNmu_hadr(1), pPart)
                case (103)
                   call fillTmunu(arrTmunuNmu_hadr(2), pPart)
                end select
             end if
             if (perEnsemble_Tmunu) call fillTmunu(arrTmunuNmu(i), pPart)
          end if

          if (do_TmunuZ) then
             if (abs(pPart%pos(3))<gridsize(3)) then
                iCh = int( (pPart%pos(3)+gridsize(3))/widthZ )+1
                call fillTmunu(arrTmunuNmu(iCh), pPart)
             end if
          end if

          if (do_cumulants) call calcCumulants1

          if (do_Mult) call calcMult
       end do
    end do

    if (do_Mult) then
       if (densMes0(101) < 0) then
          densMes0 = densMes
          densBar0 = densBar
       end if
    end if

    !
    ! produce output:
    !

    if (timestep>=0) then ! this is a real timestep
       call doOutputTimestep
    else
       call doOutputFinal
    end if



  contains

    !**************************************************************************
    subroutine calcCumulants1

      real, dimension(3) :: scalePos
      real :: mVol
      integer :: iVol

      ! scale the position to the box extends:
      scalePos(1:3) = abs(pPart%pos(1:3))/gridsize(1:3)

      ! the minimal box volume with all coordinates included:
      mVol = maxval(scalePos)**3

      ! select minimal volume bin:
      iVol = int(mVol*10)+1

      if ((iVol<1).or.(iVol>10)) then
         write(*,*) 'mVol,iVol: ',mVol,iVol
         write(*,*) 'ooops!!!!!!!!!'
      end if

      if (pPart%charge /= 0) then
         countQ(iVol:10) = countQ(iVol:10) + pPart%charge
         countNch(iVol:10) = countNch(iVol:10) + 1.0
      end if

      if (pPart%ID==101) then
         if (pPart%charge==-1) then
            countNpi(iVol:10,1) = countNpi(iVol:10,1) + 1.0
         else if (pPart%charge==1) then
            countNpi(iVol:10,2) = countNpi(iVol:10,2) + 1.0
         end if
      end if

    end subroutine calcCumulants1

    !**************************************************************************
    subroutine calcCumulants2

      integer :: ii
      real :: h

      do ii=1,10
         call AveragerAdd(AveQ(ii), countQ(ii))     ! Q = sum_i q_i
         call AveragerAdd(AveNch(ii), countNch(ii)) ! Nch = sum_i (q_i/=0?1:0)
         call AveragerAdd(AveNpi(ii,1), countNpi(ii,1)) ! = N(pi-)
         call AveragerAdd(AveNpi(ii,2), countNpi(ii,2)) ! = N(pi+)

         h = countNch(ii)
         if (h > 0.0) then
            h = countQ(ii)/h
         else
            h = 10000.0
         end if
         call AveragerAdd(AveF(ii), h ) ! = Q/Nch

         h = sum(countNpi(ii,:))
         if (h > 0.0) then
            h = (countNpi(ii,2)-countNpi(ii,1))/h
         else
            h = 10000.0
         end if
         call AveragerAdd(AveFpi(ii), h ) ! = (N(pi+)-N(pi-))/(N(pi+)+N(pi+))

         h = countNpi(ii,2)
         if (h > 0.0) then
            h = countNpi(ii,1)/h
         else
            h = 10000.0
         end if
         call AveragerAdd(AveR(ii), h ) ! = N(pi+)/N(pi-)

         write(7000+ii,*) timestep,CountNpi(ii,1)
         write(7100+ii,*) timestep,CountNpi(ii,2)

         flush(7000+ii)
         flush(7100+ii)

      end do

      do ii=1,9
         call AveragerWrite(5000+ii,timestep*1.0,AveNpi(ii,1))
         call AveragerWrite(5100+ii,timestep*1.0,AveNpi(ii,2))
         call AveragerWrite(5200+ii,timestep*1.0,AveFpi(ii))
         call AveragerWrite(5300+ii,timestep*1.0,AveR(ii))
         flush(5000+ii)
         flush(5100+ii)
         flush(5200+ii)
         flush(5300+ii)
      end do

      if (mod(timestep,10)==0) then
         rewind(6001)
         rewind(6002)
         rewind(6003)
         rewind(6004)
         do ii=1,10
            call AveragerWrite(6001,ii*0.1,AveNpi(ii,1))
            call AveragerWrite(6002,ii*0.1,AveNpi(ii,2))
            call AveragerWrite(6003,ii*0.1,AveFpi(ii))
            call AveragerWrite(6004,ii*0.1,AveR(ii))
         end do
         flush(6001)
         flush(6002)
         flush(6003)
         flush(6004)
      end if

    end subroutine calcCumulants2

    !**************************************************************************
    subroutine calcMult
      use IdTable, only: isMeson,isBaryon

      if (isMeson(pPart%ID)) then
         densMes(pPart%ID) = densMes(pPart%ID) + mulfak
      else if (isBaryon(pPart%ID)) then
         if (pPart%anti) then
            densBar(2,pPart%ID) = densBar(2,pPart%ID) + mulfak
         else
            densBar(1,pPart%ID) = densBar(1,pPart%ID) + mulfak
         end if
      end if
    end subroutine calcMult

    !**************************************************************************
    subroutine doOutputTimestep

      real :: h

      if (do_P) then
         if (mod(timestep,5)==1) then
            do iSet=1,nSet
               if (.not.useSet(iSet)) cycle

               call WriteHistMP(hMP_pSet(iSet), &
                    file='p_Set'//achar(48+iSet)//'_'//intTochar4(timestep)//'.dat', &
                    add=1e-20, mul=mulfak, iColumn=1)

               call WriteHistMP(hMP_ESet(iSet), &
                    file='E_Set'//achar(48+iSet)//'_'//intTochar4(timestep)//'.dat', &
                    add=1e-20, mul=mulfak, iColumn=1)

            end do

            ! call WriteHist(hMassRho, file='massRho_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
            ! call WriteHistMC(hMassDelta, file='massDelta_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
            ! call WriteParticleVector('parts_'//intTochar(timestep),realPart)
            ! call WriteHistMC(hMCMomPion, file='MomPion_'//intTochar4(timestep)//'.dat', add=1e-20, mul=mulfak)
         end if
      end if

!!$      if (mod(timestep,5)==1) then
!!$         call WriteHist(hMassRho, file='massRho.dat', add=1e-20, mul=mulfak)
!!$      end if
!!$
!!$      write(1040,'(1p,300e13.5)') timestep*1.0,(hMassRho%yVal(i,2),i=1,ubound(hMassRho%yVal,dim=1))
!!$      flush(1040)

      do iSet=1,nSet
         if (.not.useSet(iSet)) cycle

         write(1010+iSet,'(i11,1P,100E12.4,0P)') timestep, &
              sum(hMP_ESet(iSet)%yVal(:,:,2),dim=2)*mulfak, &
              0. ! 0 for historical reasons
         flush(1010+iSet)
      end do


      if (do_Tmunu) then

         if (iand(selectTmunuFormat,1)==1) then
            write(1021,'(i11,1P,100E16.8,0P)') timestep, &
              & TmunuNmu%Tmunu(:)*mulfak, &
              & TmunuNmu%Nmu(:)*mulfak, &
                 & TmunuNmu%Jmu(:)*mulfak, &
                 & TmunuNmu%B, TmunuNmu%S
            flush(1021)
         end if
         if (iand(selectTmunuFormat,2)==2) then
            write(1031) timestep, &
                 & TmunuNmu%Tmunu(:)*mulfak, &
                 & TmunuNmu%Nmu(:)*mulfak, &
                 & TmunuNmu%Jmu(:)*mulfak, &
                 & TmunuNmu%B, TmunuNmu%S
            flush(1031)
         end if

         if (do_Tmunu_sub) then
            if (iand(selectTmunuFormat,1)==1) then
               write(1022,'(i11,1P,100E16.8,0P)') timestep, &
                    & TmunuNmuSub%Tmunu(:)*mulfak, &
                    & TmunuNmuSub%Nmu(:)*mulfak, &
                    & TmunuNmuSub%Jmu(:)*mulfak, &
                    & TmunuNmuSub%B, TmunuNmuSub%S
               flush(1022)
            end if
            if (iand(selectTmunuFormat,2)==2) then
               write(1032) timestep, &
                    & TmunuNmuSub%Tmunu(:)*mulfak, &
                    & TmunuNmuSub%Nmu(:)*mulfak, &
                    & TmunuNmuSub%Jmu(:)*mulfak, &
                    & TmunuNmuSub%B, TmunuNmuSub%S
               flush(1032)
            end if
         end if


         if (do_Tmunu_pirho) then
            if (iand(selectTmunuFormat,1)==1) then
               write(1023,'(i11,1P,100E16.8,0P)') timestep, &
                 & arrTmunuNmu_hadr(1)%Tmunu(:)*mulfak, &
                 & arrTmunuNmu_hadr(1)%Nmu(:)*mulfak, &
                    & arrTmunuNmu_hadr(1)%Jmu(:)*mulfak, 0d0, 0d0
               flush(1023)
            end if
            if (iand(selectTmunuFormat,2)==2) then
               write(1033) timestep, &
                    & arrTmunuNmu_hadr(1)%Tmunu(:)*mulfak, &
                    & arrTmunuNmu_hadr(1)%Nmu(:)*mulfak, &
                    & arrTmunuNmu_hadr(1)%Jmu(:)*mulfak, 0d0, 0d0
               flush(1033)
            end if

            if (iand(selectTmunuFormat,1)==1) then
               write(1024,'(i11,1P,100E16.8,0P)') timestep, &
                 & arrTmunuNmu_hadr(2)%Tmunu(:)*mulfak, &
                 & arrTmunuNmu_hadr(2)%Nmu(:)*mulfak, &
                    & arrTmunuNmu_hadr(2)%Jmu(:)*mulfak, 0d0, 0d0
               flush(1024)
            end if
            if (iand(selectTmunuFormat,2)==2) then
               write(1034) timestep, &
                    & arrTmunuNmu_hadr(2)%Tmunu(:)*mulfak, &
                    & arrTmunuNmu_hadr(2)%Nmu(:)*mulfak, &
                    & arrTmunuNmu_hadr(2)%Jmu(:)*mulfak, 0d0, 0d0
               flush(1034)
            end if
         end if

         if (perEnsemble_Tmunu) then
            ! since we print all information for every ensemble, we must not divide by nEns here
            do i=1,nEns
               open(123,file="BoxAnalysis_Tmunu."//intToChar4(i)//".dat", status="old",position='append')
               write(123,'(i11,1P,100E12.4,0P)') timestep, &
                    & arrTmunuNmu(i)%Tmunu(:)/boxVol, &
                    & arrTmunuNmu(i)%Nmu(:)/boxVol, &
                    & arrTmunuNmu(i)%Jmu(:)/boxVol
               close(123)
            end do
         end if
      end if

      if (do_TmunuZ) then
         open(123,file="BoxAnalysis_TmunuZ."//intToChar4(timestep)//".dat", status="unknown")
         write(123,'(A)') headTmunu
         h = mulfak * (2*gridsize(3)/widthZ)
         do i=1,size(arrTmunuNmu)
            write(123,'(e12.4,1P,100E12.4,0P)') &
                 & -gridsize(3)+(real(i)-0.5)*widthZ, &
                 & arrTmunuNmu(i)%Tmunu(:)*h, &
                 & arrTmunuNmu(i)%Nmu(:)*h, &
                 & arrTmunuNmu(i)%Jmu(:)*h, 0d0, 0d0
         end do
         close(123)
      end if

      if (do_velrel) then
         call CalcAverageVelRel(realPart,timestep,nEns,nPart)
      end if

      if (do_posrel) then
         call CalcAveragePosRel(realPart,timestep,nEns,nPart)
      end if

      if (do_cumulants) call calcCumulants2

      if (do_Mult) then

         write(1025,'(i11,1P,500E12.4,0P)') timestep, &
              densBar(1,:)+densBar(2,:), &
              densMes(:)
         flush(1025)

         write(1026,'(i11,1P,500E12.4,0P)') timestep, &
              densBar(1,:)+densBar(2,:) - (densBar0(1,:)+densBar0(2,:)), &
              densMes(:) - densMes0(:)
         flush(1026)

         densMes = 0
         densBar = 0

      end if

      call cR_Write(1)
    end subroutine doOutputTimestep

    !**************************************************************************
    subroutine doOutputFinal

      use twoBodyStatistics, only: rate
      use insertion, only: DumpPartVec

      type(particle), dimension(1:2) :: dummy

      call rate(dummy,dummy,0.,.true.)

      if (do_cumulants) call calcCumulants2

      nCallsFinal = nCallsFinal+1

      close(1021)
      close(1022)
      close(1023)
      close(1024)
      close(1025)
      close(1026)
      close(1031)
      close(1032)
      close(1033)
      close(1034)

      do iSet=1,nSet
         if (.not.useSet(iSet)) cycle

         close(1010+iSet)

         open(123,file="BoxAnalysis_Final_Mult_Set"//achar(48+iSet)//".dat",&
              status="unknown")
         call WriteHistMP_Names(iSet,123)
         write(123,'(i11,1P,100E12.4,0P)') timestep, &
              sum(hMP_ESet(iSet)%yVal(:,:,2),dim=2)*mulfak, &
              0. ! 0 for historical reasons
         close(123)

         call WriteHistMP(hMP_pSet(iSet), &
              file='p_Set'//achar(48+iSet)//'_final.dat', &
              add=1e-20, mul=mulfak, iColumn=1)

         call WriteHistMP(hMP_ESet(iSet), &
              file='E_Set'//achar(48+iSet)//'_final.dat', &
              add=1e-20, mul=mulfak, iColumn=1)

         !===== Sum up different runs: =====

         call sumHistMP(hMP_ESet_final(iSet),hMP_ESet(iSet))
         call sumHistMP(hMP_pSet_final(iSet),hMP_pSet(iSet))

         open(123,file="BoxAnalysis_FinalSum_Mult_Set"//achar(48+iSet)//".dat",&
              status="unknown")
         call WriteHistMP_Names(iSet,123)
         write(123,'(i11,1P,100E12.4,0P)') timestep, &
              sum(hMP_ESet_final(iSet)%yVal(:,:,2),dim=2)*mulfak/nCallsFinal, &
              0. ! 0 for historical reasons
         close(123)

         call WriteHistMP(hMP_pSet_final(iSet), &
              file='p_Set'//achar(48+iSet)//'_finalSum.dat', &
              add=1e-20, mul=mulfak/nCallsFinal, iColumn=1)

         call WriteHistMP(hMP_ESet_final(iSet), &
              file='E_Set'//achar(48+iSet)//'_finalSum.dat', &
              add=1e-20, mul=mulfak/nCallsFinal, iColumn=1)

      end do

      if (do_DumpPartVec) then
         write(*,*) 'Dumping particle vector to file "PartVec.dat" ...'
         call DumpPartVec(realPart, 'PartVec.dat')
         write(*,*) 'Dumping particle vector to file "PartVec.dat" done.'
      end if


    end subroutine doOutputFinal

  end subroutine DoBoxAnalysisTime

  !****************************************************************************
  !****s* BoxAnalysis/CalcAverageVelRel
  ! NAME
  ! subroutine CalcAverageVelRel(realPart,timestep,nEns,nPart)
  ! PURPOSE
  ! calculate the average v_rel of all particles with each other. Due to
  ! speed reasons, it may be a good idea to correlate only particles in the
  ! same ensemble, but this s only a approximation.
  !****************************************************************************
  subroutine CalcAverageVelRel(realPart,timestep,nEns,nPart)
    use particleDefinition
!    use lorentzTrafo, only: eval_sigmaBoost

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep
    integer, intent(in) :: nEns,nPart

    integer :: iEns1,iEns2, iPart1,iPart2
    type(particle), POINTER :: pPart1, pPart2
    real :: sum0,sum1,velrel,m1,m2,s
    real :: ptot(0:3)
    real :: vrel, velMol
    real, dimension(1:3) :: vrel_vector


    ! due to speed reasons, I only correlate particles in the same ensemble

!    write(*,*) 'calculating velrel....'
    sum0 = 0.0
    sum1 = 0.0
    do iEns1=1,nEns
       do iPart1=1,nPart
          pPart1 => realPart(iEns1,iPart1)
          if (pPart1%Id <  0) exit
          if (pPart1%Id <= 0) cycle
          m1 = pPart1%mass**2

          iEns2 = iEns1
          do iPart2=iPart1+1,nPart

             pPart2 => realPart(iEns2,iPart2)
             if (pPart2%Id <  0) exit
             if (pPart2%Id <= 0) cycle
             m2 = pPart2%mass**2

             ptot = pPart1%mom + pPart2%mom
             s = ptot(0)**2-sum(ptot(1:3)**2)
             ! this is v_Moller:
             velMol = sqrt( max(0.0,(s-m1-m2)**2/4-m1*m2) )/( pPart1%mom(0)*pPart2%mom(0) )

             ! this is the 'real' v_rel:
             velrel = sqrt( max(0.0,(s-m1-m2)**2-4*m1*m2) )/( s-m1-m2 )

             vrel_vector=pPart1%vel-pPart2%vel
             vrel = sqrt(Dot_product(vrel_vector,vrel_vector))

!!$             write(*,*) velMol,velrel,vrel,vrel*eval_sigmaBoost(pPart1%mom,pPart2%mom), &
!!$                  (velMol-vrel*eval_sigmaBoost(pPart1%mom,pPart2%mom))/velMol, &
!!$                  (velrel-vrel*eval_sigmaBoost(pPart1%mom,pPart2%mom))/velrel


             sum0 = sum0 + 1
             sum1 = sum1 + velrel

          end do
       end do
    end do

    write(*,'(A,f8.5,f15.0,i13)') 'Average velrel: ',sum1/sum0, sum0, timestep

  end subroutine CalcAverageVelRel

  !****************************************************************************
  !****s* BoxAnalysis/CalcAveragePosRel
  ! NAME
  ! subroutine CalcAverageVelPos(realPart,timestep,nEns,nPart)
  ! PURPOSE
  ! calculate the average relative distance of particles in the same ensemble.
  !****************************************************************************
  subroutine CalcAveragePosRel(realPart,timestep,nEns,nPart)
    use particleDefinition
    use densityModule, only: gridsize
    use output, only: intToChar4

    type(particle),dimension(:,:),intent(in), target :: realPart
    integer, intent(in) :: timestep
    integer, intent(in) :: nEns,nPart

    integer :: iEns1,iEns2, iPart1,iPart2
    type(particle), POINTER :: pPart1, pPart2
    real, dimension(3) :: d2,h2
    real :: sum0, sum1, dd

    type(histogram) :: hPosRel

    if (mod(timestep,10).ne.1) return


    call CreateHist(hPosRel, "relative distance", 0.0, 10.0, 0.02)

!    write(*,*) 'calculating posrel....'

    sum0 = 0.0
    sum1 = 0.0
    do iEns1=1,nEns
       do iPart1=1,nPart
          pPart1 => realPart(iEns1,iPart1)
          if (pPart1%Id <  0) exit
          if (pPart1%Id <= 0) cycle

          iEns2 = iEns1
          do iPart2=iPart1+1,nPart

             pPart2 => realPart(iEns2,iPart2)
             if (pPart2%Id <  0) exit
             if (pPart2%Id <= 0) cycle

             d2 = (pPart1%pos(1:3) - pPart2%pos(1:3))**2

             h2 = (pPart1%pos(1:3) - pPart2%pos(1:3) - 2*gridsize(1:3))**2
             where (h2 < d2)
                d2 = h2
             end where
             h2 = (pPart1%pos(1:3) - pPart2%pos(1:3) + 2*gridsize(1:3))**2
             where (h2 < d2)
                d2 = h2
             end where

             dd = sqrt(sum(d2))


             sum0 = sum0 + 1
             sum1 = sum1 + dd

             call AddHist(hPosRel, dd, 1.0/dd**2)

          end do
       end do
    end do

    write(*,'(A,f8.5,f15.0,i13)') 'Average posrel: ',sum1/sum0, sum0, timestep

    call WriteHist(hPosRel, file='PosRel_'//intTochar4(timestep)//'.dat', &
         add=1e-20, mul=1./size(realPart,dim=1))

  end subroutine CalcAveragePosRel


  !**************************************************************************
  !**************************************************************************
  pure logical function isInSubBox(Part)

    use particleDefinition
    use densityModule, only: gridsize

    type(particle), intent(in) :: Part

    real, dimension(3) :: scalePos

    ! scale the position to the box extends:
    scalePos(1:3) = abs(Part%pos(1:3))/gridsize(1:3)

    isInSubBox = (maxval(scalePos)**3 < factorSubBoxVolume)
  end function isInSubBox


end module BoxAnalysis
