program HeavyIonMixRuns

  use callstack, only: traceback

  implicit none

  character(len=:), allocatable :: Dirs(:)
  integer :: nDirs = 0

  integer :: timestep

  ! ===== Build directory list =====

  call BuildDirs
  if (nDirs==0) then
     call traceback('List of directories is empty. STOP!')
  end if

  ! ===== Mix Tmunu output =====
  do timestep=0, 1000, 20
     call mixTmunu(timestep)
  end do

  call mixTmunu(   5)


contains
  !****************************************************************************
  subroutine mixTmunu(timestep)

    use TmunuDefinition
    use output, only: intTochar, intTochar4

    integer, intent(in) :: timestep

    logical, parameter :: readBinary = .true.
    logical, parameter :: writeBinary = .true.

    integer, parameter :: nBin = 100
    integer, parameter :: nArr = 2

    integer :: iDir, iBin, iArr, ios, nRead
    real :: x, mulFak

    type(tTmunuNmu), save :: tTmunuNmu0 ! used to reset the array !

    type(tTmunuNmu), dimension(:), allocatable, save :: ArrX, ArrY, ArrZ
    type(tTmunuNmu) :: tT
    real, dimension(:), allocatable, save :: xVal,yVal,zVal

    logical, save :: isFirst = .true.

    if (isFirst) then 
       allocate(ArrX(nBin))
       allocate(ArrY(nBin))
       allocate(ArrZ(nBin))
       allocate(xVal(nBin))
       allocate(yVal(nBin))
       allocate(zVal(nBin))

       isFirst = .false.
    end if


    do iArr=1,nArr
       if (readBinary) then
          if (.not.TryFile('Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat.bin')) cycle
       else
          if (.not.TryFile('Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat')) cycle
       end if

       ArrX = tTmunuNmu0
       ArrY = tTmunuNmu0
       ArrZ = tTmunuNmu0
       xVal = 0.
       yVal = 0.
       zVal = 0.

       nRead = 0

       do iDir=1,nDirs
          if (readBinary) then
             open(123,file=trim(Dirs(iDir))//'/Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat.bin', &
                  status="old", form="unformatted",iostat=ios)

             write(*,*) 'reading: ',trim(Dirs(iDir))//'/Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat.bin'
          else
             open(123,file=trim(Dirs(iDir))//'/Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat', &
                  status="old",iostat=ios)
          end if

          if (ios.ne.0) then
             write(*,*) 'error opening file ',trim(Dirs(iDir)),iArr,timestep
             cycle
          end if

          nRead = nRead+1

          do iBin=1,nBin
             if (readBinary) then
                read(123,iostat=ios) x, &
                     tT%Tmunu(:), tT%Nmu(:), tT%Jmu(:), tT%B, tT%S
             else
                read(123,*,iostat=ios) x, &
                     tT%Tmunu(:), tT%Nmu(:), tT%Jmu(:), tT%B, tT%S
             end if
             if (ios.ne.0) cycle

             xVal(iBin) = x
             ArrX(iBin) = ArrX(iBin) + tT
          end do

          do iBin=1,nBin
             if (readBinary) then
                read(123,iostat=ios) x, &
                     tT%Tmunu(:), tT%Nmu(:), tT%Jmu(:), tT%B, tT%S
             else
                read(123,*,iostat=ios) x, &
                     tT%Tmunu(:), tT%Nmu(:), tT%Jmu(:), tT%B, tT%S
             end if
             if (ios.ne.0) cycle

             yVal(iBin) = x
             ArrY(iBin) = ArrY(iBin) + tT
          end do

          do iBin=1,nBin
             if (readBinary) then
                read(123,iostat=ios) x, &
                     tT%Tmunu(:), tT%Nmu(:), tT%Jmu(:), tT%B, tT%S
             else
                read(123,*,iostat=ios) x, &
                     tT%Tmunu(:), tT%Nmu(:), tT%Jmu(:), tT%B, tT%S
             end if
             if (ios.ne.0) cycle

             zVal(iBin) = x
             ArrZ(iBin) = ArrZ(iBin) + tT
          end do

          close(123)

       end do

       mulFak = 1.0/nRead


       if (writeBinary) then
          open(123,file='Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat.bin', &
               status="unknown", form="unformatted")

          do iBin=1,nBin
             write(123) xVal(iBin), &
                  & ArrX(iBin)%Tmunu(:)*mulfak, &
                  & ArrX(iBin)%Nmu(:)*mulfak, &
                  & ArrX(iBin)%Jmu(:)*mulfak, &
                  & ArrX(iBin)%B*mulfak, ArrX(iBin)%S*mulfak
          end do

          ! how to separate data blocks in binary files????

          do iBin=1,nBin
             write(123) yVal(iBin), &
                  & ArrY(iBin)%Tmunu(:)*mulfak, &
                  & ArrY(iBin)%Nmu(:)*mulfak, &
                  & ArrY(iBin)%Jmu(:)*mulfak, &
                  & ArrY(iBin)%B*mulfak, ArrY(iBin)%S*mulfak
          end do

          do iBin=1,nBin
             write(123) zVal(iBin), &
                  & ArrZ(iBin)%Tmunu(:)*mulfak, &
                  & ArrZ(iBin)%Nmu(:)*mulfak, &
                  & ArrZ(iBin)%Jmu(:)*mulfak, &
                  & ArrZ(iBin)%B*mulfak, ArrZ(iBin)%S*mulfak
          end do

          close(123)
       else
          open(123,file='Tmunu_'//Achar(48+iArr)//'_'//intTochar4(timestep)//'.dat', status="unknown")
          write(123,'(A)') headTmunu

          do iBin=1,nBin
             write(123,'(f11.4,1P,100E14.6,0P)') xVal(iBin), &
                  & ArrX(iBin)%Tmunu(:)*mulfak, &
                  & ArrX(iBin)%Nmu(:)*mulfak, &
                  & ArrX(iBin)%Jmu(:)*mulfak, &
                  & ArrX(iBin)%B*mulfak, ArrX(iBin)%S*mulfak
          end do
          write(123,*)
          write(123,*)

          do iBin=1,nBin
             write(123,'(f11.4,1P,100E14.6,0P)') yVal(iBin), &
                  & ArrY(iBin)%Tmunu(:)*mulfak, &
                  & ArrY(iBin)%Nmu(:)*mulfak, &
                  & ArrY(iBin)%Jmu(:)*mulfak, &
                  & ArrY(iBin)%B*mulfak, ArrY(iBin)%S*mulfak
          end do
          write(123,*)
          write(123,*)

          do iBin=1,nBin
             write(123,'(f11.4,1P,100E14.6,0P)') zVal(iBin), &
                  & ArrZ(iBin)%Tmunu(:)*mulfak, &
                  & ArrZ(iBin)%Nmu(:)*mulfak, &
                  & ArrZ(iBin)%Jmu(:)*mulfak, &
                  & ArrZ(iBin)%B*mulfak, ArrZ(iBin)%S*mulfak
          end do

          close(123)
       end if


    end do

  end subroutine mixTmunu


  !****************************************************************************
  logical function TryFile(file)
    implicit none
    character*(*),  intent(in)          :: file

    integer :: ios

    write(*,*) 'Trying: ',trim(Dirs(1))//"/"//trim(file)

    open(9932,file=trim(Dirs(1))//"/"//trim(file),status='OLD',IOSTAT=ios)
    close(9932)
    TryFile = (ios.eq.0)

  end function TryFile

  !****************************************************************************
  subroutine BuildDirs()
    implicit none

    integer :: ios
    character(len=1000) :: line
    integer :: nLen, i

    open(9932,file='directories.txt',status='OLD',IOSTAT=ios)
    if (ios/=0) then
       write(*,*) "WARNING: file 'directories.txt' does not exists."
       write(*,*) "Reading directories '0/' ... '9/' instead."

       nDirs=10
       allocate(character(1) :: Dirs(nDirs))
       do i=1,10
          Dirs(i)(1:1) = Achar(i-1+48)
       end do

!!$       do i=1,nDirs
!!$          write(*,*) i,': >',trim(Dirs(i)),'<'
!!$       end do

       return
    end if

    write(*,*) "reading 'directories.txt'..."
    nDirs = 0
    nLen = 0
    do
       read (9932,'(A)',iostat=ios) line
       if (ios/=0) exit
       nDirs = nDirs+1
       nLen = max(nLen, len(trim(line)))
    end do
    write(*,*) 'nDirs = ',nDirs,'   nLen=',nLen

    allocate(character(nLen) :: Dirs(nDirs))

    rewind(9932)
    nDirs = 0
    do
       read (9932,'(A)',iostat=ios) line
       if (ios/=0) exit
       nDirs = nDirs+1
       nLen = len(trim(line))
       if (line(nLen:nLen)=='/') nLen=nLen-1
       Dirs(nDirs) = trim(line(:nLen))
    end do
    write(*,*) "reading 'directories.txt' done"

!!$    do i=1,nDirs
!!$       write(*,*) i,': >',trim(Dirs(i)),'<'
!!$    end do


  end subroutine BuildDirs
  !****************************************************************************


end program HeavyIonMixRuns
