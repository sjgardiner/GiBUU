!******************************************************************************
!****p* /PrintHAFT
! NAME
! program PrintHAFT
! PURPOSE
! Read a HAFT acceptance file and print (some) values to stdout
!
!******************************************************************************
program PrintHAFT

  implicit none

  integer, parameter :: nids  = 14

  type :: AcceptanceMatrix
     real(4), dimension(:), allocatable :: matrix     ! the actual matrix
     integer(4) :: xdim = 0, ydim = 0, zdim = 0       ! dimensions
     real(4) :: pmin = 0., thmin = 0., phmin = 0., &  ! limits
          pmax = 0., thmax = 0., phmax = 0.
     real :: dp = 0., dth = 0., dph = 0.              ! step sizes
  end type AcceptanceMatrix

  ! acceptance matrices for e+, e-, pi+, pi-, K+, K- and p
  type(AcceptanceMatrix), dimension(nids), save :: acc

  type :: ResParameterTable
     real(4), dimension(:,:), allocatable :: tab  ! the actual table
     integer(4) :: xdim, ydim                     ! dimensions
     real(4) :: pmin, pmax, thmin, thmax          ! limits
     real :: dp, dth                              ! step sizes
     logical :: resflag = .false.                 ! parameters loaded?
  end type ResParameterTable

  ! resolution parameter tables for e+ and e-
  type(ResParameterTable), dimension(nids), save :: par




  character(len=1000)           :: arg

  if (command_argument_count() /= 1) then
     write(*,*) 'Yoy have to give the filename as CLI argument. STOP!'
     stop
  end if

  call get_command_argument(1, arg)
  call doRead(trim(arg))

  call writeResolution(par(3),103)
  call writeResolution(par(2),102)

contains

  !****************************************************************************
  !****************************************************************************
  subroutine doRead(fName)
    character(*),intent(in) :: fname

    integer, parameter :: runit = 77
    character(len=80) :: comment
    real(4) :: sigpA(3), sigpB(3), sigth, sigph, XX0  ! resolution parameters
    integer(4) :: pid
    integer(4) :: ntab
    integer :: i, bins, bytes


    !    write(*,*) 'Reading file "',fname,'"...'

    write(*,'(A)') '======================================'
    write(*,'(A)') fname



    open(unit=runit,file=fname,access='stream',status='old',err=99)
    bytes=1
    read(runit,pos=bytes,err=100) comment
    bytes = bytes + 80
    write(*,'(a80)') comment
    write(*,'(A)') '--------------------------------------'
    read(runit,pos=bytes,err=100) sigpA(1:3), sigpB(1:3), sigth, sigph, XX0
    bytes = bytes + 9*4

    write(*,*) sigpA
    write(*,*) sigpB
    write(*,*) sigth
    write(*,*) sigph
    write(*,*) XX0
    write(*,'(A)') '--------------------------------------'


    do

       read(runit,pos=bytes,end=50,err=100) pid  ! break out if EOF reached
       bytes = bytes + 4

!       write(*,*) pid
!       exit

       if (pid>=0) then

          ! read acceptance matrix
          if (pid<1 .or. pid>nids) then
             write(6,*) 'acceptance not yet supported for PID ', pid, ' File = ',trim(fname)
             stop
          end if

          read(runit,pos=bytes,err=100) acc(pid)%xdim, acc(pid)%ydim, acc(pid)%zdim
          bytes = bytes + 3*4

          bins = acc(pid)%xdim * acc(pid)%ydim * acc(pid)%zdim

          read(runit,pos=bytes,err=100) acc(pid)%pmin, acc(pid)%pmax,   &
               acc(pid)%thmin, acc(pid)%thmax, &
               acc(pid)%phmin, acc(pid)%phmax
          bytes = bytes + 6*4

          allocate(acc(pid)%matrix(bins), source=0._4)
          read(runit,pos=bytes,err=100) acc(pid)%matrix(1:bins)
          bytes = bytes + bins*4

          write(6,*) 'Acceptance matrix for ID= ', pid
          write(6,*) 'dims= ',acc(pid)%xdim, ' ', acc(pid)%ydim, ' ', acc(pid)%zdim
          write(6,*) 'lims= ',acc(pid)%pmin, ' ', acc(pid)%pmax, ' ', acc(pid)%thmin, &
               ' ', acc(pid)%thmax, ' ', acc(pid)%phmin, ' ', acc(pid)%phmax
          write(6,*) 'size of matrix :', bins
          write(6,*) '--------------------------------------'
          acc(pid)%dp = (acc(pid)%pmax-acc(pid)%pmin)/real(acc(pid)%xdim)
          acc(pid)%dth = (acc(pid)%thmax-acc(pid)%thmin)/real(acc(pid)%ydim)
          acc(pid)%dph = (acc(pid)%phmax-acc(pid)%phmin)/real(acc(pid)%zdim)

       else

          ! read resolution parameters
          pid = -pid
          if (pid<1 .or. pid>nids) then
             write(6,*) 'resolution not yet supported for PID ', pid, ' File = ',trim(fname)
            stop
         end if

         read(runit,pos=bytes,err=100) par(pid)%xdim, par(pid)%ydim
         bytes = bytes + 2*4

         bins = par(pid)%xdim*par(pid)%ydim

         read(runit,pos=bytes,err=100) par(pid)%pmin, par(pid)%pmax, &
              par(pid)%thmin, par(pid)%thmax
         bytes = bytes + 4*4

         read(runit,pos=bytes,err=100) ntab ! nb. of parameter tables
         bytes = bytes + 4

         allocate(par(pid)%tab(ntab,bins), source=0._4)
         do i=1,ntab
            read(runit,pos=bytes,err=100) par(pid)%tab(i,1:bins)
            bytes = bytes + bins*4
         end do

         par(pid)%resflag = .true.

         write(6,*) 'Resolution tables for ID= ', pid
         write(6,*) 'dims= ',par(pid)%xdim, ' ', par(pid)%ydim
         write(6,*) 'lims= ',par(pid)%pmin, ' ', par(pid)%pmax, ' ', &
              par(pid)%thmin, ' ', par(pid)%thmax
         write(6,*) 'size of parameter tables :', ntab, ' x', bins
         write(6,*) '--------------------------------------'
         par(pid)%dp  = (par(pid)%pmax-par(pid)%pmin)   / real(par(pid)%xdim)
         par(pid)%dth = (par(pid)%thmax-par(pid)%thmin) / real(par(pid)%ydim)

      end if

    end do

50  close(runit)
!    readHAFTmatrix = bytes-1 ! return number of bytes read
    return

    ! Error opening or reading
99  close(runit)
    write(6,*) 'Open error on unit ', runit, ' File = ',trim(fname)
!    readHAFTMatrix = -1
    return

100 close(runit)
    write(6,*) 'Read error on unit ', runit, ' File = ',trim(fname)
!    readHAFTmatrix = -1
    return

  end subroutine doRead

  !****************************************************************************
  !****************************************************************************
  subroutine writeDummy()

    integer, parameter :: runit = 77
    character(len=80) :: comment
    real(4) :: sigpA(3), sigpB(3), sigth, sigph, XX0  ! resolution parameters



  end subroutine writeDummy


  !****************************************************************************
  !****************************************************************************
  subroutine writeResolution(par,iF)
    type(ResParameterTable), intent(in) :: par
    integer, intent(in) :: iF

    integer :: ix,iy, ilin, itab
    real :: x,y
    real, dimension(5) :: z

    write(*,*) 'writeResolution...',par%xdim,par%ydim

    do ix = 1,par%xdim
       x = par%pmin + (ix-0.5)*par%dp
       do iy = 1,par%ydim
          y = par%thmin + (iy-0.5)*par%dth
          ilin = ix+par%xdim*(iy-1)    ! linearized index
          do itab=1,5
             z(itab) = par%tab(itab,ilin)
          end do
          write(iF,*) ix,iy,x,y,z
       end do
    end do
  end subroutine writeResolution

end program PrintHAFT
