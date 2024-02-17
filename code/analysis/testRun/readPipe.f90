program readPipe
  ! This program tries to check the output into the named pipe (fifo)
  ! "FinalEvents.pipe"
  ! It reads the output and does some dummy output.

  use CallStack, only: TraceBack

  implicit none

  logical :: ex
  integer :: ios
  character(len=*), parameter :: filename ="FinalEvents.pipe"

  integer :: run, iPart, ID, charge, history, prod_id
  real :: perweight, enu
  real, dimension(1:3) :: pos
  real, dimension(0:3) :: mom


  inquire(file=filename,exist=ex)

  if (.not.ex) then
     write(*,*) 'ERROR: Pipe file "',trim(filename),'" does not exist.'
     write(*,*) 'Generate it with mkfifo.'
     call TRACEBACK('stop')
  end if
  open(47, file=filename, access="stream",action="read", &
       form='formatted', iostat=ios, status='old')
  if (ios /= 0) then
     write(*,*) 'ERROR: Opening ',trim(filename),' failed.'
     call TRACEBACK('stop')
  end if

  do
     read(47,*,iostat=ios) run, iPart, ID, charge, perweight, &
          pos(:), mom(:), history, prod_id, enu

     if (is_iostat_end(ios)) then
        ! End of file.
        exit
     else if (ios /= 0) then
        write(*,*) '!!! Input error !!!'
        cycle
     end if

     write(*,*) run, iPart, ID ! do some dummy output
  end do

  close(47)

end program readPipe
