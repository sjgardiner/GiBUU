program readFinalEvents

  use hist

  implicit none

  integer, parameter :: nPtot = 100

  integer, dimension(nPtot) :: run, iPart, ID, charge, history, prod_id
  real, dimension(nPtot) :: perweight, enu
  real, dimension(nPtot,3) :: pos
  real, dimension(nPtot,0:3) :: mom

  integer :: iP, nP, nRun
  real :: addFak, mulFak
  integer :: ios
  logical :: doIt

  character(len=300) :: line

  type(histogram):: hPmu

  call createHist(hPmu, "P_mu", 0.,10.,0.1)

  read (*,'(A)',iostat=ios) line

  nP = 0
  doIt = .false.
  do
     nP = nP+1
     read(*,*,iostat=ios) run(nP), iPart(nP), ID(nP), charge(nP), perweight(nP), &
          pos(nP,:), mom(nP,:), history(nP), prod_id(nP), enu(nP)

     if (ios /= 0) doIt = .true.
     if ((nP>1).and.( (run(nP)/=run(nP-1)) .or. (iPart(nP)/=iPart(nP-1)) )) doIt = .true.

     if (doIt) then
        ! The event is stored in the entries (1:nP-1)...
!!$        do iP=1,nP-1
!!$           write(*,*) run(iP), iPart(iP), ID(iP), charge(iP), perweight(iP)
!!$        end do
!!$        write(*,*)

        if (accept()) then
           call AddHist(hPmu, mom(1,0), perweight(1))
        end if

        nRun = run(1)
        nP = 0
        doIt = .false.

        if (ios /= 0) exit
        backspace(5) ! undo the last read
     end if

  enddo

  addFak = 1e-20
  mulFak = 1.0/nRun

  call writeHist(hPmu, add=addFak,mul=mulFak, file="h.dat")


contains

  logical function accept()

    accept = .false.

    do iP=1,nP-1
       if ((ID(iP)>100) .and. (ID(iP)<200)) return ! failure
    end do

    accept = .true.

  end function accept


end program readFinalEvents
