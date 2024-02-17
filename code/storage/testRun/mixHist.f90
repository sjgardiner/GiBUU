program MixHist

  ! this program uses Fortran 2003 features for the cli

  use hist

  integer :: iArg, nHist
  character(len=100) :: Arg, fNameOut
  type(histogram) :: H, Htmp
  real :: add,mul
  logical :: flagOK

  if (command_argument_count() < 2) then
     call get_command_argument(0, Arg)
     write(*,*) "Error: not enough arguments."
     write(*,*) "Usage:"
     write(*,*) trim(Arg)," outFile.dat inFile1.bin [inFile2.bin [inFile3.bin ...]]"
     stop
  end if

  ! iArg=0 returns the name of the program
  ! iArg=1 = name of output file

  call get_command_argument(1, fNameOut)
  if (len_trim(fNameOut) == 0) stop

  iArg = 2
  do
     call get_command_argument(iArg, Arg)
     if (len_trim(Arg) == 0) exit

     write(*,*) iArg,trim(Arg)

     write(*,*) "Reading file: '",trim(Arg),"' ..."
     if (iArg==2) then
        call FetchHist(H,trim(Arg),add=add,mul=mul,flagOK=flagOK)
        if (.not.flagOK) then
           write(*,*) 'ERROR while reading. Stop!'
           stop
        end if
     else
        call FetchHist(Htmp,trim(Arg),flagOK=flagOK)
        if (.not.flagOK) then
           write(*,*) 'ERROR while reading. Stop!'
           stop
        end if

        call sumHist(H,Htmp)

!        H%xExtreme(1) = min(H%xExtreme(1),Htmp%xExtreme(1))
!        H%xExtreme(2) = max(H%xExtreme(2),Htmp%xExtreme(2))
!        H%yVal = H%yVal + Htmp%yVal

     end if

     iArg = iArg+1

  end do
  nHist = iArg - 2

  write(*,*) 'nHist=',nHist

end program MixHist
