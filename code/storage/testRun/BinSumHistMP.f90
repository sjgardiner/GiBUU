! This is a program to sum different HistMP (written in binary mode)
program BinSumHistMP

  use histMP

  implicit none

  character(len=1000) :: arg, arg0, dir, file
  integer :: nDirs, ios, iC
  logical :: flagOK
  type(histogramMP) :: hSum, hTmp
  real :: addFak, mulFak

  if (command_argument_count() < 1) then
     write(*,*) 'you have to give (at least) the name of the file.'
     stop
  end if

  if (command_argument_count() == 2) then
     call get_command_argument(2, arg)
     write(*,*) 'list of directories provided by user: ',trim(arg)
     open(31,FILE=trim(arg),action="read",iostat=ios)
     if (ios /= 0) then
        write(*,*) 'input not a vaild file.'
        stop
     end if
  else
     call system('ls -d ./*/ > dirs.txt 2>/dev/null')
     open(31,FILE='dirs.txt',action="read")
  end if

  call get_command_argument(1, arg)
  if (index(arg,".bin")==0) then
     arg0 = arg
     arg = trim(arg)//".bin"
  else
     arg0 = arg(1:index(arg,".bin")-1)
  end if

  write(*,*) 'processing: ',trim(arg)


  nDirs = 0
  do
     read(31,FMT='(a)',iostat=ios) dir
     if (ios/=0) exit

     ! remove comments:
     iC = index(dir, "#")
     if (iC > 0) dir(iC:) = " "

     ! skip lines which do not contain '/' (including empty lines):
     iC = index(dir, "/")
     if (iC == 0) cycle

     write(*,*) 'reading: ',trim(dir)//trim(arg)
     if (nDirs==0) then
        call fetchHistMP(hSum,file=trim(dir)//trim(arg), &
             add=addFak,mul=mulFak, flagOK=flagOK)
        write(*,*) 'add,mul = ',addFak,mulFak
     else
        call fetchHistMP(hTmp,file=trim(dir)//trim(arg), flagOK=flagOK)
     end if
     if (.not.flagOK) cycle

     if (nDirs>0) then
        call sumHistMP(hSum, hTmp)
     end if

     nDirs = nDirs+1

  end do

  write(*,*) 'count = ',nDirs

  call writeHistMP(hSum, file=trim(arg0), add=addFak,mul=mulFak/nDirs)

end program BinSumHistMP
