! This is a program to sum different HistMC (written in binary mode)
program BinSumHistMC_avg

  use histMC_avg

  implicit none

  character(len=1000) :: arg, arg0, dir, file
  integer :: nDirs, ios
  logical :: flagOK
  type(histogramMC_avg) :: hSum, hTmp
  real :: addFak, mulFak

  if (command_argument_count() /= 1) then
     write(*,*) 'you have to give the name of the file.'
     stop
  end if

  call get_command_argument(1, arg)

  if (index(arg,".bin")==0) then
     arg0 = arg
     arg = trim(arg)//".bin"
  else
     arg0 = arg(1:index(arg,".bin")-1)
  end if

  write(*,*) 'processing:',trim(arg)

  call system('ls -d ./*/ > dirs.txt 2>/dev/null')
  open(31,FILE='dirs.txt',action="read")

  nDirs = 0
  do
     read(31,FMT='(a)',iostat=ios) dir
     if (ios/=0) exit

     write(*,*) 'reading: ',trim(dir)//trim(arg)
     if (nDirs==0) then
        call fetchHistMC_avg(hSum,file=trim(dir)//trim(arg), &
             add=addFak,mul=mulFak, flagOK=flagOK)
        write(*,*) 'add,mul = ',addFak,mulFak
     else
        call fetchHistMC_avg(hTmp,file=trim(dir)//trim(arg), flagOK=flagOK)
     end if
     if (.not.flagOK) cycle

     if (nDirs>0) then
        call sumHistMC_avg(hSum, hTmp)
     end if

     nDirs = nDirs+1

  end do

  write(*,*) 'count = ',nDirs

  call writeHistMC_avg(hSum, file=trim(arg0)) !, add=addFak,mul=mulFak/nDirs)

end program BinSumHistMC_avg
