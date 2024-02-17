module aux

contains
  
    !*************************************************************************
    function intTochar(zahl)
      !wandelt integer in character um

      implicit none
      integer, intent(in)          :: zahl
      character(3) :: intToChar

      integer einser,zehner,hunderter

      hunderter=Mod(zahl/100,10)
      zehner=Mod(zahl/10,10)
      einser=Mod(zahl,10)

      intToChar=Achar(hunderter+48)//Achar(zehner+48)//Achar(einser+48)

!      if (zahl < 9) then
!         intToChar=Achar(einser+48)
!      else if (zahl >= 9 .and. zahl < 99) then
!         intToChar=Achar(zehner+48)//Achar(einser+48)
!      else
!         intToChar=Achar(hunderter+48)//Achar(zehner+48)//Achar(einser+48)
!      endif

    end function intTochar

  function intTochar4(zahl)
    !wandelt integer in character um

    implicit none
    integer, intent(in)          :: zahl
    character(4) :: intToChar4

    integer tausender

    tausender=Mod(zahl/1000,10)
    intToChar4=Achar(tausender+48)//intToChar(zahl-1000*tausender)

  end function intTochar4

end module aux


! Averaging over runs of the output in file dens.dat 
program densAver

!  use output, only: intTochar, intToChar4

  use aux
  
  implicit none

  integer, parameter ::  numTimeSteps        = 300
  real, dimension(1:numTimeSteps) ::  time, dens, temp, mub, Qzz
  real :: dens_inp, rhoz_antibar, p_inp(0:3), baryon_number, charge_number, strangeness, temp_inp, mub_inp
  real :: Qzz_inp
  integer :: idir,ndir,i,ndir1
  logical :: exist_file
  character(100) :: Directory,file
  
  ndir=70

  dens=0.
  temp=0.
  mub=0.
  Qzz=0.
  ndir1=0
  
  file='dens.dat'

  do idir=1,ndir

        write(6,*)' idir: ', idir       

        if(idir.le.999) then
           Directory=trim(intToChar(idir))
        else
           Directory=trim(intToChar4(idir))
        end if

        inquire(file=trim(Directory)//'/'//trim(file),exist=exist_file)
        if(.not.exist_file) then
           write(6,*)' file ', trim(file), ' does not exist'       
           cycle
        end if

        ndir1=ndir1+1
        
        open(1,file=trim(Directory)//'/'//trim(file),status='old',action='read')
        read(1,*)

        do i=1,numTimeSteps
           read(1,5) time(i), dens_inp, rhoz_antibar, p_inp(:), baryon_number, charge_number, strangeness,&
                   & temp_inp, mub_inp, Qzz_inp
5           format(1P,50(1x,e13.6))
            dens(i)=dens(i)+dens_inp
            temp(i)=temp(i)+temp_inp
            mub(i)=mub(i)+mub_inp
            Qzz(i)=Qzz(i)+Qzz_inp
        end do

        close(1)
        
  end do

  
  dens=dens/float(ndir1)
  temp=temp/float(ndir1)
  mub=mub/float(ndir1)
  Qzz=Qzz/float(ndir1)

  open(2,file='densAver.dat',status='unknown')
  write(2,*)'# t, fm/c:     dens_bar, fm^-3:      temp, GeV:     mub, GeV:   Qzz, (GeV/c)^2: '
  do i=1,numTimeSteps  
     write(2,5) time(i), dens(i), temp(i), mub(i), Qzz(i)
  end do
  close(2)

end program densAver
