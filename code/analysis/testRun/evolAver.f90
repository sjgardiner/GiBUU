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


! Averaging over runs of the output in file evol.dat 
program evolAver

!  use output, only: intTochar, intToChar4

  use aux
  
  implicit none

  integer, parameter ::  numTimeSteps        = 300
  real, dimension(1:numTimeSteps) ::  time,mul_nucleon,mul_delta,mul_1535,mul_1520,mul_otherRes,&
                                    & mul_pion,mul_eta,mul_rho  
  real :: mul_nucleon_inp,mul_delta_inp,mul_1535_inp,mul_1520_inp,mul_otherRes_inp,&
        & mul_pion_inp,mul_eta_inp,mul_rho_inp
  integer :: idir,ndir,i,ndir1
  logical :: exist_file
  character(100) :: Directory,file
  
  ndir=70
  
  mul_nucleon=0.
  mul_delta=0.
  mul_1535=0.
  mul_1520=0.
  mul_otherRes=0.
  mul_pion=0.
  mul_eta=0.
  mul_rho=0.
  ndir1=0
  
  file='evol.dat'
  
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
        do i=1,33
           read(1,*)
        end do
        
        do i=1,numTimeSteps
           read(1,5) time(i), mul_nucleon_inp,mul_delta_inp,mul_1535_inp,mul_1520_inp,mul_otherRes_inp,&
                   & mul_pion_inp,mul_eta_inp,mul_rho_inp
5          format(1P,50(1x,e13.6))
           mul_nucleon(i)=mul_nucleon(i)+mul_nucleon_inp
           mul_delta(i)=mul_delta(i)+mul_delta_inp
           mul_1535(i)=mul_1535(i)+mul_1535_inp
           mul_1520(i)=mul_1520(i)+mul_1520_inp
           mul_otherRes(i)=mul_otherRes(i)+mul_otherRes_inp
           mul_pion(i)=mul_pion(i)+mul_pion_inp
           mul_eta(i)=mul_eta(i)+mul_eta_inp
           mul_rho(i)=mul_rho(i)+mul_rho_inp
        end do

        close(1)

  end do
     
  mul_nucleon=mul_nucleon/float(ndir1)
  mul_delta=mul_delta/float(ndir1)
  mul_1535=mul_1535/float(ndir1)
  mul_1520=mul_1520/float(ndir1)
  mul_otherRes=mul_otherRes/float(ndir1)
  mul_pion=mul_pion/float(ndir1)
  mul_eta=mul_eta/float(ndir1)
  mul_rho=mul_rho/float(ndir1)
  
  open(2,file='evolAver.dat',status='unknown')
  write(2,*)'# t, fm/c:   nuc:   delta:   N*(1535):   N*(1520):   otherNonstrRes:   pion:  eta:  rho:' 
  do i=1,numTimeSteps  
     write(2,5) time(i), mul_nucleon(i),mul_delta(i),mul_1535(i),mul_1520(i),mul_otherRes(i),&
                   & mul_pion(i),mul_eta(i),mul_rho(i)
  end do
  close(2)

end program evolAver
