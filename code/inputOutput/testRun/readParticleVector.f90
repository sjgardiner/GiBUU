!******************************************************************************
!****p* /readParticleVector
! NAME
! program readParticleVector
! PURPOSE
! Illustrate how to read data files produced with the subroutine
! 'WriteParticleVector', like e.g. "RealParticles_Init_ALL.dat",
! by a selfcontained, GiBUU independend code.
!
! This program is a stubb, you have to fill it with your own analysis
!******************************************************************************
program readParticleVector

  implicit none

  character(*), parameter :: fName = "RealParticles_Init_ALL.dat"

  integer :: ios
  character(10) :: B


!  character*(*), parameter :: FormatMom = '4e11.3'
  character*(*), parameter :: FormatMom = '4e15.7'

!!$  character(250), parameter :: FormatPart2 = &
!!$       & '(i6,"|",i6,"[",i9,"]:",i3,i3,l2," ",f7.3,1P," |",'//FormatMom// &
!!$       & ',"|",3e11.3,"| ",0P,f5.3," ",i9,"[",2i9,"] ",'&
!!$       & //'1P,e11.3," time:",0P,3f9.3,l2," parents:",3i4)'

!!$  character(250), parameter :: FormatPart2 = &
!!$       & '(i6,A1,i6,A1,i9,A2,i3,i3,l2,A1,f7.3,A2,'//FormatMom//'&
!!$       & ,A1,3e11.3,A1,i9,A1,2i9)'

  character(250), parameter :: FormatPart2 = &
       & '(i6,A1,i6,A1,i9,A2,i3,i3,l2,A1,f7.3,A2,'//FormatMom//'&
       & ,A1,3e11.3,A1,0P,f6.3,i9,A1)'




!!$" |",'//FormatMom// &
!!$       & ',"|",3e11.3,"| ",0P,f5.3," ",i9,"[",2i9,"] ",'&
!!$       & //'1P,e11.3," time:",0P,3f9.3,l2," parents:",3i4)'




  integer :: iEns,iPart,number,ID,Q,firstevent,event(2)
  logical :: anti
  real :: mass, mom(0:3), pos(1:3), scaleXS

  write(*,*) "Reading the file '",fName,"': ..."
  open(131,file=fName,status='old',iostat=ios)
  if (ios/=0) then
     write(*,*) "Error while opening the file."
     stop
  endif

  do
     read(131,FormatPart2,iostat=ios) &
          iEns,B,iPart,B,number,B,ID,Q,anti,B,mass,B,mom,B,pos,B,scaleXS,firstevent,B,event(1)


!!$               iEns,B,iPart,B,number,B,ID,Q,anti,B,mass,B,mom,B,pos,B,scaleXS,&
!!$          B,firstEvent,B,event
     if (ios > 0) cycle
     if (ios < 0) exit

     write(*,*) iEns,iPart,number,ID,Q,anti,mass,mom,pos,scaleXS,firstevent,&
          event,">",B,"<"

  end do

  close(131)

end program readParticleVector
