!******************************************************************************
!****p* /readParticleVectorBinary
! NAME
! program readParticleVectorBinary
! PURPOSE
! Illustrate how to read binary data files produced with the subroutine
! 'WriteParticleVector', like e.g. "RealParticles_Init.bin",
! by a selfcontained, GiBUU independend code.
!
! This program is a stubb, you have to fill it with your own analysis
!******************************************************************************
program readParticleVectorBinary

  implicit none

  character(*), parameter :: fName = "RealParticles_Init.bin"

  integer :: ios

  integer :: nEns,nPart,iEns,iPart,number,ID,charge,firstevent,event(2),history
  logical :: anti,pert,in_Formation
  real :: mass, mom(0:3), pos(1:3), vel(1:3), scaleXS, &
       lastCollisionTime, productionTime,formationTime, perWeight, &
       offshellParameter

  write(*,*) "Reading the file '",fName,"': ..."
  open(131,file=fName,status='old',form='UNFORMATTED',iostat=ios)
  if (ios/=0) then
     write(*,*) "Error while opening the file."
     stop
  endif

  read(131) nEns

  do iEns=1,nEns
     read(131) nPart
     do iPart=1,nPart

        read(131,iostat=ios) &
             pos,mom,vel,mass,lastCollisionTime,productionTime,formationTime,&
             perWeight,scaleXS,offshellParameter,ID,number,charge,&
             event,firstEvent,history,anti,pert,in_Formation

        if (ios > 0) cycle
        if (ios < 0) exit


        ! here you can do your own stuff...

!        write(*,*) iEns,iPart,number,ID,charge,anti,mass,mom,pos,scaleXS,&
!             firstevent,event

     end do
  end do

  close(131)



end program readParticleVectorBinary
