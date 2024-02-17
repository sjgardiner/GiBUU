program main
use IDtable
use particleDefinition
use particleProperties
use energycalc
implicit none

type(particle) :: teilchen
real, dimension(1:3) :: betaToLRF
integer :: i
real :: vorher

call init_Database

betaToLRF=0.

teilchen%ID=pion
teilchen%mass=meson(pion)%mass
teilchen%pert=.false.

!!$Do i=0,200
!!$teilchen%mom(1:3)=(/0.,0.,float(i)*0.01/)
!!$teilchen%mom(0)=freeEnergy(teilchen)
!!$vorher=teilchen%mom(0)
!!$Print * , 'vorher', vorher
!!$call energyDetermination(teilchen, betaToLRF)
!!$Print * , 'nachher',teilchen%mom
!!$write(10,*) teilchen%mom(3),freeEnergy(teilchen), teilchen%mom(0)
!!$End do


betaToLRF=(/0.98,0.,0./)


teilchen%charge=-1
teilchen%ID=101
teilchen%mass=meson(pion)%mass
teilchen%pert=.true.
teilchen%mom(0:3)=(/ 112270.098891883 ,       79430.0891287182,   22526.8094651811   ,   76078.7677444666 /)
vorher=teilchen%mom(0)
Print * , 'vorher', vorher
call energyDetermination(teilchen, betaToLRF)
Print * , 'nachher',teilchen%mom
write(10,*) teilchen%mom(3),freeEnergy(teilchen), teilchen%mom(0)
teilchen%charge=0
teilchen%ID=101
teilchen%mass=meson(pion)%mass
teilchen%pert=.true.
teilchen%mom(0:3)=(/ 112270.098891883 ,       79430.0891287182,   22526.8094651811   ,   76078.7677444666 /)
vorher=teilchen%mom(0)
Print * , 'vorher', vorher
call energyDetermination(teilchen, betaToLRF)
Print * , 'nachher',teilchen%mom
write(10,*) teilchen%mom(3),freeEnergy(teilchen), teilchen%mom(0)
teilchen%charge=1
teilchen%ID=101
teilchen%mass=meson(pion)%mass
teilchen%pert=.true.
teilchen%mom(0:3)=(/ 112270.098891883 ,       79430.0891287182,   22526.8094651811   ,   76078.7677444666 /)
vorher=teilchen%mom(0)
Print * , 'vorher', vorher
call energyDetermination(teilchen, betaToLRF)
Print * , 'nachher',teilchen%mom
write(10,*) teilchen%mom(3),freeEnergy(teilchen), teilchen%mom(0)

end program main


