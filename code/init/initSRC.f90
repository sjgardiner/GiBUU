!******************************************************************************
!****m* /initSRC
! NAME
! module initSRC
!
! PURPOSE
! Contains a routine which samples a short-range-correlated pn pair
! within the target nucleus
!******************************************************************************
module initSRC

  implicit none
  private

  integer :: ibin
  real, dimension(:), allocatable :: dNSRC
  real, parameter :: dr3=0.2
  integer, parameter :: nr=1000
 
  integer, parameter :: imode=2   ! 1 -- local choice of positions (r1=r2), 2 -- nonlocal choice
  logical, parameter :: flagProb=.false.  ! turn-on the calculation of SRC probability 
  
  logical, parameter :: flagHist=.false.

  integer, save :: Ntrials=0, Nn
  real, save :: Prob=0.
  
  public :: DoInit_pnSRC, writeHistoSRC, writeProbSRC
  
contains

  !****s* initSRC/DoInit_pnSRC
  ! NAME
  ! subroutine DoInit_pnSRC(r1,r2)  
  !
  ! PURPOSE
  ! Samples positions of the proton and neutron of the SRC pair
  !
  ! OUTPUT
  ! * real, dimension(3), intent(out) :: r1, r2    ! p and n position vectors (fm)
  !****************************************************************************
  subroutine DoInit_pnSRC(r1,r2)
    use nucleusDefinition
    use dichtedefinition
    use nucleus, only: getTarget
    use densityStatic, only: staticDensity
    use random, only: rn
    real, dimension(3), intent(out) :: r1, r2    ! p and n position vectors (fm)

    type(tNucleus), save, pointer :: nuc
    real, save :: maxDist,maxP,maxN,maxDist2
    logical, save :: flagIni=.true.
    
    type(dichte) :: density
    integer :: i
    real :: rsq,r12,Pacc
    
    if(flagIni) then
       if(flagHist) then
          allocate(dNSRC(nr))
          dNSRC=0.
       end if
       nuc => getTarget()
       maxDist = nuc%MaxDist
       maxP = nuc%MaxDens(1)
       maxN = nuc%MaxDens(2)
       maxDist2=maxDist**2       
       Nn=nuc%mass-nuc%charge
       flagIni=.false.
    end if
        
    select case (imode)

    case(1)  ! local choice
    
       do
          do 
               r1(1)=(1.-2.*rn())*maxDist
               r1(2)=(1.-2.*rn())*maxDist
               r1(3)=(1.-2.*rn())*maxDist
               rsq=r1(1)**2+r1(2)**2+r1(3)**2
               if (rsq .le. maxDist2) exit
          end do
         
          density=staticDensity(r1,nuc)

          if(rn()*maxP*maxN.le.density%proton(0)*density%neutron(0)) exit

       end do
    
       r2=r1
     
       if(flagHist) then
          ibin=int(rsq**1.5/dr3)+1
          if(ibin.ge.1.and.ibin.le.nr) dNSRC(ibin)=dNSRC(ibin)+1.
       end if

    case(2)  ! nonlocal choice
  
       i=0     
       do
          do
             do 
                  r1(1)=(1.-2.*rn())*maxDist
                  r1(2)=(1.-2.*rn())*maxDist
                  r1(3)=(1.-2.*rn())*maxDist
                  rsq=r1(1)**2+r1(2)**2+r1(3)**2
                  if (rsq .le. maxDist2) exit
             end do   
         
             density=staticDensity(r1,nuc)

             if(rn()*maxP.le.density%proton(0)) exit

          end do

          do
             do 
                  r2(1)=(1.-2.*rn())*maxDist
                  r2(2)=(1.-2.*rn())*maxDist
                  r2(3)=(1.-2.*rn())*maxDist
                  rsq=r2(1)**2+r2(2)**2+r2(3)**2
                  if (rsq .le. maxDist2) exit
             end do   
         
             density=staticDensity(r2,nuc)

             if(rn()*maxN.le.density%neutron(0)) exit

          end do

          r12=sqrt((r1(1)-r2(1))**2+(r1(2)-r2(2))**2+(r1(3)-r2(3))**2)

          Pacc=CorFun(r12)

          if(.not.flagProb) then
             if(rn()*1.1.le.Pacc) exit
          else
             Prob=Prob+Pacc
             Ntrials=Ntrials+1
             i=i+1
             if(i.eq.100) exit
          end if
          
       end do

    end select
    
  end subroutine DoInit_pnSRC


  subroutine writeHistoSRC

    if(flagHist) then
       open(1,file='dNSRCdr3.dat',status='unknown')
       write(*,*)'# r:     dNSRC/dr3:'
       do ibin=1,nr
          write(1,*) (dr3*(ibin-0.5))**0.333333,dNSRC(ibin)/dr3
       end do
       write(1,*)'# Norma : ', sum(dNSRC(:))
       close(1)
    end if

  end subroutine writeHistoSRC

  
  subroutine writeProbSRC

    if(flagProb) then
       write(*,*)' Nn, Ntrials : ', Nn, Ntrials
       write(*,*)' Probability of a proton to belong to SRC : ', float(Nn)*Prob/float(Ntrials)
    end if

  end subroutine writeProbSRC
  
    
  real function CorFun(r)

    real, intent(in) :: r  ! distance between proton and neutron, fm
    
    ! Sharp cutoff:
    !real, parameter :: r0=0.92 ! fm  - Pb208, fitted to P_pn=21 %
    real, parameter :: r0=0.97 ! fm  - C12, fitted to P_pn=18 %  
    if(r.le.r0) then
       CorFun=1.
    else
       CorFun=0.
    end if
    
  end function CorFun

 
end module initSRC

