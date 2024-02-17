! This program calculates the total cross section as function of W and Q2.
! It also does a splitting into different channels:
! * W<2GeV: ...
! * W>2GeV: ...
! The program does a Rosenbluth separation, i.e. it calculates the cross
! section for two distinct values of epsilon and prints it for the
! extrapolation eps=0 (==sigma_T) and eps=1 (==sigma_T+sigma_L).
!
! This does the same as code/init/ElectronGenerator/testRun/plotAllXS.f90,
! but here with the neutrino init.


program plotAllXS

  use constants
  use inputGeneral
  use particleProperties, only: initParticleProperties
  use particleDefinition
  use Coll_gammaN
  use PythiaSpecFunc
  use hadronFormation, only : forceInitFormation
  use Coll_Pythia
  use CollTools

  use eN_eventDefinition
  use eN_event
!  use eventGenerator_eN_lowEnergy
  !  use eventGenerator_eN_HiEnergy
  use neutrinoXsection
  use neutrinoSigma

  use ParamEP
  use output
  use minkowski, only: abs4

  use AnaEventDefinition
  use AnaEvent, only : event_init,event_clear,event_add
  use MultiplicityAnalysis

  use collisionTerm, only:  ForceDecays
  use insertion, only: GarbageCollection

  use hist

  implicit none

  integer :: iW,iQ2,iN,iN1,iN2, i,k

  integer :: iWmin,iWmax,iWdelta
  integer :: iQ2min,iQ2max,iQ2delta
  integer :: iNmax
  logical :: verbose

  real    :: W, Q2
  real :: nu,Wfree,eps,fT1,fT2
  Type(electronNucleon_event), save :: eNev1,eNev2, eNev
  type(particle), dimension(1,50) :: realPart
  type(particle),dimension(1,1:50),TARGET  :: finalState, finalStateBAK
  logical :: flagOK
  real :: s1,s2,x1,x2

  real, dimension(37) :: sig
  real :: XS, XS_sum1, XS_sum2, R1
  real,dimension(0:4) :: XS_Arr, XS_Arr_sum1, XS_Arr_sum2
  real,dimension(8) :: XS_Arr_low, XS_Arr_low_sum1, XS_Arr_low_sum2

  real, parameter :: eps1=0.05, eps2=0.99

  type(particle) :: TargetNuc
  integer :: chargeTarget

  logical :: doMultAna = .true.
!  logical :: doMultAna = .false.
  character*(20) :: Prefix_MultAna
  type(particle), POINTER :: pPart
  type(tAnaEvent)            :: tEv
  integer :: nPart0, nPion
  real, dimension(2,2,-1:2) :: nPionArr
  real, dimension(0:3) :: MomSum1, MomSum2

  type(histogram) :: hW1,hW2
  logical :: flag1,flag2
  real :: W1, W2

  NAMELIST /datatable/ iWmin,iWmax,iWdelta,iQ2min,iQ2max,iQ2delta,iNmax,&
       & verbose, chargeTarget


  call readinputGeneral
  call initParticleProperties
  call forceInitFormation
!  call InitParticleProperties

!  call ReadHiGammaNucleus
  call get_sigma_namelist(0)

  call SetSomeDefaults_PY
  call event_INIT(tEv)

  iWmin   = 110
  iWmax   = 199
  iWdelta =   2

  iQ2min   =   0
  iQ2max   = 200
  iQ2delta = 100

  iNmax = 10000

  verbose = .false.

  chargeTarget = 1

  call Write_ReadingInput('datatable',0)
  rewind(5)
  read(5,nml=datatable)
  write(*,'(A,3i7)') '  iW   = ',iWmin,iWmax,iWdelta
  write(*,'(A,3i7)') '  iQ2  = ',iQ2min,iQ2max,iQ2delta
  write(*,'(A,3i7)') '  iN   = ',1,iNmax,1
  write(*,*)
  write(*,*) '  verbose=',verbose
  write(*,*)
  write(*,*) '  charge of target: ',chargeTarget

  call Write_ReadingInput('datatable',1)
  write(*,*)
  write(*,*) 'eps:',eps1,eps2

  call setToDefault(TargetNuc)
  TargetNuc%ID = 1
  TargetNuc%charge = chargeTarget
  TargetNuc%mass = 0.938
  TargetNuc%mom = (/0.938, 0.0, 0.0, 0.0 /)
!  TargetNuc%mom = (/0.1, 0.0, 0.0, 0.0 /)
  TargetNuc%pos = 9999999d0

  realPart(1,1)%ID = 0

  call eNev_SetProcess(eNev1, 1,1)  ! set to EM and electron
  call eNev_SetProcess(eNev2, 1,1)  ! set to EM and electron

  call createHist(hW1, "hadron W eps1", 0.9, 3.5, 0.01)
  call createHist(hW2, "hadron W eps2", 0.9, 3.5, 0.01)

  do iQ2=iQ2min,iQ2max,iQ2delta
     Q2 = iQ2*0.01
     if (Q2<0.1) Q2=0.1

     do iW=iWmin,iWmax,iWdelta
        W = iW*0.01

        call eNev_init_eWQ(eNev1,eps1,W,Q2,flagOK)
        call eNev_init_Target(eNev1,TargetNuc,flagOK)

        call eNev_init_eWQ(eNev2,eps2,W,Q2,flagOK)
        call eNev_init_Target(eNev2,TargetNuc,flagOK)



!        call write_electronNucleon_event(eNev1,.FALSE.,.FALSE.)
        call eNeV_GetKinV(eNev1, nu,Q2,W,Wfree,eps,fT1)
        fT1 = fT1/ ( 1e-3* pi/(eNev1%lepton_out%mom(0)*eNev1%lepton_in%mom(0)))

        s1=abs4(eNev1%nucleon%mom+eNev1%lepton_in%mom)
        x1=eNeV_Get_LightX(eNev1)
!        write(*,*) 'nu :',nu
!        write(*,*) 'xB :',Q2/(2*0.938*nu)
!        write(*,*) 'sqrt(s) :',s1
!        write(*,*) 'x  :',x1

!        call write_electronNucleon_event(eNev2,.FALSE.,.FALSE.)
        call eNeV_GetKinV(eNev2, nu,Q2,W,Wfree,eps,fT2)
        fT2 = fT2/ ( 1e-3* pi/(eNev2%lepton_out%mom(0)*eNev2%lepton_in%mom(0)))
        s2=abs4(eNev2%nucleon%mom+eNev2%lepton_in%mom)
        x2=eNeV_Get_LightX(eNev2)
!        write(*,*) 'nu :',nu
!        write(*,*) 'xB :',Q2/(2*0.938*nu)
!        write(*,*) 'sqrt(s) :',s2
!        write(*,*) 'x  :',x2

        write(111,'(2f7.3,1P,12e13.4)') W,Q2, nu,Q2/(2*0.938*nu),&
             &s1,x1,(1-x1)*(s1**2-2*0.938**2),&
             &s2,x2,(1-x2)*(s2**2-2*0.938**2)

!        stop

        XS_sum1=0.0
        XS_sum2=0.0

        XS_Arr_sum1=0.0
        XS_Arr_sum2=0.0

        XS_Arr_low_sum1=0.0
        XS_Arr_low_sum2=0.0

        iN1=0
        iN2=0

        write(*,'(2f7.3,1P,12e13.4)') W,Q2

        if (doMultAna) then
           write(Prefix_MultAna,'(2f8.4)') W,Q2
           call Multiplicity_Reset
           nPionArr = 0.0
        endif


        do iN=1,iNmax
           if (verbose) write(*,*) '=======iN =',iN

           !===============================================================
           ! eps1
           !===============================================================

           sig = 0.
           do k=2,37 ! no QE!
              eNev = eNev1
              call XsecdCosthetadElepton(eNev,k,finalState(1,:),sig(k))

              if ((doMultAna).and.(k==34)) finalStateBAK = finalState
           end do

           sig = sig/fT1
           sig = sig/(2*pi) ! ???

           XS_Arr_low = 0.
           XS_Arr_low(2) = sum(sig(2:31))
           XS_Arr_low(3) = sum(sig(32:33))
           XS_Arr_low(4) = sig(37)
           XS_Arr_low(5) = sig(34)

           XS_Arr = 0.0

           XS = sum(sig)

           flagOK = (sum(abs(sig))>0.)

           if (flagOK) then
              iN1=iN1+1
              XS_sum1=XS_sum1+XS
              XS_Arr_sum1=XS_Arr_sum1+XS_Arr
              XS_Arr_low_sum1=XS_Arr_low_sum1+XS_Arr_low

              if ((doMultAna).and.(sig(34).gt.0.0)) then

!                 write(*,*) "Sig=",sig(34)

                 !...Force all particles to decay...
!                 call ForceDecays(finalstateBAK,realPart, 0.)
                 call GarbageCollection(finalstateBAK)

!                 call writeParticle(6,99,finalstateBAK(1,:))

                 finalstateBAK%perweight = 1.0
                 nPart0 = 0
                 call event_CLEAR(tEv)
                 MomSum1 = 0
                 do i=1,size(finalstateBAK,dim=2)
                    if (finalstateBAK(1,i)%ID.lt.0) exit
                    if (finalstateBAK(1,i)%ID.eq.0) cycle
                    pPart => finalstateBAK(1,i)
                    call event_ADD(tEv,pPart)
                    nPart0 = nPart0+1
                    MomSum1 = MomSum1 + pPart%mom
                 end do
                 call Multiplicity_AddEvent(tEv)

                 W1 = abs4(MomSum1,flag1)
                 if (flag1) call addHist(hW1, W1, sig(34))

                 nPion = sum(tEv%numberParticles(1,:))
!                 write(*,*) 'nPart0,nPion=',nPart0,tEv%numberParticles(1,:)

                 if (nPart0==1) then
                    nPionArr(1,1,-1) = nPionArr(1,1,-1) + sig(34)
                    if (nPion<3) nPionArr(1,1,nPion) = nPionArr(1,1,nPion) + sig(34)
                 else
                    nPionArr(1,2,-1) = nPionArr(1,2,-1) + sig(34)
                    if (nPion<3) nPionArr(1,2,nPion) = nPionArr(1,2,nPion) + sig(34)
                 end if

              end if


           end if

           !===============================================================
           ! eps2
           !===============================================================

           sig = 0.
           do k=2,37 ! no QE!
              eNev = eNev2
              call XsecdCosthetadElepton(eNev,k,finalState(1,:),sig(k))

              if ((doMultAna).and.(k==34)) finalStateBAK = finalState
           end do
           sig = sig/fT2
           sig = sig/(2*pi) ! ???

           XS_Arr_low = 0.
           XS_Arr_low(2) = sum(sig(2:31))
           XS_Arr_low(3) = sum(sig(32:33))
           XS_Arr_low(4) = sig(37)
           XS_Arr_low(5) = sig(34)

           XS_Arr = 0.0

           XS = sum(sig)

           flagOK = (sum(abs(sig))>0.)

           if (flagOK) then
              iN2=iN2+1
              XS_sum2=XS_sum2+XS
              XS_Arr_sum2=XS_Arr_sum2+XS_Arr
              XS_Arr_low_sum2=XS_Arr_low_sum2+XS_Arr_low

              if ((doMultAna).and.(sig(34).gt.0.0)) then

!                 write(*,*) "Sig=",sig(34)

                 !...Force all particles to decay...
!                 call ForceDecays(finalstateBAK,realPart, 0.)
                 call GarbageCollection(finalstateBAK)

!                 call writeParticle(6,99,finalstateBAK(1,:))

                 finalstateBAK%perweight = 1.0
                 nPart0 = 0
                 call event_CLEAR(tEv)
                 MomSum2 = 0
                 do i=1,size(finalstateBAK,dim=2)
                    if (finalstateBAK(1,i)%ID.lt.0) exit
                    if (finalstateBAK(1,i)%ID.eq.0) cycle
                    pPart => finalstateBAK(1,i)
                    call event_ADD(tEv,pPart)
                    nPart0 = nPart0+1
                    MomSum2 = MomSum2 + pPart%mom
                 end do
                 call Multiplicity_AddEvent(tEv)

                 W2 = abs4(MomSum2,flag2)
                 if (flag2) call addHist(hW2, W2, sig(34))

!                 write(*,*) "W,sum: ",W,W1,W2

                 nPion = sum(tEv%numberParticles(1,:))
!                 write(*,*) 'nPart0,nPion=',nPart0,tEv%numberParticles(1,:)

                 if (nPart0==1) then
                    nPionArr(2,1,-1) = nPionArr(2,1,-1) + sig(34)
                    if (nPion<3) nPionArr(2,1,nPion) = nPionArr(2,1,nPion) + sig(34)
                 else
                    nPionArr(2,2,-1) = nPionArr(2,2,-1) + sig(34)
                    if (nPion<3) nPionArr(2,2,nPion) = nPionArr(2,2,nPion) + sig(34)
                 end if

              end if
           end if

        end do

        if (iN1.gt.0) then
           XS_sum1=XS_sum1/iN1
           XS_Arr_sum1=XS_Arr_sum1/iN1
           XS_Arr_low_sum1=XS_Arr_low_sum1/iN1
           nPionArr(1,:,:) = nPionArr(1,:,:)/iN1
        end if

        if (iN2.gt.0) then
           XS_sum2=XS_sum2/iN2
           XS_Arr_sum2=XS_Arr_sum2/iN2
           XS_Arr_low_sum2=XS_Arr_low_sum2/iN2
           nPionArr(2,:,:) = nPionArr(2,:,:)/iN2
        end if

        write(121,'(2f7.3,1P,99e13.4)') W,Q2, XS_sum1, XS_sum2, &
             XS_Arr_sum1, XS_Arr_sum2,&
             XS_Arr_low_sum1, XS_Arr_low_sum2


        write(122,'(2f7.3,1P,99e13.4)') W,Q2, (XS_sum1-XS_sum2)/(eps1-eps2), &
             & (eps2*XS_sum1-eps1*XS_sum2)/(eps2-eps1), &
             & (XS_Arr_sum1-XS_Arr_sum2)/(eps1-eps2), &
             & (eps2*XS_Arr_sum1-eps1*XS_Arr_sum2)/(eps2-eps1), &
             & (XS_Arr_low_sum1-XS_Arr_low_sum2)/(eps1-eps2), &
             & (eps2*XS_Arr_low_sum1-eps1*XS_Arr_low_sum2)/(eps2-eps1)


        call flush(121)
        call flush(122)

        call CalcParamEP(W,Q2,eps1, XS_sum1)
        call CalcParamEP(W,Q2,eps2, XS_sum2)

        write(221,'(2f7.3,1P,12e13.4)') W,Q2, XS_sum1, XS_sum2
        write(222,'(2f7.3,1P,12e13.4)') W,Q2, (XS_sum1-XS_sum2)/(eps1-eps2),&
             & (eps2*XS_sum1-eps1*XS_sum2)/(eps2-eps1)

        call flush(221)
        call flush(222)

!!$        if (isHigh) then

        XS_sum1 = 0.
        XS_sum2 = 0.
        if (W>1.75) then
           call CalcParamEP_ALLM(W,Q2,XS_sum1)
           call CalcParamEP_ALLM97(W,Q2,XS_sum2)
        endif
        call CalcParamEP_R1990(W,Q2,R1)

        write(301,'(2f7.3,1P,12e13.4)') W,Q2, Q2/(2*0.938*nu), XS_sum1, XS_sum2, R1, Q2/nu**2
        call flush(301)

!        if (doMultAna) call Multiplicity_Write(Prefix_MultAna)

        if (doMultAna) then
           write(401,'(2f7.3,1P,12e13.4)') W,Q2,nPionArr(1,1,:),nPionArr(1,2,:)
           call flush(401)
           write(402,'(2f7.3,1P,12e13.4)') W,Q2,nPionArr(2,1,:),nPionArr(2,2,:)
           call flush(402)
        end if

        call writeHist(hW1, add=1e-20, mul=1./iN2, file="hadrW1.dat")
        call writeHist(hW2, add=1e-20, mul=1./iN2, file="hadrW2.dat")

        stop


     end do




!!$     write(111,*)
!!$     write(111,*)
!!$
!!$     write(121,*)
!!$     write(121,*)
!!$     write(122,*)
!!$     write(122,*)
!!$
!!$     write(221,*)
!!$     write(221,*)
!!$     write(222,*)
!!$     write(222,*)
!!$
!!$     write(301,*)
!!$     write(301,*)

  end do

  write(*,*) 'Done.'


end program plotAllXS
