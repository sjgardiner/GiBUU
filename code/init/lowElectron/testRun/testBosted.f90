!******************************************************************************
! program to print out the Bosted
!
! run it with:
! ./testBosted.x < jobBosted
!******************************************************************************
program testBosted
  use ParamEP
  use inputGeneral
  use particleProperties
  use output
  use eventprops

  implicit none
  
  ! Kinematics:  
  real :: W,W2,Q2,nu,eps,Gamma,xbj
  
  ! Structure functions:
  real :: F1n,F2n,F1p,F2p,FLp,FLp_res,FLp_nr,FLn,FLn_res,FLn_nr
  real :: F1p_res,F1p_nr,F2p_res,F2p_nr,F1n_res,F1n_nr,F2n_res,F2n_nr
  real :: F1d,F1d_res,F1d_nr,F2d,F2d_res,F2d_nr,F1_res,F1_nr,F2_res,F2_nr
  real, dimension(2,7) :: F1pres,F1nres,F1dres,sigresp,sigresn,sigresd 
  
  ! sigmaL/sigmaT ratios:
  real :: rp,rp_res,rp_nr,rn,rn_res,rn_nr,rd
  
  ! Cross sections:
  real :: XSp,XSpT,XSpL,XSpT_res,XSpT_nr,XSpL_res,XSpL_nr,XSpTL_nr,XSpTL_res
  real :: XSp_res,Xsp_nr,Xsn_res,Xsn_nr,XSpTres,XSpLres,XSnTres,XSnTL_nr,XSnTL_res
  real :: XSn,XSnT,XSnL,XSnT_res,XSnT_nr,XSnL_res,XSnL_nr,XSpTL,XSnTL,XS,XSd,XSd_res,XSd_nr
  real, dimension(7):: XSnTLrescontrib,XSpTLrescontrib,sigmarescontrib
  real :: sigma,sigma_nr,sigma_res  
  
  real :: Ein,theta,thetadeg
  real, save :: pi = 3.14159265
  real, save :: M = 0.938
  
  integer :: i,ios,k   
  
  
  NAMELIST /epkin/ Ein,thetadeg

  call readinputGeneral
  call initParticleProperties
  
  rewind(5,IOSTAT=ios)
  read(5,nml=epkin,IOSTAT=ios)
  call Write_ReadingInput("epkin",0,ios) 
!  Ein = 3.245
!  thetadeg = 26.98
  
  theta = thetadeg * pi/180. 
  write(118,"('#',1X,'Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(119,"('#',1X,'Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(120,"('#',1X,'Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(121,"('#',1X,'Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(122,"('#',1X,'Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(123,"('#',1X,'proton  ','Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(124,"('#',1X,'neutron  ','Ein =',f6.3,3X,'theta(deg) =',f6.2,3X,'theta(rad) =',f6.4,/)") Ein, thetadeg, theta
  write(118,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',10X,'XS')")
  write(119,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',10X,'XSp',10X,'XSp_res',8X,'XSp_nr',8X,'XSn',10X,'XSn_res',8X,'XSn_nr',8X,'XSd',8X,'XSd_res',8X,'XSd_nr')") 
  write(120,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',12X,'XSp',12X,'XSpTL',10X,'XSn',10X,'XSnTL',10X, &
		 & 'XSpT',10X,'XSpT_res',10X,'XSpT_nr',8X,'XSpL',9X,'XSpL_res',8X,'XSpL_nr',8X,   &
		 & 'XSnT',10X,'XSnT_res',8X,'XSnT_nr',8X,'XSnL',10X,'XSnL_res',8X,'XSnL_nr')")	
		 
  write(121,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',12X,'F1n',10X,'F1n_res',10X,'F1n_nr',10X,'F2n',10X,'F2n_res',8X,'F2n_nr',10X,'FLn',10X,'FLn_res',8X,'FLn_nr')")

  write(122,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',12X,'Rp',12X,'Rp_res',10X,'Rp_nr',10X,'Rn',12X,'Rn_res',8X,'Rn_nr')")

  write(123,"('#',4X,'1: W',5X,'2: nu',5X,'3: Q2',5X,'4: xbj',5X,'5: eps',8X,'6: Gamma',6X,'7: XSpTL',4X,'8: XSpTL_res',3X,'9: XSpTL_nr',3X,'10: XSpTLrescontrib')")

  write(124,"('#',4X,'1: W',5X,'2: nu',5X,'3: Q2',5X,'4: xbj',5X,'5: eps',8X,'6: Gamma',6X,'7: XSnTL',4X,'8: XSnTL_res',3X,'9: XSnTL_nr',3X,'10: XSnTLrescontrib')")

  write(125,"('#',4X,'W',9X,'Q2',8X,'nu',8X,'eps',8X,'Gamma',12X,'F1p',10X,'F1p_res',10X,'F1p_nr',10X,'F2p',10X,'F2p_res',8X,'F2p_nr',10X,'FLp',10X,'FLp_res',8X,'FLp_nr')")
  
  Q2=1
  eps=0.5
  
  call set_evprops(1,(/0.938,0.,0.,0./),(/0.,0.,0./))
  call CalcParamEP(1.5,Q2,eps, XS) ! dummy in order to read input
  
 
  
   
  do i=110,300
! do i=120,120
     W = i*0.01
     
     nu = Ein - Ef(W,Ein,theta)
     Q2 = Qsq(W,Ein,theta) 
     eps = (1 + 2*qvec2(W,Ein,theta)/Q2 * tan(theta/2.)**2)**(-1)
	 xbj = Q2/(2*M*nu)
     
     call CalcParamEP(W,Q2,eps, XS) 
	 !call eNBosted
	 !The returned cross section XS is
     !    XS = \sigma^* = \sigma_T+\epsilon\sigma_L
     !             = \frac{1}{\Gamma} \frac{d\sigma}{dE' d\Omega}
     
! Now multiply sum of longitudinal and transverse cross section with flux factor Gamma
! to get dd X-section      
     
     Gamma = (nu - Q2/(2*0.93822))/(2*pi**2*Q2) * (Ein - nu)/Ein &
           & * 1./(1 - eps) * 1./137.    ! Flux of virtual photons
     XS = XS * Gamma
                
! XS is in mubar/(sr * GeV)        
     write(118,'(4f10.4,2e14.4)') W,Q2,nu,eps, Gamma, XS

write(*,*)	  
write(*,*) 'first proton'
      call set_evprops(1,(/0.938,0.,0.,0./),(/0.,0.,0./))	  
	  w2 = W**2
	  call F1F2IN09(1.D0, 1.D0, q2, w2, F1pres,F1p, F2p,rp,F1p_res,F1p_nr,F2p_res,F2p_nr)
	  ! F1pres contains the contributions of individual resonances
	  
	  ! Now calculate dd X-sections for total, resonance and nonresonance contributions of 
	  ! structure functions.
	  
	  XSp = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2p/nu + 2.*tan(theta/2.)**2*F1p/M) * 0.3894e3  ! factor 0.3894 e3 to give X-sect in microbarn
	  XSp_res = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2p_res/nu + 2.*tan(theta/2.)**2*F1p_res/M) * 0.3894e3
      XSp_nr = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2p_nr/nu + 2.*tan(theta/2.)**2*F1p_nr/M) * 0.3894e3		
write(*,*) 'XSp=',XSp,'F1p=',F1p,'F1p_res=',F1p_res,'F1p_nr=',F1p_nr,'F2p=',F2p,'F2p_res=',F2p_res,'F2p_nr=',F2p_nr	 
    
	FLp =     (1+Q2/nu**2)*F2p - 2.*xbj*F1p        ! longitudinal structure function
	FLp_res = (1+Q2/nu**2)*F2p_res - 2.*xbj*F1p_res 
	FLp_nr =  (1+Q2/nu**2)*F2p_nr -  2.*xbj*F1p_nr
	
	XSpT =     4.*pi**2/Q2 * 1./137. * 2.*xbj/(1. - xbj) * F1p * 0.3894e3	    ! transverse X-section sigma_T
	XSpT_res = 4.*pi**2/Q2 * 1./137. * 2.*xbj/(1. - xbj) * F1p_res * 0.3894e3	! res contribution to sigma_T
	XSpT_nr =  4.*pi**2/Q2 * 1./137. * 2.*xbj/(1. - xbj) * F1p_nr * 0.3894e3	! nonres contribution to sigma_T
	
	XSpL =     4.*pi**2/Q2 * 1./137. /(1.-xbj) * FLp * 0.3894e3	        ! longitudinal X-section sigma_L	
	XSpL_res = 4.*pi**2/Q2 * 1./137. /(1.-xbj) * FLp_res * 0.3894e3	    ! same for ressigma_L
	XSpL_nr =  4.*pi**2/Q2 * 1./137. /(1.-xbj) * FLp_nr * 0.3894e3	
	
	rp_res = XSpL_res/XSpT_res             ! L/T ratio for resonance sigma
	rp_nr = XSpL_nr/XSpT_nr
	
	XSpTL =    Gamma * (XSpT + eps*XSpL)      ! measured dd X-section as sum of T and L 
	XSpTL_res =    Gamma * (XSpT_res + eps*XSpL_res)      ! res contrib to dd X-section as sum of T and L
	XSpTL_nr = Gamma * (XSpT_nr + eps*XSpL_nr)      ! background dd X-section as sum of T and L
	
	XSpTres = 0.0
	XSpLres = 0.0

    sigresp(1,:) = 2.*xbj/(1. - xbj) * 4.*pi**2/Q2 *1./137.* F1pres(1,:) * 0.3894e3
	sigresp(2,:) = sigresp(1,:) * rp_res
	! now check of resonance contributions by summing over res contributions:
	do k=1,7
       XSpTres = XSpTres + sigresp(1,k)
	   XSpLres = XSpLres + sigresp(2,k)
	enddo
	
	XSpTLrescontrib = Gamma * (sigresp(1,:) + eps*sigresp(2,:))  !contribution of individual resonances to dd X-sect
	
write(*,*) '153 proton   ','XSpTres=',XSpTres,'XSpT_res=',XSpT_res,'XSpLres=',XSpLres,'XSpL_res=',XSpL_res
	
	
    write(123,'(4f10.4,17e14.4)') W,nu,Q2,xbj,eps,Gamma,XSpTL,XSpTL_res,XSpTL_nr,XSpTLrescontrib
	!write (123,*)'proton   ','Q2=',q2,'W2=',W2,'XSpTL=',XSpTL,'XSpTL_res=',XSpTL_res,'XSpTL_nr=',XSpTL_nr,'XSpT_res=',XSpT_res,'XspTres=',XSpTres,'XSpTLres=',XSpTLres


    call eNucleon(+1,Ein,theta,nu,sigmarescontrib,sigma_res,sigma_nr,sigma,F1_res,F1_nr,F2_res,F2_nr)		
	
write (*,*) 'proton..164..', 'Ein=',Ein,'theta=',theta,'Q2=',Q2,'nu=',nu,'sigma=',sigma,'sigma_res=',sigma_res,'sigma_nr=',sigma_nr, 'sigmares=',sigmarescontrib
write(*,*)	
write(*,*) 'now neutron'
	  call set_evprops(0,(/0.938,0.,0.,0./),(/0.,0.,0./))
  
	  call F1F2IN09(0.D0, 1.D0, q2, w2, F1nres,F1n, F2n,rn,F1n_res,F1n_nr,F2n_res,F2n_nr)
write(*,*) 'F1nres(1,1)=',F1nres(1,1)
	  XSn = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2n/nu + 2.*tan(theta/2.)**2*F1n/M) * 0.3894e3
	  XSn_res = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2n_res/nu + 2.*tan(theta/2.)**2*F1n_res/M) * 0.3894e3
      XSn_nr = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2n_nr/nu + 2.*tan(theta/2.)**2*F1n_nr/M) * 0.3894e3		
write(*,*) 'Xsn=',Xsn,'F1n=',F1n,'F1n_res=',F1n_res,'F1n_nr=',F1n_nr,'F2n=',F2n,'F2n_res=',F2n_res,'F2n_nr=',F2n_nr
	
	FLn =     (1+Q2/nu**2)*F2n - 2.*xbj*F1n 
	FLn_res = (1+Q2/nu**2)*F2n_res - 2.*xbj*F1n_res
	FLn_nr =  (1+Q2/nu**2)*F2n_nr - 2.*xbj*F1n_nr
	
	XSnT =     4.*pi**2/Q2 * 1./137. * 2.*xbj/(1. - xbj) *  F1n * 0.3894e3	
	XSnT_res = 4.*pi**2/Q2 * 1./137. * 2.*xbj/(1. - xbj) * F1n_res * 0.3894e3	
	XSnT_nr =  4.*pi**2/Q2 * 1./137. * 2.*xbj/(1. - xbj) * F1n_nr * 0.3894e3	
	
	XSnL =     4.*pi**2/Q2 * 1./137. /(1.-xbj) * FLn * 0.3894e3	
	XSnL_res = 4.*pi**2/Q2 * 1./137. /(1.-xbj) * FLn_res * 0.3894e3	
	XSnL_nr =  4.*pi**2/Q2 * 1./137. /(1.-xbj) * FLn_nr * 0.3894e3	
  
    XSnTL = Gamma * (XSnT + eps*XSnL)
	XSnTL_nr = Gamma * (XSnT_nr + eps*XSnL_nr)      ! background dd X-section as sum of T and L
	XSnTL_res =    Gamma * (XSnT_res + eps*XSnL_res)      ! dd X-section as sum of T and L
	
	rn_res = XSnL_res/XSnT_res
	rn_nr = XSnL_nr/XSnT_nr
	
	XSnTres =  0.0

    sigresn(1,:) = 2.*xbj/(1. - xbj) * 4.*pi**2/Q2 *1./137.* F1nres(1,:) * 0.3894e3 ! T
	sigresn(2,:) = sigresn(1,:) * rn                                                ! L
	! now check of resonance contributions by summing over res contributions:
	do k=1,7
       XSnTres = XSnTres + sigresn(1,k)
	enddo
	XSnTLrescontrib = Gamma * (sigresn(1,:) + eps*sigresn(2,:))     !contribution of individual resonances to dd X-sect
	
	write (124,'(4f10.4,17e14.4)') W,nu,Q2,xbj,eps,Gamma,XSnTL,XSnTL_res,XSnTL_nr,XSnTLrescontrib
	!write (124,*)'neutron   ','Q2=',q2,'W2=',W2,'XSnTL=',XSnTL,'XSnTL_res=',XSnTL_res,'XSnTL_nr=',XSnTL_nr,'XSnT_res=',XSnT_res,'XSnTres=',XSnTres,'XSnTLres=',XSnTLres
      
	
write(*,*)
write(*,*) 'now deuteron'	  
	  call F1F2IN09(1.D0, 1.D0, q2, w2, F1dres,F1d, F2d,rd,F1d_res,F1d_nr,F2d_res,F2d_nr)
	  XSd = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2d/nu + 2.*tan(theta/2.)**2*F1d/M) * 0.3894e3
	  XSd_res = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2d_res/nu + 2.*tan(theta/2.)**2*F1d_res/M) * 0.3894e3
      XSd_nr = 4.* (1./137.)**2 * (Ein - nu)**2/Q2**2 * cos(theta/2.)**2 &
	      & * (F2d_nr/nu + 2.*tan(theta/2.)**2*F1d_nr/M) * 0.3894e3			  
   
write(*,*) 'F1d=',F1d,'F1d_res=',F1d_res,'Fd_nr=',F1d_nr,'F2d=',F2d,'F2d_res=',F2d_res,'F2d_nr=',F2d_nr
	  
	  
	   write(119,'(4f10.4,10e14.4)') W,Q2,nu,eps, Gamma, XSp,XSp_res,XSp_nr,XSn,XSn_res,XSn_nr,XSd,XSd_res,XSd_nr
	   
	   write(120,'(4f10.4,17e15.4)') W,Q2,nu,eps,Gamma,XSp,XSpTL,XSn,XSnTL,XSpT,XSpT_res,XSpT_nr,XspL,XSpL_res,XSpL_nr,XSnT,XSnT_res,XsnT_nr,XSnL,XSnL_res,XSnL_nr
	   
	   write(121,'(4f10.4,10e15.4)') W,Q2,nu,eps,Gamma,F1n,F1n_res,F1n_nr,F2n,F2n_res,F2n_nr,FLn,FLn_res,FLn_nr
	   
	   write(125,'(4f10.4,10e15.4)') W,Q2,nu,eps,Gamma,F1p,F1p_res,F1p_nr,F2p,F2p_res,F2p_nr,FLp,FLp_res,FLp_nr
	   
	   write(122, '(4f10.4,7e15.4)') W,Q2,nu,eps,Gamma,rp,rp_res,rp_nr,rn,rn_res,rn_nr
  end do
  
CONTAINS

REAL FUNCTION Qsq(W,Ein,theta)
Real W, Ein,theta

Qsq = 4 * Ein *Ef(W,Ein,theta) * sin(theta/2)**2

END FUNCTION Qsq

  
end program testBosted
