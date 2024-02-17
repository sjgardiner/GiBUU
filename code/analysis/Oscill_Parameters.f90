!******************************************************************************
!****m* /Oscill_Parameters
! NAME
! module Oscill_Parameters
!
! PURPOSE
! * This module contains the neutrino mixing angles and squared mass differences
! * Parameters taken from PDG (2017)
! * Also subroutine Posc:
! * Approx. Oscillation Formula from M. Freund, Phys.Rev. D64 (2001) 053003,
! * as used by Lindner et al, (2017), arXiv:1709.10252 
!*****************************************************************************
module Oscill_Parameters
   PUBLIC
    
    !! deltaM2_21/deltaM2_31 << 1   used by Freund
      
    ! all values for normal mass hierarchy

    real, parameter :: sin2_theta12 = 0.307,         & 
                     & sin2_theta23 = 0.51,          & ! sin^2(theta23)
                     & sin2_theta13 = 0.0210                               
                                
                                          
    real, parameter :: deltaM2_21 = 7.53e-5,         &
                     & deltaM2_32 = 2.45e-3                               
 !  deltaM2 given in eV^2
    
 !  matter density parameter, combined with Fermi-constant: GFNE = GF * Ne
    real, parameter :: GFNE = 7.01e-14   ! eV, value from Lindner et al.                     
                     
contains


  !****************************************************************************
  !****s* neutrinoAnalysis/Posc
  ! NAME
  ! subroutine Posc(pid,L,Enu,sin_deltaCP,Posc_mue,Posc_mutau,Posc_mumu)
  ! PURPOSE
  ! Calculate Oscillation Probability
  ! input: 
  ! pid (=+1 for neutrino, =-1 for antineutrino)
  ! L : baseline length in km, Enu: neutrino energy in GeV, 
  ! sin_DeltaCP = sin(delta_CP)
  ! output:
  ! Posc_mue gives probability for numu -> nue
  ! Posc_mutau gives probability for numu -> nutau 
  ! Posc_mumu gives probability for mu survival numu -> numu
  !****************************************************************************
 subroutine Posc(pid,L,Enu,sin_deltaCP,Posc_mue,Posc_mutau,Posc_mumu)  
    
    implicit none
     
    real, intent(in) :: pid,L,Enu,sin_deltaCP
    real, intent(out):: Posc_mue 
    real, intent(out), optional :: Posc_mutau,Posc_mumu
    
    real :: deltaM2_31, sin_2theta12, sin_2theta13,sin_2theta23,     &
       & s13,c13,s23,c23,s12,c12
       
    real ::  JCP,Delta,alpha,JCP_cotdeltaCP,A,cos_deltaCP  
     
    !! L in km and E in GeV
    !! numerical factor 5.068 in oscillation expressions from conversion of
    !! units in deltaM2 * L/E
            
    deltaM2_31 = deltaM2_32 + deltaM2_21
     
    ! sin_2theta = sin(2 * theta)
    
    sin_2theta12 = 2 * sqrt(sin2_theta12)*sqrt(1 - sin2_theta12)
    sin_2theta13 = 2 * sqrt(sin2_theta13)*sqrt(1 - sin2_theta13)
    sin_2theta23 = 2 * sqrt(sin2_theta23)*sqrt(1 - sin2_theta23)                                                            
    
    s13 = sqrt(sin2_theta13)
    c13  = sqrt(1 - sin2_theta13)                               
    s23  = sqrt(sin2_theta23)                               
    c23  = sqrt(1 - sin2_theta23) 
    s12 = sqrt(sin2_theta12)
    c12 = sqrt(1 - sin2_theta12)    
    
!*********************************************************************************  
   
    JCP = pid * 1./8. * sin_deltaCP * sin_2theta12*sin_2theta13*sin_2theta23*c13
    Delta = deltaM2_31 * L/(4.*Enu) * 5.068 ! 5.068 factor from dimension conversion
    alpha = deltaM2_21/deltaM2_31     
    
    ! pid = -1 for antineutrinos, +1 for neutrinos
    
    A = pid * 2. * sqrt(2.) * GFNE * Enu/deltaM2_32 * 1.e9         ! matter effects 
    
    cos_deltaCP = sqrt(1. - sin_deltaCP**2)        
    
    JCP_cotdeltaCP = 1./8. * cos_deltaCP * sin_2theta12*sin_2theta13*sin_2theta23*c13  
    
    Posc_mue = 4. * s13**2 * c13**2 * s23**2 * (sin((1. - A)*Delta)/(1. - A))**2  &
             & - 8 * alpha * JCP                                                  &
             & * sin(Delta) * sin(A*Delta)/A * sin((1. - A) * Delta)/(1. - A)     &
             & + 8 * alpha * (JCP_cotdeltaCP)                                     &
             & * cos(Delta) * sin(A*Delta)/A                                      &
             & * sin((1. - A)*Delta)/(1. - A)                                     &
             & + 4 * alpha**2 * s12**2 * c12**2 * c23**2 * (sin(A*Delta)/A)**2
 ! For antiparticles this gives probability for antinumu -> antinue 
   
    if(present(Posc_mutau))  then         
       Posc_mutau = sin_2theta23**2 * c13**4 * (sin(deltaM2_32 * L/(4*Enu)*5.068))**2
       Posc_mumu  = 1. - Posc_mutau  - Posc_mue    !muon survival probability
    end if            
         
  end subroutine Posc                                                      
    
end module Oscill_Parameters    