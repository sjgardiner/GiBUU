!******************************************************************************
!****m* /AZN
! NAME
! module AZN
!
! PURPOSE
! provides mass number A, charge Z and neutrino number N of target nucleus
!
!******************************************************************************

module AZN
implicit none 
real  :: Atarget,Ntarget,Ztarget

contains

subroutine AZNsub(Atarget,Ntarget,Ztarget)
    use nucleus, only: getTarget
    use nucleusdefinition
    
    type(tnucleus), pointer :: targetNuc
    
    real :: Atarget,Ntarget,Ztarget
    
    targetNuc => getTarget()
    Atarget = targetNuc%mass
    Ztarget = targetNuc%charge
    Ntarget = Atarget - Ztarget 
end subroutine AZNsub    

end module AZN