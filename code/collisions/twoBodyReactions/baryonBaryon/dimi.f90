!******************************************************************************
!****m* /dimi
! NAME
! module dimi
! PURPOSE
! Includes the N N <-> Delta N cross sections according to the OPE model by
! Dmitriev et al.
!
! See:
! * Diploma Thesis Effenberger, chapter 3.3.1
! * Dmitriev/Sushkov/Gaarde, Nucl. Phys. A 459 (1986) 503-524
!******************************************************************************
module dimi
  use constants, only: mN, mPi
  implicit none
  private

  public :: dimiSigma
  public :: dimiIntegrated
  public :: massNN_NDelta
  public :: mnnnd2
  public :: dimidm
  public :: setMaxSqrts

  real, parameter :: dMass = 0.001

  real, parameter :: dSrts = 0.001
  real, parameter :: srtsMin = 2.*mN + mPi
  integer, save :: maxPoints_srts = 4000     ! this will be overwritten by 'setMaxSqrts', according to the high energy thresholds in master_2body

contains

  !****************************************************************************
  !****f* dimi/massNN_NDelta
  ! NAME
  ! function massNN_NDelta(srts) result (massDelta)
  ! PURPOSE
  ! Evaluates masses in N N -> N Delta process using dimidm.
  ! This is done by a Monte-Carlo process utilizing dsigma/dm.
  !****************************************************************************
  function massNN_NDelta(srts) result (massDelta)
    use constants, only: mN
    use random, only: rn
    use particleProperties, only: hadron
    use idTable, only: Delta

    real, intent(in) :: srts
    real :: massDelta

    integer :: iMass
    real :: sum_dsdm, x, dsdm, dsdm_Integrated, sigma

    iMass=0
    sum_dsdm=0.
    x=rn()
    do
       massDelta=hadron(Delta)%minMass+float(iMass)*dMass
       call dimidm(dsdm,dsdm_Integrated,sigma,srts,massDelta)
       sum_dsdm=sum_dsdm+dsdm*dMass
       if (sum_dsdm >= x*dsdm_Integrated) exit
       if (massDelta > 6.) then
          write(*,*) 'Problem in massNN_NDelta. massDelta > 6 GeV'
          write(*,*) 'mass of Delta=', massDelta,srts
          write(*,*) 'sqrt(s)=', srts
          write(*,*) sum_dsdm, dsdm_Integrated,x
          stop
       end if
       iMass=iMass+1
    end do

    if (iMass==0) then
       massDelta= hadron(Delta)%minMass &
            + rn()*(srts - hadron(Delta)%minMass - mN)
    else
       massDelta= massDelta + (rn()-0.5)*dMass
       if (massDelta > srts-mN) then
          massDelta = hadron(Delta)%minMass &
               + (float(iMass)-0.5)*dMass &
               + rn() * (srts - mN - hadron(Delta)%minMass - (float(iMass)-0.5)*dMass)
          if (massDelta > srts-mN) then
             write(*,*) 'problems in massNN_NDelta:',srts,massDelta,iMass
             stop
          end if
       end if
    end if

  end function massNN_NDelta


  !****************************************************************************
  !****f* dimi/dimiSigma
  ! NAME
  ! function dimiSigma(srts, massDelta) result (sigma)
  ! PURPOSE
  ! Returns the cross section for Delta++ n -> pp scattering.
  ! NOTES
  ! The return value is not a cross section, but has to be divided by p_in
  !****************************************************************************
  function dimiSigma(srts, massDelta) result (sigma)
    real, intent(in)  :: srts, massDelta
    real :: sigma
    real :: dummy1, dummy2
    call dimidm(dummy1, dummy2, sigma, srts, massDelta)
  end function dimiSigma


  !****************************************************************************
  !****f* dimi/dimiIntegrated
  ! NAME
  ! function dimiIntegrated(srts) result (dsdm_Integrated)
  ! PURPOSE
  ! Returns the cross section for p p -> Delta++ n scattering, integrated
  ! over the Delta mass.
  ! NOTES
  ! The return value is not a cross section, but has to be divided by p_in
  !****************************************************************************
  function dimiIntegrated(srts) result (dsdm_Integrated)
    real, intent(in) ::  srts
    real :: dsdm_Integrated
    real, parameter :: dummy1=1.2
    real :: dummy2, dummy3
    call dimidm(dummy2, dsdm_Integrated, dummy3, srts, dummy1)
  end function dimiIntegrated


  !****************************************************************************
  !****s*  dimi/setMaxSqrts
  ! NAME
  ! subroutine setMaxSqrts(srtsMax)
  ! PURPOSE
  ! set the upper bound of the tabulation
  !****************************************************************************
  subroutine setMaxSqrts(srtsMax)
    real, intent(in) :: srtsMax
    maxPoints_srts = ceiling((srtsMax - srtsMin) / dSrts)
  end subroutine setMaxSqrts


  !****************************************************************************
  !****s* dimi/dimidm
  ! NAME
  ! subroutine dimidm(dsdm, dsdmIntegrated, sigma, srts, mdel)
  ! PURPOSE
  ! this subroutine calculates the cross sections for p p <-> Delta++ n
  !
  ! Stores cross sections at first call to a field and then it only returns
  ! those field values.
  ! INPUTS
  ! * real :: srts -- sqrt(s) in the process
  ! * real :: mdel -- mass of delta
  ! OUTPUT
  ! * real :: sigma -- total cross section (for reaction with incoming delta)
  ! * real :: dsdm  --  mass differential cross section as function of mdel
  !   and srts (for reaction with outgoing delta)
  ! * real :: dsdm_Integrated -- mass differential cross section integrated
  !   over mass as function of srts (for reaction with outgoing delta)
  ! NOTES
  ! The return values are not cross sections, but have to be divided by p_in
  !****************************************************************************
  subroutine dimidm(dsdm, dsdm_Integrated, sigma, srts, mdel)

    use particleProperties, only: hadron
    use idTable, only: Delta

    real, intent(out) :: dsdm, dsdm_Integrated, sigma
    real, intent(in) :: srts, mdel

    logical,save       :: initFlag=.true.
    integer, parameter :: maxPoints_mass=4000

    real, allocatable, dimension(:,:), save :: dsdm_save, sigma_save
    real, allocatable, dimension(:)  , save :: dsdm_integrated_save
    real :: mass, dMass_new, srts_down, mass_down, wSrts, wMass
    integer :: iSrts, iMass

    if (initFlag) call init

    iSrts = floor((srts-srtsMin)/dSrts)
    iMass = floor((mdel-hadron(Delta)%minMass)/dMass)

    if (iSrts<0 .or. iMass<0) then
       dsdm = 0.
       dsdm_integrated = 0.
       sigma = 0.
    else if (iSrts>maxPoints_srts-1 .or. iMass>maxPoints_mass-1) then
       ! no fast possibility to get cross section
       ! give warning to tell user that system slows down
       write(*,*) 'Warning!!!!!! dimidm is out of bounds.'
       write(*,*) ' Need to calculate it. System therefore slow. If this happens too often, then change bounds.'
       write(*,*) 'SquareRoot(s)= ',srts, iSrts, maxPoints_srts
       write(*,*) 'Mass of delta= ',mdel, iMass, maxPoints_mass
       dsdm_integrated=0.

       if (iMass > maxPoints_mass) then
          ! Increase mass step size
          dMass_new=dMass*10.
          iMass=Max(0, NINT((mdel-hadron(Delta)%minMass)/dMass_new) )
          if (iMass.gt.maxPoints_mass) then
             write(*,*) 'Critical error.stop!'
          end if
          write(*,*) 'New indizes for mass= ',iMass, maxPoints_mass
       else
          dMass_new=dMass
       end if

       do iMass=1,maxPoints_mass
          mass = hadron(Delta)%minMass + Real(iMass) * dMass_new
          call calculate(dsdm,sigma,srts,mass)
          dsdm_integrated = dsdm_integrated+dsdm * dMass_new
       end do

       call calculate(dsdm,sigma,srts,mdel)
    else
       srts_down = srtsMin + float(iSrts) * dSrts
       mass_down = hadron(Delta)%minMass + float(iMass) * dMass
       wSrts = (srts-srts_down)/dSrts
       wMass = (mdel-mass_down)/dMass

       dsdm = dsdm_save(iSrts  , iMass  ) * (1.-wSrts) * (1.-wMass) &
            + dsdm_save(iSrts  , iMass+1) * (1.-wSrts) * wMass      &
            + dsdm_save(iSrts+1, iMass  ) * wSrts      * (1.-wMass) &
            + dsdm_save(iSrts+1, iMass+1) * wSrts      * wMass

       dsdm_integrated = dsdm_integrated_save(iSrts  ) * (1.-wSrts) &
                       + dsdm_integrated_save(iSrts+1) * wSrts

       sigma = sigma_save(iSrts  , iMass  ) * (1.-wSrts) * (1.-wMass) &
             + sigma_save(iSrts  , iMass+1) * (1.-wSrts) * wMass      &
             + sigma_save(iSrts+1, iMass  ) * wSrts      * (1.-wMass) &
             + sigma_save(iSrts+1, iMass+1) * wSrts      * wMass

    end if

  contains

    subroutine init
      use output, only: Write_InitStatus

      integer :: iSrts, iMass
      real :: wurzelS

      call Write_InitStatus('Dimitriev Xsections',0)

      allocate(dsdm_save(0:maxPoints_srts, 0:maxPoints_mass))
      allocate(sigma_save(0:maxPoints_srts, 0:maxPoints_mass))
      allocate(dsdm_integrated_save(0:maxPoints_srts))

      do iSrts=0,maxPoints_srts
        do iMass=0,maxPoints_mass
          mass    = hadron(Delta)%minMass + Real(iMass) * dMass
          wurzelS = srtsMin               + Real(iSrts) * dSrts
          call calculate(dsdm,sigma,wurzelS,mass)
          dsdm_save(iSrts, iMass)=dsdm
          sigma_save(iSrts, iMass)=sigma
        end do
        ! Integrate dsDm_save over the mass index
        dsdm_integrated_save(iSrts) = Sum(dsdm_save(iSrts, : )) * dMass
      end do
      initFlag=.false.

      call Write_InitStatus('Dimitriev Xsections',1)

    end subroutine init


    subroutine calculate(dsdm, sigma, srts, mdel)
      use constants, only: pi, mN, mPi, GeVSquared_times_mb
      use baryonWidth, only: fullWidthBaryon

      real, intent(out) :: dsdm, sigma
      real, intent(in)  :: srts, mdel

      real :: t,pi2,pf2,g,f,cost,u,E1,E2,E3,E4,melem2,dsdodm,dsdo,fac
      integer :: i
      integer, parameter :: nt = 100

      if (srts<2*mN+mPi .or. srts<mdel+mN) then
         dsdm = 0.0
         sigma = 0.0
         return
      end if

      pi2 = srts**2/4.-mN**2
      pf2 = max((srts**2-(mdel+mN)**2)*(srts**2-(mdel-mN)**2)/(4.*srts**2),0.)

      E1 = sqrt(mN**2 + pi2)
      E2 = E1
      E3 = sqrt(mN**2+pf2)
      E4 = sqrt(mdel**2+pf2)

      if (mdel<mN+mPi) then
         f = 0.
      else
         g=FullWidthBaryon(delta,mdel)
         !         g=0.120
         f=1./pi*mdel*g/((mdel**2-hadron(Delta)%mass**2)**2+(mdel*g)**2)
      end if

      dsdm=0.
      sigma=0.
      do i=-nt,nt

         cost=float(i)/float(nt)

         t= mdel**2 + mN**2 - 2.*E4*E2 + 2.*cost*sqrt(pi2*pf2)
         u= 2.*mN**2 - 2.*E3*E2 - 2.*cost*sqrt(pi2*pf2)

         !****    Relativistic matrix element^2, with some factors
         melem2 = 1./(32.*pi*srts**2) * mnnnd2(t,u,mdel) / GeVSquared_times_mb

         !****    p p -> n Delta++ :
         dsdodm = melem2 * sqrt(pf2) * 2.*mdel*f

         !****    n Delta++ -> p p, factor 1/4 because of spin and identical particles :
         dsdo = melem2 * sqrt(pi2) / 4.

         fac=1.
         if (abs(i).eq.nt) fac=0.5
         dsdm=dsdm+dsdodm/float(nt)*fac
         sigma=sigma+dsdo/float(nt)*fac

      end do

    end subroutine calculate

  end subroutine dimidm



  !****************************************************************************
  !****f* dimi/mnnnd2
  ! NAME
  ! real function mnnnd2(t,u,mdel)
  !
  ! PURPOSE
  ! return Matrix element squared and averaged over initial and summed over
  ! final spins for p_1 p_2 -> n_3 delta_4^++
  ! (within Dmitriev and Sushkov model with effective nucleon and delta masses).
  !
  ! INPUTS
  ! * real :: t, u  -- Mandelstam variables
  ! * real :: mdel  -- current mass of delta
  !
  !****************************************************************************
  real function mnnnd2(t, u, mdel)

    use constants, only: mN
    use particleProperties, only: hadron
    use idTable, only: Delta
    use twoBodyTools, only: pCM_sqr

    real, intent(in) :: t, u, mdel

    real :: zt,zu,m1,m2,m12,ft,fu,pt,pu,dmass2,mdel2

    real, parameter :: fps    = 2.202     ! coupling constant pi-N-N vertex
    real, parameter :: fp     = 1.008     ! coupling constant pi-N-Delta vertex
    real, parameter :: lamda2 = 0.63**2   ! cutoff parameter in the form factor [in GeV^2]
    real, parameter :: c2     = 0.2**2    ! kappa^2 [in GeV^2]
    real, parameter :: mpi    = 0.14      ! pion mass in GeV, slightly different from GiBUU value. Don't change !!!
    real, parameter :: mpi2 = mpi**2, mn2 = mn**2
    real, parameter :: gp = fp * 2. * mn/mpi
    logical, parameter :: debugflag = .false.

    dmass2 = hadron(Delta)%mass**2
    mdel2  = mdel**2

    zt = (pCM_sqr(dmass2,mn2,t) + c2) / (pCM_sqr(mdel2, mn2,t) + c2)  ! decay form factor, eq. (12) in Dimitriev paper
    zu = (pCM_sqr(dmass2,mn2,u) + c2) / (pCM_sqr(mdel2 ,mn2,u) + c2)  !   "

    if (debugflag) write(*,*) 'zt,zu:',zt,zu

    m1  = t * (t-(mdel-mn)**2) * ((mdel+mn)**2-t)**2 / (3.*mdel2)   ! matrix element, eq. (7)
    m2  = u * (u-(mdel-mn)**2) * ((mdel+mn)**2-u)**2 / (3.*mdel2)   !  "

    ! eq. (8), note: there is a typo in the Dimitriev paper!
    ! eq. (3.26) in Effenberger Diploma, with minus sign from isospin factor, table (3.3)
    m12 = 1./(2.*mdel2) * &
          (          (t*u+(mdel2-mn2)*(t+u)-mdel**4+mn**4) * (t*u+mn*(mdel+mn)*(mdel2-mn2)) &
           - 1./3. * (t*u-(mdel+mn)**2*(t+u)+(mdel+mn)**4) * (t*u-mn*(mdel-mn)*(mdel2-mn2)) )

    if (debugflag) write(*,*) ' m1,m2,m12:',m1,m2,m12

    ft = (lamda2-mpi2)/(lamda2-t)   ! vertex form factor, eq. (4)
    fu = (lamda2-mpi2)/(lamda2-u)   !   "
    pt = 1./(t-mpi2)
    pu = 1./(u-mpi2)

    if (debugflag) write(*,*) ' ft,fu,pt,pu:',ft,fu,pt,pu

    mnnnd2 = (fps*gp/mpi)**2 * &
             ( ft**4*zt*pt**2*m1 + fu**4*zu*pu**2*m2 + (ft*fu)**2*sqrt(zt*zu)*pt*pu*m12 )  ! eq. (7),(8)

    if (debugflag) write(*,*) ' mnnnd2:', mnnnd2

  end function mnnnd2


end module dimi
