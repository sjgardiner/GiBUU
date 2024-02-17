!******************************************************************************
!****m* /helicityAmplitudes
! NAME
! module helicityAmplitudes
! PURPOSE
! * Provides the helicity amplitudes of the MAID analysis
!   -> http://www.kph.uni-mainz.de/MAID/maid2003/maid2003.html
! * Modified version of a source code provided by L.Tiator
! * Citation: D.Drechsel, S.S.Kamalov, L.Tiator; Nucl. Phys. A645 (1999) 145-174
! * Also holds the definitions of MAID 2005 and MAID 2007
!******************************************************************************
module helicityAmplitudes

  private

  public :: get_helicityAmplitudes,writeOutFORMFACTOR
  public :: HP33_MAID05  ! Delta(1232) form factor

  public :: setScale1535


  real, parameter :: X3P33=1
  real, parameter :: X1P33=1
  real, parameter :: XSP33=1
  real, parameter :: X1S31=1
  real, parameter :: XSS31=1
  real, parameter :: X1D33=1
  real, parameter :: X3D33=1
  real, parameter :: XSD33=1
  real, parameter :: X1P11p=1
  real, parameter :: XSP11p=1
  real, parameter :: X1S11p=1
  real, parameter :: XSS11p=1
  real, parameter :: X1P11n=1
  real, parameter :: XSP11n=1
  real, parameter :: X1S11n=1
  real, parameter :: XSS11n=1
  real, parameter :: X1D13p=1
  real, parameter :: X3D13p=1
  real, parameter :: XSD13p=1
  real, parameter :: X1F15p=1
  real, parameter :: X3F15p=1
  real, parameter :: XSF15p=1
  real, parameter :: X1D13n=1
  real, parameter :: X3D13n=1
  real, parameter :: XSD13n=1
  real, parameter :: X1F15n=1
  real, parameter :: X3F15n=1
  real, parameter :: XSF15n=1
  real, parameter :: X1S2p=1
  real, parameter :: XSS2p=1
  real, parameter :: X1S2n=1
  real, parameter :: XSS2n=1
  real, parameter :: X1P13p=1
  real, parameter :: X3P13p=1
  real, parameter :: XSP13p=1
  real, parameter :: X1P31=1
  real, parameter :: XSP31=1
  real, parameter :: X1P13n=1
  real, parameter :: X3P13n=1
  real, parameter :: XSP13n=1
  real, parameter :: X1D15p=1
  real, parameter :: X3D15p=1
  real, parameter :: XSD15p=1
  real, parameter :: X1D15n=1
  real, parameter :: X3D15n=1
  real, parameter :: XSD15n=1
  real, parameter :: X1F35=1
  real, parameter :: X3F35=1
  real, parameter :: XSF35=1
  real, parameter :: X1F37=1
  real, parameter :: X3F37=1
  real, parameter :: XSF37=1

  logical, save :: scale1535 = .false. ! [GiBUU]
  ! only allowed for photoproduction, cf. setScale1535 below!!!

contains

  !****************************************************************************
  !****s* helicityAmplitudes/setScale1535
  ! NAME
  ! subroutine setScale1535(flag)
  !
  ! PURPOSE
  ! set a flag, whether the helicity amplitudes of the N*(1535) should be
  ! rescaled or not. Rescaling is necesarry in GiBUU for getting the
  ! photoproduction of eta mesons (quite) right.
  !
  ! Since only the A_1/2 amplitudes are rescaled and not the S_1/2 (which get
  ! importance for Q^2>0), this rescaling works only for photoproduction
  ! and not for electroproduction.
  !
  ! INPUTS
  ! * logical :: flag --- value the internal flag should get
  !****************************************************************************
  subroutine setScale1535(flag)
    logical, intent(in) :: flag

    scale1535 = flag
  end subroutine setScale1535


  !****************************************************************************

  subroutine writeOutFORMFACTOR()
    implicit none
    integer :: ig,iso
    real(8) :: Q2G
    do  IG=1,3
       Q2G=(IG-1)*0.5
       do  ISO=1,2
          CALL HEL_OUT_MAID03(ISO,Q2G)
       end do
       write(6,*)
    end do
    STOP
  end subroutine writeOutFORMFACTOR



  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !
  ! MAID 2003
  !
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  subroutine HP33_MAID03(Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    Fq=exp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)
    !c *****************************************************************
    Phi_R=0.
    A10=-140.250*X1P33
    A30=-265.437*X3P33
    S10=  27.577*XSP33

    A1= A10*(1.+ 0.0214486*Q2G)*Fq
    A3= A30*(1.- 0.0065431*Q2G)*Fq
    S1= S10*(1.+ 0.0166667*Q2G**3)*Fq

    AE= (1./2.)*(A3/sqrt(3.) -A1)
    AM=-(1./2.)*(Sqrt(3.0)*A3+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1P33, A1,X3P33, A3, XSP33, S1
123 format( 5x,'P33(1232):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HP33_MAID03

  subroutine HP11_MAID03(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)

    PI=3.1415926536D0
    Phi_R = -15.48
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    !c ***************************************************************
    A10 = -74.585/COSR*X1P11p
    S10 = 32.995/COSR*XSP11p
    A1= A10*(1.- 0.95*Q2G-0.95*Q2G**4)*exp(-1.55*Q2G)
    S1= S10*(1.+ 7.001*Q2G)*exp(-2.1*Q2G)
    GO TO 20
    !c ***************************************
10  A10 = 49.926/COSR*X1P11n
    S10 = -40./COSR*XSP11n
    A1= A10*(1.+ 0.7*Q2G)*exp(-1.77*Q2G)
    S1=S10*(1.+ 2.97843*Q2G)*exp(-1.55*Q2G)
    !c **************************************************
20  AM=A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123)  Phi_R, X1P11p, A1, XSP11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123 format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HP11_MAID03

  subroutine HD13_MAID03(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 32.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-25.7796/COSR*X1D13p
    A30=140.439/COSR*X3D13p
    S10=-27.14/COSR*XSD13p

    A1= A10*(1.+6.2193*Q2G)*exp(-1.08576*Q2G)
    A3= A30*(1.+2.2357*Q2G)*exp(-2.54217*Q2G)
    S1= S10*(1.+ 8.429*Q2G)*exp(-3.4534*Q2G)
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 19.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-81.02/COSR*X1D13n
    A30=-140.6/COSR*X3D13n
    S10= 12.85/COSR*XSD13n

    A1= A10*(1.-2.31468*Q2G)*exp(-1.55*Q2G)
    A3= A30*(1.+1.0903*Q2G)*exp(-1.75*Q2G)

    S1= S10*(1.+7.13437*Q2G)*exp(-1.57380*Q2G)
    !c ***************************************************
20  AE=-(1./2.)*(sqrt(3.)*A3+A1)
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123)  Phi_R, X1D13p, A1, X3D13p, A3, XSD13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123)  Phi_R, X1D13n, A1, X3D13n, A3, XSD13n, S1
123 format( 5x,'D13(1520):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HD13_MAID03

  subroutine HD33_MAID03(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)

    PI=3.1415926536D0
    Phi_R = 61.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=65.5267/COSR*X1D33
    A30=103.0903/COSR*X3D33
    S10=1./COSR*XSD33
    !c *************************************************

    A1= A10*(1.+3.83921*Q2G)*exp(-1.77207*Q2G)
    A3= A30*(1.+1.9722*Q2G)*exp(-2.2*Q2G)
    S1= S10/(1.+Q2G/0.71)**2

    AE=-(1./2.)*(A3*sqrt(3.)+A1)
    AM=-(1./2.)*(A3/sqrt(3.)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1D33, A1, X3D33, A3, XSD33, S1
123 format( 5x,'D33(1700):',1x,F6.2,3(1x,2(F8.3,1x)))

    Return
  end subroutine HD33_MAID03

  subroutine HS11f_MAID03(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !    COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
    !    COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=8.193
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=72.482/COSR*X1S11p
    S10=-15.70/COSR*XSS11p

    A1=A10*(1.+1.34470*Q2G)*exp(-0.75879*Q2G)
    S1=S10*(1.+2.8261*Q2G)*exp(-0.73735*Q2G)

    GO TO 20
    !c **************************************************

10  A10=-41.94/COSR*X1S11n
    S10=28.18/COSR*XSS11n

    A1=A10*(1.+4.33783*Q2G)*exp(-1.68723*Q2G)
    S1=S10*(1.+0.35874*Q2G)*exp(-1.55*Q2G)

    !c ****************************************************
20  AE=-A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123)  Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1535):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11f_MAID03

  subroutine HS11s_MAID03(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !   COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=6.961
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=32.0913/COSR*X1S11p
    S10=-3.489/COSR*XSS11p

    A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
    S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

    GO TO 20
    !c **************************************************

10  A10=26.32/COSR*X1S11n
    S10=10./COSR*XSS11n

    A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
    S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

    !c ********************************************
20  AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11s_MAID03


  subroutine HS31_MAID03(Q2G,PhiR,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !    COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=22.54
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    A10 = 65.2101/COSR*X1S31
    S10 =1./COSR*XSS31

    A1=A10*(1.+1.5*Q2G)*exp(-3.*Q2G)
    S1=S10/(1.+Q2G/0.71)**2
    AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1S31, A1, XSS31, S1
123 format( 5x,'S31(1620):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))
    return
  end subroutine HS31_MAID03


  subroutine HF15_MAID03(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !   COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
    !  COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
    !c ***************************************************************
    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 10.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-24.2131/COSR*X1F15p
    A30=132.2624/COSR*X3F15p
    S10=-25.65/COSR*XSF15p

    A1= A10*(1.+3.72353*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+1.53671*Q2G)*exp(-2.22357*Q2G)
    S1= S10*(1.+5.5987*Q2G)*exp(-1.55*Q2G)
    GO TO 20
    !c ***************************************************************
10  CONTINUE
    Phi_R = 15.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=24.49/COSR*X1F15n
    A30=-33.70/COSR*X3F15n
    S10=1./COSR*XSF15n

    A1= A10*(1.+0.*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+4.*Q2G)*exp(-1.75*Q2G)
    S1= S10/(1.+Q2G/0.71)**2
    !c *******************************************************
20  AE=-(1./3.)*( A3*sqrt(2.)+A1)
    AM=-(1./3.)*( A3/sqrt(2.)-A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123) Phi_R,  X1F15p, A1, X3F15p, A3, XSF15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123) Phi_R,  X1F15n, A1, X3F15n, A3, XSF15n, S1
123 format( 5x,'F15(1680):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF15_MAID03

  subroutine HP31_MAID03(Q2G,PhiR,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    ! COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31

    PI=3.1415926536D0
    Phi_R=35.
    PhiR = Phi_R*PI/180.
    COSR=COS(PhiR)
    A10=34.03/COSR*X1P31
    S10=1./COSR*XSP31

    A1= A10/(1.+Q2G/0.71)**2
    S1= S10/(1.+Q2G/0.71)**2

    AM= A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1P31, A1, XSP31, S1
123 format( 5x,'P31(1910):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    Return
  end subroutine HP31_MAID03

  subroutine HD15_MAID03(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !    COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 20.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=21.397/COSR*X1D15p
    A30=22.600/COSR*X3D15p
    S10=1./COSR*XSD15p

    A1=A10/(1.+Q2G/0.71)**2
    A3= A30/(1.+Q2G/0.71)**2
    S1=S10/(1.+Q2G/0.71)**2
    GO TO 20
    !c ********************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=-61.45/COSR*X1D15n
    A30=-74.31/COSR*X3D15n
    S10=-1./COSR*XSD15n

    A1=A10/(1.+Q2G/0.71)**2
    A3= A30/(1.+Q2G/0.71)**2
    S1=S10/(1.+Q2G/0.71)**2


    !c *************************************************
20  AE= (1./3.)*(A3/sqrt(2.)-A1)
    AM=-(1./3.)*(A3*sqrt(2.)+A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3)  write(6,123) Phi_R,  X1D15p, A1, X3D15p, A3, XSD15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4)  write(6,123) Phi_R,  X1D15n, A1, X3D15n, A3, XSD15n, S1
123 format( 5x,'D15(1675):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HD15_MAID03


  subroutine HF35_MAID03(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !   COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R =40.
    PhiR =Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=11.092/COSR*X1F35
    A30=-16.63/COSR*X3F35
    S10=1./COSR*XSF35

    A1= A10/(1.+Q2G/0.71)**2
    A3= A30/(1.+Q2G/0.71)**2
    S1= S10/(1.+Q2G/0.71)**2

    AE=-(1./3.)*(A3*sqrt(2.)+A1)
    AM=-(1./3.)*(A3/sqrt(2.)-A1)
    AS=-(1./3.)*S1*sqrt(2.)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F35, A1, X3F35, A3, XSF35, S1
123 format( 5x,'F35(1905):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF35_MAID03

  subroutine HF37_MAID03(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !    COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R = 30.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-67.55/COSR*X1F37
    A30=-87.21/COSR*X3F37
    S10=1./COSR*XSF37

    A1=A10/(1.+Q2G/0.71)**2
    A3=A30/(1.+Q2G/0.71)**2
    S1=S10/(1.+Q2G/0.71)**2

    AE= (1./4.)*(A3*3./sqrt(15.)-A1)
    AM=-(1./4.)*(A3*5./sqrt(15.)+A1)
    AS=-(1./4.)*Sqrt(2.)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F37, A1, X3F37, A3, XSF37, S1
123 format( 5x,'F37(1950):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF37_MAID03

  subroutine HP13_MAID03(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !   COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
    !  COMMON/parPPn/ X1P13n,X3P13n,XSP13n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=54.89/COSR*X1P13p
    A30=-31.69/COSR*X3P13p
    S10=-53.03/COSR*XSP13p

    A1= A10*(1.+4.24137*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+3.9974*Q2G)*EXP(-1.55*Q2G)
    S1= S10*(1.+2.45819*Q2G)*EXP(-1.55*Q2G)
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=16.68/COSR*X1P13n
    A30=-74.93/COSR*X3P13n
    S10=-1./COSR*XSP13n

    A1= A10*(1.+4.24137*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+0.9974*Q2G)*EXP(-1.55*Q2G)
    S1= S10/(1.+Q2G/0.71)**2

    !c ****************************************************

20  AE= (1./2.)*(A3/sqrt(3.)-A1)
    AM=-(1./2.)*(A3*sqrt(3.)+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1P13p, A1, X3P13p, A3, XSP13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R, X1P13n, A1, X3P13n, A3, XSP13n, S1
123 format( 5x,'P13(1720):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HP13_MAID03


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !
  ! MAID 2005
  !
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  subroutine HP33_MAID05(Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    ! ***************************************************************
    Fq=exp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)
    ! *****************************************************************
    Phi_R=0.
    !C ******  fit parameters changed to XE and XM, 19 April 2005
    !C ***     changes onlyn inside of HP33
    !C ***     outside XMP33 is still calles X3P33 and XEP33 -> X1P33
    XMP33=X3P33
    XEP33=X1P33

    A10=-137.445 *X1P33
    A30=-260.128 *X3P33
    S10=  27.577 *XSP33

    A1= A10*(1.+ 0.0214486*Q2G)*Fq
    A3= A30*(1.- 0.0065431*Q2G)*Fq
    S1= S10*(1.+ 0.0166667*Q2G**3)*Fq

    AE= (1./2.)*(A3/sqrt(3.) -A1)*XEP33
    AM=-(1./2.)*(Sqrt(3.0)*A3+A1)*XMP33
    AS=-(1./2.)*Sqrt(2.0)*S1

    A1=-(1./2.)*(3*AE+AM)
    A3= (1./2.)*Sqrt(3.)*(AE-AM)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1P33, A1,X3P33, A3, XSP33, S1
123 format( 5x,'P33(1232):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HP33_MAID05

  subroutine HP11_MAID05(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    PI=3.1415926536D0
    Phi_R = -15.48
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    !c ***************************************************************
    A10 = -59.134/COSR  *X1P11p
    S10 = 32.995/COSR   *XSP11p
    A1= A10*(1.- 0.95*Q2G -1.007*Q2G**4)*exp(-1.51*Q2G)
    S1= S10*(1.+ 7.001*Q2G)*exp(-2.1*Q2G)
    GO TO 20
    !c ***************************************
10  A10 = 52.137/COSR *X1P11n
    S10 = -40./COSR   *XSP11n
    A1= A10*(1.+0.9450*Q2G)*exp(-1.77*Q2G)
    S1=S10*(1.+ 2.97843*Q2G)*exp(-1.55*Q2G)
    !c **************************************************
20  AM=A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123)  Phi_R, X1P11p, A1, XSP11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123 format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HP11_MAID05

  subroutine HD13_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
    !       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
    !c ***************************************************************
    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    XED13p=1
    XMD13p=1

    Phi_R = 32.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-23.2016/COSR *X1D13p
    A30=136.2258/COSR *X3D13p
    S10=-27.14/COSR   *XSD13p

    A1= A10*(1.+6.84123*Q2G)*exp(-1.08576*Q2G)
    A3= A30*(1.+2.2357*Q2G)*exp(-2.54217*Q2G)
    S1= S10*(1.+6.06888*Q2G)*exp(-3.4534*Q2G)
    AE=-(1./2.)*(sqrt(3.)*A3+A1)   * XED13p
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)  * XMD13p
    AS=-(1./2.)*Sqrt(2.0)*S1
    A1=(3*AM-AE)/2.0
    A3=-Sqrt(3.0)*(AE+AM)/2.0
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 19.
    PhiR=Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-72.362/COSR   *X1D13n
    A30=-145.620/COSR  *X3D13n
    S10= 12.85/COSR    *XSD13n

    A1= A10*(1.-1.20549*Q2G)*exp(-1.55*Q2G)
    A3= A30*(1.+0.24466*Q2G)*exp(-1.75*Q2G)
    S1= S10*(1.+7.13437*Q2G)*exp(-1.57380*Q2G)
    AE=-(1./2.)*(sqrt(3.)*A3+A1)
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1
    !c ***************************************************
20  CONTINUE

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123)  Phi_R, X1D13p, A1, X3D13p, A3, XSD13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123)  Phi_R, X1D13n, A1, X3D13n, A3, XSD13n, S1
123 format( 5x,'D13(1520):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HD13_MAID05

  subroutine HD33_MAID05(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
    !c*****************************************************************
    PI=3.1415926536D0
    Phi_R = 61.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=109.631/COSR *X1D33
    A30=101.914/COSR *X3D33
    S10=1./COSR *XSD33
    !c *************************************************

    A1= A10*(1.+1.906628*Q2G)*exp(-1.77207*Q2G)
    A3= A30*(1.+1.9722*Q2G)*exp(-2.2*Q2G)
    S1= S10*exp(-2.0*Q2G)

    AE=-(1./2.)*(A3*sqrt(3.)+A1)
    AM=-(1./2.)*(A3/sqrt(3.)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1D33, A1, X3D33, A3, XSD33, S1
123 format( 5x,'D33(1700):',1x,F6.2,3(1x,2(F8.3,1x)))

    Return
  end subroutine HD33_MAID05

  subroutine HS11f_MAID05(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
    !COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=8.193
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=65.751/COSR   *X1S11p
    S10=-15.70/COSR   *XSS11p

    A1=A10*(1.+1.61364*Q2G)*exp(-0.75879*Q2G)
    S1=S10*(1.+2.8261*Q2G)*exp(-0.73735*Q2G)

    GO TO 20
    !c **************************************************

10  A10=-50.148/COSR *X1S11n
    S10=28.18/COSR   *XSS11n

    A1=A10*(1.+2.86297*Q2G)*exp(-1.68723*Q2G)
    S1=S10*(1.+0.35874*Q2G)*exp(-1.55*Q2G)

    !c ****************************************************
20  AE=-A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123)  Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1535):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11f_MAID05

  subroutine HS11s_MAID05(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=6.961
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=33.0210/COSR *X1S11p
    S10=-3.489/COSR  *XSS11p

    A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
    S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

    GO TO 20
    !c **************************************************

10  A10=9.186/COSR *X1S11n
    S10=10./COSR   *XSS11n

    A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
    S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

    !c ********************************************
20  AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11s_MAID05

  subroutine HS31_MAID05(Q2G,PhiR,AE,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
    !c *****************************************************************
    PI=3.1415926536D0
    Phi_R=22.54
    PhiR=Phi_R*PI/180.
    COSR=COS(PhiR)

    A10 = 60.6258/COSR*X1S31
    S10 =1./COSR*XSS31

    A1=A10*(1.+1.5*Q2G)*exp(-3.*Q2G)
    S1=S10*exp(-2.*Q2G)

    AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1S31, A1, XSS31, S1
123 format( 5x,'S31(1620):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))
    return
  end subroutine HS31_MAID05

  subroutine HF15_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !c ***************************************************************
    !       COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
    !       COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n
    !c ***************************************************************
    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 10.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-24.7443/COSR   *X1F15p
    A30=132.2624/COSR   *X3F15p
    S10=-25.65/COSR     *XSF15p

    A1= A10*(1.+3.72353*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+1.352305*Q2G)*exp(-2.22357*Q2G)
    S1= S10*(1.+4.47896*Q2G)*exp(-1.55*Q2G)
    GO TO 20
    !c ***************************************************************
10  CONTINUE
    Phi_R = 15.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=26.94/COSR   *X1F15n
    A30=-37.07/COSR  *X3F15n
    S10=1./COSR*XSF15n

    A1= A10*(1.+0.001*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+3.*Q2G)*exp(-1.75*Q2G)
    S1= S10 *exp(-1.55*Q2G)
    !c *******************************************************
20  AE=-(1./3.)*( A3*sqrt(2.)+A1)
    AM=-(1./3.)*( A3/sqrt(2.)-A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R,  X1F15p, A1, X3F15p, A3, XSF15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1F15n, A1, X3F15n, A3, XSF15n, S1
123 format( 5x,'F15(1680):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF15_MAID05

  subroutine HP31_MAID05(Q2G,PhiR,AM,AS,A1,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31

    PI=3.1415926536D0
    Phi_R=35.
    PhiR = Phi_R*PI/180.
    COSR=COS(PhiR)
    A10=14.786/COSR*X1P31
    S10=1./COSR*XSP31

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AM= A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1P31, A1, XSP31, S1
123 format( 5x,'P31(1910):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    Return
  end subroutine HP31_MAID05

  subroutine HD15_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 20.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=14.356/COSR  *X1D15p
    A30=20.322/COSR  *X3D15p
    S10=1./COSR*XSD15p


    A1= A10*(1.+0.1*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.1*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)
    GO TO 20
    !c ********************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=-61.738/COSR *X1D15n
    A30=-83.868/COSR *X3D15n
    S10=-1./COSR*XSD15n

    A1= A10*(1.+0.01*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.01*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.01*Q2G)*exp(-2.*Q2G)


    !c *************************************************
20  AE= (1./3.)*(A3/sqrt(2.)-A1)
    AM=-(1./3.)*(A3*sqrt(2.)+A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R,  X1D15p, A1, X3D15p, A3, XSD15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R,  X1D15n, A1, X3D15n, A3, XSD15n, S1
123 format( 5x,'D15(1675):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HD15_MAID05

  subroutine HF35_MAID05(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R =40.
    PhiR =Phi_R*Pi/180.
    COSR=COS(PhiR)
    A10=13.934/COSR*X1F35
    A30=-21.427/COSR*X3F35
    S10=1./COSR*XSF35

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AE=-(1./3.)*(A3*sqrt(2.)+A1)
    AM=-(1./3.)*(A3/sqrt(2.)-A1)
    AS=-(1./3.)*S1*sqrt(2.)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F35, A1, X3F35, A3, XSF35, S1
123 format( 5x,'F35(1905):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF35_MAID05

  subroutine HF37_MAID05(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    PI=3.1415926536D0
    Phi_R = 30.
    PhiR = Phi_R*Pi/180.
    COSR=COS(PhiR)

    A10=-81.06/COSR*X1F37
    A30=-104.65/COSR*X3F37
    S10=1./COSR*XSF37

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AE= (1./4.)*(A3*3./sqrt(15.)-A1)
    AM=-(1./4.)*(A3*5./sqrt(15.)+A1)
    AS=-(1./4.)*Sqrt(2.)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F37, A1, X3F37, A3, XSF37, S1
123 format( 5x,'F37(1950):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF37_MAID05

  subroutine HP13_MAID05(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
    IMPLICIT real(8) (A-H,O-Z)
    !COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
    !COMMON/parPPn/ X1P13n,X3P13n,XSP13n

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=73.002/COSR   *X1P13p
    A30=-11.465/COSR  *X3P13p
    S10=-53.03/COSR   *XSP13p

    A1= A10*(1.+1.891651*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+15.9896*Q2G)*EXP(-1.55*Q2G)
    S1= S10*(1.+2.45819*Q2G)*EXP(-1.55*Q2G)
    GO TO 20
    !c ***************************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = COS(PhiR)

    A10=-2.904/COSR  *X1P13n
    A30=-30.972/COSR *X3P13n
    S10=-1./COSR   *XSP13n

    A1= A10*(1.+12.72411*Q2G)*EXP(-1.55*Q2G)
    A3= A30*(1.+4.987*Q2G)*EXP(-1.55*Q2G)
    S1= S10*EXP(-1.55*Q2G)
    !c ****************************************************

20  AE= (1./2.)*(A3/sqrt(3.)-A1)
    AM=-(1./2.)*(A3*sqrt(3.)+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,123) Phi_R, X1P13p, A1, X3P13p, A3, XSP13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,123) Phi_R, X1P13n, A1, X3P13n, A3, XSP13n, S1
123 format( 5x,'P13(1720):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HP13_MAID05

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************
  !
  ! MAID 2007
  !
  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

  subroutine HP33_MAID07(Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/param3/ XMP33,XEP33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33
!!$    COMMON /Mass/ ami, amf, amPION
!!$    COMMON /KinVar/ OmegL, Q2, W, wGacm, akcm

    implicit none
    real*8 :: Q2G,qcm,qcm0,AE,AM,AS,A1,A3,S1
    integer :: IPRN

    real*8 :: Fq,Phi_R, AE0,AM0, AE1,AM1, A11,A31
    real*8 :: XMP33,XEP33
    real*8 :: hX1P33,hX3P33
    real*8 :: CNORM,ami,W,Q2pt,Db,Q2MEV,qcm_D,tau,AS0

    ! ***************************************************************
    Fq=dexp(-0.21*Q2G)/(1+Q2G/0.71)**2*(qcm/qcm0)
    ! *****************************************************************
    Phi_R=0.

    XMP33=X3P33 ! [GiBUU] cf. MAID2005
    XEP33=X1P33 ! [GiBUU]

    AE0 = -6.36992
    AM0 = 299.880

    AE1= AE0*(1.-0.0205657*Q2G)*Fq*dexp(0.05*Q2G)
    !      AM= AM0*(1.-0.012*Q2G)*Fq
    AM1= AM0*(1.+0.0095*Q2G)*Fq*DEXP(-0.02*Q2G)
    A11=-(1./2.)*(3*AE1+AM1)
    A31= (1./2.)*Sqrt(3.)*(AE1-AM1)

    AE=AE1 * XEP33
    AM=AM1 * XMP33
    A1=-(1./2.)*(3*AE+AM)
    A3= (1./2.)*Sqrt(3.)*(AE-AM)

    hX1P33=A1/A11
    hX3P33=A3/A31
    ! ***********  parametrization based on Siegert theorem *******************
    CNORM=1.
!    if (IPRN.NE.1) GO TO 111 ! [GiBUU] probably a bug!
    ami=938.2723
    W = 1232.
    CNORM=1000.
!111 Q2pt= -(W-ami)**2*1.D-6 ! [GiBUU]
    Q2pt= -(W-ami)**2*1.D-6
    ! --------  our new parametrization -------------
    Db=4.9
    Q2MEV=Q2G*1.E+6
    qcm_D=(1232**2-ami**2)/2./1232.
    tau=Q2MEV/(4.*ami*ami)
    AS0=12.403
    S1=sqrt(2.)*AS0*(1.+0.12*Q2G)/(1+Db*tau)*Fq*DEXP(-0.02*Q2G) &
         *qcm/qcm_D * CNORM* XSP33

    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  hX1P33, A1,hX3P33, A3, XSP33, S1
123 format( 5x,'P33(1232):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HP33_MAID07

  subroutine HP11_MAID07(ISO,PhiR,Q2G,AM,AS,A1,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
!!$    COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n

    implicit none
    real*8 :: PhiR,Q2G,AM,AS,A1,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,S10
    real*8 :: a11,a12,a14,a1x, s11,s12,s14,s1x

    PI=3.1415926536D0
    Phi_R = -15.48
    PhiR=Phi_R*Pi/180.
    COSR=DCOS(PhiR)
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    ! ****************  2007 ****************************************
    !      A10 = -59.134/COSR  *X1P11p
    !      S10 = 4.02/COSR     *XSP11p
    !
    !      A1=A10*(1.-1.221691*Q2G-0.55*Q2G**4)*exp(-1.51*Q2G)
    !      S1= S10*(1.+ 41.001*Q2G+1.5*Q2G**4)*exp(-1.75*Q2G)
    !
    ! ****************  2008 ****************************************
    A10 = -61.360
    S10 =   4.171
    a11= 0.87139
    a12=-3.51582
    a14=-0.15832
    a1x=1.36169
    s11=40.0
    s12=0.0
    s14=1.50
    s1x=1.75
    A1=A10*(1.+a11*Q2G+a12*Q2G**2+a14*Q2G**4)*exp(-a1x*Q2G)
    S1=S10*(1.+s11*Q2G+s12*Q2G**2+s14*Q2G**4)*exp(-s1x*Q2G)
    GO TO 20
    ! ***************************************
10  A10 = 52.137/COSR *X1P11n
    S10 = -40./COSR   *XSP11n
    A1= A10*(1.+0.9450*Q2G)*exp(-1.77*Q2G)
    S1=S10*(1.+ 2.97843*Q2G)*exp(-1.55*Q2G)
    ! **************************************************
20  AM=A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123)  Phi_R, X1P11p, A1, XSP11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123)  Phi_R, X1P11n, A1, XSP11n, S1
123 format( 5x,'P11(1440):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HP11_MAID07

  subroutine HD13_MAID07(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
!!$    COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10
    real*8 :: a11,a12,a14,a1x, a31,a32,a34,a3x, s11,s12,s14,s1x

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 32.
    PhiR=Phi_R*Pi/180.
    COSR=DCOS(PhiR)
    ! ****************  2007 ****************************************
    !        A10=-23.2016/COSR *X1D13p
    !        A30=136.2258/COSR *X3D13p
    !        S10=-53.9019/COSR *XSD13p
    !      A1= A10*(1.+7.7698*Q2G)*exp(-1.08576*Q2G)
    !      A3= A30*(1.+0.69263*Q2G)*exp(-2.104*Q2G)
    !      S1= S10*(1.+4.19237*Q2G)*exp(-3.4*Q2G)
    !
    ! ****************  2008 ****************************************
    A10 =-27.359
    A30 =160.635
    S10 =-63.50
    a11= 8.5800
    a12=-0.25216
    a14= 0.35735
    a1x=1.200
    a31=-0.81961
    a32= 0.54119
    a34=-0.01552
    a3x=1.05625
    s11=4.190
    s12=0.0
    s14=0.0
    s1x=3.40
    A1=A10*(1.+a11*Q2G+a12*Q2G**2+a14*Q2G**4)*exp(-a1x*Q2G)
    A3=A30*(1.+a31*Q2G+a32*Q2G**2+a34*Q2G**4)*exp(-a3x*Q2G)
    S1=S10*(1.+s11*Q2G+s12*Q2G**2+s14*Q2G**4)*exp(-s1x*Q2G)
    GO TO 20
    ! ***************************************************
10  CONTINUE
    Phi_R = 19.
    PhiR=Phi_R*Pi/180.
    COSR=DCOS(PhiR)

    A10=-72.362/COSR   *X1D13n
    A30=-145.620/COSR  *X3D13n
    S10= 12.85/COSR    *XSD13n

    A1= A10*(1.-0.533924*Q2G)*exp(-1.55*Q2G)
    A3= A30*(1.+0.578587*Q2G)*exp(-1.75*Q2G)
    S1= S10*(1.+15.74199*Q2G)*exp(-1.5738*Q2G)
    ! ***************************************************
20  AE=-(1./2.)*(sqrt(3.)*A3+A1)
    AM=-(1./2.)*(A3/Sqrt(3.0)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123)  Phi_R, X1D13p, A1, X3D13p, A3, XSD13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123)  Phi_R, X1D13n, A1, X3D13n, A3, XSD13n, S1
123 format( 5x,'D13(1520):',1x,F6.2,3(1x,2(F8.3,1x)))
    return
  end subroutine HD13_MAID07

  subroutine HD33_MAID07(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10
    real*8 :: a11,a12,a14,a1x, a31,a32,a34,a3x, s11,s12,s14,s1x


    PI=3.1415926536D0
    Phi_R = 61.
    PhiR = Phi_R*Pi/180.
    COSR=DCOS(PhiR)

    ! ****************  2007 ****************************************
    !        A10=109.631/COSR *X1D33
    !        A30=101.914/COSR *X3D33
    !	   S10=1./COSR *XSD33
    !      A1= A10*(1.+1.906628*Q2G)*exp(-1.77207*Q2G)
    !      A3= A30*(1.+1.9722*Q2G)*exp(-2.2*Q2G)
    !      S1= S10*exp(-2.0*Q2G)
    !
    ! ****************  2008 ****************************************
    A10 =226.132
    A30 =210.214
    S10 =2.063
    a11= 1.910
    a12= 0.0
    a14= 0.0
    a1x=1.770
    a31= 0.87946
    a32= 1.70620
    a34= 0.0
    a3x=2.02291
    s11=0.0
    s12=0.0
    s14=0.0
    s1x=2.00
    A1=A10*(1.+a11*Q2G+a12*Q2G**2+a14*Q2G**4)*exp(-a1x*Q2G)
    A3=A30*(1.+a31*Q2G+a32*Q2G**2+a34*Q2G**4)*exp(-a3x*Q2G)
    S1=S10*(1.+s11*Q2G+s12*Q2G**2+s14*Q2G**4)*exp(-s1x*Q2G)
    ! *************************************************

    AE=-(1./2.)*(A3*sqrt(3.)+A1)
    AM=-(1./2.)*(A3/sqrt(3.)-A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1D33, A1, X3D33, A3, XSD33, S1
123 format( 5x,'D33(1700):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HD33_MAID07

  subroutine HS11f_MAID07(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/par11p/ X1P11p,XSP11p,X1S11p,XSS11p
!!$    COMMON/par11n/ X1P11n,XSP11n,X1S11n,XSS11n

    implicit none
    real*8 :: PhiR,Q2G,AE,AS,A1,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,S10

    PI=3.1415926536D0
    Phi_R=8.193
    PhiR=Phi_R*PI/180.
    COSR=DCOS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=65.751/COSR   *X1S11p
    S10=-2.0/COSR     *XSS11p

    if (scale1535) A10=A10*1.60 ! [GiBUU]

    A1=A10*(1.+1.6083226*Q2G)*exp(-0.70*Q2G)
    S1=S10*(1.+23.90148*Q2G)*exp(-0.81*Q2G)

    GO TO 20
    ! **************************************************

10  A10=-50.148/COSR *X1S11n
    S10=28.18/COSR   *XSS11n

    if (scale1535) A10=A10*1.70 ! [GiBUU]

    A1=A10*(1.+4.746117*Q2G)*exp(-1.68723*Q2G)
    S1=S10*(1.+0.35874*Q2G)*exp(-1.55*Q2G)

    ! ****************************************************
20  AE=-A1
    AS=-sqrt(2.)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123)  Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1535):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11f_MAID07

  subroutine HS11s_MAID07(ISO,PhiR,Q2G,AE,AS,A1,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/secS11/ X1S11p,XSS11p, X1S11n,XSS11n

    implicit none
    real*8 :: PhiR,Q2G,AE,AS,A1,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,S10

    PI=3.1415926536D0
    Phi_R=6.961
    PhiR=Phi_R*PI/180.
    COSR=DCOS(PhiR)

    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10

    A10=33.0210/COSR *X1S11p
    S10=-3.489/COSR  *XSS11p

    A1=A10*(1.+1.45359*Q2G)*exp(-0.6167*Q2G)
    S1=S10*(1.+2.878*Q2G)*exp(-0.75879*Q2G)

    GO TO 20
    ! **************************************************

10  A10=9.186/COSR *X1S11n
    S10=10./COSR   *XSS11n

    A1=A10*(1.+0.13305*Q2G)*exp(-1.55*Q2G)
    S1=S10*(1.-0.5*Q2G)*exp(-1.55*Q2G)

    ! ********************************************
20  AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123) Phi_R, X1S11p, A1, XSS11p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123) Phi_R,  X1S11n, A1, XSS11n, S1
123 format( 5x,'S11(1650):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    return
  end subroutine HS11s_MAID07


  subroutine HS31_MAID07(Q2G,PhiR,AE,AS,A1,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/param3/ X3P33,X1P33,XSP33,X1S31,XSS31,X1D33,X3D33,XSD33

    implicit none
    real*8 :: PhiR,Q2G,AE,AS,A1,S1
    integer :: IPRN !,ISO

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,S10

    PI=3.1415926536D0
    Phi_R=22.54
    PhiR=Phi_R*PI/180.
    COSR=DCOS(PhiR)

    A10 = 60.6258/COSR *X1S31
    S10 =15./COSR      *XSS31

    A1=A10*(1.+1.858125*Q2G)*exp(-2.5*Q2G)
    S1=S10*(1.+2.82996*Q2G)*exp(-2.*Q2G)

    AE=-A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R,  X1S31, A1, XSS31, S1
123 format( 5x,'S31(1620):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))
    return
  end subroutine HS31_MAID07


  subroutine HF15_MAID07(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parDFp/ X1D13p,X3D13p,XSD13p,X1F15p,X3F15p,XSF15p
!!$    COMMON/parDFn/ X1D13n,X3D13n,XSD13n,X1F15n,X3F15n,XSF15n

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10
    real*8 :: a11,a12,a14,a1x, a31,a32,a34,a3x, s11,s12,s14,s1x

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 10.
    PhiR = Phi_R*Pi/180.
    COSR=DCOS(PhiR)

    ! ****************  2007 ****************************************
    !        A10=-24.7443/COSR  *X1F15p
    !        A30=132.2624/COSR  *X3F15p
    !        S10=-43.2904/COSR  *XSF15p
    !      A1= A10*(1.+3.978924*Q2G )*exp(-1.2*Q2G)
    !      A3= A30*(1.+0.996276*Q2G )*exp(-2.22357*Q2G)
    !      S1= S10*(1.+3.138554*Q2G )*exp(-1.58*Q2G)
    !
    ! ****************  2008 ****************************************
    A10 =-25.126
    A30 =134.303
    S10 =-43.958
    a11= 3.780
    a12=-0.29242
    a14= 0.07961
    a1x=1.250
    a31= 1.01594
    a32= 0.22228
    a34= 0.23652
    a3x=2.41329
    s11=3.78338
    s12=0.0
    s14=0.0
    s1x=1.84897
    A1=A10*(1.+a11*Q2G+a12*Q2G**2+a14*Q2G**4)*exp(-a1x*Q2G)
    A3=A30*(1.+a31*Q2G+a32*Q2G**2+a34*Q2G**4)*exp(-a3x*Q2G)
    S1=S10*(1.+s11*Q2G+s12*Q2G**2+s14*Q2G**4)*exp(-s1x*Q2G)
    GO TO 20
    ! ***************************************************************
10  CONTINUE
    Phi_R = 15.
    PhiR = Phi_R*Pi/180.
    COSR=DCOS(PhiR)

    A10=26.94/COSR   *X1F15n
    A30=-37.07/COSR  *X3F15n
    S10=1./COSR*XSF15n

    A1= A10*(1.+0.001*Q2G)*exp(-1.2*Q2G)
    A3= A30*(1.+4.09308*Q2G)*exp(-1.75*Q2G)
    S1= S10 *exp(-1.55*Q2G)
    ! *******************************************************
20  AE=-(1./3.)*( A3*sqrt(2.)+A1)
    AM=-(1./3.)*( A3/sqrt(2.)-A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123) Phi_R,  X1F15p, A1, X3F15p, A3, XSF15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123) Phi_R,  X1F15n, A1, X3F15n, A3, XSF15n, S1
123 format( 5x,'F15(1680):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF15_MAID07

  subroutine HP31_MAID07(Q2G,PhiR,AM,AS,A1,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31

    implicit none
    real*8 :: PhiR,Q2G,AM,AS,A1,S1
    integer :: IPRN !,ISO

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,S10

    PI=3.1415926536D0
    Phi_R=35.
    PhiR = Phi_R*PI/180.
    COSR=DCOS(PhiR)
    A10=14.786/COSR*X1P31
    S10=1./COSR*XSP31

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AM= A1
    AS=-Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1P31, A1, XSP31, S1
123 format( 5x,'P31(1910):',1x,F6.2,1x,2(F8.3,1x),20x,2(F8.3,1x))

    Return
  end subroutine HP31_MAID07

  subroutine HD15_MAID07(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parD15/ X1D15p,X3D15p,XSD15p,X1D15n,X3D15n,XSD15n

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10
    real*8 :: a11,a12,a14,a1x, a31,a32,a34,a3x, s11,s12,s14,s1x

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 20.
    PhiR = Phi_R*Pi/180.
    COSR=DCOS(PhiR)
    !
    ! ****************  2007 ****************************************
    !      A10=14.356/COSR  *X1D15p
    !      A30=20.322/COSR  *X3D15p
    !      S10=1./COSR*XSD15p
    !      A1= A10*(1.+0.1*Q2G)*exp(-2.*Q2G)
    !      A3= A30*(1.+0.1*Q2G)*exp(-2.*Q2G)
    !      S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)
    !
    ! ****************  2008 ****************************************
    A10 =15.277
    A30 =21.626
    S10 =1.064
    a11= 1.55916
    a12=-1.26017
    a14= 0.0
    a1x=0.79729
    a31= 0.20020
    a32= 1.62340
    a34= 0.0
    a3x=0.83560
    s11=0.0
    s12=0.0
    s14=0.0
    s1x=2.00
    A1=A10*(1.+a11*Q2G+a12*Q2G**2+a14*Q2G**4)*exp(-a1x*Q2G)
    A3=A30*(1.+a31*Q2G+a32*Q2G**2+a34*Q2G**4)*exp(-a3x*Q2G)
    S1=S10*(1.+s11*Q2G+s12*Q2G**2+s14*Q2G**4)*exp(-s1x*Q2G)
    GO TO 20
    ! ********************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*Pi/180.
    COSR=DCOS(PhiR)
    A10=-61.738/COSR *X1D15n
    A30=-83.868/COSR *X3D15n
    S10=-1./COSR   *XSD15n

    A1= A10*(1.+0.01*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.01*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.01*Q2G)*exp(-2.*Q2G)


    ! *************************************************
20  AE= (1./3.)*(A3/sqrt(2.)-A1)
    AM=-(1./3.)*(A3*sqrt(2.)+A1)
    AS=-(1./3.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123) Phi_R,  X1D15p, A1, X3D15p, A3, XSD15p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123) Phi_R,  X1D15n, A1, X3D15n, A3, XSD15n, S1
123 format( 5x,'D15(1675):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HD15_MAID07


  subroutine HF35_MAID07(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: IPRN !,ISO

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10

    PI=3.1415926536D0
    Phi_R =40.
    PhiR =Phi_R*Pi/180.
    COSR=DCOS(PhiR)
    A10=13.934/COSR*X1F35
    A30=-21.427/COSR*X3F35
    S10=1./COSR*XSF35

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AE=-(1./3.)*(A3*sqrt(2.)+A1)
    AM=-(1./3.)*(A3/sqrt(2.)-A1)
    AS=-(1./3.)*S1*sqrt(2.)

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F35, A1, X3F35, A3, XSF35, S1
123 format( 5x,'F35(1905):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF35_MAID07

  subroutine HF37_MAID07(Q2G,PhiR,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parF3/ X1F35,X3F35,XSF35,X1F37,X3F37,XSF37

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: IPRN !,ISO

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10

    PI=3.1415926536D0
    Phi_R = 30.
    PhiR = Phi_R*Pi/180.
    COSR=DCOS(PhiR)

    A10=-81.06/COSR*X1F37
    A30=-104.65/COSR*X3F37
    S10=1./COSR*XSF37

    A1= A10*(1.+0.*Q2G)*exp(-2.*Q2G)
    A3= A30*(1.+0.*Q2G)*exp(-2.*Q2G)
    S1= S10*(1.+0.*Q2G)*exp(-2.*Q2G)

    AE= (1./4.)*(A3*3./sqrt(15.)-A1)
    AM=-(1./4.)*(A3*5./sqrt(15.)+A1)
    AS=-(1./4.)*Sqrt(2.)*S1

    if (IPRN.EQ.0) return
    write(6,123) Phi_R, X1F37, A1, X3F37, A3, XSF37, S1
123 format( 5x,'F37(1950):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HF37_MAID07

  subroutine HP13_MAID07(ISO,PhiR,Q2G,AE,AM,AS,A1,A3,S1,IPRN)
!!$    IMPLICIT REAL*8 (A-H,O-Z)
!!$    COMMON/parPPp/ X1P13p,X3P13p,XSP13p,X1P31,XSP31
!!$    COMMON/parPPn/ X1P13n,X3P13n,XSP13n

    implicit none
    real*8 :: PhiR,Q2G,AE,AM,AS,A1,A3,S1
    integer :: ISO,IPRN

    real*8 :: PI,Phi_R,COSR
    real*8 :: A10,A30,S10
    real*8 :: a11,a12,a14,a1x, a31,a32,a34,a3x, s11,s12,s14,s1x

    PI=3.1415926536D0
    if (ISO.EQ.2 .OR. ISO.EQ.4) GO TO 10
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = DCOS(PhiR)
    !
    ! ****************  2007 ****************************************
    !        A10=73.002/COSR   *X1P13p
    !        A30=-11.465/COSR  *X3P13p
    !	  S10=-53.03/COSR   *XSP13p
    !      A1= A10*(1.+1.891651*Q2G)*DEXP(-1.55*Q2G)
    !      A3= A30*(1.+15.9896*Q2G)*DEXP(-1.55*Q2G)
    !      S1= S10*(1.+2.45819*Q2G)*DEXP(-1.55*Q2G)
    !
    ! ****************  2008 ****************************************
    A10 = 73.002
    A30 =-11.465
    S10 =-53.030
    a11= 2.30411
    a12=-0.97720
    a14= 0.0
    a1x= 0.76990
    a31= 8.21035
    a32=12.49687
    a34= 0.0
    a3x= 1.02324
    s11= 2.460
    s12= 0.0
    s14= 0.0
    s1x= 1.550
    A1=A10*(1.+a11*Q2G+a12*Q2G**2+a14*Q2G**4)*exp(-a1x*Q2G)
    A3=A30*(1.+a31*Q2G+a32*Q2G**2+a34*Q2G**4)*exp(-a3x*Q2G)
    S1=S10*(1.+s11*Q2G+s12*Q2G**2+s14*Q2G**4)*exp(-s1x*Q2G)
    GO TO 20
    ! ***************************************************
10  CONTINUE
    Phi_R = 0.
    PhiR = Phi_R*PI/180.
    COSR = DCOS(PhiR)

    A10=-2.904/COSR  *X1P13n
    A30=-30.972/COSR *X3P13n
    S10=-1./COSR   *XSP13n

    A1= A10*(1.+12.72411*Q2G)*DEXP(-1.55*Q2G)
    A3= A30*(1.+4.987*Q2G)*DEXP(-1.55*Q2G)
    S1= S10*DEXP(-1.55*Q2G)
    ! ****************************************************

20  AE= (1./2.)*(A3/sqrt(3.)-A1)
    AM=-(1./2.)*(A3*sqrt(3.)+A1)
    AS=-(1./2.)*Sqrt(2.0)*S1

    if (IPRN.EQ.0) return
    if (ISO.EQ.1.OR.ISO.EQ.3) &
         write(6,123) Phi_R, X1P13p, A1, X3P13p, A3, XSP13p, S1
    if (ISO.EQ.2.OR.ISO.EQ.4) &
         write(6,123) Phi_R, X1P13n, A1, X3P13n, A3, XSP13n, S1
123 format( 5x,'P13(1720):',1x,F6.2,3(1x,2(F8.3,1x)))

    return
  end subroutine HP13_MAID07


  !****************************************************************************
  !****************************************************************************
  !****************************************************************************


  SUBROUTINE HEL_OUT_MAID03(ISO,Q2G)

    ! Q2G -- Q^2 in GeV
    ! ISO -- 1=proton
    !        2=neutron

    IMPLICIT real(8)(A-H,O-Z)
    DOUBLE PRECISION mi, kgcm0
    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,2) Q2G
2   FORMAT(5x,'proton e.m. helicity amplitudes at Q^2=',F7.3,2x,'in units 10^-3/sqrt(GeV)')
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,3) Q2G
3   FORMAT(5x,'neutron e.m. helicity amplitudes at Q^2=',F7.3,2x,'in units 10^-3/sqrt(GeV)')
    write(6,4)
4   format(17X,'Phi_R',6x,'X1',5x,'A1/2',8x,'X3',5x,'A3/2',8x,'XS', 5x,'S1/2')

    W0=1.232
    mi=0.9382723
    kgcm0 = (W0*W0-mi*mi)/2./W0
    egcm = (W0*W0-Q2G-mi*mi)/2./W0
    qcm=sqrt(egcm**2+Q2G)

    CALL HP33_MAID03( Q2G,qcm,real(kgcm0,8),DE,DM,DS,A1,A3,S1,1)
    CALL HP11_MAID03( ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    CALL HD13_MAID03( ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HS11f_MAID03(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    CALL HS31_MAID03( Q2G,PhiR,DM,DS,A1,S1,1)
    CALL HS11s_MAID03(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    CALL HD15_MAID03( ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HF15_MAID03( ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HD33_MAID03( Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL HP13_MAID03( ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL HF35_MAID03( Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL HP31_MAID03( Q2G,PhiR,DM,DS,A1,S1,1)
    CALL HF37_MAID03( Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)

    RETURN
  end subroutine HEL_OUT_MAID03

  SUBROUTINE HEL_OUT_MAID07(ISO,Q2G)
    IMPLICIT REAL*8(A-H,O-Z)
    REAL*8 mp, mdelta, kgcm0
    CHARACTER*10 Name
    CHARACTER*1 nuc(2)
    DATA nuc/'p','n'/
    DATA IFST/0/
    SAVE IFST
    mdelta=1.232
    mp=0.9382723
    kgcm0 = (mdelta**2-mp**2)/2./mdelta
    ! **************************************************************
    if (IFST.EQ.1) GO TO 999
    Name='P33(1232) '
    WP33=1.232
    CALL HP33_MAID07(0.D0,kgcm0,kgcm0,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WP33,GM0P33,GE0P33,GC0P33)
    write(7,90) Name,nuc(iso),A1,A3,S1,GM0P33,GE0P33,GC0P33
    !
    Name='P11(1440) '
    WP11=1.440
    CALL HP11_MAID07(ISO,PhiR,0.D0,DM,DS,A1,S1,1)
    DE=0
    A3=0
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WP11,GM0P11,GE0P11,GC0P11)
    write(8,90) Name,nuc(iso),A1,A3,S1,GM0P11,GE0P11,GC0P11
    GE0P11=1
    !
    Name='D13(1520) '
    WD13=1.520
    CALL HD13_MAID07(ISO,PhiR,0.D0,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WD13,GM0D13,GE0D13,GC0D13)
    write(9,90) Name,nuc(iso),A1,A3,S1,GM0D13,GE0D13,GC0D13
    !
    Name='S11(1535) '
    WS11f=1.535
    CALL HS11f_MAID07(ISO,PhiR,0.D0,DE,DS,A1,S1,1)
    DM=0
    A3=0
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WS11f,GM0S11f,GE0S11f,GC0S11f)
    write(10,90) Name,nuc(iso),A1,A3,S1,GM0S11f,GE0S11f,GC0S11f
    GM0S11f=1
    !
    Name='S11(1650) '
    WS11s=1.650
    CALL HS11s_MAID07(ISO,PhiR,0.D0,DE,DS,A1,S1,1)
    DM=0
    A3=0
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WS11s,GM0S11s,GE0S11s,GC0S11s)
    write(11,90) Name,nuc(iso),A1,A3,S1,GM0S11s,GE0S11s,GC0S11s
    GM0S11s=1
    !
    Name='S31(1620) '
    WS31=1.620
    CALL HS31_MAID07(0.D0,PhiR,DE,DS,A1,S1,1)
    DM=0
    A3=0
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WS31,GM0S31,GE0S31,GC0S31)
    write(12,90) Name,nuc(iso),A1,A3,S1,GM0S31,GE0S31,GC0S31
    GM0S31=1
    !
    Name='D15(1675) '
    WD15=1.675
    CALL HD15_MAID07(ISO,PhiR,0.D0,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WD15,GM0D15,GE0D15,GC0D15)
    write(13,90) Name,nuc(iso),A1,A3,S1,GM0D15,GE0D15,GC0D15
    !
    Name='F15(1680) '
    WF15=1.680
    CALL HF15_MAID07(ISO,PhiR,0.D0,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WF15,GM0F15,GE0F15,GC0F15)
    write(14,90) Name,nuc(iso),A1,A3,S1,GM0F15,GE0F15,GC0F15
    !
    Name='D33(1700) '
    !	WD33=1.740
    WD33=1.700
    CALL HD33_MAID07(0.D0,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WD33,GM0D33,GE0D33,GC0D33)
    write(15,90) Name,nuc(iso),A1,A3,S1,GM0D33,GE0D33,GC0D33
    !
    Name='P13(1720) '
    WP13=1.720
    CALL HP13_MAID07(ISO,PhiR,0.D0,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WP13,GM0P13,GE0P13,GC0P13)
    write(16,90) Name,nuc(iso),A1,A3,S1,GM0P13,GE0P13,GC0P13
    !
    Name='F35(1905) '
    WF35=1.905
    CALL HF35_MAID07(0.D0,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WF35,GM0F35,GE0F35,GC0F35)
    write(17,90) Name,nuc(iso),A1,A3,S1,GM0F35,GE0F35,GC0F35
    !
    Name='P31(1910) '
    WP31=1.910
    CALL HP31_MAID07(0.D0,PhiR,DM,DS,A1,S1,1)
    DE=0
    A3=0
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WP31,GM0P31,GE0P31,GC0P31)
    write(18,90) Name,nuc(iso),A1,A3,S1,GM0P31,GE0P31,GC0P31
    GE0P31=1
    !
    Name='F37(1950) '
    WF37=1.950
    CALL HF37_MAID07(0.D0,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(0.D0,DE,DM,DS,WF37,GM0F37,GE0F37,GC0F37)
    write(19,90) Name,nuc(iso),A1,A3,S1,GM0F37,GE0F37,GC0F37

90  format(2x,A10,A1,3x, &
         'A1/2(0) =',F8.2,',  A3/2(0) =',F8.2,',  S1/2(0) =',F8.3, &
         ', Gm*(0) =',F8.4,',  Ge*(0) =',F8.4,',  Gc*(0) =',F8.4,//, &
         ' Q2(GeV)   A1/2(Q2)    A3/2(Q2)   S1/2(Q2)    Gm*(Q2)   ', &
         '  Ge*(Q2)     Gc*(Q2)     Gm*(norm) Ge*(norm)', &
         ' Gc*(norm)    Gm*/GD  Ge*/GD  Gc*/GD  ' )
    IFST=1
999 CONTINUE

    if (ISO.EQ.1.OR.ISO.EQ.3) write(6,2) Q2G
2   FORMAT(5x,'proton e.m. helicity amplitudes at Q^2=',F7.3,2x, &
         'in units 10^-3/sqrt(GeV)')
    if (ISO.EQ.2.OR.ISO.EQ.4) write(6,3) Q2G
3   FORMAT(5x,'neutron e.m. helicity amplitudes at Q^2=',F7.3,2x, &
         'in units 10^-3/sqrt(GeV)')
    write(6,4)
4   format(17X,'Phi_R',6x,'X1',5x,'A1/2',8x,'X3',5x,'A3/2',8x,'XS', &
         5x,'S1/2')

    egcm = (mdelta**2-Q2G-mp**2)/2./mdelta
    qcm=dsqrt(egcm**2+Q2G)
    GD=1.D0/(1.D0+Q2G/0.71D0)**2

    CALL HP33_MAID07(Q2G,qcm,kgcm0,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WP33,GMP33,GEP33,GCP33)
    write(7,100) Q2G,A1,A3,S1,GMP33,GEP33,GCP33, &
         GMP33/GM0P33,GEP33/GE0P33,GCP33/GC0P33, &
         GMP33/GM0P33/GD,GEP33/GE0P33/GD,GCP33/GC0P33/GD

    CALL HP11_MAID07(ISO,PhiR,Q2G,DM,DS,A1,S1,1)
    A3=0
    DE=0
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WP11,GMP11,GEP11,GCP11)
    write(8,100) Q2G,A1,A3,S1,GMP11,GEP11,GCP11, &
         GMP11/GM0P11,GEP11/GE0P11,GCP11/GC0P11, &
         GMP11/GM0P11/GD,GEP11/GE0P11/GD,GCP11/GC0P11/GD

    CALL HD13_MAID07(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WD13,GMD13,GED13,GCD13)
    write(9,100) Q2G,A1,A3,S1,GMD13,GED13,GCD13, &
         GMD13/GM0D13,GED13/GE0D13,GCD13/GC0D13, &
         GMD13/GM0D13/GD,GED13/GE0D13/GD,GCD13/GC0D13/GD

    CALL HS11f_MAID07(ISO,PhiR,Q2G,DE,DS,A1,S1,1)
    A3=0
    DM=0
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WS11f,GMS11f,GES11f,GCS11f)
    write(10,100) Q2G,A1,A3,S1,GMS11f,GES11f,GCS11f, &
         GMS11f/GM0S11f,GES11f/GE0S11f,GCS11f/GC0S11f, &
         GMS11f/GM0S11f/GD,GES11f/GE0S11f/GD,GCS11f/GC0S11f/GD

    CALL HS11s_MAID07(ISO,PhiR,Q2G,DE,DS,A1,S1,1)
    A3=0
    DM=0
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WS11s,GMS11s,GES11s,GCS11s)
    write(11,100) Q2G,A1,A3,S1,GMS11s,GES11s,GCS11s, &
         GMS11s/GM0S11s,GES11s/GE0S11s,GCS11s/GC0S11s, &
         GMS11s/GM0S11s/GD,GES11s/GE0S11s/GD,GCS11s/GC0S11s/GD

    CALL HS31_MAID07(Q2G,PhiR,DE,DS,A1,S1,1)
    A3=0
    DM=0
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WS31,GMS31,GES31,GCS31)
    write(12,100) Q2G,A1,A3,S1,GMS31,GES31,GCS31, &
         GMS31/GM0S31,GES31/GE0S31,GCS31/GC0S31, &
         GMS31/GM0S31/GD,GES31/GE0S31/GD,GCS31/GC0S31/GD

    CALL HD15_MAID07(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WD15,GMD15,GED15,GCD15)
    write(13,100) Q2G,A1,A3,S1,GMD15,GED15,GCD15, &
         GMD15/GM0D15,GED15/GE0D15,GCD15/GC0D15, &
         GMD15/GM0D15/GD,GED15/GE0D15/GD,GCD15/GC0D15/GD

    CALL HF15_MAID07(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WF15,GMF15,GEF15,GCF15)
    write(14,100) Q2G,A1,A3,S1,GMF15,GEF15,GCF15, &
         GMF15/GM0F15,GEF15/GE0F15,GCF15/GC0F15, &
         GMF15/GM0F15/GD,GEF15/GE0F15/GD,GCF15/GC0F15/GD

    CALL HD33_MAID07(Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WD33,GMD33,GED33,GCD33)
    write(15,100) Q2G,A1,A3,S1,GMD33,GED33,GCD33, &
         GMD33/GM0D33,GED33/GE0D33,GCD33/GC0D33, &
         GMD33/GM0D33/GD,GED33/GE0D33/GD,GCD33/GC0D33/GD

    CALL HP13_MAID07(ISO,PhiR,Q2G,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WP13,GMP13,GEP13,GCP13)
    write(16,100) Q2G,A1,A3,S1,GMP13,GEP13,GCP13, &
         GMP13/GM0P13,GEP13/GE0P13,GCP13/GC0P13, &
         GMP13/GM0P13/GD,GEP13/GE0P13/GD,GCP13/GC0P13/GD

    CALL HF35_MAID07(Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WF35,GMF35,GEF35,GCF35)
    write(17,100) Q2G,A1,A3,S1,GMF35,GEF35,GCF35, &
         GMF35/GM0F35,GEF35/GE0F35,GCF35/GC0F35, &
         GMF35/GM0F35/GD,GEF35/GE0F35/GD,GCF35/GC0F35/GD

    CALL HP31_MAID07(Q2G,PhiR,DM,DS,A1,S1,1)
    A3=0
    DE=0
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WP31,GMP31,GEP31,GCP31)
    write(18,100) Q2G,A1,A3,S1,GMP31,GEP31,GCP31, &
         GMP31/GM0P31,GEP31/GE0P31,GCP31/GC0P31, &
         GMP31/GM0P31/GD,GEP31/GE0P31/GD,GCP31/GC0P31/GD

    CALL HF37_MAID07(Q2G,PhiR,DE,DM,DS,A1,A3,S1,1)
    CALL FORMSTAR_MAID07(Q2G,DE,DM,DS,WF37,GMF37,GEF37,GCF37)
    write(19,100) Q2G,A1,A3,S1,GMF37,GEF37,GCF37, &
         GMF37/GM0F37,GEF37/GE0F37,GCF37/GC0F37, &
         GMF37/GM0F37/GD,GEF37/GE0F37/GD,GCF37/GC0F37/GD


100 FORMAT(F7.4,6E12.4,2x,3F10.6,2x,3F8.4)
    RETURN
  end subroutine HEL_OUT_MAID07

  subroutine FORMSTAR_MAID07(Q2G,DE,DM,DS,XMSTAR,GMSTAR,GESTAR,GCSTAR)
    IMPLICIT REAL*8 (A-H,O-Z)
    REAL*8 kWcm,kgcm
    PI=3.1415926536D0
    XMp=0.9382723D0
    kWcm = (XMstar**2-XMp**2)/2./XMstar
    egcm = (XMstar**2-Q2G-XMp**2)/2./XMstar
    kgcm=dsqrt(egcm**2+Q2G)
    cstar=Sqrt((kWcm*XMp**3)/(4*PI/137.D0*XMstar*kgcm**2))
    ! all signs (zeta phases) are set to +1 here
    GMstar=2*cstar*DM*1.0D-3
    GEstar=2*cstar*DE*1.0D-3
    GCstar=2*cstar*2*XMstar/kgcm*DS*1.0D-3
    RETURN
  end subroutine FORMSTAR_MAID07


  !****************************************************************************
  !****s* helicityAmplitudes/get_helicityAmplitudes
  ! NAME
  ! subroutine get_helicityAmplitudes(targetCharge,resonanceID,QSquared,A1,A3,S1,MAID_version)
  ! PURPOSE
  ! Interface for the MAID helicity amplitudes.
  ! Returns the helicity amplitudes for a given resonance.
  !
  ! INPUTS
  ! * real, intent(in) :: QSquared        ! QSquared in GeV**2
  ! * integer, intent(in) :: targetCharge ! Charge of target nucleon
  ! * integer, intent(in) :: resonanceID  ! ID of the resonance according to GiBUU
  ! * integer,optional, intent(in) :: MAID_version ! =2003,=2005,=2007 (default)
  ! OUTPUT
  ! * real, intent(out):: A1,A3,S1 ! In units of GeV^(-1/2)
  !
  ! NOTES
  ! * In this routine we are converting from real to real(8) since the
  !   MAID routines are made for real(8) input and output.
  !****************************************************************************
  subroutine get_helicityAmplitudes(targetCharge,resonanceID,QSquared_IN,A1_OUT,A3_OUT,S1_OUT,MAID_version_in)
    use IDTABLE
    use callstack, only: traceback
    implicit none

    real, intent(in) :: QSquared_IN
    integer, intent(in) :: targetCharge
    integer, intent(in) :: resonanceID
    real, intent(out):: A1_OUT,A3_OUT,S1_OUT
    integer,optional, intent(in) :: MAID_version_in

    real(8) :: A1,A3,S1
    real(8) :: QSquared
    integer :: iso
    integer, parameter :: printInfos = 0     ! 0= no printout, 1= printout
    real(8) :: egcm,qcm,kgcm0,DE,DM,DS,phiR
    real(8), parameter :: W0=1.232
    real(8), parameter :: mi=0.9382723

    integer :: MAID_version ! 2003,2005,2007

    MAID_version = 2007
    if (present(MAID_version_in)) then
       select case (MAID_version_in)
       case (2003,2005,2007)
          MAID_version = MAID_version_in
       case default
          call traceback("wrong MAID version")
       end select
    end if


    !Checks
    if (.not.(targetCharge.eq.0.or.targetCharge.eq.1)) then
       write(*,*) 'Wrong nucleon charge', targetCharge
       call traceback()
    end if

    QSquared=QSquared_In

    A1=0.
    A3=0.
    S1=0.

    iso=2-targetCharge

    select case (MAID_version)
    case (2003)

       select case (resonanceID)
       case (Delta)
          kgcm0 = (W0*W0-mi*mi)/2./W0
          egcm = (W0*W0-QSquared-mi*mi)/2./W0
          qcm=sqrt(egcm**2+QSquared)
          CALL HP33_MAID03( QSquared,qcm,kgcm0,DE,DM,DS,A1,A3,S1,printInfos)
       case (P11_1440)
          CALL HP11_MAID03( ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D13_1520)
          CALL HD13_MAID03( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (S11_1535)
          CALL HS11f_MAID03(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (S31_1620)
          CALL HS31_MAID03( QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (S11_1650)
          CALL HS11s_MAID03(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D15_1675)
          CALL HD15_MAID03( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F15_1680)
          CALL HF15_MAID03( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (D33_1700)
          CALL HD33_MAID03( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P13_1720)
          CALL HP13_MAID03( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F35_1905)
          CALL HF35_MAID03( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P31_1910)
          CALL HP31_MAID03( QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (F37_1950)
          CALL HF37_MAID03( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       end select

    case (2005)

       select case (resonanceID)
       case (Delta)
          kgcm0 = (W0*W0-mi*mi)/2./W0
          egcm = (W0*W0-QSquared-mi*mi)/2./W0
          qcm=sqrt(egcm**2+QSquared)
          CALL HP33_MAID05( QSquared,qcm,kgcm0,DE,DM,DS,A1,A3,S1,printInfos)
       case (P11_1440)
          CALL HP11_MAID05( ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D13_1520)
          CALL HD13_MAID05( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (S11_1535)
          CALL HS11f_MAID05( ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (S31_1620)
          CALL HS31_MAID05( QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (S11_1650)
          CALL HS11s_MAID05(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D15_1675)
          CALL HD15_MAID05( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F15_1680)
          CALL HF15_MAID05( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (D33_1700)
          CALL HD33_MAID05( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P13_1720)
          CALL HP13_MAID05( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F35_1905)
          CALL HF35_MAID05( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P31_1910)
          CALL HP31_MAID05( QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (F37_1950)
          CALL HF37_MAID05( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       end select

    case (2007)

       select case (resonanceID)
       case (Delta)
          kgcm0 = (W0*W0-mi*mi)/2./W0
          egcm = (W0*W0-QSquared-mi*mi)/2./W0
          qcm=sqrt(egcm**2+QSquared)
          CALL HP33_MAID07( QSquared,qcm,kgcm0,DE,DM,DS,A1,A3,S1,printInfos)
       case (P11_1440)
          CALL HP11_MAID07( ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D13_1520)
          CALL HD13_MAID07( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (S11_1535)
          CALL HS11f_MAID07( ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (S31_1620)
          CALL HS31_MAID07( QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (S11_1650)
          CALL HS11s_MAID07(ISO,PhiR,QSquared,DM,DS,A1,S1,printInfos)
       case (D15_1675)
          CALL HD15_MAID07( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F15_1680)
          CALL HF15_MAID07( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (D33_1700)
          CALL HD33_MAID07( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P13_1720)
          CALL HP13_MAID07( ISO,PhiR,QSquared,DE,DM,DS,A1,A3,S1,printInfos)
       case (F35_1905)
          CALL HF35_MAID07( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       case (P31_1910)
          CALL HP31_MAID07( QSquared,PhiR,DM,DS,A1,S1,printInfos)
       case (F37_1950)
          CALL HF37_MAID07( QSquared,PhiR,DE,DM,DS,A1,A3,S1,printInfos)
       end select

    end select


    ! Converting to GeV^(-1/2)
    A1_OUT=A1/1000.
    A3_OUT=A3/1000.
    S1_OUT=S1/1000.

    !we use different sign in S_1/2 amplitude, thus:
    S1_OUT=-S1_OUT


  end subroutine get_helicityAmplitudes
end module helicityAmplitudes
