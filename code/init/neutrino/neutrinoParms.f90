module neutrinoParms

  use CALLSTACK, only: TRACEBACK

  implicit none
  private

  !****************************************************************************
  !****g* neutrinoParms/twopiBG_on
  ! SOURCE
  real, save, public :: twopiBG_on = 1.5 ! GeV
  ! PURPOSE
  ! W for transition from 1piBG to 2piBG for BC parametrization of nr X-section
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/Wtrans
  ! SOURCE
  real, save, public :: Wtrans = 2.7 ! GeV
  ! PURPOSE
  ! W for transition from Bosted-Christy Parametrization to PYTHIA DIS
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/NormpiBG
  ! SOURCE
  real, save, public :: NormpiBG = 1.0
  ! PURPOSE
  ! overall normalization factor for pi BG and Bloom-Gilman X-section,
  ! only relevant for neutrinos
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/normRes
  ! SOURCE0
  ! PURPOSE
  real, save, public :: NormRes = 1.0
  ! overall normalization factor for neutrino-induced resonance contributions
  ! beyond the Delta
  ! only relevant for neutrinos
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/normBC
  ! SOURCE
  real, save, public :: normBC = 1.0
  ! PURPOSE
  ! overall normalization factor for neutrino-induced Christy-Bosted
  ! contributions between 2 GeV and DIS onset
  ! only relevant for neutrinos
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/Q2cut
  ! SOURCE
  real, save, public ::  Q2cut = 4.5 ! GeV^2
  ! PURPOSE
  ! threshold value
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/A00
  ! SOURCE
  real, save, public ::  A00 = -0.43
  real, save, public ::  A10 = -0.97
  real, save, public ::  A01 = 0.87
  real, save, public ::  A20 = 8.01
  real, save, public ::  A11 = 4.59
  real, save, public ::  A02 = 3.18
  real, save, public ::  A30 = -10.42
  real, save, public ::  A21 = 2.04
  real, save, public ::  A12 = -9.14
  real, save, public ::  A03 = -15.15

  real, save, public ::  W00 = 1.55
  real, save, public ::  x00 = 0.35
  ! PURPOSE
  ! coefficients in Taylor expansion of attenuation function Att(W,x)
  !****************************************************************************

  ! Now parameters used pimarily in 2p2h:

  !****************************************************************************
  !****g* neutrinoParms/ME_Version
  ! PURPOSE
  ! indicate the type of matrix element parametrisation
  !
  ! SOURCE
  integer, save, public :: ME_Version = 6
  ! NOTES
  ! possible values:
  ! * 1: const ME_Norm_XX  ! const for CC  fitted to MiniBooNE is 1.8e-6
  ! * 2: constant transverse and decreasing with Enu
  ! * 3: "Dipole transverse" transverse, fall with Q2 as 4-th power
  ! * 4: MEC from E. Christy (8/2015), with parametrization for longitudinal
  ! * 5: MEC from Bosted arXiV:1203.2262, with parametrization for longitudinal
  ! * 6: MEC additional parametrization, with parametrization for longitudinal
  !   not yet implemented
  !
  ! remarks:
  ! * case 1 is model-I in Lalakulich,Gallmeister,Mosel PRC86(2012)014614
  ! * case 2 is model-II from Lalakulich,Gallmeister,Mosel PRC86(2012)014614
  ! * case 3 gives a good description of MiniBooNE data with MA ~ 1.5 GeV
  !****************************************************************************


  ! The following are all tunable strength parameters for 2p2h hadronic
  ! structure functions. Default is no tuning, i.e. all parameters = 1.
  ! except for ME_Long, for which default is =0 (no longitudinal component)



  !****************************************************************************
  !****g* neutrinoParms/ME_Norm_QE
  ! PURPOSE
  ! Overall strength of 2p2h matrix element with 2N out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_Norm_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 gives the coded strength
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_Norm_Delta
  ! PURPOSE
  ! Overall strength of 2p2h matrix element with NDelta out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_Norm_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_Mass_QE
  ! PURPOSE
  ! Cutoff-mass in some parametrizations of 2p2h matrix element for NN out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_Mass_QE    = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_Mass_Delta
  ! PURPOSE
  ! Cutoff-mass in some parametrizations of matrix element for NDelta out
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_Mass_Delta = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value == 1 is a dummy value
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_Transversity
  ! PURPOSE
  ! Parametrisation of structure functions
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_Transversity = (/1.0, 1.0, 1.0/)
  ! NOTES
  ! The value = 1 chooses structure function W2 so that 2p2h is pure transverse
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_LONG
  ! PURPOSE
  ! Parametrization of structure functions
  !
  ! for (EM,CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_LONG = (/0.0, 0.0, 0.0/)
  ! NOTES
  ! The value = 0 turns any additional longitudinal contribution
  ! to structure funct. W2 off
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_W3
  ! PURPOSE
  ! Overall strength factor for structure function W3
  !
  ! only for (CC,NC)
  ! SOURCE
  real,dimension(1:3), save, public :: ME_W3 = (/0.0, 1.0, 1.0/)
  ! NOTES
  ! overall strength parameter for structure function W3
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/ME_ODW
  ! PURPOSE
  ! switch for choosing the connection between structure functions
  ! W1(electron) and W1(neutrino) and W3(neutrino):
  ! * 1: for expressions from Martini et al.
  ! * 2: for expressions from O'Connell et al.
  !
  ! only for (CC,NC)
  ! SOURCE
  integer, save, public :: ME_ODW = 2
  ! NOTES
  ! * O'Connell et al: PR C6 (1972) 719
  ! * Martini et al: PR C80 (2009) 065501
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/inmedW
  ! PURPOSE
  ! Controls which inv mass W is used in parametrization of 2p2h W1
  !
  ! SOURCE
  integer, save, public :: inmedW = 1
  ! NOTES
  ! * 1: W = static inv. mass in 2p2h parametrization of W1
  ! * 2: W = inv mass for Fermi moving nucleons in potential
  ! * 3: W = inv mass for Fermi moving nucleons without potential
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/T
  ! PURPOSE
  ! target isospin, affects only neutrino 2p2h structure function
  !
  ! SOURCE
  real, save, public ::  T = 1
  ! NOTES
  ! * T = 0, 1 , ...
  ! * T = 99 gives T = abs((N-Z)/2)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/Adep
  ! PURPOSE
  ! Switch for A-dependence of 2p2h structure function
  !
  ! SOURCE
  integer, save, public ::  Adep = 2
  ! NOTES
  ! * 1: A-dependence for zero-range force (Mosel, Gallmeister, 2016)
  ! * 2: linear A-dependence, normalized to C12
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/new_eN
  ! SOURCE
  logical, save, public :: new_eN = .true.
  ! PURPOSE
  ! switch to choose new electron-Nucleon X-sections
  ! with Bosted-Christy parametrization
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoParms/new_eNres
  ! SOURCE
  logical, save, public :: new_eNres = .true.
  ! PURPOSE
  ! switch to choose new electron-Nucleon X-sections for resonances
  ! with Bosted-Christy parametrization
  !
  ! NOTES
  ! * only used if new_eN = .T.
  ! * currently, .T. should be used for electrons and .F. for neutrinos
  ! * .T. and .F. can be set as the default value using 'setDefault_Res'
  !   giving the process_ID (=EM,CC,NC and anti) as argument
  ! * any value given in the jobcard will override this default value
  !****************************************************************************

  public :: readInput_neutrino
  public :: readInput_2p2h
  public :: setDefault_Res
  public :: VAfact

  integer, save :: store_process_ID = -99

contains

  !****************************************************************************
  !****s* neutrinoParms/readInput_neutrino
  ! NAME
  ! subroutine readInput_neutrino
  !****************************************************************************
  subroutine readInput_neutrino
    use output, only: Write_InitStatus,Write_ReadingInput
    use AZN

    integer :: ios
    logical, save:: initflag = .true.

    !************************************************************************
    !****n* neutrinoParms/nl_Neutrino2piBack
    ! NAME
    ! NAMELIST /nl_Neutrino2piBack/
    ! PURPOSE
    ! Parameters for transition from RES to DIS:
    ! * Wtrans
    ! * NormpiBG
    ! * normRES
    ! * normBC
    ! Includes parameters of neutrino 2 pion background:
    ! * Q2cut
    ! * A00
    ! * A10
    ! * A01
    ! * A02
    ! * A11
    ! * A20
    ! * A03
    ! * A12
    ! * A21
    ! * A30
    ! * W00
    ! * x00
    ! * twopiBG_on
    !************************************************************************
    NAMELIST /nl_neutrino2piBack/ Wtrans, &
         NormpiBG,normRES,normBC,Q2cut, &
         A00,A10,A01,A02,A11,A20,A03,A12,A21,A30,W00,x00,twopiBG_on

    !************************************************************************
    !****n* neutrinoParms/nl_neweN
    ! NAME
    ! NAMELIST /nl_neweN/
    ! PURPOSE
    ! Parameters for new eN implementation
    ! * new_eN
    ! * new_eNres
    ! * ME_ODW
    ! * T
    !************************************************************************
    NAMELIST /nl_neweN/ &
         new_eN,new_eNres,ME_ODW,T
    !************************************************************************

    if (.not.initFlag) return

    call Write_InitStatus('neutrinoParms',0)

    call Write_ReadingInput('nl_neweN',0)
    rewind(5)
    read(5,nml=nl_neweN,IOSTAT=ios)
    call Write_ReadingInput('nl_neweN',0,ios)
    write(*,*) 'new_eN   =',new_eN
    write(*,*) 'new_eNres=',new_eNres
    call Write_ReadingInput('nl_neweN',1)

    call Write_ReadingInput('nl_neutrino2piBack',0)
    rewind(5)
    read(5,nml=nl_neutrino2piBack,IOSTAT=ios)
    call Write_ReadingInput('nl_neutrino2piBack',0,ios)

    if (new_eN .and. abs(store_process_ID) == 1) then
       if (abs(NormRes*NormpiBG*NormBC - 1.0) > 1e-5) then

          write(*,*) 'Attention: if new_eN and EM process, Norm has to be 1.0!'
          write(*,*) 'Values in jobcard overwritten by defaults!'

          NormRES  = 1.0
          NormpiBG = 1.0
          NormBC   = 1.0
       end if
    end if

    call AZNsub(Atarget,Ntarget,Ztarget)
    if(T == 99) T = abs((Ntarget - Ztarget)/2.0)
    write(*,*) 'ME_ODW=', ME_ODW
    write(*,*) 'T2pi=', T
    write(*,*) 'Wtrans  =',Wtrans
    write(*,*) 'NormpiBG=',NormpiBG
    write(*,*) 'NormRES =',NormRES
    write(*,*) 'NormBC  =',NormBC
    write(*,fmt='(12(F7.2))') Q2cut,A00,A10,A01,A20,A11,A02,A30,A21,A12,A03,T
    call Write_ReadingInput('nl_neutrino2piBack',1)

    call Write_InitStatus('neutrinoParms',1)

    initFlag=.false.
  end subroutine readInput_neutrino

  !****************************************************************************
  !****s* neutrinoParms/readInput_2p2h
  ! NAME
  ! subroutine readInput_2p2h
  !****************************************************************************
  subroutine readInput_2p2h

    use output, only: Write_ReadingInput
    use AZN

    integer :: ios, i

    logical, save:: initflag = .true.

    !**************************************************************************
    !****n* neutrinoParms/Lepton2p2h
    ! NAME
    ! NAMELIST /Lepton2p2h/
    ! PURPOSE
    ! Includes parameters for 2p2h events:
    ! * ME_Version
    ! * ME_Norm_QE
    ! * ME_Norm_Delta
    ! * ME_Mass_QE
    ! * ME_Mass_Delta
    ! * ME_Transversity
    ! * ME_LONG
    ! * ME_W3
    ! * inmedW
    ! * Adep
    !**************************************************************************
    NAMELIST /lepton2p2h/ ME_Version, &
         ME_Norm_QE, ME_Norm_Delta, &
         ME_Mass_QE, ME_Mass_Delta,ME_Transversity,ME_LONG, &
         ME_W3,inmedW,Adep

    if (.not.initFlag) return

    call Write_ReadingInput('lepton2p2h',0)
    rewind(5)
    read(5,nml=lepton2p2h,IOSTAT=ios)
    call Write_ReadingInput("lepton2p2h",0,ios)

    select case (ME_Version)     ! for 2p2h

    case (1)
       write(*,'(A)') 'ME1  const =  4.0e-6 * ME_Norm_XX'

    case (2)
       write(*,'(A)') 'ME2  4.8e4 * 0.635 / Enu^2 * ME_Norm_XX  &
            &  in transverse part only, decreasing with Enu'

    case (3)
       write(*,'(A)') 'ME3  dipole parameterization of the "form factor" &
            & ME=8e4*(1+Q2/MA2)^{-4},transverse part only'

    case (4)
       write(*,'(A)') 'ME 4, parametrization of structure functions W1,W2,W3,&
            & W1 from E. Christy'
       !       write(*,*) 'ME_W3   =', ME_W3

    case (5)
       write(*,'(A)') 'parametrization of structure functions W1,W2,W3,&
            & W1 from Bosted'

    case (6)
       write(*,'(A)') 'parametrization of structure functions W1,W2,W3,&
            & W1 new parametrization'

    case default
       write(*,*) 'ME_Version = ',ME_Version
       call TRACEBACK('wrong value for ME_Version')

    end select

   
    write(*,*) 'inmedW=',inmedW,'   T=',T
    write(*,*) 'Adep=',Adep

    write(*,'(A)') 'parameters ( N N final state ):  [i=EM,CC,NC]'
    do i=1,3
       write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
            & ME_Norm_QE(i),ME_Mass_QE(i),ME_Transversity(i),ME_LONG(i)
    end do
    write(*,'(A)') 'parameters ( N Delta final state ):'
    do i=1,3
       write(*,'("   A=",ES13.5,"  M=",ES13.5)') &
            & ME_Norm_Delta(i),ME_Mass_Delta(i)
    end do

    call Write_ReadingInput('lepton2p2h',1)

    initFlag=.false.

  end subroutine readInput_2p2h

  !****************************************************************************
  !****s* neutrinoParms/setDefault_Res
  ! NAME
  ! subroutine setDefault_Res(process_ID)
  ! PURPOSE
  ! select the default resonance treatment for 'new_eN' according the used
  ! process.
  !****************************************************************************
  subroutine setDefault_Res(process_ID)
    use leptonicID

    integer, intent(in) :: process_ID

    select case(process_ID)
    case (EM,antiEM)
       new_eNres = .true.
    case default
       new_eNres = .false.
    end select

    store_process_ID = process_ID

  end subroutine setDefault_Res
  
  
  subroutine VAfact(ME_ODW,qvecsq,nu,kinfact)
  use constants, only: mN
  integer, intent (in) :: ME_ODW
  real, intent (in) :: qvecsq,nu
  real, intent (out) :: kinfact
  
   select case (ME_ODW)
    case (1)       ! Martini et al
       kinfact = nu**2/qvecsq
    case (2)        ! O'Connell et al.
       kinfact = qvecsq/(2*mN)**2
    case (3) 
       kinfact = (qvecsq - nu**2)/(2*mN)**2
    case default
       write(*,*) 'kinfact error in RT'
    end select
    
   end subroutine VAfact

end module neutrinoParms
