!******************************************************************************
!****p* /HiLeptonMixRuns
! NAME
! program HiLeptonMixRuns
! PURPOSE
! Read histograms from different files to accumulate better statistics.
!
! This program reads histograms created by the HiLepton analysis
! from different directories, adds them up,
! and prints them out again in the calling directory.
!
! The directories to be walked through are read in from the file
! 'directories.txt', which in every line gives on single directory
! (trailing slashes are not required). If the code does not find this
! file, it will run over the directories '0/', '1/' ... '9/'. Thus it
! resembles the bahaviour of the old 'HiLeptonMix10Runs' program.
!
! You could get the same default behavior by generating the input file via
!    ls -d ? > directories.txt
!
!******************************************************************************
program HiLeptonMixRuns

  use hist
  use hist2D
  use histMP
  use histMC
  use output

  implicit none

  type(histogram)   :: HHH,hLep_Q2,hLep_nu
  type(histogramMP) :: hMP
  type(histogram2D) :: h2D
  type(histogramMC) :: hMC
  type(histogram)   :: HHH1,HHH2,HHH3

  real, save :: nLeptons
  character(3), dimension(-1:1), parameter :: piName  = (/'pi-','pi0','pi+'/)
  character(2), dimension(0:1),  parameter :: nucName = (/'N0','N+'/)
  integer :: i, iNu,iQ2,iZH, ii
  logical :: lDum, isFirst, isOK
  character(500) :: Buf

  character(len=:), allocatable :: Dirs(:)
  integer :: nDirs = 0

  ! ===== Build directory list =====

  call BuildDirs
  if (nDirs==0) then
     write(*,*) 'List of directories is empty. STOP!'
     stop
  end if

  ! ===== Lepton Kinematics =====


  call FetchHist_multi(hLep_Q2,"HiLep.lep.Q2.kinematics.dat")
  call WriteHist(hLep_Q2,101,add=1e-20,file="HiLep.lep.Q2.kinematics.dat",dump=.true.)
  call WriteHist(hLep_Q2,101,add=1e-20,DoAve=.true.,&
       &file="HiLep.lep.nuAve_Q2.kinematics.dat")

  if (.not.hLep_Q2%initialized) then
     write(*,*) 'ERROR: no data found. STOP!'
     stop
  end if

  nLeptons = SUM(hLep_Q2%yVal(:,1),dim=1)
  write(*,*) nLeptons

  call FetchHist_multi(hLep_nu,"HiLep.lep.nu.kinematics.dat")
  call WriteHist(hLep_nu,101,add=1e-20,file="HiLep.lep.nu.kinematics.dat")
  call WriteHist(hLep_nu,101,add=1e-20,DoAve=.true.,&
       &file="HiLep.lep.Q2Ave_nu.kinematics.dat")

  ! to be continued...


  ! ===== nu-spectra =====

  call FetchHistMP_multi(hMP,"HiLep.nu.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,H2=hLep_nu,iColumn=3,&
       &file='HiLep.nu.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,H2=hLep_nu,iColumn=1,&
       &file='HiLep.nu.idH.Acc.dat') ! Acc

  call FetchHistMP_multi(hMP,"HiLep.AvePT2.nu.dat")
  call WriteHistMP(hMP,101,DoAve=.true.,&
       &file='HiLep.AvePT2.nu.dat')


  ! to be continued...


  ! ===== Q2-spectra =====

  call FetchHistMP_multi(hMP,"HiLep.Q2.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,H2=hLep_Q2,iColumn=3,&
       &file='HiLep.Q2.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,H2=hLep_Q2,iColumn=1,&
       &file='HiLep.Q2.idH.Acc.dat') ! Acc

  call FetchHistMP_multi(hMP,"HiLep.AvePT2.Q2.dat")
  call WriteHistMP(hMP,101,DoAve=.true.,&
       &file='HiLep.AvePT2.Q2.dat')


  ! to be continued...

  ! ===== DoClassifyFirst =====

  if (TryFile("HiLep.zHclass.000.dat.bin")) then

     do i=-1,10 ! loop maybe larger than necessary
        ii = i
        if (i==-1) ii=999

        call FetchHistMP_multi(hMP,"HiLep.zHclass."//trim(intToChar(ii))//".dat",isOK)
        if (isOK) then
           call WriteHistMP(hMP,add=1e-20,mul=1./NLeptons,iColumn=1,&
                file='HiLep.zHclass.'//trim(intToChar(ii))//'.dat')
        end if

        call FetchHistMP_multi(hMP,"HiLep.pT2class."//trim(intToChar(ii))//".dat",isOK)
        if (isOK) then
           call WriteHistMP(hMP,add=1e-20,DoAve=.true.,&
                file='HiLep.pT2class.'//trim(intToChar(ii))//'.dat')
        end if

        call FetchHistMP_multi(hMP,"HiLep.zH_Generation."//trim(intToChar(ii))//".dat",isOK)
        if (isOK) then
           call WriteHistMP(hMP,add=1e-20,mul=1./NLeptons,iColumn=1,&
                file='HiLep.zH_Generation.'//trim(intToChar(ii))//'.dat')
        end if

        call FetchHistMP_multi(hMP,"HiLep.pT2_Generation."//trim(intToChar(ii))//".dat",isOK)
        if (isOK) then
           call WriteHistMP(hMP,add=1e-20,mul=1./NLeptons,iColumn=1,&
                file='HiLep.pT2_Generation.'//trim(intToChar(ii))//'.dat')
        end if

     end do
  end if

  ! ===== zH-spectra =====

  call FetchHistMP_multi(hMP,"HiLep.zH.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,mul=1./NLeptons,iColumn=3,&
       &file='HiLep.zH.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,mul=1./NLeptons,iColumn=1,&
       &file='HiLep.zH.idH.Acc.dat') ! Acc

  call FetchHistMP_multi(hMP,"HiLep.AvePT2.zH.dat")
  call WriteHistMP(hMP,101,DoAve=.true.,&
       &file='HiLep.AvePT2.zH.dat')


  ! to be continued...

  ! ===== pT2-spectra =====

  call FetchHistMP_multi(hMP,"HiLep.pT2.idH.Acc.dat")
  call WriteHistMP(hMP,101,add=1e-20,mul=1./NLeptons,iColumn=3,&
       &file='HiLep.pT2.idH.noAcc.dat') ! noAcc
  call WriteHistMP(hMP,102,add=1e-20,mul=1./NLeptons,iColumn=1,&
       &file='HiLep.pT2.idH.Acc.dat') ! Acc

  ! =====================

  do i=-1,1
     call FetchHist2D_multi(h2D,'HiLep.pT2Ave_Plane.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140,DoAve=.true.,maxval=0.0,&
          &file='HiLep.pT2Ave_Plane.'//piName(i)//'.dat')

     call FetchHist2D_multi(h2D,'HiLep.pT2_zH.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          &file='HiLep.pT2_zH.'//piName(i)//'.dat')

     call FetchHist2D_multi(h2D,'HiLep.pT2_nu.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140, H2=hLep_nu, add=1e-20,&
          &file='HiLep.pT2_nu.'//piName(i)//'.dat')

     call FetchHist2D_multi(h2D,'HiLep.pT2_Q2.'//piName(i)//'.dat')
     call WriteHist2D_Gnuplot(h2D,140, H2=hLep_Q2, add=1e-20,&
          &file='HiLep.pT2_Q2.'//piName(i)//'.dat')
  end do

  ! =====================
  ! Binning_5:
  ! =====================

  if (TryFile("HiLep.JLAB5_R.nu_1_1.dat.bin")) then

     ! the upper limits of 5 are dummy. TryFile(...) checks for existance
     do iZH=1,5
        isFirst = .true.
        do iQ2=1,5
           Buf = 'HiLep.JLAB5_R.nu_'//Achar(iQ2+48)//'_'//Achar(iZH+48)//'.dat'
           if (TryFile(Buf)) then
              call FetchHist_multi(HHH,Buf)
              call WriteHist(HHH,add=1e-20,mul=1./nLeptons,file=trim(Buf))
              if (isFirst) then
                 call CopyHist(HHH1, HHH)
                 isFirst = .false.
              else
                 call SumHist(HHH1, HHH)
              end if
           end if
        end do
        if (.not.isFirst) then
           call WriteHist(HHH1,add=1e-20,mul=1./nLeptons,&
                file='HiLep.JLAB5_R.nu_all_'//Achar(iZH+48)//'.dat')
        end if
     end do

     do iNu=1,5
        do iZH=1,5
           Buf = 'HiLep.JLAB5_R.Q2_'//Achar(iNu+48)//'_'//Achar(iZH+48)//'.dat'
           if (TryFile(Buf)) then
              call FetchHist_multi(HHH,Buf)
              call WriteHist(HHH,add=1e-20,mul=1./nLeptons,file=trim(Buf))
           end if
        end do
     end do

     do iNu=1,5
        isFirst = .true.
        do iQ2=1,5
           Buf = 'HiLep.JLAB5_R.zH_'//Achar(iQ2+48)//'_'//Achar(iNu+48)//'.dat'
           if (TryFile(Buf)) then
              call FetchHist_multi(HHH,Buf)
              call WriteHist(HHH,add=1e-20,mul=1./nLeptons,file=trim(Buf))
              if (isFirst) then
                 call CopyHist(HHH1, HHH)
                 isFirst = .false.
              else
                 call SumHist(HHH1, HHH)
              end if
           end if
        end do
        if (.not.isFirst) then
           call WriteHist(HHH1,add=1e-20,mul=1./nLeptons,&
                file='HiLep.JLAB5_R.zH_all_'//Achar(iNu+48)//'.dat')
        end if
     end do

     do iNu=1,5
        do iQ2=1,5
           do iZH=1,10
              Buf = 'HiLep.JLAB5_R.pT2_'//Achar(iNu+48)//'_'//Achar(iQ2+48)//'_'//Achar(iZH+48)//'.dat'
              if (TryFile(Buf)) then
                 call FetchHist_multi(HHH,Buf)
                 call WriteHist(HHH,add=1e-20,mul=1./nLeptons,file=trim(Buf))
              end if
           end do
        end do
     end do

  end if


  ! =====================
  ! Binning_11:
  ! =====================

  if (TryFile("HiLep.dN_id.nu.zH.001.dat.bin")) then
     do i=1,3
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.nu.zH.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.nu.zH.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.Q2.zH.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.Q2.zH.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.pT2.zH.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.pT2.zH.'//trim(intToChar(i))//'.dat') ! Acc

        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.zH.nu.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.zH.nu.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.Q2.nu.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.Q2.nu.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.pT2.nu.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.pT2.nu.'//trim(intToChar(i))//'.dat') ! Acc
     end do
     do i=1,2
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.nu.pT2.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.nu.pT2.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.zH.pT2.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.zH.pT2.'//trim(intToChar(i))//'.dat') ! Acc
        call FetchHistMP_multi(hMP,&
             &file='HiLep.dN_id.Q2.pT2.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP,141,add=1e-20,mul=1./NLeptons,iColumn=1,&
             &file='HiLep.dN_id.Q2.pT2.'//trim(intToChar(i))//'.dat') ! Acc
     end do

  end if

!!$  call FetchHist2D_multi(h2D,'HiLep.CollHistPT2.1.dat')
!!$  call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
!!$       & file='HiLep.CollHistPT2.1.dat') ! Acc
!!$
!!$  call FetchHist2D_multi(h2D,'HiLep.CollHistPT2.2.dat')
!!$  call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
!!$       & file='HiLep.CollHistPT2.2.dat') ! Acc

  if (TryFile("HiLep.nleadPT.N.000.dat.bin")) then
     do i=0,3
        call FetchHistMP_multi(hMP,'HiLep.nleadPT.N.'//trim(intToChar(i))//'.dat')
        call WriteHistMP(hMP, 56, add=1e-20,mul=1./NLeptons,iColumn=1,&
             & file='HiLep.nleadPT.N.'//trim(intToChar(i))//'.dat') ! noAcc
        call WriteHistMP(hMP, 56 ,mul=1./NLeptons,DoAve=.true.,&
             & file='HiLep.nleadPT.PT2.'//trim(intToChar(i))//'.dat') ! noAcc
     end do
  end if

  if (TryFile("HiLep.Rho0_NuQ2.dat.bin")) then
     call FetchHist2D_multi(h2D,file='HiLep.Rho0_NuQ2.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_NuQ2.dat',dump=.true.) ! Acc
     call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_NuQ2.AveMis.dat',dump=.false.) ! Acc

     call IntegrateHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
          & file="HiLep.Rho0_NuQ2.Int2.dat",dump=.true.)
     call AverageHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,&
          & file="HiLep.Rho0_NuQ2.Ave2.dat",dump=.true.)

  end if

  do i=0,4
     if (TryFile("HiLep.Rho0_NuQ2_proc"//trim(intToChar(i))//".dat.bin")) then
        call FetchHist2D_multi(h2D,file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.dat')
        call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
        call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.AveMis.dat',dump=.false.) ! Acc

        call IntegrateHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.Int2.dat',dump=.true.) ! Acc
        call AverageHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,&
             & file='HiLep.Rho0_NuQ2_proc'//trim(intToChar(i))//'.Ave2.dat',dump=.true.) ! Acc
     end if
  end do

  if (TryFile("HiLep.Rho0_WQ2.dat.bin")) then
     call FetchHist2D_multi(h2D,file='HiLep.Rho0_WQ2.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_WQ2.dat',dump=.true.) ! Acc
     call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_WQ2.AveMis.dat',dump=.false.) ! Acc

     call IntegrateHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
          & file="HiLep.Rho0_WQ2.Int2.dat",dump=.true.)
     call AverageHist2D(h2D,HHH,2)
     call WriteHist(HHH,140,&
          & file="HiLep.Rho0_WQ2.Ave2.dat",dump=.true.)

  end if

  do i=0,4
     if (TryFile("HiLep.Rho0_WQ2_proc"//trim(intToChar(i))//".dat.bin")) then
        call FetchHist2D_multi(h2D,file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.dat')
        call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.dat',dump=.true.) ! Acc
        call WriteHist2D_Gnuplot(h2D,140, DoAve=.true., mul=1./nLeptons, add=1e-20,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.AveMis.dat',dump=.false.) ! Acc

        call IntegrateHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,add=1e-20,mul=1./nLeptons,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Int2.dat',dump=.true.) ! Acc
        call AverageHist2D(h2D,HHH,2)
        call WriteHist(HHH,140,&
             & file='HiLep.Rho0_WQ2_proc'//trim(intToChar(i))//'.Ave2.dat',dump=.true.) ! Acc
     end if
  end do


  if (TryFile("HiLep.Rho0_NuQ2_VMD.dat.bin")) then
     call FetchHist2D_multi(h2D,file='HiLep.Rho0_NuQ2_VMD.dat')
     call WriteHist2D_Gnuplot(h2D,140, mul=1./nLeptons, add=1e-20,&
          & file='HiLep.Rho0_NuQ2_VMD.dat',dump=.true.) ! Acc
  end if

  if (TryFile("HiLep.Rho0_MV.dat.bin")) then
     call FetchHist_multi(hhh,file='HiLep.Rho0_MV.dat')
     call WriteHist(hhh,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.Rho0_MV.dat',dump=.true.)

     call FetchHist_multi(hhh,file='HiLep.Rho0_MX.dat')
     call WriteHist(hhh,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.Rho0_MX.dat',dump=.true.)

     call FetchHist_multi(hhh,file='HiLep.Rho0_DE.dat')
     call WriteHist(hhh,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.Rho0_DE.dat',dump=.true.)
  end if

  if (TryFile("HiLep.Shadowing.dat.bin")) then
     call FetchHist2D_multi(h2D,file='HiLep.Shadowing.dat')
     call WriteHist2D_Gnuplot(h2D,140, DoAve=.true.,MaxVal=-1.0,&
          & file='HiLep.Shadowing.dat',dump=.true.) ! Acc
  end if

  if (TryFile("HiLep.nuQ2.Flux.dat")) then
     call FetchHist2D_multi(h2D,file='HiLep.nuQ2.Flux.dat')
     call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.nuQ2.Flux.dat',dump=.true.)
     call FetchHist2D_multi(h2D,file='HiLep.nuQ2.FluxW.dat')
     call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
          & file='HiLep.nuQ2.FluxW.dat',dump=.true.)
  end if

!!$  do i=0,4
!!$     if (TryFile("0/HiLep.nuQ2."//intTochar(i)//".dat")) then
!!$        call FetchHist2D_multi(h2D,file='HiLep.nuQ2.'//intTochar(i)//'.dat')
!!$        call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
!!$             & file='HiLep.nuQ2.'//intTochar(i)//'.dat',dump=.true.)
!!$     end if
!!$  end do

  if (TryFile("HiLep.Brooks.dat")) then
     call FetchHistMC_multi(hMC,file='HiLep.Brooks.dat')
     call WriteHistMC(hMC,'HiLep.Brooks.dat',add=1e-20,mul=1./nLeptons)
  end if

  if (TryFile("HiLep.NuQ2planeXS.SYS.dat")) then
     do i=0,9
        call FetchHist2D_multi(h2D,file='HiLep.NuQ2planeXS.'//trim(intToChar(i))//'.dat')
        call WriteHist2D_Gnuplot(h2D,140,add=1e-20,mul=1./nLeptons,&
             & file='HiLep.NuQ2planeXS.'//trim(intToChar(i))//'.dat',dump=.true.)
     end do
  end if

  if (TryFile("HiLep.CentralN_b.N0.dat")) then
     call FetchHist_multi(HHH1,file='HiLep.CentralN_b.C_CN_b.dat')
     call FetchHist_multi(HHH2,file='HiLep.CentralN_b.C_CN_bT.dat')
     call FetchHist_multi(HHH3,file='HiLep.CentralN_b.C_CN_bZ.dat')

     do i=0,1
        call FetchHist2D_multi(h2D,file='HiLep.CentralN_b.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_b.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH1,&
               & file='HiLep.CentralN_b.norm.'//nucName(i)//'.dat')

        call FetchHist2D_multi(h2D,file='HiLep.CentralN_bT.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_bT.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH2,&
               & file='HiLep.CentralN_bT.norm.'//nucName(i)//'.dat')

        call FetchHist2D_multi(h2D,file='HiLep.CentralN_bZ.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_bZ.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH3,&
               & file='HiLep.CentralN_bZ.norm.'//nucName(i)//'.dat')

        call FetchHist2D_multi(h2D,file='HiLep.CentralN_bZ_Ekin.'//nucName(i)//'.dat')

        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN_bZ_Ekin.'//nucName(i)//'.dat',dump=.true.)
        call WriteHist2D_Gnuplot(H2D, MaxVal=0.0,H2=HHH3,&
               & file='HiLep.CentralN_bZ_Ekin.norm.'//nucName(i)//'.dat')

        call FetchHist2D_multi(h2D,file='HiLep.CentralN.pTz.'//nucName(i)//'.dat')
        call WriteHist2D_Gnuplot(H2D, mul=1./nLeptons, add=1e-20,&
               & file='HiLep.CentralN.pTz.'//nucName(i)//'.dat',dump=.true.)

     end do
  end if

  if (TryFile("HiLep.EIC_zH.idH.Acc.dat")) then
     call FetchHistMP_multi(hMP,"HiLep.EIC_zH.idH.Acc.dat")
     call WriteHistMP(hMP,add=1e-20,mul=1./nLeptons,iColumn=1,&
          &file='HiLep.EIC_zH.idH.Acc.dat',dump=.true.)

     call FetchHistMP_multi(hMP,"HiLep.EIC_nu.idH.Acc.dat")
     call WriteHistMP(hMP,add=1e-20,H2=hLep_nu,iColumn=1,&
          &file='HiLep.EIC_nu.idH.Acc.dat',dump=.true.)

     call FetchHistMP_multi(hMP,"HiLep.EIC_Q2.idH.Acc.dat")
     call WriteHistMP(hMP,add=1e-20,H2=hLep_Q2,iColumn=1,&
          &file='HiLep.EIC_Q2.idH.Acc.dat',dump=.true.)

     do i=1,5
        call FetchHistMP_multi(hMP,'HiLep.EIC_zH.idH.Q2_'//Achar(i+48)//'.Acc.dat')
        call WriteHistMP(hMP,add=1e-20,mul=1./nLeptons,iColumn=1,&
             &file='HiLep.EIC_zH.idH.Q2_'//Achar(i+48)//'.Acc.dat',dump=.true.)

        call FetchHistMP_multi(hMP,'HiLep.EIC_nu.idH.Q2_'//Achar(i+48)//'.Acc.dat')
        call WriteHistMP(hMP,add=1e-20,H2=hLep_nu,iColumn=1,&
             &file='HiLep.EIC_nu.idH.Q2_'//Achar(i+48)//'.Acc.dat',dump=.true.)
     end do

  end if

contains

  !****************************************************************************
  subroutine FetchHist_multi(H,file,flagOK)
    use output

    implicit none
    type(histogram), intent(inout)      :: H
    character*(*), intent(in)           :: file
    logical, intent(out), OPTIONAL      :: flagOK

    type(histogram) :: H1
    integer :: i
    logical :: flag

    if (present(flagOK)) flagOK = .false.

    if (.not.TryFile(trim(file)//".bin")) return

    call FetchHist(H,trim(Dirs(1))//"/"//trim(file)//".bin",flagOK=flag)
    if (.not.flag) return

    do i=2,nDirs
       call FetchHist(H1,trim(Dirs(i))//"/"//trim(file)//".bin",flagOK=flag)
       if (.not.flag) cycle
       call sumHist(H,H1)
    end do

    if (present(flagOK)) flagOK = .true.

  end subroutine FetchHist_multi


  !****************************************************************************
  subroutine FetchHistMP_multi(H,file,flagOK)
    use output

    implicit none
    type(histogramMP), intent(inout)    :: H
    character*(*), intent(in)           :: file
    logical, intent(out), OPTIONAL      :: flagOK

    type(histogramMP) :: H1
    integer :: i
    logical :: flag

    if (present(flagOK)) flagOK = .false.

    if (.not.TryFile(trim(file)//".bin")) return

    call FetchHistMP(H,trim(Dirs(1))//"/"//trim(file)//".bin",flagOK=flag)
    if (.not.flag) return

    do i=2,nDirs
       call FetchHistMP(H1,trim(Dirs(i))//"/"//trim(file)//".bin",flagOK=flag)
       if (.not.flag) cycle
       call sumHistMP(H,H1)
    end do

    if (present(flagOK)) flagOK = .true.

  end subroutine FetchHistMP_multi

  !****************************************************************************
  subroutine FetchHistMC_multi(H,file,flagOK)
    use output

    implicit none
    type(histogramMC), intent(inout)    :: H
    character*(*), intent(in)           :: file
    logical, intent(out), OPTIONAL      :: flagOK

    type(histogramMC) :: H1
    integer :: i
    logical :: flag

    if (present(flagOK)) flagOK = .false.

    if (.not.TryFile(trim(file)//".bin")) return

    call FetchHistMC(H,trim(Dirs(1))//"/"//trim(file)//".bin",flagOK=flag)
    if (.not.flag) return

    do i=2,nDirs
       call FetchHistMC(H1,trim(Dirs(i))//"/"//trim(file)//".bin",flagOK=flag)
       if (.not.flag) cycle
       call sumHistMC(H,H1)
    end do

    if (present(flagOK)) flagOK = .true.

  end subroutine FetchHistMC_multi

  !****************************************************************************
  subroutine FetchHist2D_multi(H,file,flagOK)
    use output

    implicit none
    type(histogram2D), intent(inout)    :: H
    character*(*), intent(in)           :: file
    logical, intent(out), OPTIONAL      :: flagOK

    type(histogram2D) :: H1
    integer :: i
    logical :: flag

    if (present(flagOK)) flagOK = .false.

    if (.not.TryFile(trim(file)//".bin")) return

    call FetchHist2D(H,trim(Dirs(1))//"/"//trim(file)//".bin",flagOK=flag)
    if (.not.flag) return

    do i=2,nDirs
       call FetchHist2D(H1,trim(Dirs(i))//"/"//trim(file)//".bin",flagOK=flag)
       if (.not.flag) cycle
       call sumHist2D(H,H1)
    end do

    if (present(flagOK)) flagOK = .true.

  end subroutine FetchHist2D_multi

  !****************************************************************************
  logical function TryFile(file)
    implicit none
    character*(*),  intent(in)          :: file

    integer :: ios

    open(9932,file=trim(Dirs(1))//"/"//trim(file),status='OLD',IOSTAT=ios)
    close(9932)
    TryFile = (ios.eq.0)

  end function TryFile

  !****************************************************************************
  subroutine BuildDirs()
    implicit none

    integer :: ios
    character(len=1000) :: line
    integer :: nLen, i

    open(9932,file='directories.txt',status='OLD',IOSTAT=ios)
    if (ios/=0) then
       write(*,*) "WARNING: file 'directories.txt' does not exists."
       write(*,*) "Reading directories '0/' ... '9/' instead."

       nDirs=10
       allocate(character(1) :: Dirs(nDirs))
       do i=1,10
          Dirs(i)(1:1) = Achar(i-1+48)
       end do

!!$       do i=1,nDirs
!!$          write(*,*) i,': >',trim(Dirs(i)),'<'
!!$       end do

       return
    end if

    write(*,*) "reading 'directories.txt'..."
    nDirs = 0
    nLen = 0
    do
       read (9932,'(A)',iostat=ios) line
       if (ios/=0) exit
       nDirs = nDirs+1
       nLen = max(nLen, len(trim(line)))
    end do
    write(*,*) 'nDirs = ',nDirs,'   nLen=',nLen

    allocate(character(nLen) :: Dirs(nDirs))

    rewind(9932)
    nDirs = 0
    do
       read (9932,'(A)',iostat=ios) line
       if (ios/=0) exit
       nDirs = nDirs+1
       nLen = len(trim(line))
       if (line(nLen:nLen)=='/') nLen=nLen-1
       Dirs(nDirs) = trim(line(:nLen))
    end do
    write(*,*) "reading 'directories.txt' done"

!!$    do i=1,nDirs
!!$       write(*,*) i,': >',trim(Dirs(i)),'<'
!!$    end do


  end subroutine BuildDirs
  !****************************************************************************

end program HiLeptonMixRuns
