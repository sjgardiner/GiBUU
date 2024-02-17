program read_Varirates

  implicit none

  character*(200) :: BUF

  integer :: ios
  integer :: val1,val2,n1,n2,i
  integer, dimension(3) :: ID_in, Q_in
  integer, dimension(20) :: ID_out, Q_out



  open(31,file="VariRate.rates.csv",action="read",iostat=ios)
  if (ios /= 0) then
     write(*,*) "input file not found."
     stop
  end if

  do
     read(31,FMT='(a)',iostat=ios) BUF
     if (ios/=0) exit

     !     read(BUF,'(4(i6,","))') val1,val2,n1,n2
     read(BUF,*) val1,val2,n1,n2

     ID_in = 0
     Q_in = 0
     ID_out = 0
     Q_out = 0

     read(BUF,*) val1,val2,n1,n2, &
          (ID_in(i), Q_in(i),i=1,n1), &
          (ID_out(i), Q_out(i),i=1,n2)

!     write(*,*) '>',trim(BUF),'<'
!     write(*,*) val1,val2,n1,n2, ID_in(1:n1), ID_out(1:n2)

     call calcLambdaSigma(.false.)

  end do

  call calcLambdaSigma(.true.)

contains

  subroutine calcLambdaSigma(doPrint)

    logical, intent(in) :: doPrint

    integer :: i,j, iL_in, iS_in, iL_out, iS_out
    logical, save :: first = .true.

    integer, dimension(3,3,2), save :: arrBB_L, arrBB_S
    integer, dimension(2), save :: arrBM_L, arrBM_S, &
         arrDec_L, arrDec_S, &
         arrChEx_S
    integer, dimension(2,2), save :: arrConv_SL
    integer :: iB1,iB2


    if (n1==2 .and. n2==2) then
       if ( ID_in(1)==ID_out(1) .and. &
            ID_in(2)==ID_out(2) .and. &
            Q_in(1)==Q_out(1) .and. &
            Q_in(2)==Q_out(2) ) return ! ==> elastic
    end if

    if (first) then

       ! B = N,Delta,R

       arrBB_L = 0 ! BB -> L + X
       arrBB_S = 0 ! BB -> S + X
       arrBM_L = 0 ! Bm -> L + X
       arrBM_S = 0 ! Bm -> S + X
       arrDec_L = 0 ! X -> L + X'
       arrDec_S = 0 ! X -> S + X'
       arrConv_SL = 0 ! S + X -> L + X'
       arrChEx_S = 0 ! S' + X -> S + X'

       first = .false.
    end if

    if (doPrint) then
       do i=1,3
          do j=i,3
             write(*,*) "Lambda, BB:", i,j, arrBB_L(i,j,:)
             write(331,*) arrBB_L(i,j,:)
          end do
       end do
       write(*,*) "Lambda, Bm:                        ",arrBM_L(:)
       write(*,*) "Lambda, Decay:                     ",arrDec_L(:)

       write(331,*) arrBM_L(:)
       write(331,*) arrDec_L(:)


       do i=1,3
          do j=i,3
             write(*,*) "Sigma, BB:", i,j, arrBB_S(i,j,:)
             write(331,*) arrBB_S(i,j,:)
          end do
       end do
       write(*,*) "Sigma, Bm:                        ", arrBM_S(:)
       write(*,*) "Sigma, Decay:                     ", arrDec_S(:)

       write(331,*) arrBM_S(:)
       write(331,*) arrDec_S(:)

       write(*,*) "Conversion Sigma0  -> Lambda0:    ", arrConv_SL(1,:)
       write(*,*) "Conversion Sigma+- -> Lambda0:    ", arrConv_SL(2,:)

       write(*,*) "Sigma, Charge exchange:           ", arrChEx_S(:)

       write(331,*) arrConv_SL(1,:)
       write(331,*) arrConv_SL(2,:)
       write(331,*) arrChEx_S(:)


       return
    end if

    iL_in = 0
    do i=1,n1
       if ((ID_in(i)==32).and.(Q_in(i)==0)) then
          iL_in = i
          exit
       end if
    end do

    iS_in = 0
    do i=1,n1
       if ((ID_in(i)==33).and.(Q_in(i)==0)) then
          iS_in = i
          exit
       end if
    end do

    iL_out = 0
    do i=1,n2
       if ((ID_out(i)==32).and.(Q_out(i)==0)) then
          iL_out = i
          exit
       end if
    end do

    iS_out = 0
    do i=1,n2
       if ((ID_out(i)==33).and.(Q_out(i)==0)) then
          iS_out = i
          exit
       end if
    end do

    if (iL_in+iL_out+ iS_in+iS_out == 0) return ! ==> no Lambda/Sigma

    if (iL_out > 0 .and. iS_in > 0) then
       arrConv_SL(1,:) = arrConv_SL(1,:) + (/val1,val2/)
       return
    end if
    if (iL_in > 0 .and. iS_out > 0) then
       arrConv_SL(1,:) = arrConv_SL(1,:) + (/val2,val1/)
       return
    end if

    if (iL_out > 0) then
       select case (n1)
       case (1)
          arrDec_L(1:2) = arrDec_L(1:2)  + (/val1,val2/)
          return
       case (2)
          select case(ID_in(1))
          case (1)
             iB1 = 1
          case (2)
             iB1 = 2
          case (3:31)
             iB1 = 3
          case default
             iB1 = 0
          end select

          select case(ID_in(2))
          case (1)
             iB2 = 1
          case (2)
             iB2 = 2
          case (3:31)
             iB2 = 3
          case default
             iB2 = 0
          end select

          if (iB1>0 .and.iB2>0) then
             arrBB_L(iB1,iB2,1:2) = arrBB_L(iB1,iB2,1:2) + (/val1,val2/)
             return
          end if

          if ( (iB1>0 .and. ID_in(2) > 100) .or. &
               (iB2>0 .and. ID_in(1) > 100) ) then
             arrBM_L(1:2) = arrBM_L(1:2)  + (/val1,val2/)
             return
          end if

       end select
    end if

    if (iS_out > 0) then
       select case (n1)
       case (1)
          arrDec_S(1:2) = arrDec_S(1:2)  + (/val1,val2/)
          return
       case (2)
          select case(ID_in(1))
          case (1)
             iB1 = 1
          case (2)
             iB1 = 2
          case (3:31)
             iB1 = 3
          case default
             iB1 = 0
          end select

          select case(ID_in(2))
          case (1)
             iB2 = 1
          case (2)
             iB2 = 2
          case (3:31)
             iB2 = 3
          case default
             iB2 = 0
          end select

          if (iB1>0 .and.iB2>0) then
             arrBB_S(iB1,iB2,1:2) = arrBB_S(iB1,iB2,1:2) + (/val1,val2/)
             return
          end if

          if ( (iB1>0 .and. ID_in(2) > 100) .or. &
               (iB2>0 .and. ID_in(1) > 100) ) then
             arrBM_S(1:2) = arrBM_S(1:2)  + (/val1,val2/)
             return
          end if
       end select
    end if

    if (iS_out >  0) then
       do i=1,n1
          if (ID_in(i)== 33 .and. Q_in(i)/=0) then
             arrChEx_S(1:2) = arrChEx_S(1:2) + (/val1,val2/)
             return
          end if
       end do
    end if

    if (iS_in >  0) then
       do i=1,n2
          if (ID_out(i)== 33 .and. Q_out(i)/=0) then
             arrChEx_S(1:2) = arrChEx_S(1:2) + (/val2,val1/)
             return
          end if
       end do
    end if

    if (iL_out > 0) then
       do i=1,n1
          if (ID_in(i)== 33 .and. Q_in(i)/=0) then
             arrConv_SL(2,1:2) = arrConv_SL(2,1:2) + (/val1,val2/)
             return
          end if
       end do
    end if

    if (iL_in > 0) then
       do i=1,n2
          if (ID_out(i)== 33 .and. Q_out(i)/=0) then
             arrConv_SL(2,1:2) = arrConv_SL(2,1:2) + (/val2,val1/)
             return
          end if
       end do
    end if

    write(*,*) iL_in,iS_in,iL_out,iS_out, ' >',trim(BUF),'<'



  end subroutine calcLambdaSigma


end program read_Varirates
