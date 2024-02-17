!******************************************************************************
!****m* /rhoWidthMedium_tables
! NAME
! module rhoWidthMedium_tables
! PURPOSE
! Implements the routines for the medium width of the rho-meson
! at finite density and temperature
!******************************************************************************
module rhoWidthMedium_tables

  use IDTable
  use bzip

  implicit none
  private

  real, save :: delta_mass,delta_p,delta_dens,delta_temp
  real, save :: temp_max, dens_max, p_max, min_mass, max_mass
  integer, save :: n_temp, n_dens, n_p, n_mass
  integer, save :: n_dens_min=100000,  n_dens_max=0
  integer, save :: n_temp_min=100000, n_temp_max=0
  integer, save :: n_p_min=100000, n_p_max=0

  character(1000) :: fileName
  type(bzFile) :: f
  character(len=400) :: buf

  real(4), save, dimension(:,:,:,:), allocatable :: widthTable

  public :: get_GammaColl_rho, readTable, writeTable

contains


  !****************************************************************************
  !****f* rhoWidthMedium_tables/get_GammaColl_rho
  ! NAME
  ! function get_GammaColl_rho(mass,p,med) result(GammaColl)
  !
  ! PURPOSE
  ! * Returns the collisonal width of rho-meson in its rest frame according to GammaColl=sigma*rho*v * gamma_Lor
  ! * An average over the Fermi distribution at finite temperature is performed
  ! * Sigma is given by the actual full collision term, with the incoming particles in the vacuum
  !
  ! INPUTS
  ! * real, intent(in)          :: mass  --- off-shell mass of meson (GeV)
  ! * real, intent(in)          :: p     --- absolute Momentum in LRF (GeV/c)
  ! * type(medium), intent(in)  :: med   --- medium information
  !
  ! OUTPUT
  ! * real :: collisional width GammaColl (GeV)
  !
  !****************************************************************************
  function get_GammaColl_rho(mass,p,med) result(GammaColl)
    use mediumDefinition
    use tabulation, only: interpolate4

    real, intent(in)               :: mass
    real, intent(in)               :: p
    type(medium), intent(in)       :: med

    real :: GammaColl, dens, temp
    real(4), dimension(1:4), save :: minX, maxX, deltaX
    real(4), dimension(1:4) :: X
!    character(1000) :: fname
    logical, save :: readTable_flag=.true., first=.true.

    dens=med%density       ! baryon density in LRF (fm^-3)
    temp=med%temperature   ! temperature (GeV)

    if (readTable_flag) then
       ! READ IN THE DATA FILES:
!       fname='GammaColl_mod_rhoN.103.dat.bz2'
!       call readTable(fname)
       call readTable
       readTable_flag=.false.
    end if

!   widthTable(index_mass,index_p,index_dens,index_temp)

    if (first) then
       ! Initialize arrays containing grid information
       minX = (/ min_mass, 0., 0., 0. /)
       maxX = (/ max_mass, p_max, dens_max,temp_max /)
       deltaX = (/ delta_mass, delta_p, delta_dens, delta_temp /)
       first = .false.
    end if

    x(1:4) = (/ mass,p,dens,temp /)
    GammaColl = interpolate4(X,minX,maxX,deltaX,widthTable)

  end function get_GammaColl_rho


  !****************************************************************************
  !****s* rhoWidthMedium_tables/readTable
  ! NAME
  ! subroutine readTable(fileName_in)
  ! INPUTS
  ! *  character(*), optional :: fileName_in --- name of input bz2-file
  !
  ! PURPOSE
  ! * Reads tabulated array widthTable(index_mass,index_p,index_dens,index_temp)
  !   from file buuinput/inMediumWidth/GammaColl.103.dat.bz2
  !   or from file 'fileName_in'
  !****************************************************************************
  subroutine readTable(fileName_in)
    use inputGeneral, only: path_To_Input
    use output, only: Write_ReadingInput

    character(*), optional, intent(in) :: fileName_in

    integer :: n_dens_min_in,  n_dens_max_in, n_temp_min_in, n_temp_max_in
    integer :: n_p_min_in, n_p_max_in
    integer :: ll,ios,index_mass,index_p,index_dens !,index_temp
    !    real :: mass,p,dens,temp
    character(5) ::raute
    logical, save :: firstCall=.true.

    call Write_ReadingInput("Collisional Width of rho",0)

    if(present(fileName_in)) then
       fileName=trim(path_to_Input)//'/inMediumWidth/'//trim(fileName_in)
    else
       fileName=trim(path_to_Input)//'/inMediumWidth/GammaColl.103.dat.bz2'
    end if

    write(*,*)' fileName: ', fileName

    f = bzOpenR(trim(fileName))

    ll = 0
    call bzReadLine(f,buf,ll)
    read(buf(1:ll),*,iostat=ios) raute, n_temp, n_dens, n_p, n_mass
    if (ios.ne.0) then
       write(*,*) 'Error in opening input file (1): ',trim(fileName)
       stop
    end if
    write(*,*)' n_temp, n_dens, n_p, n_mass : ', n_temp, n_dens, n_p, n_mass

    ll = 0
    call bzReadLine(f,buf,ll)
    read(buf(1:ll),*,iostat=ios) raute, temp_max, dens_max, p_max, min_mass, max_mass
    if (ios.ne.0) then
       write(*,*) 'Error in opening input file (2): ',trim(fileName)
       stop
    end if
    write(*,*)' temp_max, dens_max, p_max, min_mass, max_mass : ',&
              & temp_max, dens_max, p_max, min_mass, max_mass

    ll = 0
    call bzReadLine(f,buf,ll)
    read(buf(1:ll),*,iostat=ios) raute, n_p_min_in, n_p_max_in, n_dens_min_in, n_dens_max_in, n_temp_min_in, n_temp_max_in
    if (ios.ne.0) then
       write(*,*) 'Error in opening input file (3): ',trim(fileName)
       stop
    end if
    write(*,*)' n_p_min_in, n_p_max_in : ', n_p_min_in, n_p_max_in
    write(*,*)' n_dens_min_in, n_dens_max_in : ', n_dens_min_in, n_dens_max_in
    write(*,*)' n_temp_min_in, n_temp_max_in : ', n_temp_min_in, n_temp_max_in

    if(firstCall) then
       allocate(widthTable(0:n_mass,0:n_p,0:n_dens,0:n_temp))
       delta_mass=(max_mass-min_mass)/float(n_mass)
       delta_p=p_max/float(n_p)
       delta_dens=dens_max/float(n_dens)
       delta_temp=temp_max/float(n_temp)
       firstCall=.false.
    end if

    ll = 0
    do index_dens=n_dens_min_in,n_dens_max_in
        do index_p=n_p_min_in,n_p_max_in
            do index_mass=0,n_mass
                call bzReadLine(f,buf,ll)
                read(buf(1:ll),*,iostat=ios) widthTable(index_mass,index_p,index_dens,n_temp_min_in:n_temp_max_in)
                if (IOS.ne.0) then
                   write(*,*) ' Error in opening input file (3): ', trim(fileName)
                   write(*,*) ' index_dens, index_p, index_mass : ', index_dens, index_p, index_mass
                   stop
                end if
            end do
        end do
    end do

    call bzCloseR(f)

    call Write_ReadingInput("Collisional Width of rho",1)

    if(n_p_min_in.lt.n_p_min) n_p_min=n_p_min_in
    if(n_p_max_in.gt.n_p_max) n_p_max=n_p_max_in

    if(n_dens_min_in.lt.n_dens_min) n_dens_min=n_dens_min_in
    if(n_dens_max_in.gt.n_dens_max) n_dens_max=n_dens_max_in

    if(n_temp_min_in.lt.n_temp_min) n_temp_min=n_temp_min_in
    if(n_temp_max_in.gt.n_temp_max) n_temp_max=n_temp_max_in

!    open(1,file='GammaColl.103.chk',status='unknown')
!    do index_temp=0,n_temp
!        temp=float(index_temp)*delta_temp
!        do index_dens=0,n_dens
!            dens=float(index_dens)*delta_dens
!            write(1,'(2(A,F15.4))')'# temp=', temp,' dens=', dens
!            write(1,'(A)')'# p:             mass:          Gamma_coll:'
!            do index_p=0,n_p
!                p=float(index_p)*delta_p
!                do index_mass=0,n_mass
!                    mass=min_mass+float(index_mass)*delta_mass
!                    write(1,'(4(2x,e13.7))') p,mass,widthTable(index_mass,index_p,index_dens,index_temp)
!                end do
!            end do
!        end do
!    end do
!    close(1)

  end subroutine readTable

  !****************************************************************************
  !****s* rhoWidthMedium_tables/writeTable
  ! NAME
  ! subroutine writeTable
  !
  ! PURPOSE
  ! * Writes tabulated array widthTable(index_mass,index_p,index_dens,index_temp)
  !   to the file GammaColl.103.dat.bz2
  !****************************************************************************
  subroutine writeTable

    use output, only: intToChar

    character(25) :: format
    integer :: index_mass,index_p,index_dens,index_temp
    real :: mass,p,dens,temp

    format='(' // intToChar(n_temp_max-n_temp_min+1) // '(1x,E13.7))'

    ! Write the table to file:
    fileName='./GammaColl.103.dat.bz2'

    f = bzOpenW(trim(fileName))
    write(buf,'(A,4I6)')'#', n_temp, n_dens, n_p, n_mass
    call bzWriteLine(f,buf)
    write(buf,'(A,5(1x,E13.7))')'#', temp_max, dens_max, p_max, min_mass, max_mass
    call bzWriteLine(f,buf)
    write(buf,'(A,6I6)')'#', n_p_min, n_p_max, n_dens_min, n_dens_max, n_temp_min, n_temp_max
    call bzWriteLine(f,buf)
    do index_dens=n_dens_min,n_dens_max
        do index_p=n_p_min,n_p_max
            do index_mass=0,n_mass
                write(buf,format) widthTable(index_mass,index_p,index_dens,n_temp_min:n_temp_max)
                call bzWriteLine(f,buf)
            end do
        end do
    end do
    call bzCloseW(f)


    open(1,file='GammaColl.103.chk',status='unknown')
    do index_temp=n_temp_min,n_temp_max
        temp=float(index_temp)*delta_temp
        do index_dens=n_dens_min,n_dens_max
            dens=float(index_dens)*delta_dens
            write(1,'(2(A,F15.4))')'# temp=', temp,' dens=', dens
            write(1,'(A)')'# p:             mass:          Gamma_coll:'
            do index_p=n_p_min,n_p_max
                p=float(index_p)*delta_p
                do index_mass=0,n_mass
                    mass=min_mass+float(index_mass)*delta_mass
                    write(1,'(4(2x,e13.7))') p,mass,widthTable(index_mass,index_p,index_dens,index_temp)
                end do
            end do
        end do
    end do
    close(1)

    deallocate(widthTable)

  end subroutine writeTable


end module rhoWidthMedium_tables
