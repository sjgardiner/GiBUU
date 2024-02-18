!******************************************************************************
!****m* /EventOutput
! NAME
! module EventOutput
!
! PURPOSE
! This module provides classes for writing event output to disk
! in different formats.
!
! Currently the following formats are supported:
! * "Les Houches" event files
! * "OSCAR 2013" event files
! * "Shanghai 2014"
! * "Root"
! * "NuHepMC"
!
! For a description of the Les Houches format, please refer to:
! * http://arxiv.org/abs/hep-ph/0609017
! * https://gibuu.hepforge.org/trac/wiki/LesHouches
!
! For a description of the OSCAR 2013 format, see:
! * http://phy.duke.edu/~jeb65/oscar2013
!
! For a description of the NuHepMC format, see:
! * https://arxiv.org/abs/2310.13211
!
! INPUTS
! (none)
!******************************************************************************
module EventOutput

  implicit none
  private

  !****************************************************************************
  !****t* EventOutput/EventOutputFile
  ! NAME
  ! type EventOutputFile
  ! PURPOSE
  ! This is an abstract base type to represent a file for event output.
  ! It is used as a common interface for LHOutputFile and OscarOutputFile.
  !
  ! SOURCE
  !
  type, abstract, public :: EventOutputFile
     integer, private :: iFile = 0  ! private file handle
  contains
    ! deferred type-bound procedures that need to be implemented in the derived classes
    procedure(open_ifc),       deferred :: open
    procedure(close_ifc),      deferred :: close
    procedure(write_EH_ifc),   deferred :: write_event_header
    procedure(write_EF_ifc),   deferred :: write_event_footer
    procedure(write_part_ifc), deferred :: write_particle
    procedure :: write_additionalInfo
    ! type-bound procedures that are implemented in the base class
    procedure :: write_real
    procedure :: write_pert
  end type
  !****************************************************************************


  ! interfaces for the deferred methods
  abstract interface
    subroutine open_ifc(this, pert, nCall, nTimeStep)
      import :: EventOutputFile
      class(EventOutputFile) :: this
      logical, intent(in) :: pert
      integer, intent(in) :: nCall, nTimeStep
    end subroutine
    subroutine close_ifc(this)
      import :: EventOutputFile
      class(EventOutputFile), intent(in) :: this
    end subroutine
    subroutine write_EH_ifc(this, nParts, nEvent, wgt, iFE)
      import :: EventOutputFile
      class(EventOutputFile) :: this
      integer, intent(in) :: nParts
      integer, intent(in) :: nEvent            ! number of current event
      real, intent(in), optional :: wgt
      integer, intent(in), optional :: iFE
    end subroutine
    subroutine write_EF_ifc(this)
      import :: EventOutputFile
      class(EventOutputFile), intent(in) :: this
    end subroutine
    subroutine write_part_ifc(this, part)
      use particleDefinition
      import :: EventOutputFile
      class(EventOutputFile) :: this
      type(particle), intent(in) :: part
    end subroutine
  end interface


  !****************************************************************************
  !****t* EventOutput/LHOutputFile
  ! NAME
  ! type LHOutputFile
  ! PURPOSE
  ! This is an extended type for event output in LesHouches format.
  ! It is derived from the base type EventOutputFile and implements its
  ! interfaces.
  !
  ! For a description of the Les Houches format, please refer to:
  ! * http://arxiv.org/abs/hep-ph/0609017
  ! * https://gibuu.hepforge.org/trac/wiki/LesHouches
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: LHOutputFile
  contains
    procedure :: open                 => LH_open
    procedure :: close                => LH_close
    procedure :: write_event_header   => LH_write_event_header
    procedure :: write_event_footer   => LH_write_event_footer
    procedure :: write_particle       => LH_write_particle
    procedure :: write_additionalInfo => LH_write_additionalInfo
  end type
  !****************************************************************************


  !****************************************************************************
  !****t* EventOutput/OscarOutputFile
  ! NAME
  ! type OscarOutputFile
  ! PURPOSE
  ! This is an extended type for event output in OSCAR 2013 format.
  ! It is derived from the base type EventOutputFile and implements its
  ! interfaces.
  !
  ! For a description of the OSCAR 2013 format, see:
  ! * http://phy.duke.edu/~jeb65/oscar2013
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: OscarOutputFile
     integer, private :: nEvent = 0   ! store the event number
     integer, private :: iFE = 0      ! store the iFE value
  contains
    procedure :: open                 => Oscar_open
    procedure :: close                => Oscar_close
    procedure :: write_event_header   => Oscar_write_event_header
    procedure :: write_event_footer   => Oscar_write_event_footer
    procedure :: write_particle       => Oscar_write_particle
  end type
  !****************************************************************************

  !****************************************************************************
  !****t* EventOutput/OscarExtOutputFile
  ! NAME
  ! type OscarExtOutputFile
  ! PURPOSE
  ! This is an extended type for event output in OSCAR 2013 format.
  ! It is derived from the base type EventOutputFile and implements its
  ! interfaces.
  !
  ! For a description of the OSCAR 2013 format, see:
  ! * http://phy.duke.edu/~jeb65/oscar2013
  !
  ! This is an extended version whith much more output columns
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: OscarExtOutputFile
     integer, private :: nEvent = 0   ! store the event number
     integer, private :: iFE = 0      ! store the iFE value
   contains
    procedure :: open                 => OscarExt_open
    procedure :: close                => OscarExt_close
    procedure :: write_event_header   => OscarExt_write_event_header
    procedure :: write_event_footer   => OscarExt_write_event_footer
    procedure :: write_particle       => OscarExt_write_particle
  end type
  !****************************************************************************


  !****************************************************************************
  !****t* EventOutput/ShanghaiOutputFile
  ! NAME
  ! type ShanghaiOutputFile
  ! PURPOSE
  ! This is an extended type for event output in Shanghai2014 format.
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: ShanghaiOutputFile
  contains
    procedure :: open                 => Shanghai_open
    procedure :: close                => Shanghai_close
    procedure :: write_event_header   => Shanghai_write_event_header
    procedure :: write_event_footer   => Shanghai_write_event_footer
    procedure :: write_particle       => Shanghai_write_particle
  end type
  !****************************************************************************

  !****************************************************************************
  !****t* EventOutput/RootOutputFile
  ! NAME
  ! type RootOutputFile
  ! PURPOSE
  ! This is an extended type for event output in Root format.
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: RootOutputFile
    real, private :: weight
  contains
    procedure :: open                 => Root_open
    procedure :: close                => Root_close
    procedure :: write_event_header   => Root_write_event_header
    procedure :: write_event_footer   => Root_write_event_footer
    procedure :: write_particle       => Root_write_particle
    procedure :: write_additionalInfo => Root_write_additionalInfo
  end type RootOutputFile
  !****************************************************************************

  !****************************************************************************
  !****t* EventOutput/NuHepMCOutputFile
  ! NAME
  ! type NuHepMCOutputFile
  ! PURPOSE
  ! This is an extended type for event output in the NuHepMC format
  ! defined in https://arxiv.org/abs/2310.13211
  !
  ! SOURCE
  !
  type, extends(EventOutputFile), public :: NuHepMCOutputFile
    real, private :: weight
    integer, private :: particle_count
  contains
    procedure :: open                 => NuHepMC_open
    procedure :: close                => NuHepMC_close
    procedure :: write_event_header   => NuHepMC_write_event_header
    procedure :: write_event_footer   => NuHepMC_write_event_footer
    procedure :: write_particle       => NuHepMC_write_particle
    procedure :: write_additionalInfo => NuHepMC_write_additionalInfo

    ! Helper functions for preparing the output
    procedure :: get_proc_ID => nuhepmc_proc_ID_from_buu
  end type NuHepMCOutputFile
  !****************************************************************************

contains

!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/write_additionalInfo
  ! NAME
  ! subroutine write_additionalInfo(this, iFE, pNode)
  ! PURPOSE
  ! Write additional info about the event, depending on eventtype.
  !****************************************************************************
  subroutine write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition, only: tParticleListNode
    class(EventOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode
    ! do nothing here by default
  end subroutine write_additionalInfo

!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/LH_open
  ! NAME
  ! subroutine LH_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output according to the
  ! "Les Houches Event Files" standard.
  !****************************************************************************
  subroutine LH_open(this, pert, nCall, nTimeStep)
    class(LHOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (pert) then
      buf1 = '.Pert.'
      this%iFile = 721
    else
      buf1 = '.Real.'
      this%iFile = 722
    end if
    write(buf2,'(I8.8)') nCall
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.lhe'

    ! open file
    open(this%iFile, file=fName, status='unknown')
    rewind(this%iFile)

    write(this%iFile,'(A)') '<LesHouchesEvents version="1.0">'
    write(this%iFile,'(A)') '<!-- File generated by GiBUU. For documentation see ' //  &
                            'https://gibuu.hepforge.org/trac/wiki/LesHouches -->'

    write(this%iFile,'(A)') '<header>'
    write(this%iFile,'(A)') '     <!-- individual XML tags may follow -->'
    write(this%iFile,'(A)') '</header>'

    write(this%iFile,'(A)') '<init>'
    write(this%iFile,'(1P,2I8,2E14.6,6I6)') 0,0, 0.,0., 0,0,0,0,0,0
    write(this%iFile,'(A)') '</init>'

  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_close
  ! NAME
  ! subroutine LH_close(this)
  ! PURPOSE
  ! Close a file after outputting Les-Houches event information.
  !****************************************************************************
  subroutine LH_close(this)
    class(LHOutputFile), intent(in) :: this

    write(this%iFile, '(A)') '</LesHouchesEvents>'
    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_event_header
  ! NAME
  ! subroutine LH_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for a Les-Houches event, including the number of particles
  ! and the event weight.
  !****************************************************************************
  subroutine LH_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(LHOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(len=15), parameter :: f1 = '(1P,2I6,4E14.6)'
    integer, parameter :: IDPRUP = 0
    real, parameter :: SCALUP = 0.0, AQEDUP = 0.0, AQCDUP = 0.0
    real :: weight

    if (present(wgt)) then
      weight = wgt
    else
      weight = 1.
    end if

    write(this%iFile,'(A)') '<event>'
    write(this%iFile,f1) nParts, IDPRUP, weight, SCALUP, AQEDUP, AQCDUP
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_event_footer
  ! NAME
  ! subroutine LH_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a Les-Houches event.
  !****************************************************************************
  subroutine LH_write_event_footer(this)
    class(LHOutputFile), intent(in) :: this
    write(this%iFile,'(A)') '</event>'
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_particle
  ! NAME
  ! subroutine LH_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in Les-Houches format.
  !****************************************************************************
  subroutine LH_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU

    class(LHOutputFile) :: this
    type(particle), intent(in) :: part

    character(len=22), parameter :: f2 = '(1P,I8,5I5,5E18.10,A6)'
    integer :: KF

    KF = KFfromBUU(part)
    write(this%iFile,f2) KF, 0, 0,0, 0,0, &
                         part%mom(1:3), part%mom(0), &
                         sqrts(part), '0. 9.'
  end subroutine


  !****************************************************************************
  !****s* EventOutput/LH_write_additionalInfo
  ! NAME
  ! subroutine LH_write_additionalInfo(this, iFE)
  ! PURPOSE
  ! Write additional info about the event, depending on eventtype.
  !
  ! This routine tries to find additional information about the event.
  ! It tries routines for different event types, which only return
  ! some information, if it was really stored.
  !
  ! The following cases are handled:
  ! * For eventtype "HiLep", the following line is added:
  !     # 14 nu Q2 eps phiLepton Eventtype
  !   (14 is the magic number of "HiLepton")
  ! * For eventtype "neutrino", the following line is added:
  !     # 5 Eventtype Weight momLepIn(0:3) momLepOut(0:3) momNuc(0:3)
  !   (5 is the magic number for neutrino events)
  ! * For eventtype "heavyIon", the following line is added:
  !     # 1 b
  !   (1 is the magic number of "heavyIon", b is the impact parameter in fm)
  ! * For eventtype "hadron", the following line is added:
  !     # 300 b
  !   (300 is the magic number of "hadron", b is the impact parameter in fm)
  !****************************************************************************
  subroutine LH_write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use inputGeneral, only: eventType
    use eventtypes, only: hiLepton, neutrino, heavyIon, hadron
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use FreezeoutAnalysis, only: getFreezeoutAnalysis_Pert
    use PIL_freezeout, only: PIL_freezeout_GET

    class(LHOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode

    real :: weight,nu,Q2,eps,phiL
    integer :: evtType, chrg_nuc
    real,dimension(0:3) :: momLepIn, momLepOut, momBos, momNuc
    type(particle), pointer :: pPart
    real, dimension(0:3) :: pos
    integer :: history
    logical :: escaped

    select case (eventType)
    case (heavyIon)
      write(this%iFile,'(A,ES13.4)') '# 1 ', b_HI
    case (hadron)
      write(this%iFile,'(A,ES13.4)') '# 300 ', b_had
    case (neutrino)
      if (.not. present(iFE)) return
      if (NeutrinoProdInfo_Get(iFE,evtType,Weight,momLepIn,momLepOut,momBos,momNuc,chrg_nuc)) &
        write(this%iFile,'(A,I5,1P,e18.10,1P,3(" ",4e18.10),0P,A)') &
           '# 5 ', evtType, Weight, momLepIn, momLepOut, momNuc
    case (hiLepton)
      if (.not. present(iFE)) return
      if (EventInfo_HiLep_Get(0,iFE,Weight,nu,Q2,eps,evtType,phi_Lepton=phiL)) &
        write(this%iFile,'(A,1P,4e13.4,0P,I8)') '# 14 ', nu,Q2,eps,phiL,evtType
    end select

    if (getFreezeoutAnalysis_Pert() .and. present(pNode)) then
      do
        if (.not. associated(pNode)) exit
        pPart => pNode%V
        if (PIL_freezeout_GET(pPart%number, history, escaped, pos)) then
          write(this%iFile,'(A,1P,4e13.4,0P,I12,L3)') '# 1001 ', pos, history, escaped
        else
          write(this%iFile,'(A,1P,4e13.4,0P,I12,L3)') '# 1002 '
        end if

        pNode => pNode%next
      end do
    end if

  end subroutine LH_write_additionalInfo


!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/Oscar_open
  ! NAME
  ! subroutine Oscar_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output according to the
  ! "OSCAR 2013" standard.
  !****************************************************************************
  subroutine Oscar_open(this, pert, nCall, nTimeStep)
    class(OscarOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=4)  :: buf

    if (pert) then
      buf = 'Pert'
      this%iFile = 723
    else
      buf = 'Real'
      this%iFile = 724
    end if
    fName = 'EventOutput.' // trim(buf) // '.oscar'

    if (nCall == 1) then
      ! open file for the first time
      open(this%iFile, file=fName, status='unknown')
      ! write header
      write(this%iFile,'(A)') '#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID'
      write(this%iFile,'(A)') '# Units: fm fm fm fm GeV GeV GeV GeV GeV none none'
      write(this%iFile,'(A)') '# File generated by GiBUU (https://gibuu.hepforge.org)'
    else
      ! append to exiting file
      open(this%iFile, file=fName, status='old', position='append')
    end if

  end subroutine


  !****************************************************************************
  !****s* EventOutput/Oscar_close
  ! NAME
  ! subroutine Oscar_close(this)
  ! PURPOSE
  ! Close a file after outputting OSCAR 2013 event information.
  !****************************************************************************
  subroutine Oscar_close(this)
    class(OscarOutputFile), intent(in) :: this

    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Oscar_write_event_header
  ! NAME
  ! subroutine Oscar_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an OSCAR 2013 event, including the number of particles
  ! and the event weight.
  !****************************************************************************
  subroutine Oscar_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(OscarOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(*), parameter :: f1 = '("# event ",I9," out",I4," weight",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," out",I4)'

    this%nEvent = nEvent
    this%iFE = -1
    if (present(iFE)) this%iFE = iFE

    if (present(wgt)) then
       write(this%iFile,f1) this%nEvent, nParts, wgt
    else
       write(this%iFile,f0) this%nEvent, nParts
    end if

  end subroutine

  !****************************************************************************
  !****s* EventOutput/Oscar_write_event_footer
  ! NAME
  ! subroutine Oscar_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes an OSCAR 2013 event.
  !****************************************************************************
  subroutine Oscar_write_event_footer(this)
    use inputGeneral, only: eventType
    use eventtypes, only: heavyIon, hadron, LoPion, HiPion
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use initPion, only: getImpact_Lo => getImpact
    use initHiPion, only: getImpact_Hi => getImpact

    class(OscarOutputFile), intent(in) :: this

    character(*), parameter :: f1 = '("# event ",I9," end 0 impact",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," end 0")'

    select case (eventType)
    case (heavyIon)
       write(this%iFile,f1) this%nEvent, b_Hi
    case (hadron)
       write(this%iFile,f1) this%nEvent, b_had
    case (LoPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Lo(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case (HiPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Hi(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case default
       write(this%iFile,f0) this%nEvent
    end select

  end subroutine

  !****************************************************************************
  !****s* EventOutput/Oscar_write_particle
  ! NAME
  ! subroutine Oscar_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in OSCAR 2013 format.
  !****************************************************************************
  subroutine Oscar_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU

    class(OscarOutputFile) :: this
    type(particle), intent(in) :: part

    character(len=22), parameter :: f = '(9ES17.9,2I9)'
    integer :: PDGcode

    PDGcode = KFfromBUU(part)
    write(this%iFile,f) part%lastCollTime, part%pos(1:3), &
                        sqrts(part), part%mom(0:3), PDGcode, part%number

  end subroutine


!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/OscarExt_open
  ! NAME
  ! subroutine OscarExt_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output according to the
  ! "OSCAR 2013" standard.
  !****************************************************************************
  subroutine OscarExt_open(this, pert, nCall, nTimeStep)
    class(OscarExtOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=4)  :: buf

    if (pert) then
      buf = 'Pert'
      this%iFile = 723
    else
      buf = 'Real'
      this%iFile = 724
    end if
    fName = 'EventOutput.' // trim(buf) // '.oscar'

    if (nCall == 1) then
      ! open file for the first time
      open(this%iFile, file=fName, status='unknown')
      ! write header
      write(this%iFile,'(A)') '#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID mass0 lastCollisionTime mother1 mother2 mother3 generation'
      write(this%iFile,'(A)') '# Units: fm fm fm fm GeV GeV GeV GeV GeV none none GeV fm none none none none'
      write(this%iFile,'(A)') '# File generated by GiBUU (https://gibuu.hepforge.org)'
    else
      ! append to exiting file
      open(this%iFile, file=fName, status='old', position='append')
    end if

  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_close
  ! NAME
  ! subroutine OscarExt_close(this)
  ! PURPOSE
  ! Close a file after outputting OSCAR 2013 event information.
  !****************************************************************************
  subroutine OscarExt_close(this)
    class(OscarExtOutputFile), intent(in) :: this

    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_write_event_header
  ! NAME
  ! subroutine OscarExt_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an OSCAR 2013 event, including the number of particles
  ! and the event weight.
  !****************************************************************************
  subroutine OscarExt_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(OscarExtOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(*), parameter :: f1 = '("# event ",I9," out",I4," weight",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," out",I4)'

    this%nEvent = nEvent
    this%iFE = -1
    if (present(iFE)) this%iFE = iFE

    if (present(wgt)) then
       write(this%iFile,f1) this%nEvent, nParts, wgt
    else
       write(this%iFile,f0) this%nEvent, nParts
    end if
  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_write_event_footer
  ! NAME
  ! subroutine OscarExt_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes an OSCAR 2013 event.
  !****************************************************************************
  subroutine OscarExt_write_event_footer(this)
    use inputGeneral, only: eventType
    use eventtypes, only: heavyIon, hadron, LoPion, HiPion
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use initPion, only: getImpact_Lo => getImpact
    use initHiPion, only: getImpact_Hi => getImpact

    class(OscarExtOutputFile), intent(in) :: this

    character(*), parameter :: f1 = '("# event ",I9," end 0 impact",E14.6)'
    character(*), parameter :: f0 = '("# event ",I9," end 0")'

    select case (eventType)
    case (heavyIon)
       write(this%iFile,f1) this%nEvent, b_Hi
    case (hadron)
       write(this%iFile,f1) this%nEvent, b_had
    case (LoPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Lo(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case (HiPion)
       if (this%iFE>0) then
          write(this%iFile,f1) this%nEvent, getImpact_Hi(this%iFE)
       else
          write(this%iFile,f0) this%nEvent
       end if
    case default
       write(this%iFile,f0) this%nEvent
    end select

  end subroutine


  !****************************************************************************
  !****s* EventOutput/OscarExt_write_particle
  ! NAME
  ! subroutine OscarExt_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in OSCAR 2013 extended format.
  !****************************************************************************
  subroutine OscarExt_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU
    use history, only: history_getParents,history_getGeneration

    class(OscarExtOutputFile) :: this
    type(particle), intent(in) :: part

    character(len=30), parameter :: f = '(9ES17.9,2I9,2ES17.9,4I9)'
    integer :: PDGcode, generation, k
    integer, dimension(1:3) :: parents

    PDGcode = KFfromBUU(part)
    parents = history_getParents(part%history)
    generation=history_getGeneration(part%history)
    do k=1,2
       if (parents(k)>200) parents(k)=200-parents(k)
    end do
    write(this%iFile,f) -99.9, part%pos(1:3), &
         sqrts(part), part%mom(0:3), PDGcode, part%number, &
         part%mass, part%lastCollTime, &
         parents(1:3), generation

  end subroutine


!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/Shanghai_open
  ! NAME
  ! subroutine Shanghai_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_open(this, pert, nCall, nTimeStep)
    class(ShanghaiOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (pert) then
      buf1 = '.Pert.'
      this%iFile = 725
    else
      buf1 = '.Real.'
      this%iFile = 726
    end if
    write(buf2,'(I8.8)') nTimeStep
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.dat'

    ! open file
    open(this%iFile, file=fName, status='unknown', position='append')

  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_close
  ! NAME
  ! subroutine Shanghai_close(this)
  ! PURPOSE
  ! Close a file after outputting event information in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_close(this)
    class(ShanghaiOutputFile), intent(in) :: this

    close(this%iFile)
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_write_event_header
  ! NAME
  ! subroutine Shanghai_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an event in Shanghai2014 format,
  ! including the number of particles and the event weight.
  !****************************************************************************
  subroutine Shanghai_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(ShanghaiOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    character(len=20) :: f

    if (present(wgt)) then
      f = '(A,I4,A,I4,A,E14.6)'
      write(this%iFile,f) '# event ', nEvent, ' out ', nParts, ' weight', wgt
    else
      f = '(A,I4,A,I4)'
      write(this%iFile,f) '# event ', nEvent, ' out ', nParts
    end if
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_write_event_footer
  ! NAME
  ! subroutine Shanghai_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a event in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_write_event_footer(this)
    class(ShanghaiOutputFile), intent(in) :: this
    ! empty, no footer!
  end subroutine


  !****************************************************************************
  !****s* EventOutput/Shanghai_write_particle
  ! NAME
  ! subroutine Shanghai_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in Shanghai2014 format.
  !****************************************************************************
  subroutine Shanghai_write_particle(this, part)
    use particleDefinition
    use minkowski, only: abs4

    class(ShanghaiOutputFile) :: this
    type(particle), intent(in) :: part

    character(len=22), parameter :: f = '(2I9,7ES17.9)'
    integer :: fact

    if (part%anti) then
      fact = -1
    else
      fact = 1
    end if

    write(this%iFile,f) fact*part%ID, part%charge, abs4(part%mom), &
                        part%pos(1:3), part%mom(1:3)

  end subroutine

!******************************************************************************
!******************************************************************************
!******************************************************************************


  !****************************************************************************
  !****s* EventOutput/Root_open
  ! NAME
  ! subroutine Root_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output in Root format.
  !****************************************************************************
  subroutine Root_open(this, pert, nCall, nTimeStep)
    class(RootOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (pert) then
       buf1 = '.Pert.'
    else
       buf1 = '.Real.'
    end if
    write(buf2,'(I8.8)') nCall
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.root'

    ! open file
    call rootinit(fName)

  end subroutine Root_open


  !****************************************************************************
  !****s* EventOutput/Root_close
  ! NAME
  ! subroutine Root_close(this)
  ! PURPOSE
  ! Close a file after outputting event information in Root format.
  !****************************************************************************
  subroutine Root_close(this)
    class(RootOutputFile), intent(in) :: this

    call rootclose()

  end subroutine Root_close


  !****************************************************************************
  !****s* EventOutput/Root_write_event_header
  ! NAME
  ! subroutine Root_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an event in Root format,
  ! including the number of particles and the event weight.
  !****************************************************************************
  subroutine Root_write_event_header(this, nParts, nEvent, wgt, iFE)
    class(RootOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    this%weight = 1.0
    if (present(wgt)) this%weight=wgt

  end subroutine Root_write_event_header


  !****************************************************************************
  !****s* EventOutput/Root_write_event_footer
  ! NAME
  ! subroutine Root_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a event in Root format.
  !****************************************************************************
  subroutine Root_write_event_footer(this)
    class(RootOutputFile), intent(in) :: this

    call rootaddevent(this%weight)

  end subroutine Root_write_event_footer


  !****************************************************************************
  !****s* EventOutput/Root_write_particle
  ! NAME
  ! subroutine Root_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in Root format.
  !****************************************************************************
  subroutine Root_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU

    class(RootOutputFile) :: this
    type(particle), intent(in) :: part

    real, dimension(1:3) :: pos

    integer :: KF

    KF = KFfromBUU(part)

    if (useProductionPos) then
       pos = getProductionPos(part)
    else
       pos = part%pos
    end if

    call rootaddparticle(KF, &
         part%mom(1),part%mom(2),part%mom(3), part%mom(0), &
         pos(1),pos(2),pos(3))

  end subroutine Root_write_particle

  !****************************************************************************
  !****s* EventOutput/Root_write_additionalInfo
  ! NAME
  ! subroutine Root_write_additionalInfo(this, iFE)
  ! PURPOSE
  ! add additional info about the event, depending on eventtype.
  !
  ! This routine tries to find additional information about the event.
  ! It tries routines for different event types, which only return
  ! some information, if it was really stored.
  !
  ! The following cases are handled:
  ! * For eventtype "HiLep", the following line is added:
  !     # 14 nu Q2 eps phiLepton Eventtype
  !   (14 is the magic number of "HiLepton")
  ! * For eventtype "neutrino", the following line is added:
  !     # 5 Eventtype Weight momLepIn(0:3) momLepOut(0:3) momNuc(0:3)
  !   (5 is the magic number for neutrino events)
  ! * For eventtype "heavyIon", the following line is added:
  !     # 1 b
  !   (1 is the magic number of "heavyIon", b is the impact parameter in fm)
  ! * For eventtype "hadron", the following line is added:
  !     # 300 b
  !   (300 is the magic number of "hadron", b is the impact parameter in fm)
  !****************************************************************************
  subroutine Root_write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use inputGeneral, only: eventType, numEnsembles, num_runs_SameEnergy
    use initNeutrino, only: process_ID, flavor_ID
    use nucleusDefinition
    use nucleus, only: getTarget
    use eventtypes, only: hiLepton, neutrino, heavyIon, hadron
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use FreezeoutAnalysis, only: getFreezeoutAnalysis_Pert
    use PIL_freezeout, only: PIL_freezeout_GET

    class(RootOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode
    type(tNucleus), pointer :: targetNuc

    ! variable have to have attribute 'save' due to rootadd... routine:
    real, save :: weight,nu,Q2,eps,phiL
    integer, save :: evtType, chrg_nuc
    real,dimension(0:3),save :: momLepIn, momLepOut, momBos, momNuc

!    integer :: history
!    logical :: escaped
    real,save :: tmp

    targetNuc => getTarget()

    select case (eventType)
    case (heavyIon)
       tmp = b_HI
       call rootadddouble(tmp,"b")
    case (hadron)
       tmp = b_had
       call rootadddouble(tmp,"b")
    case (neutrino)
      if (.not. present(iFE)) return
      if (NeutrinoProdInfo_Get(iFE,evtType,Weight,momLepIn,momLepOut,momBos,momNuc,chrg_nuc)) then
         call rootaddint(numEnsembles, "numEnsembles")
         call rootaddint(num_runs_SameEnergy, "numRuns")
         call rootaddint(targetNuc%mass, "nucleus_A")
         call rootaddint(targetNuc%charge, "nucleus_Z")
         call rootaddint(flavor_ID, "flavor_ID")
         call rootaddint(process_ID, "process_ID")
         call rootaddint(evtType, "evType")
         call rootadddouble(momLepIn(0), "lepIn_E")
         call rootadddouble(momLepIn(1), "lepIn_Px")
         call rootadddouble(momLepIn(2), "lepIn_Py")
         call rootadddouble(momLepIn(3), "lepIn_Pz")
         call rootadddouble(momLepOut(0), "lepOut_E")
         call rootadddouble(momLepOut(1), "lepOut_Px")
         call rootadddouble(momLepOut(2), "lepOut_Py")
         call rootadddouble(momLepOut(3), "lepOut_Pz")
         call rootadddouble(momNuc(0), "nuc_E")
         call rootadddouble(momNuc(1), "nuc_Px")
         call rootadddouble(momNuc(2), "nuc_Py")
         call rootadddouble(momNuc(3), "nuc_Pz")
         call rootaddint(chrg_nuc, "nuc_charge")

      end if
    case (hiLepton)
       if (.not. present(iFE)) return
       if (EventInfo_HiLep_Get(0,iFE,Weight,nu,Q2,eps,evtType,phi_Lepton=phiL)) then
          call rootaddint(evtType, "evType")
          call rootadddouble(nu, "nu")
          call rootadddouble(Q2, "Q2")
          call rootadddouble(eps, "eps")
          call rootadddouble(phiL, "phiL")

       end if
    end select

  end subroutine Root_write_additionalInfo



!******************************************************************************
!******************************************************************************
!******************************************************************************

  !****************************************************************************
  !****f* EventOutput/nuhepmc_proc_ID_from_buu
  ! NAME
  ! function nuhepmc_proc_ID_from_buu(this, buu_prod_id)
  ! PURPOSE
  ! Convert integer codes describing a GiBUU event into a NuHepMC proc_ID code
  !****************************************************************************
  function nuhepmc_proc_ID_from_buu(this, buu_prod_id) result(nuhepmc_proc_ID)
    use initNeutrino, only: process_ID !, flavor_ID

    class(NuHepMCOutputFile) :: this
    integer, intent(in) :: buu_prod_id !(1=N, 2=Delta, ..., 34=DIS, ...)
    integer :: nuhepmc_proc_ID, abs_process_ID
    integer, parameter :: EM = 1, CC = 2, NC = 3 ! positive process_ID values

    ! NuHepMC procID codes for CC events based on E.C.1
    integer, parameter :: QEL_ID = 200, MEC_ID = 300, RES_ID = 400
    integer, parameter :: SIS_ID = 500, DIS_ID = 600
    integer, parameter :: NC_OFFSET = 50

    ! Negative process_ID values in GiBUU represent antilepton channels.
    ! We take the absolute value here to get a general treatment.
    abs_process_ID = abs(process_ID)
    if ( abs(process_ID) .eq. EM ) then
      write(*,*) 'Output in NuHepMC format currently not implemented for EM processes. STOP!'
      stop
    end if

    ! Default to zero ("unknown"), then try to assign a different code
    nuhepmc_proc_ID = 0
    if (buu_prod_id .eq. 1) then
      nuhepmc_proc_ID = QEL_ID ! 200 = Quasielastic
    else if (buu_prod_id .ge. 2 .and. buu_prod_id .le. 31) then
      nuhepmc_proc_ID = RES_ID ! 400 = Delta RES
      if (buu_prod_id .gt. 2) then
        nuhepmc_proc_ID = nuhepmc_proc_ID + 1 ! 401 = Other RES
      end if
    else if (buu_prod_id .ge. 32 .and. buu_prod_id .le. 33) then
      nuhepmc_proc_ID = SIS_ID ! 500 = 1 pi neutron-background
      if (buu_prod_id .gt. 32) then
        nuhepmc_proc_ID = nuhepmc_proc_ID + 1 ! 501 = 1 pi proton-background
      end if
    else if (buu_prod_id .eq. 34) then
      nuhepmc_proc_ID = DIS_ID ! 600 = DIS
    else if (buu_prod_id .ge. 35 .and. buu_prod_id .le. 36) then
      nuhepmc_proc_ID = MEC_ID ! 300 = 2p2h QE
      if (buu_prod_id .gt. 35) then
        nuhepmc_proc_ID = nuhepmc_proc_ID + 1 ! 301 = 2p2h Delta
      end if
    else if (buu_prod_id .eq. 37) then
      nuhepmc_proc_ID = SIS_ID + 2 ! 502 = 2 pi background
    end if

    if (nuhepmc_proc_ID .ne. 0 .and. abs(process_ID) .eq. NC) then
      ! Add 50 to the codes above for NC channels
      nuhepmc_proc_ID = nuhepmc_proc_ID + NC_OFFSET
    end if

  end function nuhepmc_proc_ID_from_buu

  !****************************************************************************
  !****s* EventOutput/NuHepMC_open
  ! NAME
  ! subroutine NuHepMC_open(this, pert, nCall, nTimeStep)
  ! PURPOSE
  ! Open a file for event-information output in NuHepMC format.
  !****************************************************************************
  subroutine NuHepMC_open(this, pert, nCall, nTimeStep)
    use eventtypes, only: neutrino !, hiLepton, heavyIon, hadron
    use inputGeneral, only: eventType, numEnsembles, num_runs_SameEnergy
    use initNeutrino, only: process_ID, flavor_ID

    class(NuHepMCOutputFile) :: this
    logical, intent(in) :: pert
    integer, intent(in) :: nCall, nTimeStep

    character(len=40) :: fName  ! file name
    character(len=6)  :: buf1
    character(len=8)  :: buf2

    if (eventType .ne. neutrino) then
       write(*,*) 'Output in NuHepMC format currently implemented only for neutrino events. STOP!'
       stop
    end if

    if (pert) then
       buf1 = '.Pert.'
    else
       buf1 = '.Real.'
    end if
    write(buf2,'(I8.8)') nCall
    fName = 'EventOutput' // trim(buf1) // trim(buf2) // '.hepmc3'

    ! open file
    open(this%iFile, file=fName, status='unknown')
    rewind(this%iFile)

    write(this%iFile,'(A)') 'HepMC::Version 3.02.05'
    write(this%iFile,'(A)') 'HepMC::Asciiv3-START_EVENT_LISTING'
    write(this%iFile,'(A)') 'W CV'
    write(this%iFile,'(A)') 'T GiBUU\|2023.0\|'

    write(this%iFile,'(A17,1X,I0)') 'A GiBUU.Ensembles', numEnsembles
    write(this%iFile,'(A12,1X,I0)') 'A GiBUU.Runs', num_runs_SameEnergy
    write(this%iFile,'(A16,1X,I0)') 'A GiBUU.FlavorID', flavor_ID
    write(this%iFile,'(A17,1X,I0)') 'A GiBUU.ProcessID', process_ID

    write(this%iFile,'(A)') 'A NuHepMC.Citations.Generator.DOI 10.1016/j.physrep.2011.12.001 10.1088/1361-6471/ab3830'
    write(this%iFile,'(A)') 'A NuHepMC.Citations.Generator.arXiv 1106.1344 1904.11506 2308.16161'

    ! TODO: fix this
    write(this%iFile,'(A)') 'A NuHepMC.Conventions E.C.1 E.C.2 E.C.4 E.C.5 G.C.1 G.C.4 G.C.5 G.C.6'
    ! TODO: add this
    write(this%iFile,'(A)') 'A NuHepMC.FluxAveragedTotalCrossSection'
    write(this%iFile,'(A)') 'A NuHepMC.ParticleStatusIDs'

    ! Lookup table for NuHepMC procID codes
    ! These are produced from GiBUU event information via the helper function
    ! nuhepmc_proc_ID_from_buu defined above
    !
    ! 0 = Unknown
    ! 200 = Quasielastic
    ! 300 = 2p2h QE
    ! 301 = 2p2h Delta
    ! 400 = Delta RES
    ! 401 = Other RES
    ! 500 = 1 pi neutron-background
    ! 501 = 1 pi proton-background
    ! 502 = 2 pi background
    ! 600 = DIS
    ! Add 50 to the codes above for NC channels. EM codes not yet assigned.
    write(this%iFile,'(A)') 'A NuHepMC.ProcessIDs 200 250 300 301 350 351 400 401 450 451 500 501 550 551 600 650'

    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[0].Description Unrecognized interaction type'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[0].Name Unknown'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[200].Description Charged-current quasielastic'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[200].Name CCQE'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[250].Description Neutral-current quasielastic'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[250].Name NCQE'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[300].Description Charged-current 2p2h quasielastic'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[300].Name CC 2p2h QE'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[350].Description Neutral-current 2p2h quasielastic'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[350].Name NC 2p2h QE'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[301].Description Charged-current 2p2h Delta production'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[301].Name CC 2p2h Delta'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[351].Description Neutral-current 2p2h Delta production'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[351].Name NC 2p2h Delta'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[400].Description Charged-current resonant Delta production'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[400].Name CC RES Delta'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[450].Description Neutral-current resonant Delta production'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[450].Name NC RES Delta'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[401].Description Charged-current other resonance production'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[401].Name CC RES other'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[451].Description Neutral-current other resonance production'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[451].Name NC RES other'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[500].Description Charged-current single-pion neutron background'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[500].Name CC Bkgd-n'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[550].Description Neutral-current single-pion neutron background'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[550].Name NC Bkgd-n'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[501].Description Charged-current single-pion proton background'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[501].Name CC Bkgd-p'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[551].Description Neutral-current single-pion proton background'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[551].Name NC Bkgd-p'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[502].Description Charged-current two-pion background'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[502].Name CC Bkgd 2pi'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[552].Description Neutral-current two-pion background'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[552].Name NC Bkgd 2pi'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[600].Description Charged-current deep inelastic scattering'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[600].Name CC DIS'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[650].Description Neutral-current deep inelastic scattering'
    write(this%iFile,'(A)') 'A NuHepMC.ProcessInfo[650].Name NC DIS'

    ! TODO: adjust the cross-section units when EM channels are added to the
    ! code for NuHepMC-format output
    write(this%iFile,'(A)') 'A NuHepMC.Units.CrossSection.TargetScale PerTargetNucleon'
    write(this%iFile,'(A)') 'A NuHepMC.Units.CrossSection.Unit 1e-38 cm2'

    write(this%iFile,'(A)') 'A NuHepMC.Version.Major 0'
    write(this%iFile,'(A)') 'A NuHepMC.Version.Minor 9'
    write(this%iFile,'(A)') 'A NuHepMC.Version.Patch 0'
    write(this%iFile,'(A)') 'A NuHepMC.VertexStatusIDs'

  end subroutine NuHepMC_open


  !****************************************************************************
  !****s* EventOutput/NuHepMC_close
  ! NAME
  ! subroutine NuHepMC_close(this)
  ! PURPOSE
  ! Close a file after outputting event information in NuHepMC format.
  !****************************************************************************
  subroutine NuHepMC_close(this)
    class(NuHepMCOutputFile), intent(in) :: this

    write(this%iFile,'(A)') 'HepMC::Asciiv3-END_EVENT_LISTING'

  end subroutine NuHepMC_close


  !****************************************************************************
  !****s* EventOutput/NuHepMC_write_event_header
  ! NAME
  ! subroutine NuHepMC_write_event_header(this, nParts, nEvent, wgt, iFE)
  ! PURPOSE
  ! Write the header for an event in NuHepMC format,
  ! including the number of particles and the event weight.
  !****************************************************************************
  subroutine NuHepMC_write_event_header(this, nParts, nEvent, wgt, iFE)
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use initNeutrino, only: process_ID, flavor_ID

    class(NuHepMCOutputFile) :: this
    integer, intent(in) :: nParts            ! number of particles
    integer, intent(in) :: nEvent            ! number of current event
    real, intent(in), optional :: wgt        ! weight of event
    integer, intent(in), optional :: iFE     ! firstFvent

    real :: per_weight, nu_mass, hit_nuc_mass, lep_mass
    integer :: prod_id_buu, chrg_nuc, my_nuhepmc_proc_id
    integer :: neutrino_pdg_code, hit_nuc_pdg_code, lep_pdg_code
    real, dimension(0:3) :: momLepIn, momLepOut, momBos, momNuc

    this%weight = 1.0
    if (present(wgt)) this%weight=wgt

    write(this%iFile,'(A2,I0,1X,I0,A2)') 'E ', nEvent, nParts + 3, ' 1'
    write(this%iFile,'(A)') 'U GEV CM'
    write(this%iFile,'(A2,E28.22)') 'W ', this%weight
    write(this%iFile,'(A)') 'A 0 LabPos 0.000000 0.000000 0.000000 0.000000'

    neutrino_pdg_code = 0
    lep_pdg_code = 0
    select case (flavor_ID)
      case (1)
        neutrino_pdg_code = 12 ! nu_e
      case (2)
        neutrino_pdg_code = 14 ! nu_mu
      case (3)
        neutrino_pdg_code = 16 ! nu_tau
    end select

    ! Antileptons have negative process_ID values
    if ( process_ID .lt. 0 ) then
      neutrino_pdg_code = -1 * neutrino_pdg_code
    end if

    ! Assign the outgoing lepton PDG code based on whether this is a CC
    ! event (process_ID = ±2) or not
    if ( abs(process_ID) .eq. 2 ) then
       if (process_ID .eq. 2) then
         lep_pdg_code = neutrino_pdg_code - 1
       else
         lep_pdg_code = neutrino_pdg_code + 1
       end if
    else
       lep_pdg_code = neutrino_pdg_code
    end if

    !Get(iFE,... ?????
    if (NeutrinoProdInfo_Get(nEvent,prod_id_buu,per_weight,momLepIn,momLepOut,momBos,momNuc,chrg_nuc)) then

       my_nuhepmc_proc_id = this%get_proc_ID(prod_id_buu)
       write(this%iFile,'(A11,I0)') 'A 0 ProcID ', my_nuhepmc_proc_id
       write(this%iFile,'(A20,E28.22E2)') 'A 0 GiBUU.PerWeight ', per_weight

       !call rootaddint(targetNuc%mass, "nucleus_A")
       !call rootaddint(targetNuc%charge, "nucleus_Z")

       ! Add particle definition for projectile
       nu_mass = sqrt(max(0., momLepIn(0)**2 - momLepIn(1)**2 - momLepIn(2)**2 - momLepIn(3)**2))
       write(this%iFile,'(A6,I0,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,A2)') 'P 1 0 ', neutrino_pdg_code, momLepIn(1), momLepIn(2), momLepIn(3), momLepIn(0), nu_mass, ' 4'

       ! Add particle definition for the struck nucleon
       hit_nuc_mass = sqrt(max(0., momNuc(0)**2 - momNuc(1)**2 - momNuc(2)**2 - momNuc(3)**2))

       hit_nuc_pdg_code = 2112 ! n
       if (chrg_nuc .gt. 0) then
         hit_nuc_pdg_code = 2212 ! p
       end if

       write(this%iFile,'(A6,I0,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,A3)') 'P 2 0 ', hit_nuc_pdg_code, momNuc(1), momNuc(2), momNuc(3), momNuc(0), hit_nuc_mass, ' 20'

       ! Add fixed primary vertex definition
       write(this%iFile,'(A)') 'V -1 1 [1,2] @ 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00'

       lep_mass = sqrt(max(0., momLepOut(0)**2 - momLepOut(1)**2 - momLepOut(2)**2 - momLepOut(3)**2))
       write(this%iFile,'(A6,I0,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,A2)') 'P 3 -1 ', lep_pdg_code, momLepOut(1), momLepOut(2), momLepOut(3), momLepOut(0), lep_mass, ' 1'

       ! Update the total particle count now that we've added the first
       ! three in this subroutine
       this%particle_count = 3

    end if

  end subroutine NuHepMC_write_event_header


  !****************************************************************************
  !****s* EventOutput/NuHepMC_write_event_footer
  ! NAME
  ! subroutine NuHepMC_write_event_footer(this)
  ! PURPOSE
  ! Write the footer that closes a event in NuHepMC format.
  !****************************************************************************
  subroutine NuHepMC_write_event_footer(this)
    class(NuHepMCOutputFile), intent(in) :: this

    ! Do nothing

  end subroutine NuHepMC_write_event_footer


  !****************************************************************************
  !****s* EventOutput/NuHepMC_write_particle
  ! NAME
  ! subroutine NuHepMC_write_particle(iFile, part)
  ! PURPOSE
  ! Write a single particle in NuHepMC format.
  !****************************************************************************
  subroutine NuHepMC_write_particle(this, part)
    use particleDefinition
    use ID_translation, only: KFfromBUU

    class(NuHepMCOutputFile) :: this
    type(particle), intent(in) :: part

    ! TODO: reintroduce this information as needed with appropriate vertices
    !real, dimension(1:3) :: pos

    integer :: KF
    real :: part_mass

    KF = KFfromBUU(part)

    !if (useProductionPos) then
    !   pos = getProductionPos(part)
    !else
    !   pos = part%pos
    !end if

    ! Increment the particle count and write out the particle in the
    ! event record
    this%particle_count = this%particle_count + 1
    part_mass = sqrt(max(0., part%mom(0)**2 - part%mom(1)**2 - part%mom(2)**2 - part%mom(3)**2))

    write(this%iFile,'(A2,I0,A4,I0,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,1X,E23.16,1X,A2)') 'P ', this%particle_count, ' -1 ', KF, part%mom(1), part%mom(2), part%mom(3), part%mom(0), part_mass, ' 1'

  end subroutine NuHepMC_write_particle

  !****************************************************************************
  !****s* EventOutput/NuHepMC_write_additionalInfo
  ! NAME
  ! subroutine NuHepMC_write_additionalInfo(this, iFE)
  ! PURPOSE
  ! add additional info about the event, depending on eventtype.
  !****************************************************************************
  subroutine NuHepMC_write_additionalInfo(this, iFE, pNode)
    use particlePointerListDefinition
    use EventInfo_HiLep, only: EventInfo_HiLep_Get
    use neutrinoProdInfo, only: NeutrinoProdInfo_Get
    use inputGeneral, only: eventType, numEnsembles, num_runs_SameEnergy
    use initNeutrino, only: process_ID, flavor_ID
    use nucleusDefinition
    use nucleus, only: getTarget
    use eventtypes, only: hiLepton, neutrino, heavyIon, hadron
    use initHeavyIon, only: b_HI => b
    use initHadron, only: b_had => b
    use FreezeoutAnalysis, only: getFreezeoutAnalysis_Pert
    use PIL_freezeout, only: PIL_freezeout_GET

    class(NuHepMCOutputFile), intent(in) :: this
    integer, intent(in), optional :: iFE
    type(tParticleListNode), pointer, optional :: pNode
    !type(tNucleus), pointer :: targetNuc

    !real :: weight,nu,Q2,eps,phiL
    !integer :: evtType, chrg_nuc
    !real,dimension(0:3) :: momLepIn, momLepOut, momBos, momNuc

    !targetNuc => getTarget()

    !select case (eventType)
    !case (neutrino)
    !  if (.not. present(iFE)) return
    !  if (NeutrinoProdInfo_Get(iFE,evtType,Weight,momLepIn,momLepOut,momBos,momNuc,chrg_nuc)) then
    !     !call rootaddint(targetNuc%mass, "nucleus_A")
    !     !call rootaddint(targetNuc%charge, "nucleus_Z")
    !  end if
    !end select

  end subroutine NuHepMC_write_additionalInfo

!******************************************************************************
!******************************************************************************
!******************************************************************************

  !****************************************************************************
  !****s* EventOutput/write_real
  ! NAME
  ! subroutine write_real(this, parts)
  ! PURPOSE
  ! Do the actual printout for real particles.
  ! NOTES
  ! For the case of real particles, one event simply corresponds to one
  ! ensemble.
  !****************************************************************************
  subroutine write_real(this, parts)
    use particleDefinition
    use IdTable, only: EOV, NOP
    use inputGeneral, only: current_run_number

    class(EventOutputFile), intent(in) :: this
    type(particle), intent(in), dimension(:,:), target :: Parts

    integer :: iEns, iPart, NUP, nEnsembles, iEvent
    integer, dimension(:), allocatable :: nParts

    allocate(nParts(1:size(Parts,dim=1)))
    nParts=0

    ! count particles per ensemble
    do iEns = 1,size(Parts,dim=1)
      do iPart = 1,size(Parts,dim=2)
        if (Parts(iEns,iPart)%ID==EOV) exit
        if (Parts(iEns,iPart)%ID==NOP) cycle
        nParts(iEns) = nParts(iEns) + 1
      end do
    end do

    ! Loop over all events and write them to file:
    nEnsembles = size(Parts,dim=1)
    do iEns = 1, nEnsembles
       NUP = nParts(iEns) ! number of particles
       if (NUP == 0) cycle

       iEvent = (current_run_number-1)*nEnsembles + iEns
       call this%write_event_header(NUP, iEvent)  ! no weight here

       do iPart = 1,size(Parts,dim=2)
         if (Parts(iEns,iPart)%ID==EOV) exit
         if (Parts(iEns,iPart)%ID==NOP) cycle
         call this%write_particle(Parts(iEns,iPart))
       end do

       call this%write_additionalInfo()

       call this%write_event_footer()
    end do

  end subroutine write_real


  !****************************************************************************
  !****s* EventOutput/ValueListAllocate
  ! NAME
  ! subroutine ValueListAllocate(n1)
  ! PURPOSE
  ! Do the allocation stuff for the Particle Info List.
  !****************************************************************************
  subroutine ValueListAllocate(ValueList, n1)
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_INIT

    type(tParticleList), allocatable :: ValueList(:)
    integer, intent(in) :: n1

    integer :: n0,i
    type(tParticleList),allocatable :: L0(:)

    if (.not.allocated(ValueList)) then
        allocate(ValueList(n1))
        do i=1,n1
          call ParticleList_INIT(ValueList(i))
        end do
        return
    end if

    n0 = size(ValueList)            ! old size

    allocate(L0(n0))
    do i=1,n0
        L0(i)%first => ValueList(i)%first
        L0(i)%last  => ValueList(i)%last
        L0(i)%nEntries = ValueList(i)%nEntries
    end do
    deallocate(ValueList)
    allocate(ValueList(n1))
    do i=1,n0
        ValueList(i)%first => L0(i)%first
        ValueList(i)%last  => L0(i)%last
        ValueList(i)%nEntries = L0(i)%nEntries
    end do
    do i=n0+1,n1
        call ParticleList_INIT(ValueList(i))
    end do
    deallocate(L0)

  end subroutine ValueListAllocate


  !****************************************************************************
  !****s* EventOutput/write_pert
  ! NAME
  ! subroutine write_pert(this, parts)
  ! PURPOSE
  ! Do the actual printout for perturbative particles.
  ! NOTES
  ! We have to sort the particles according their "firstevent" field.
  ! Therefore we allocate an array of "tParticleList". Unfortunately we can
  ! not use the "firstevent" entry directly as array index, since this
  ! is not starting with 1 and continously increasing for all kind
  ! of eventtypes. Therefore we (ab)use the module "PILIndex", which
  ! implements methods of "indexing". (We do not use the possibility of
  ! reallocating as provided by the module "PILIndex".)
  !****************************************************************************
  subroutine write_pert(this, parts)
    use particleDefinition
    use particlePointerListDefinition
    use particlePointerList, only: ParticleList_APPEND, ParticleList_CLEAR
    use PILIndex, only: tIndexList, PILIndex_PUT
    use inputGeneral, only: current_run_number

    class(EventOutputFile), intent(in) :: this
    type(particle), intent(in), dimension(:,:), target :: Parts

    type(tIndexList), save :: IndexList
    type(tParticleList), allocatable, save :: ValueList(:)

    integer :: i,iEns,iPart, iFE,iiFE, NUP, iEvent

    type(particle), pointer :: pPart
    type(tParticleListNode),Pointer  :: pNode

    ! Clean up the arrays:
    if (allocated(ValueList)) then
       do i=1,size(ValueList)
          call ParticleList_CLEAR(ValueList(i))
       end do
    end if
    IndexList%nEntry = 0

    ! Loop over all particles and group them according their first event:
    do iEns = 1,size(Parts,dim=1)
       PartLoop:do iPart = 1,size(Parts,dim=2)
          pPart => Parts(iEns,iPart)
          if (pPart%ID <= 0) cycle PartLoop

          iFE = pPart%firstEvent
          if (iFE.eq.0) cycle PartLoop ! particle did not interact !
          iiFE = PILIndex_PUT(IndexList, iFE, 'LesHouches')
          if (iiFE>0) then
             call ParticleList_APPEND(ValueList(iiFE), pPart)
          else
             call ValueListAllocate(ValueList, size(IndexList%PartNumber))
             call ParticleList_APPEND(ValueList(-iiFE), pPart)
          end if

       end do PartLoop
    end do

    ! Loop over all events and write them to file:
    if (.not.allocated(ValueList)) return

    do iiFE=1,size(ValueList)
       NUP = ValueList(iiFE)%nEntries ! number of particles
       if (NUP == 0) cycle

       pNode => ValueList(iiFE)%first
       iFE = pNode%V%firstEvent
       iEvent = (current_run_number-1)*size(ValueList) + iiFE
       call this%write_event_header(NUP, iEvent, pNode%V%perweight, iFE)

       do
          if (.not. associated(pNode)) exit
          pPart => pNode%V
          call this%write_particle(pPart)
          pNode => pNode%next
       end do

       call this%write_additionalInfo(iFE, ValueList(iiFE)%first)

       call this%write_event_footer()
    end do

  end subroutine write_pert


end module EventOutput
