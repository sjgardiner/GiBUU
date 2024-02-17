!******************************************************************************
!****m* /eventtypes
! NAME
! module eventtypes
! PURPOSE
! Declare constants which define the different event classes:
! * 000 -- Elementary collisions (hadron+hadron, real particles)
! * 001 -- Heavy-Ion collisions (nucleus+nucleus, real particles)
! * 002 -- Pion-induced reactions (low energy, perturbative particles)
! * 003 -- Real-photon-induced reactions
! * 004 -- Lepton (i.e. virtual photon) induced reactions (low energy)
! * 005 -- Neutrino-induced reactions
! * 012 -- Pion- and nucleon-induced reactios (high energy, perturbative particles)
! * 014 -- Lepton (i.e. virtual photon) induced reactions (high energy)
! * 022 -- Transport of an external hadronic source from a data file
! * 031 -- Nucleons in a [box of nucleons] (continuous boundary conditions)
! * 032 -- Pions in a [box of nucleons] (continuous boundary conditions)
! * 033 -- Deltas in a [box of nucleons] (continuous boundary conditions)
! * 041 -- Box of particles
! * 100 -- Groundstate calculation
! * 200 -- Simple transport of a given particle
! * 300 -- Hadron-induced reactions (hadron+nucleus, real particles)
!******************************************************************************
module eventtypes

  implicit none

  integer, parameter :: elementary             = 0
  integer, parameter :: HeavyIon               = 1

  integer, parameter :: LoPion                 = 2
  integer, parameter :: RealPhoton             = 3
  integer, parameter :: LoLepton               = 4

  integer, parameter :: Neutrino               = 5

  integer, parameter :: HiPion                 = 12
  integer, parameter :: HiLepton               = 14

  integer, parameter :: ExternalSource         = 22

  integer, parameter :: InABox                 = 31
  integer, parameter :: InABox_pion            = 32
  integer, parameter :: InABox_delta           = 33

  integer, parameter :: Box                    = 41

  integer, parameter :: groundState            = 100
  integer, parameter :: transportGivenParticle = 200
  integer, parameter :: hadron                 = 300

contains

  character(22) function cEventType(iEventType)
    integer, intent(in) :: iEventType

    character(22), dimension(0:16), parameter :: cName = (/ &
         'elementary            ', & ! 0
         'HeavyIon              ', & ! 1
         'LoPion                ', & ! 2
         'RealPhoton            ', & ! 3
         'LoLepton              ', & ! 4
         'Neutrino              ', & ! 5
         'HiPion                ', & ! 12   6
         'HiLepton              ', & ! 14   7
         'ExternalSource        ', & ! 22   8
         'InABox                ', & ! 31   9
         'InABox_pion           ', & ! 32  10
         'InABox_delta          ', & ! 33  11
         'Box                   ', & ! 41  12
         'groundState           ', & ! 100 13
         'transportGivenParticle', & ! 200 14
         'hadron                ', & ! 300 15
         '*** unknown ***       ' /) ! xxx 16

    select case(iEventType)
    case (0:5)
       cEventType = cName(iEventType)
    case (12)
       cEventType = cName(6)
    case (14)
       cEventType = cName(7)
    case (22)
       cEventType = cName(8)
    case (31:33)
       cEventType = cName(iEventType-22)
    case (41)
       cEventType = cName(12)
    case (100)
       cEventType = cName(13)
    case (200)
       cEventType = cName(14)
    case (300)
       cEventType = cName(15)
    case default
       cEventType = cName(16)
    end select

  end function cEventType

end module eventtypes
