!******************************************************************************
!****m* /neutrinoAnalysis
! NAME
! module neutrinoAnalysis
!
! PURPOSE
! * This module does the analysis of the output of neutrino induced processes.
!******************************************************************************
module neutrinoAnalysis

  use AnaEvent
  use hist
  use hist2D
  use initNeutrino, only: max_Hist, includeHist, K2Hist, &
      get_init_namelist,OscLength,isOsc,process_ID
  use AnaEventDefinition
  use CALLSTACK
  use Oscill_Parameters
  use initNeutrino, only: cost_min,cost_max,delta_cost, &
       Elept_min,Elept_max,delta_Elept
  use constants, only: mN
  use PreEvListDefinition

  implicit none
  private

  Public:: neutrino_Analyze, cleanUp



  !****************************************************************************
  !****g* neutrinoAnalysis/detailed_diff_output
  ! SOURCE
  logical, save :: detailed_diff_output = .false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forPion
  ! SOURCE
  logical, save :: forPion =.true.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forEta
  ! SOURCE
  logical, save :: forEta =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forKaon
  ! SOURCE
  logical, save :: forKaon =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forKaonBar
  ! SOURCE
  logical, save :: forKaonBar =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDmeson
  ! SOURCE
  logical, save :: forDmeson =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDbar
  ! SOURCE
  logical, save :: forDbar =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDs_plus
  ! SOURCE
  logical, save :: forDs_plus =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forDs_minus
  ! SOURCE
  logical, save :: forDs_minus =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forNucleon
  ! SOURCE
  logical, save :: forNucleon =.true.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forLambda
  ! SOURCE
  logical, save :: forLambda =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forSigmaResonance
  ! SOURCE
  logical, save :: forSigmaResonance =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forXi
  ! SOURCE
  logical, save :: forXi =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/forOmegaResonance
  ! SOURCE
  logical, save :: forOmegaResonance =.false.
  ! PURPOSE
  ! If .true. then also the detailed output of differential cross sections is
  ! produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/include_W_dist
  ! SOURCE
  !
  logical, save :: include_W_dist=.false.
  ! PURPOSE
  ! If .true. then the invariant mass distributions for events with 1 pion and
  ! 1 nucleon in the final state are produced
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dW_Npi
  ! SOURCE
  !
  real, save  :: dW_Npi=0.02
  ! PURPOSE
  ! for dsigma/d(InvariantMass);
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmin_Npi
  ! SOURCE
  !
  real, save  :: Wmin_Npi=1.08
  ! PURPOSE
  ! for dsigma/d(InvariantMass);
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmax_Npi
  ! SOURCE
  !
  real, save  :: Wmax_Npi=1.6
  ! PURPOSE
  ! for dsigma/d(InvariantMass);
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dW_mupi
  ! SOURCE
  !
  real, save  :: dW_mupi=0.04
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmin_mupi
  ! SOURCE
  !
  real, save  :: Wmin_mupi=0.24
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmax_mupi
  ! SOURCE
  !
  real, save  :: Wmax_mupi=1.2
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dW_muN
  ! SOURCE
  !
  real, save  :: dW_muN=0.04
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmin_muN
  ! SOURCE
  !
  real, save  :: Wmin_muN=1.04
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Wmax_muN
  ! SOURCE
  !
  real, save  :: Wmax_muN=2.12
  ! PURPOSE
  ! only work if include_W_dist is .true.
  ! set the min, max and steps for various W-distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/numin
  ! SOURCE
  !
  real, save  :: numin=0.
  ! PURPOSE
  ! for calorimetric analysis: values for transferred energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/numax
  ! SOURCE
  !
  real, save  :: numax=10.0
  ! PURPOSE
  ! for calorimetric analysis: values for transferred energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/nubin
  ! SOURCE
  !
  real, save  :: nubin=0.1
  ! PURPOSE
  ! for calorimetric analysis: values for transferred energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Enumin
  ! SOURCE
  !
  real, save  :: Enumin=0.
  ! PURPOSE
  ! for calorimetric analysis: values for neutrino energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Enumax
  ! SOURCE
  !
  real, save  :: Enumax=10.0
  ! PURPOSE
  ! for calorimetric analysis: values for neutrino energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Enubin
  ! SOURCE
  !
  real, save  :: Enubin=0.1
  ! PURPOSE
  ! for calorimetric analysis: values for neutrino energy;
  ! only work if calorimetric_analysis is .true.
  ! set the min, max and bins for nu distributions
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/EkinMin
  ! SOURCE
  real, save :: EkinMin=0.
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! minimal kinetic energy for dsigma/dEkin for hadrons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/EkinMax
  ! SOURCE
  real, save :: EkinMax=2.
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Maximal kinetic energy for dsigma/dEkin for hadrons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/dEkin
  ! SOURCE
  real, save :: dEkin=0.01
  ! PURPOSE
  ! if detailed_diff_output is TRUE:
  ! Delta(eKin) for dsigma/dEKin  for hadrons
  !****************************************************************************

   !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_nucleon
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_nucleon=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  ! lower detection threshold for nucleon kinetic energies
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_chargedpion
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_chargedpion=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_neutralpion
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_neutralpion=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/kineticEnergyDetectionThreshold_lepton
  ! SOURCE
  !
  real, save :: kineticEnergyDetectionThreshold_lepton=0.0
  ! PURPOSE
  ! kineticEnergyDetectionThreshold
  ! only lepton kinetic energies above this threshold can be detected
  ! This cut affects *all* events, not just the outgoing lepton!
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_nucleon
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_nucleon=180.0
  ! PURPOSE
  ! nucleon angles up to this value can be detected
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_chargedpion
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_chargedpion=180.0
  ! PURPOSE
  ! charged pion angles up to this value can be detected
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_neutralpion
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_neutralpion=180.0
  ! PURPOSE
  ! neutral pion angles angles up to this value can be detected
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/AngleUpperDetectionThresholdDegrees_lepton
  ! SOURCE
  !
  real, save :: AngleUpperDetectionThresholdDegrees_lepton=180.0
  ! PURPOSE
  ! lepton angles up to this value can be detected
  ! This cut affects *all* events, not just the outgoing lepton!
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/radialScale
  ! SOURCE
  !
  real, save :: radialScale=1.5
  ! PURPOSE
  ! If radial position of nucleon < radialScale*target radius,
  ! then the nucleon is assumed to be bound
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/applyCuts
  ! SOURCE
  !
  integer, save :: applyCuts = 0
  ! PURPOSE
  ! This parameter encodes 'binary', which cuts should be applied:
  ! * 1: lepton_acceptance
  ! * 2: isBound
  ! * 4: isBelowThreshold
  !
  ! Instead of having three indipendent flags (with values=0 or 1) as
  ! e.g. labelled 'doLepton', 'doIsBound', 'doBelowThr', applyCuts combines
  ! them formally into one number as
  !   applyCuts = 1*doLepton + 2*doIsBound + 4*doBelowThr
  ! So by setting any number between 0 and 7, one can individually switch on
  ! and off each of these cuts.
  !
  ! 'lepton_acceptance' uses the input parameters:
  ! * kineticEnergyDetectionThreshold_lepton (for all kind of outgoing leptons)
  ! * AngleUpperDetectionThresholdDegrees_lepton
  !
  ! 'isBound' tests, whether kinetic energy plus potential is <0
  !
  ! 'isBelowThreshold' uses the input parameters:
  ! * kineticEnergyDetectionThreshold_lepton (only for muons)
  ! * AngleUpperDetectionThresholdDegrees_lepton (only for muons)
  ! * kineticEnergyDetectionThreshold_nucleon (only for nucleons)
  ! * AngleUpperDetectionThresholdDegrees_nucleon (only for nucleons)
  ! * kineticEnergyDetectionThreshold_chargedpion (only for charged pion)
  ! * AngleUpperDetectionThresholdDegrees_chargedpion (only for charged pion)
  ! * kineticEnergyDetectionThreshold_neutralpion (only for neutral pions)
  ! * AngleUpperDetectionThresholdDegrees_neutralpion (only for neutral pions)
  !
  ! Some examples:
  ! * To generate full inclusive output, set the value applyCuts=0
  ! * To generate output where bound nucleons are dropped, set applyCuts=2
  ! * To generate output with specific experimental cuts for the outgoing
  !   hadrons, set applyCuts=4 or applyCuts=6 and set the corresponding
  !   threshold parameters accordingly.
  ! * If in the experiment also cuts on the outgoing lepton are used, set
  !   applyCuts=7 and set the corresponding threshold parameters accordingly.
  !
  !
  ! NOTES
  ! These cuts affects the output into the file "FinalEvents.dat".
  !
  ! The cut 'lepton_acceptance' (de-)selects the full event, while the other
  ! two cuts only decide whether a specific particle is accepted or not.
  !
  ! The kinetic energy of a bound nucleon is < 0. Therefore using the
  ! default value kineticEnergyDetectionThreshold_nucleon=0.0 also tests,
  ! whether the particle is bound or not. Set the parameter to a large negative
  ! value to become ineffective.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Fissum_analysis
  ! SOURCE
  logical, save ::  Fissum_analysis=.false.
  ! PURPOSE
  ! do analysis with cuts as needed for Fig 25 in
  ! Fissum et al, PRC 70, 034606 (2004)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/ZeroPion_analysis
  ! SOURCE
  logical, save ::  ZeroPion_analysis=.false.
  ! PURPOSE
  ! produce output of xsec for various final states with 0 pions and 2 pions
  ! see file see sigma_0pions.dat  for the list of the final states
  !
  ! see files neutrino_0pions.dat,  neutrino_0pions_QE.dat,
  ! neutrino_0pions_Delta.dat, ... for output
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/calorimetric_analysis
  ! SOURCE
  logical, save ::  calorimetric_analysis=.false.
  ! PURPOSE
  ! do calorimetric energy-transfer and neutrino-energy reconstruction
  ! (for each QE, Delta, ...)  as in the MINOS experiment
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/reconstruct_neutrino_energy
  ! SOURCE
  logical, save ::  reconstruct_neutrino_energy=.false.
  ! PURPOSE
  ! reconstruct neutrino energy for final state in "specificEvent_analysis"
  ! NOTES
  ! .true. must be combined with specificEvent_analysis=.true. and
  ! at least one specific event .true.
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/specificEvent_analysis
  ! SOURCE
  logical, save ::  specificEvent_analysis=.false.
  ! PURPOSE
  ! do analysis for specific final states
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/no_pi
  ! SOURCE
  logical, save ::  no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=1,
  ! no_pi (for example, for QE-like MiniBooNE)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/p_Xn_no_pi
  ! SOURCE
  logical, save ::  p_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=2
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/piplus
  ! SOURCE
  logical, save ::  piplus=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=3, 1 pi+ X nucleons
  ! mesons of other flavor
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/piplus_MULTI
  ! SOURCE
  logical, save ::  piplus_MULTI=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=4
  ! >=1 pi+  X other pions (incl pi+) X nucleons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pi0
  ! SOURCE
  logical, save ::  pi0=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=5,
  ! 1 pi0 X nucleons, plus mesons of other flavor
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pi0_MULTI
  ! SOURCE
  logical, save ::  pi0_MULTI=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=6,
  ! >=1 pi0  X other pions X nucleons, (pi0 K2K)
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/piminus
  ! SOURCE
  logical, save ::  piminus=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=7
  ! 1 pi-  X other pions X nucleons
  !****************************************************************************

  !***************************************************************************
  !****g* neutrinoAnalysis/piminus_MULTI
  ! SOURCE
  logical, save ::  piminus_MULTI=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=8
  ! >=1 pi-  X other pions X nucleons
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pp_no_pi
  ! SOURCE
  logical, save ::  pp_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states:  specificEvent=9
  ! 2 protons, X neutrons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pn_no_pi
  ! SOURCE
  logical, save ::  pn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=10
  ! 1 neutron, 1 proton, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/nn_no_pi
  ! SOURCE
  logical, save ::  nn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=11
  ! 2 neutrons, X protons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pp_Xn_no_pi
  ! SOURCE
  logical, save ::  pp_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=12
  ! 2 protons, X neutrons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/nn_Xp_no_pi
  ! SOURCE
  logical, save ::  nn_Xp_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=13
  ! 2 neutrons, X protons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/ppp_Xn_no_pi
  ! SOURCE
  logical, save ::  ppp_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=14
  ! 3 protons, X neutrons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pppp_Xn_no_pi
  ! SOURCE
  logical, save ::  pppp_Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=15
  ! 4 protons, X neutrons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/p_no_pi
  ! SOURCE
  logical, save ::  p_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=16
  ! 1 proton, 0 neutron, 0 pion
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/n_no_pi
  ! SOURCE
  logical, save ::  n_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=17
  ! 1 neutron, 0 proton, 0 pion
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Xn_no_pi
  ! SOURCE
  logical, save ::  Xn_no_pi=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=18,
  ! 0 proton, X neutrons, 0 pions
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/excl_hadron
  ! SOURCE
  logical, save :: excl_hadron=.false.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=19,20,21
  ! exclusive 1 pion, no other pions or other mesons of different flavor
  ! There could be still other mesons which are heavier than the D,
  ! Such events (very rare at DUNE energies) could be counted as exclusive
  ! single-meson cross section.
  ! This could be cured by extending the list of stable mesons
  !****************************************************************************

  logical, save :: excl_pi0=.false.
  logical, save :: excl_piplus=.false.
  logical, save :: excl_piminus=.false.


  !****************************************************************************
  !****g* neutrinoAnalysis/full_incl
  ! SOURCE
  logical, save :: full_incl=.true.
  ! PURPOSE
  ! do analysis for specific final states: specificEvent=22
  ! fully inclusive event, all hadrons in final state
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/QEp
  ! SOURCE
  logical ::  QEp=.false.
  ! PURPOSE
  ! if .true,
  ! do analysis for specific analysis for QE-like event with 1 mu, 0 pi, X p
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/pcut
  ! SOURCE
  real, save ::  pcut = 0.0
  ! PURPOSE
  ! determines lower acceptance cut for outgoing protons
  !
  ! values can be changed in the namelist neutrinoAnalysis via the variable
  ! kineticEnergyDetectionThreshold_nucleon
  !****************************************************************************

  integer, parameter ::  max_SpeEvent=22
  ! maximum number of special Events for which detailed analysis is being done
  ! special events in file  'includeSpeEvent'

  !****************************************************************************
  !****g* neutrinoAnalysis/binsizeQ2
  ! SOURCE
  real, save    ::  binsizeQ2=0.01
  ! PURPOSE
  ! do analysis for specific final states:
  ! binning for reconstruction of Q2 and Enu
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/binsizeEnu
  ! SOURCE
  real, save    ::  binsizeEnu=0.02
  ! PURPOSE
  ! do analysis for specific final states:
  ! binning for reconstruction of Q2 and Enu
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/maxQ2
  ! SOURCE
  real, save    ::  maxQ2=5.0
  ! PURPOSE
  ! do analysis for specific final states:
  ! max values for reconstruction of Q2 and Enu
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/maxEnu
  ! SOURCE
  real, save    ::  maxEnu=5.0
  ! PURPOSE
  ! do analysis for specific final states:
  ! max values for reconstruction of Q2 and Enu
  !****************************************************************************


  !****************************************************************************
  !****g* neutrinoAnalysis/outputEvents
  ! SOURCE
  !
  logical, save  :: outputEvents = .false.
  ! PURPOSE
  ! If .true. then all events are written to the file 'FinalEvents.dat'.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/Xsection_analysis
  ! SOURCE
  !
  logical, save  :: Xsection_analysis = .false.
  ! PURPOSE
  ! If .true. then files "..._total_Xsection_..."  and "..._dSigmadEkin_..."
  ! are printed.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/doPipe
  ! SOURCE
  !
  logical, save  :: doPipe = .false.
  ! PURPOSE
  ! If .true. then events are not written to the file 'FinalEvents.dat', but
  ! insted written into a named pipe (fifo) with the name fileNamePipe.
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/fileNamePipe
  ! SOURCE
  character(len=1000), save :: fileNamePipe = ""
  ! PURPOSE
  ! name of the pipe to be used
  !****************************************************************************

  !****************************************************************************
  !****g* neutrinoAnalysis/DoOutChannels
  ! SOURCE
  !
  logical, save :: DoOutChannels = .false.
  ! PURPOSE
  ! switch on/off: reporting of all final state channels
  !****************************************************************************


  logical, save :: initflag=.true.


  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists,dE_hists_QE, &
       & dE_hists_Delta, &
       & dE_hists_highRes,dE_hists_gen0,dE_hists_1piBG,dE_hists_2piBG,dE_hists_DIS,dE_hists_2p2h
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_Multi,dE_hists_QE_Multi, &
       & dE_hists_Delta_Multi,dE_hists_highRes_Multi,dE_hists_gen0_Multi,dE_hists_1piBG_Multi,&
       & dE_hists_2piBG_Multi,dE_hists_DIS_Multi,dE_hists_2p2h_Multi
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_1X,dE_hists_QE_1X, &
       & dE_hists_Delta_1X,dE_hists_highRes_1X,dE_hists_gen0_1X,dE_hists_1piBG_1X, &
       & dE_hists_2piBG_1X,dE_hists_DIS_1X, dE_hists_2p2h_1X
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_2X,dE_hists_QE_2X, &
       & dE_hists_Delta_2X,dE_hists_highRes_2X,dE_hists_gen0_2X,dE_hists_1piBG_2X, &
       & dE_hists_2piBG_2X,dE_hists_DIS_2X,dE_hists_2p2h_2X

  ! used for ZeroPion_analysis=.true.  for kinetic energy distributions of nucleons
  ! in events with 0 pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_0pions, &
       & dE_hists_QE_0pions, &
       & dE_hists_Delta_0pions,dE_hists_highRes_0pions, dE_hists_2p2h_0pions, &
       & dE_hists_1piBG_0pions, dE_hists_2piBG_0pions, dE_hists_DIS_0pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_Multi_0pions, &
       & dE_hists_QE_Multi_0pions,dE_hists_Delta_Multi_0pions, dE_hists_highRes_Multi_0pions, &
       & dE_hists_2p2h_Multi_0pions,dE_hists_1piBG_Multi_0pions, dE_hists_2piBG_Multi_0pions, &
       & dE_hists_DIS_Multi_0pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_1X_0pions,&
       & dE_hists_QE_1X_0pions, dE_hists_Delta_1X_0pions, &
       & dE_hists_highRes_1X_0pions, dE_hists_2p2h_1X_0pions, dE_hists_1piBG_1x_0pions, &
       & dE_hists_2piBG_1x_0pions,dE_hists_DIS_1x_0pions
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dE_hists_2X_0pions, &
       & dE_hists_QE_2X_0pions, dE_hists_Delta_2X_0pions, &
       & dE_hists_highRes_2X_0pions, dE_hists_2p2h_2X_0pions, dE_hists_1piBG_2x_0pions, &
       & dE_hists_2piBG_2x_0pions,dE_hists_DIS_2x_0pions


  ! used for include_W_dist=.true.
  type(histogram),save,dimension(-1:1) :: dW_Npi_hists, dW_mupi_hists, dW_muN_hists

  ! used for calorimetric_analysis=.true.
  type(histogram2D), dimension(0:max_Hist), save :: Ehad_versusNu_hist, Enurestored_versusEnu
  type(histogram), dimension(0:max_Hist), save :: dSigdNu_hist, dSigdEhad_hist, dSigdEnu,&
      & dSigdEnurestored


  ! used for specificEvent_analysis=.true.
  type(histogram),save, dimension(1:max_SpeEvent,0:max_Hist) :: dEnu_hist, dElepton_hist, &
      & dcoslepton_hist, dQ2lepton_hist, dQ2plepton_hist

  type(histogram2D), save, dimension(1:max_SpeEvent,0:max_Hist) :: dSigmaMC_EprimeCost_0pi,dSigmaMC_pLpT_0pi


  type(histogram),save,dimension(1:numStableParts,-2:2) :: dTheta_hists,  &
       & dPhi_hists, dTheta_hists_QE,dPhi_hists_QE, &
       & dTheta_hists_Delta,dPhi_hists_Delta,dTheta_hists_highRes,  &
       & dPhi_hists_highRes, dTheta_hists_gen0,dPhi_hists_gen0, &
       & dTheta_hists_1piBG,dPhi_hists_1piBG,dTheta_hists_2piBG,dPhi_hists_2piBG

  type(histogram2D),save,dimension(1:numStableParts,-2:2) :: dOmega_hists, &
       & dOmega_hists_QE,dOmega_hists_Delta, &
       & dOmega_hists_highRes,dOmega_hists_gen0,dOmega_hists_1piBG,  &
       & dOmega_hists_2piBG

  type(histogram),save,dimension(1:numStableParts,-2:2) :: dEcostheta_hists, &
       & dEcostheta_hists_QE,dEcostheta_hists_Delta, &
       & dEcostheta_hists_highRes, dEcostheta_hists_gen0,dEcostheta_hists_1piBG, &
       & dEcostheta_hists_2piBG,dEcostheta_hists_DIS
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dEcostheta_hists_MULTI, &
       & dEcostheta_hists_QE_MULTI, &
       & dEcostheta_hists_Delta_MULTI, dEcostheta_hists_highRes_MULTI, &
       & dEcostheta_hists_gen0_MULTI, &
       & dEcostheta_hists_1piBG_MULTI,dEcostheta_hists_2piBG_MULTI,  &
       & dEcostheta_hists_DIS_MULTI
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dcostheta_hists, &
       & dcostheta_hists_QE, dcostheta_hists_Delta, &
       & dcostheta_hists_highRes, dcostheta_hists_gen0,dcostheta_hists_1piBG, &
       & dcostheta_hists_2piBG,dcostheta_hists_DIS
  type(histogram),save,dimension(1:numStableParts,-2:2) :: dcostheta_hists_MULTI, &
       & dcostheta_hists_QE_MULTI, &
       & dcostheta_hists_Delta_MULTI, dcostheta_hists_highRes_MULTI, &
       & dcostheta_hists_gen0_MULTI,&
       & dcostheta_hists_1piBG_MULTI, dcostheta_hists_2piBG_MULTI, &
       & dcostheta_hists_DIS_MULTI

  type(histogram),save :: hNucleonVacuumMass

  logical :: initHists
  logical, dimension(1:max_speEvent), save :: includeSpeEvent

  type(tPreEvListEntry), save :: preEv
  type(tPreEvList), save :: ListPreEv


contains





  !****************************************************************************
  !****s* neutrinoAnalysis/readinput
  ! NAME
  ! subroutine readinput
  ! INPUTS
  ! NONE
  ! OUTPUT
  ! NONE
  ! PURPOSE
  ! This subroutine reads the namelist "neutrinoAnalysis".
  ! Only called once to initialize the module.
  !****************************************************************************
  subroutine readinput
    use output
    use inputGeneral, only: ExpandPath

    integer :: IOS

    !**************************************************************************
    !****n* neutrinoAnalysis/NeutrinoAnalysis
    ! NAME
    ! NAMELIST /NeutrinoAnalysis/
    ! PURPOSE
    ! This namelist includes:
    ! * detailed_diff_output
    ! * include_W_dist
    ! * kineticEnergyDetectionThreshold_lepton
    ! * AngleUpperDetectionThresholdDegrees_lepton
    ! * kineticEnergyDetectionThreshold_nucleon
    ! * AngleUpperDetectionThresholdDegrees_nucleon
    ! * kineticEnergyDetectionThreshold_chargedpion
    ! * AngleUpperDetectionThresholdDegrees_chargedpion
    ! * kineticEnergyDetectionThreshold_neutralpion
    ! * AngleUpperDetectionThresholdDegrees_neutralpion
    ! * applyCuts
    ! * Fissum_analysis
    ! * ZeroPion_analysis
    ! * calorimetric_analysis
    ! * radialScale
    ! * reconstruct_neutrino_energy
    ! * outputEvents
    ! * specificEvent_analysis
    ! * Xsection_analysis
    ! * doPipe
    ! * fileNamePipe
    ! * DoOutChannels
    !**************************************************************************
    NAMELIST /neutrinoAnalysis/ &
         detailed_diff_output, include_W_dist, &
         kineticEnergyDetectionThreshold_lepton, &
         AngleUpperDetectionThresholdDegrees_lepton, &
         kineticEnergyDetectionThreshold_nucleon, &
         AngleUpperDetectionThresholdDegrees_nucleon, &
         kineticEnergyDetectionThreshold_chargedpion, &
         AngleUpperDetectionThresholdDegrees_chargedpion,&
         kineticEnergyDetectionThreshold_neutralpion, &
         AngleUpperDetectionThresholdDegrees_neutralpion,&
         applyCuts, Fissum_analysis,            &
         ZeroPion_analysis, calorimetric_analysis, &
         radialScale, reconstruct_neutrino_energy, &
         outputEvents, specificEvent_analysis, &
         Xsection_analysis, &
         doPipe, fileNamePipe, &
         DoOutChannels

    !**************************************************************************
    !****n* neutrinoAnalysis/W_distributions
    ! NAME
    ! NAMELIST /W_distributions/
    ! PURPOSE
    ! This Namelist includes:
    ! * dW_Npi
    ! * Wmin_Npi
    ! * Wmax_Npi
    ! * dW_mupi
    ! * Wmin_mupi
    ! * Wmax_mupi
    ! * dW_muN
    ! * Wmin_muN
    ! * Wmax_muN
    !**************************************************************************
    NAMELIST /W_distributions/  dW_Npi,  Wmin_Npi,  Wmax_Npi, &
             dW_mupi, Wmin_mupi, Wmax_mupi, &
             dW_muN,  Wmin_muN,  Wmax_muN

    !**************************************************************************
    !****n* neutrinoAnalysis/nl_calorimetric_analysis
    ! NAME
    ! NAMELIST /nl_calorimetric_analysis/
    ! PURPOSE
    ! This Namelist includes:
    ! * numin
    ! * numax
    ! * nubin
    ! * Enumin
    ! * Enumax
    ! * Enubin
    !**************************************************************************
    NAMELIST /nl_calorimetric_analysis/ numin,numax,nubin,Enumin,Enumax,Enubin

    !**************************************************************************
    !****n* neutrinoAnalysis/nl_specificEvent
    ! NAME
    ! NAMELIST /nl_specificEvent/
    ! PURPOSE
    ! This namelist includes:
    ! * no_pi
    ! * p_Xn_no_pi
    ! * piplus
    ! * piplus_MULTI
    ! * pi0
    ! * pi0_MULTI
    ! * piminus
    ! * piminus_MULTI
    ! * pp_no_pi
    ! * pn_no_pi
    ! * nn_no_pi
    ! * pp_Xn_no_pi
    ! * nn_Xp_no_pi
    ! * ppp_Xn_no_pi
    ! * pppp_Xn_no_pi
    ! * p_no_pi
    ! * n_no_pi
    ! * Xn_no_pi
    ! * excl_pi0
    ! * excl_piplus
    ! * excl_piminus
    ! * full_incl
    ! * binsizeQ2
    ! * binsizeEnu
    ! * maxQ2
    ! * maxEnu
    ! * excl_hadron
    ! * QEp
    !**************************************************************************
    NAMELIST /nl_specificEvent/ no_pi, p_Xn_no_pi, &
         piplus, piplus_MULTI, pi0, pi0_MULTI, piminus, piminus_MULTI, &
         pp_no_pi, pn_no_pi, nn_no_pi, pp_Xn_no_pi, nn_Xp_no_pi, &
         ppp_Xn_no_pi, pppp_Xn_no_pi, p_no_pi, n_no_pi, Xn_no_pi, &
         excl_pi0, excl_piplus, excl_piminus, full_incl, &
         binsizeQ2, binsizeEnu, maxQ2, maxEnu, excl_hadron, QEp

    !**************************************************************************
    !****n* neutrinoAnalysis/detailed_diff
    ! NAME
    ! NAMELIST /detailed_diff/
    ! PURPOSE
    ! This namelist includes:
    ! * EkinMin
    ! * EkinMax
    ! * dEkin
    ! * forPion
    ! * forEta
    ! * forKaon
    ! * forKaonBar
    ! * forDmeson
    ! * forDbar
    ! * forDs_plus
    ! * forDs_minus
    ! * forNucleon
    ! * forLambda
    ! * forSigmaResonance
    ! * forXi
    ! * forOmegaResonance
    !**************************************************************************
    NAMELIST /detailed_diff/ ekinMin,ekinMax, dEkin,                         &
         &  forpion, foreta, forkaon, forkaonBar, forDmeson, forDbar,        &
         &  forDs_plus, forDs_minus,                                         &
         &  fornucleon, forLambda, forSigmaResonance,forXi,forOmegaResonance

    call Write_ReadingInput('neutrinoAnalysis',0)
    rewind(5)
    read(5,nml=neutrinoAnalysis,IOSTAT=IOS)
    call Write_ReadingInput('neutrinoAnalysis',0,IOS)
    write(*,*) '  detailed_diff_output  =', detailed_diff_output
    write(*,*) '  include_W_dist        =', include_W_dist
    write(*,*) '  calorimetric_analysis =', calorimetric_analysis
    write(*,*) '  specificEvent_analysis=', specificEvent_analysis
    write(*,*) '  reconstruct_neutrino_energy=', reconstruct_neutrino_energy
    write(*,*) '  DoOutChannels         = ',DoOutChannels
    if (doPipe) then
       if (len_trim(fileNamePipe)>0) then
          if (index(fileNamePipe,"/")>0) then
             call ExpandPath(fileNamePipe)
             fileNamePipe = trim(fileNamePipe)
          end if
       else
          fileNamePipe = 'FinalEvents.pipe'
       end if
       write(*,*)
       write(*,*) 'using named pipe "',trim(fileNamePipe),'"'
       write(*,*)
    end if

    write(*,nml=neutrinoAnalysis)

    call Write_ReadingInput('neutrinoAnalysis',1)

    if (detailed_diff_output) then
       call Write_ReadingInput('detailed_diff',0)
       rewind(5)
       read(5,nml=detailed_diff,IOSTAT=IOS)
       call Write_ReadingInput('detailed_diff',0,IOS)
       write(*,'(3(A,g14.5))') '  ekinMin= ',ekinMin, &
            &'  ekinMax= ',ekinMax, &
            &'  dEkin= ',dEkin


       write(*,'(10(A,L3,/))') &
            &'    pion=    ',forpion, &
            &'    Nucleon= ',fornucleon, &
            &'    Lambda=  ',forLambda, &
            &'    Sigma=   ',forSigmaResonance, &
            &'    Xi=      ',forXi, &
            &'    Omega    ',forOmegaResonance, &
            &'    eta=     ',foreta,&
            &'    kaon=    ',forkaon,&
            &'    kaonBar= ',forkaonBar,&
            &'    Dmeson=  ',forDmeson,&
            &'    Dbar=    ',forDbar,&
            &'    Ds+=     ',forDs_plus,&
            &'    Ds-=     ',forDs_minus

       call Write_ReadingInput('detailed_diff',1)
    end if

    if (include_W_dist) then
       call Write_ReadingInput('W_distributions',0)
       rewind(5)
       read(5,nml=W_distributions,IOSTAT=IOS)
       call Write_ReadingInput('W_distributions',0,IOS)
       write(*,'(3(A,g14.5,/))') &
            & '  dW_Npi = ',dW_Npi,&
            & '  dW_mupi= ',dW_mupi, &
            & '  dW_muN = ',dW_muN
       call Write_ReadingInput('W_distributions',1)
    end if

    if (calorimetric_analysis) then
       call Write_ReadingInput('nl_calorimetric_analysis',0)
       rewind(5)
       read(5,nml=nl_calorimetric_analysis,IOSTAT=IOS)
       call Write_ReadingInput('nl_calorimetric_analysis',0,IOS)
       write(*,'(3(A,g14.5,/))') &
            & '   numin= ',numin, &
            & '   numax= ',numax, &
            & '   nubin= ',nubin
       write(*,'(3(A,g14.5,/))') &
            & '   Enumin= ',Enumin, &
            & '   Enumax= ',Enumax, &
            & '   Enubin= ',Enubin
       call Write_ReadingInput('nl_calorimetric_analysis',1)
    end if

    if (specificEvent_analysis) then
       call Write_ReadingInput('nl_specificEvent',0)
       rewind(5)
       read(5,nml=nl_specificEvent,IOSTAT=IOS)
       call Write_ReadingInput('nl_specificEvent',0,IOS)
       write(*,'(A,/,20(A,L3,/))') '  specificEvent_analysis = .true. :', &
            & '    no_pi=        ', no_pi,&
            & '    p_Xn_no_pi=   ', p_Xn_no_pi, &
            & '    piplus=       ', piplus,&
            & '    piplus_MULTI= ', piplus_MULTI, &
            & '    pi0=          ', pi0, &
            & '    pi0_MULTI=    ', pi0_MULTI, &
            & '    piminus=      ', piminus, &
            & '    piminus_MULTI=', piminus_MULTI, &
            & '    pp_no_pi=     ', pp_no_pi, &
            & '    pn_no_pi=     ', pn_no_pi, &
            & '    nn_no_pi=     ', nn_no_pi, &
            & '    pp_Xn_no_pi=  ', pp_Xn_no_pi, &
            & '    nn_Xp_no_pi=  ', nn_Xp_no_pi, &
            & '    ppp_Xn_no_pi= ', ppp_Xn_no_pi, &
            & '    pppp_Xn_no_pi=', pppp_Xn_no_pi, &
            & '    p_no_pi=      ', p_no_pi, &
            & '    n_no_pi=      ', n_no_pi, &
            & '    Xn_no_pi=     ', Xn_no_pi, &
            & '    excl_pi0=     ', excl_pi0, &
            & '    excl_piplus=  ', excl_piplus, &
            & '    excl_piminus= ', excl_piminus, &
            & '    full_incl=    ', full_incl, &
            & '    excl_hadron=  ', excl_hadron, &
            & '    QEp           ', QEp
       if(QEp) then
          pcut = sqrt(kineticEnergyDetectionThreshold_nucleon**2  &
               &  + 2*mN*kineticEnergyDetectionThreshold_nucleon)
          write (*,*) 'pcut = ', pcut
       end if

 ! The switch excl_hadron selects truly exclusive 1-meson events
 ! ie. exactly 1 meson with fixed charged and no other mesons of any kind
 ! switch QEp selects QE-like events with 1 mu, 0 pi, X p


       call Write_ReadingInput('nl_specificEvent',1)
    end if

    call set_Exclusive(excl_hadron)
    if (excl_hadron) then
       excl_pi0 = .true.
       excl_piminus = .true.
       excl_piplus = .true.
       piplus = .false.
       pi0 = .false.
       piminus = .false.
    end if

    write(*,*) 'neutrinoAnalysis :',QEp,pcut
    call set_QElike(QEp,pcut)

    includeSpeEvent = (/no_pi,p_Xn_no_pi,piplus,piplus_MULTI,pi0,pi0_MULTI, &
         & piminus, piminus_MULTI,pp_no_pi, pn_no_pi, nn_no_pi, pp_Xn_no_pi,&
         & nn_Xp_no_pi, ppp_Xn_no_pi, pppp_Xn_no_pi,p_no_pi, n_no_pi, Xn_no_pi,&
         & excl_pi0,excl_piplus,excl_piminus,full_incl/)
    ! These specific events are used only in calculations of cross sections
    ! that involve lepton variables

  end subroutine readinput





  subroutine cleanUp

    call RemoveHist(dE_hists(:,:))
    call RemoveHist(dE_hists_QE(:,:))
    call RemoveHist(dE_hists_Delta(:,:))
    call RemoveHist(dE_hists_highRes(:,:))
    call RemoveHist(dE_hists_1piBG(:,:))
    call RemoveHist(dE_hists_2piBG(:,:))
    call RemoveHist(dE_hists_DIS(:,:))
    call RemoveHist(dE_hists_gen0(:,:))
    call RemoveHist(dE_hists_Multi(:,:))
    call RemoveHist(dE_hists_QE_Multi(:,:))
    call RemoveHist(dE_hists_Delta_Multi(:,:))
    call RemoveHist(dE_hists_highRes_Multi(:,:))
    call RemoveHist(dE_hists_1piBG_Multi(:,:))
    call RemoveHist(dE_hists_DIS_Multi(:,:))
    call RemoveHist(dE_hists_gen0_Multi(:,:))

    call RemoveHist(dE_hists_1X(:,:))
    call RemoveHist(dE_hists_QE_1X(:,:))
    call RemoveHist(dE_hists_Delta_1X(:,:))
    call RemoveHist(dE_hists_highRes_1X(:,:))
    call RemoveHist(dE_hists_1piBG_1X(:,:))
    call RemoveHist(dE_hists_2piBG_1X(:,:))
    call RemoveHist(dE_hists_DIS_1X(:,:))
    call RemoveHist(dE_hists_gen0_1X(:,:))
    call RemoveHist(dE_hists_2X(:,:))
    call RemoveHist(dE_hists_QE_2X(:,:))
    call RemoveHist(dE_hists_Delta_2X(:,:))
    call RemoveHist(dE_hists_highRes_2X(:,:))
    call RemoveHist(dE_hists_1piBG_2X(:,:))
    call RemoveHist(dE_hists_DIS_2X(:,:))
    call RemoveHist(dE_hists_gen0_2X(:,:))

    call RemoveHist(dTheta_hists(:,:))
    call RemoveHist(dPhi_hists(:,:))
    call RemoveHist(dTheta_hists_QE(:,:))
    call RemoveHist(dPhi_hists_QE(:,:))
    call RemoveHist(dTheta_hists_Delta(:,:))
    call RemoveHist(dPhi_hists_Delta(:,:))
    call RemoveHist(dTheta_hists_highRes(:,:))
    call RemoveHist(dPhi_hists_highRes(:,:))
    call RemoveHist(dTheta_hists_1piBG(:,:))
    call RemoveHist(dPhi_hists_1piBG(:,:))
    call RemoveHist(dTheta_hists_2piBG(:,:))
    call RemoveHist(dPhi_hists_2piBG(:,:))
    call RemoveHist(dTheta_hists_gen0(:,:))
    call RemoveHist(dPhi_hists_gen0(:,:))
    call RemoveHist2d(dOmega_hists(:,:))
    call RemoveHist2d(dOmega_hists_QE(:,:))
    call RemoveHist2d(dOmega_hists_Delta(:,:))
    call RemoveHist2d(dOmega_hists_highRes(:,:))
    call RemoveHist2d(dOmega_hists_1piBG(:,:))
    call RemoveHist2d(dOmega_hists_gen0(:,:))

    call RemoveHist(dEcostheta_hists(:,:))
    call RemoveHist(dEcostheta_hists_QE(:,:))
    call RemoveHist(dEcostheta_hists_Delta(:,:))
    call RemoveHist(dEcostheta_hists_highRes(:,:))
    call RemoveHist(dEcostheta_hists_1piBG(:,:))
    call RemoveHist(dEcostheta_hists_2piBG(:,:))
    call RemoveHist(dEcostheta_hists_DIS(:,:))
    call RemoveHist(dEcostheta_hists_gen0(:,:))
    call RemoveHist(dEcostheta_hists_MULTI(:,:))
    call RemoveHist(dEcostheta_hists_QE_MULTI(:,:))
    call RemoveHist(dEcostheta_hists_Delta_MULTI(:,:))
    call RemoveHist(dEcostheta_hists_highRes_MULTI(:,:))
    call RemoveHist(dEcostheta_hists_1piBG_MULTI(:,:))
    call RemoveHist(dEcostheta_hists_DIS_MULTI(:,:))
    call RemoveHist(dEcostheta_hists_gen0_MULTI(:,:))
    call RemoveHist(dcostheta_hists(:,:))
    call RemoveHist(dcostheta_hists_QE(:,:))
    call RemoveHist(dcostheta_hists_Delta(:,:))
    call RemoveHist(dcostheta_hists_highRes(:,:))
    call RemoveHist(dcostheta_hists_1piBG(:,:))
    call RemoveHist(dcostheta_hists_2piBG(:,:))
    call RemoveHist(dcostheta_hists_DIS(:,:))
    call RemoveHist(dcostheta_hists_gen0(:,:))
    call RemoveHist(dcostheta_hists_MULTI(:,:))
    call RemoveHist(dcostheta_hists_QE_MULTI(:,:))
    call RemoveHist(dcostheta_hists_Delta_MULTI(:,:))
    call RemoveHist(dcostheta_hists_highRes_MULTI(:,:))
    call RemoveHist(dcostheta_hists_1piBG_MULTI(:,:))
    call RemoveHist(dcostheta_hists_2piBG_MULTI(:,:))
    call RemoveHist(dcostheta_hists_DIS_MULTI(:,:))
    call RemoveHist(dcostheta_hists_gen0_MULTI(:,:))

    call RemoveHist(dE_hists_0pions(:,:))
    call RemoveHist(dE_hists_QE_0pions(:,:))
    call RemoveHist(dE_hists_Delta_0pions(:,:))
    call RemoveHist(dE_hists_highRes_0pions(:,:))
    call RemoveHist(dE_hists_Multi_0pions(:,:))
    call RemoveHist(dE_hists_QE_Multi_0pions(:,:))
    call RemoveHist(dE_hists_Delta_Multi_0pions(:,:))
    call RemoveHist(dE_hists_highRes_Multi_0pions(:,:))
    call RemoveHist(dE_hists_1X_0pions(:,:))
    call RemoveHist(dE_hists_QE_1X_0pions(:,:))
    call RemoveHist(dE_hists_Delta_1X_0pions(:,:))
    call RemoveHist(dE_hists_highRes_1X_0pions(:,:))
    call RemoveHist(dE_hists_2X_0pions(:,:))
    call RemoveHist(dE_hists_QE_2X_0pions(:,:))
    call RemoveHist(dE_hists_Delta_2X_0pions(:,:))
    call RemoveHist(dE_hists_highRes_2X_0pions(:,:))

    call RemoveHist(dW_Npi_hists(:))
    call RemoveHist(dW_mupi_hists(:))
    call RemoveHist(dW_muN_hists(:))

    call RemoveHist(dEnu_hist(:,:))
    call RemoveHist(dElepton_hist(:,:))
    call RemoveHist(dcoslepton_hist(:,:))
    call RemoveHist(dQ2lepton_hist(:,:))
    call RemoveHist(dQ2plepton_hist(:,:))
    call RemoveHist2D(dSigmaMC_EprimeCost_0pi(:,:))
    call RemoveHist2D(dSigmaMC_pLpT_0pi(:,:))

    if (doPipe) then
       close(47) ! close the pipe
    end if

  end subroutine cleanUp






  !****************************************************************************
  !****s* neutrinoAnalysis/neutrino_Analyze
  ! NAME
  ! subroutine neutrino_Analyze(Particles,finalFlag,num_runs_sameEnergy)
  ! INPUTS
  ! * type(particle), intent(in),dimension(:,:)  :: Particles
  !   -- Particles which shall be analyzed
  ! * logical, intent(in) :: finalFlag -- if .true. than the final output
  !   for a series of calls will be done
  ! * integer             :: num_runs_sameEnergy
  ! NOTES
  ! * This subroutine produces output for neutrino-nucleus scattering.
  !****************************************************************************
  subroutine neutrino_Analyze(Particles,finalFlag,num_runs_sameEnergy)
    use initNeutrino, only: getFirstEventRange, getNeutrinoInfo, nuEXP
    use particleDefinition
    use AnaEventDefinition
    use IDTable, only: nucleon,pion
    use rotation, only: get_Phi_Theta
    use degRad_conversion, only: radian
    use output, only: intToChar, intToChar4
    use history, only: history_getGeneration
    use neutrino_IDTable
    use neutrinoInfoStorage, only:  neutrinoInfoStorage_clear
    use neutrinoProdInfo, only: neutrinoProdInfo_Get!, neutrinoProdInfo_clear
    use expNeutrinofluxes, only: CCQE_recQs, CCQE_recQs_Delta, CCQE_recEnergy, &
         & CCQE_recEnergy_Delta,K2K_recEnergy, K2K_recQs
    use initNeutrino, only: includeQE, includeDELTA, includeRES, include1pi, &
         & includeDIS, include2p2hQE, include2p2hDelta, include2pi
    use minkowski, only: abs4Sq
    use vector, only: absVec
    use constants, only: pi
    use MultiplicityAnalysis, only: Multiplicity_Reset, Multiplicity_AddEvent, &
        & Multiplicity_Write
    use ZeroPionAnalysis, only: event_sigma_0pions, event_dsigma_de_0pions
    use PreEvList, only: CreateSortedPreEvent, PreEvList_CLEAR, PreEvList_INIT,&
         PreEvList_INSERT, PreEvList_Print


    type(particle), intent(in),dimension(:,:) ,target :: Particles
    logical, intent(in) :: finalFlag
    integer, intent(in) :: num_runs_sameEnergy

    ! Local variables:
    integer, dimension (1:2) :: firstEvents
    type(tAnaEvent), Allocatable, dimension(:) :: events,events_QE,    &
         events_Delta,events_highRES, &
         events_1piBG,events_DIS,events_2p2h, events_2piBG, &
         events0 ! like events, but without struck nucleon

    type(tAnaEvent), Allocatable, dimension(:) :: events_gen0,    &
         events_gen1,events_gen2, &
         events_gen3ormore ! A list of all events
    type(particle), POINTER :: Part
    integer :: i,j,first

    type(particle), Allocatable, dimension(:),target :: lepton, struckNuc

    integer, parameter :: max_generation=3
    integer :: generation, prod_id

    real  :: dPhi   ! Delta(phi) for dsigma/dOmega
    real  :: dTheta ! Delta(theta) for dsigma/dOmega
    real :: theta, phi
    real :: raiseFlagVariable

    integer,save :: numberOfCalls=0
    integer,save :: numberOfFinals=0

    logical, dimension(1:numStableParts) :: printflags

    ! In these hists we save the information
    ! which the "AnaEvent" subroutines are returning
    real, dimension(1:dimSigma,1:2),save :: sigma, sigma_QE,sigma_Delta,  &
       & sigma_highRES,sigma_1piBG, &
       & sigma_2piBG,sigma_DIS,sigma_2p2h
    real, dimension(1:dimSigma,1:2),save :: sigma_gen0,sigma_gen1,sigma_gen2, &
       & sigma_gen3ormore

    ! In these hists we save the information
    ! which the "ZeroPionAnalysis" subroutines are returning:
    !  in particular channels  with "k" protons and "n" neutrons
    real, dimension(1:dimSigma,1:2),save :: sigma_0pions, sigma_QE_0pions, &
       & sigma_Delta_0pions, &
       & sigma_highRES_0pions, sigma_1piBG_0pions, &
       & sigma_2piBG_0pions, sigma_DIS_0pions, sigma_2p2h_0pions

    real :: tkin
    integer :: ntk
    real, parameter :: EcosthetaMin=0.
    real, parameter :: EcosthetaMax=2.
    real, parameter :: dEcostheta=0.01

    real,dimension(-1:1) :: tsigmapion=0.
    real,dimension(-1:1,0:200) :: tksigmapion=0.
    real,dimension(-1:1) :: tsigmanucleon=0.          ! in DIS antiprotons are also possible
                                                      ! in final states
    real,dimension(-1:1,0:200) :: tksigmanucleon=0.   ! in DIS antiproton are also possible

    real,dimension(-1:1),save :: sum_tsigmapion=0.
    real,dimension(-1:1,0:200),save :: sum_tksigmapion=0.
    real,dimension(-1:1),save :: sum_tsigmanucleon=0. ! in DIS antiproton are also possible
                                                      ! in final states
    real,dimension(-1:1,0:200),save :: sum_tksigmanucleon=0.  ! in DIS antiproton are also
                                                              ! possible in final states

    real,dimension(0:200) :: Emiss_Fissum=0.
    real,dimension(0:200),save ::  sum_Emiss_Fissum=0.

    integer :: numberofneutrons_tot=0
    integer :: numberofprotons_tot=0
    integer :: numberofANTIprotons_tot=0
    integer :: numberofneutrons_out=0
    integer :: numberofprotons_out=0
    integer :: numberofANTIprotons_out=0
    integer,dimension(-1:1,0:max_generation) :: numberofnucleons_out
    ! in DIS antiproton are also possible in final states

    integer,save :: sum_numberofneutrons_tot=0
    integer,save :: sum_numberofprotons_tot=0
    integer,save :: sum_numberofANTIprotons_tot=0
    integer,save :: sum_numberofneutrons_out=0
    integer,save :: sum_numberofprotons_out=0
    integer,save :: sum_numberofANTIprotons_out=0
    integer,save,dimension(-1:1,0:max_generation) :: sum_numberofnucleons_out

    ! reconstruction of kinematics : as in MiniBooNN, K2K:
    type(histogram),save   :: H_Q2_real(1:max_speEvent,0:max_Hist), &
       & H_Q2_rec(1:max_speEvent,0:max_Hist)
    type(histogram2D),save :: H_Q2_rec_versus_real(1:max_speEvent,0:max_Hist)
    type(histogram),save   :: H_enu_real(1:max_speEvent,0:max_Hist), &
       & H_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram2D),save :: H_enu_rec_versus_real(1:max_speEvent,0:max_Hist)

    real :: Enureal, Enurec, Q2real, Q2rec, dummy, perweight
    real, dimension(0:3) :: lepIn_mom, lep_mom, boson_mom
    ! momenta of the ingoing and outgoing lepton and intermediate boson

    integer :: m, iHist

    character*(10) :: Prefix_MultAna
    character(100) :: filename
    character(13) :: filename1
    character(13) :: string

    real :: L, Posc_mumu,Posc_mue,Posc_mue_max,Posc_mue_antimax
    ! used for oscillation analysis

    type(histogram),save :: Oscmumu_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmumu_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram),save :: Oscmuemax_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmuemax_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram),save :: Oscmue_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmue_enu_rec(1:max_speEvent,0:max_Hist)
    type(histogram),save :: Oscmueantimax_enu_real(1:max_speEvent,0:max_Hist), &
       & Oscmueantimax_enu_rec(1:max_speEvent,0:max_Hist)

    dPhi=radian(10.)   ! Delta(phi) for dsigma/dOmega
    dTheta=radian(10.) ! Delta(theta) for dsigma/dOmega

    if (initflag) then
       call readinput
       initflag=.false.
       numberOfCalls=0
       numberOfFinals=0

       sum_numberofneutrons_tot=0
       sum_numberofprotons_tot=0
       sum_numberofANTIprotons_tot=0
       sum_numberofneutrons_out=0
       sum_numberofprotons_out=0
       sum_numberofANTIprotons_out=0
       sum_numberofnucleons_out=0
       sum_tsigmapion=0.
       sum_tksigmapion=0.
       sum_tsigmanucleon=0.
       sum_tksigmanucleon=0.
       sum_Emiss_Fissum=0.



       !  Now switches for printout of cross sections for various hadrons
       !
       !  indices in printflags(i) refer to ordering in field particleIds
       !  in module AnaEvent

       printFlags=.false.
       printFlags(1:13) = (/ forpion, foreta, forkaon, forkaonBar, &
            forDmeson, forDbar, forDs_plus, forDs_minus, fornucleon, &
            forLambda, forSigmaResonance, forXi, forOmegaResonance /)
       call set_particleIDs_flag(printFlags)


       if (Xsection_analysis) then

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows the cross sections for  preselected final states
       !
       ! Units:
       ! * For process_ID=CC and NC the units 10^{-38} cm^2 for integrated xsec
       !   (10^{-38)cm^2/GeV for dsigma/dElepton,
       !   10^{-38)cm^2/GeV^2 for dsigma/dQ^2, and so on)
       ! * For process_ID=EM the units are nanobarns=10^{-33}cm^2
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
       !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
       ! * #2-#119: see description in AnaEvent.f90, subroutine event_sigma
       !   OR in the output file sigma.dat
       !   The description of some columns is given below
       !   In each channel the outgoing lepton is presupposed
       !   (unless explicitely stated otherwise)
       !
       ! Some columns in detail:
       ! * #2: 1 pi-, no other pions of any charge and anything else
       ! * #3: 1 pi0, no other pions of any charge and anything else
       ! * #4: 1 pi+, no other pions of any charge and anything else
       ! * #5: 1 eta  and anything else
       ! * #14: 1 neutron, no other nucleons and anything else
       ! * #15: 1 proton,  no other nucleons and anything else
       ! * #68: 1 nucleon and 1 pion, no other pions or nucleons, anything else
       ! * #68: 1 proton (no other nucleons) and 0 pions and anything else
       !   (QE-like in Cherenkov detector)
       ! * #69: 0 pions and  anything else
       ! * #70: at least  1 pi-, any number of pi0 and/or pi+ , anything else,
       !   each event is counted once
       ! * #71: at least  1 pi0, any number of pi- and/or pi+ , anything else,
       !   each event is counted once
       ! * #72: at least  1 pi+, any number of pi- and/or pi0 , anything else,
       !   each event is counted once
       ! * #91: at least  1 pi0, any number of pi- and/or pi+ , anything else,
       !   each pi0 is counted
       ! * #96: no nucleons, anything else
       ! * #97: 5 or more nucleons, anything else
       ! * #120-#239:  errors to Columns #2-#119
       !***********************************************************************

       open(10,File='Neutrino_total_Xsection_multiplicities.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_QE.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_QE.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat
       ! but for QE events (=the first interaction was quasielastic or
       ! elastic scattering)
       !***********************************************************************
       if (includeQE) then
       open(10,File='Neutrino_total_Xsection_multiplicities_QE.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/neutrino_total_Xsection_multiplicities_Delta.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_Delta.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat
       ! but for Delta events (=the first interaction was production of the
       ! Delta resonance)
       !***********************************************************************
       if (includeDELTA) then
       open(10,File='Neutrino_total_Xsection_multiplicities_Delta.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_highRES.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_highRES.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat
       ! but for highRES events (=the first interaction was production any
       ! resonance beyond Delta)
       !***********************************************************************
       if (includeRES) then
       open(10,File='Neutrino_total_Xsection_multiplicities_highRES.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_1piBG.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_1piBG.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat
       ! but for 1piBG events (=the first interaction was background
       ! production of 1-pion final state)
       !***********************************************************************
       if (include1pi) then
       open(10,File='Neutrino_total_Xsection_multiplicities_1piBG.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_DIS.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_DIS.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat
       ! but for DIS events (=the first interaction was DIS)
       !***********************************************************************
       if (includeDIS) then
       open(10,File='Neutrino_total_Xsection_multiplicities_DIS.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_2p2h.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_2p2h.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat
       ! but for 2particle-2hole events
       !***********************************************************************
       if (include2p2hQE .or. include2p2hDelta) then
       open(10,File='Neutrino_total_Xsection_multiplicities_2p2h.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)
       end if


       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_gen0.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_gen0.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 0th generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='Neutrino_total_Xsection_multiplicities_gen0.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_gen1.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_gen1.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 1st generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='Neutrino_total_Xsection_multiplicities_gen1.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_gen2.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_gen2.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 2nd generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='Neutrino_total_Xsection_multiplicities_gen2.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection_multiplicities_gen3ormore.dat
       ! NAME
       ! file Neutrino_total_Xsection_multiplicities_gen3ormore.dat
       ! PURPOSE
       ! The same as Neutrino_total_Xsection_multiplicities.dat for particles
       ! of the 3rd-or-more generation
       ! The definition of generation is given in history.f90, description
       ! of the module "history"
       !***********************************************************************
       open(10,File='Neutrino_total_Xsection_multiplicities_gen3ormore.dat')
       write(10,*)'# total Xsections for defined finalstates (see sigma.dat)'
       write(10,*)'# next 120 columns: error'
       write(10,*)'# order of Xsections as in sigma.dat'
       close(10)

       open(10,File='Neutrino_total_Xsection.dat')
       write(10,*)'# Total xsec for outgoing pions and nucleons:&
          & each pion/nucleon contribute to the xsec'
       write(10,*)'# 1: raiseVariable  2: sigma(pi-)    3:pi0    4:pi+  &
          & 5:antiproton 6:neutron, 7:proton'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_total_Xsection.dat
       ! NAME
       ! file Neutrino_total_Xsection.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows the total cross sections after
       ! final state interactions for pions, protons and neutrons
       ! The cross sections are semi-inclusive, i.e. for pi+, they may contain
       ! events with only one pi+, but also events with pi+pi0.
       !
       ! Units:
       ! * For process_ID=CC and NC: 10^{-38} cm^2/GeV for integrated
       !   xsec and so on ..
       ! * For process_ID=EM: nanobarns=10^{-33}cm^2/ GeV and so on.
       ! * All cross sections are given per nucleon (1/A)
       ! * The various columns in this file are labeld by the particle species
       !
       !***********************************************************************


       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_kinetic_energy_Xsection.dat
       ! NAME
       ! file Neutrino_kinetic_energy_Xsection.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows the kinetic energy differential cross sections after
       ! final state interactions for pions, protons and neutrons
       !
       ! Units:
       ! * For process_ID=CC and NC: 10^{-38} cm^2/GeV for integrated
       !   xsec and so on ..
       ! * For process_ID=EM: nanobarns=10^{-33}cm^2/ GeV and so on.
       ! * All cross sections are given per nucleon (1/A)
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
       !   Elepton for nuXsectionMode=2=dSigmadQsdElepton and so on)
       ! * #2: kinetic energy in GeV
       ! * #3: xsec for events with at least one pi-  in the final state
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_pi_charge_-1_MULTI.dat)
       ! * #4: xsec for events with at least one pi0
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_pi_charge_+0_MULTI.dat)
       ! * #5: xsec for events with at least one pi+
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_pi_charge_+1_MULTI.dat)
       ! * #6: xsec for events with at least one neutron
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_N_charge_+0_MULTI.dat)
       ! * #7: xsec for events with at least one proton
       !   (coincides with column 2 in
       !   diff_XXX_dSigma_dEkin_N_charge_+1_MULTI.dat)
       !
       ! HERE  XXX=000,001, ... is the count of the raise variable
       !***********************************************************************
       open(10,File='Neutrino_kinetic_energy_Xsection.dat')
       write(10,*)'# order of entries: 1:raise variable  2:ekin, 3:sig pi-, &
          & 4:sig pi0, 5:sig pi+, 6:sig neut, 7:sig prot'
       close(10)

       !***********************************************************************
       !****o* neutrinoAnalysis/Neutrino_Xsection_numbers.dat
       ! NAME
       ! file Neutrino_Xsection_numbers.dat
       ! PURPOSE
       ! The file is produced in the runs with eventtype=5=neutrino .
       !
       ! The file shows number of nucleons produced in one run before and
       ! after final state interactions
       ! (should be devided by numEnsembles to obtain the average multiplicity
       ! per target nucleus and divided further by target_A to obtain the
       ! average multiplicity per target nucleon)
       !
       ! Columns:
       ! * #1: variable which was raised
       !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
       !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
       ! * #2:  all the protons produced  before final state interactions (FSI)
       ! * #3: protons that made it out of the nucleus after FSI
       ! * #4: protons produced in generation 0 and made it out of the nucleus
       !   after FSI
       ! * #5: protons produced in generation 1 and made it out of the nucleus
       !   after FSI
       ! * #6: protons produced in generation 2 and made it out of the nucleus
       !   after FSI
       ! * #7: protons produced in generation 3-or-more and made it out of the
       !   nucleus after FSI
       ! * #8: all the neutrons produced before FSI
       ! * #9:  neutrons that made it out of the nucleus after FSI
       ! * #10: neutrons produced in generation 0 and made it out of the
       !   nucleus after FSI
       ! * #11: neutrons produced in generation 1 and made it out of the
       !   nucleus after FSI
       ! * #12: neutrons produced in generation 2 and made it out of the
       !   nucleus after FSI
       ! * #13: neutrons produced in generation 3-or-more and made it out of
       !   the nucleus after FSI
       ! * #14: all the antiprotons produced before FSI
       ! * #15: antiprotons that made it out of the nucleus after FSI
       ! * #16: antiprotons produced in generation 0 and made it out of the
       !   nucleus after FSI
       ! * #17: antiprotons produced in generation 1 and made it out of the
       !   nucleus after FSI
       ! * #18: antiprotons produced in generation 2 and made it out of the
       !   nucleus after FSI
       ! * #19: antiprotons produced in generation 3-or-more and made it out
       !   of the nucleus after FSI
       !***********************************************************************
       open(10,File='Neutrino_Xsection_numbers.dat')
       write(10,*) '#1:raiseFlagVariable, #2:p_tot, #3:p_out, #4:p_out_gen0, &
            & #5:p_out_gen1, &
            & #6:p_out_gen2, #7:p_out_gen3ormore,&
            & #8:n_tot, #9:n_out, #10:n_out_gen0, #11:n_out_gen1 ,#12:n_out_gen2,&
            & #13:n_out_gen3ormore, &
            & #14:barp_tot, #15:barp_out, #16:barp_out_gen0, #17:barp_out_gen1, &
            & #18:barp_out_gen2,#19:bar p_out_gen3ormore'
       close(10)

       end if ! Xsection_analysis

       if (Fissum_analysis) then
          open(10,File='Neutrino_Emiss_spectrum.dat')
          write(10,*)'# order of entries: ekin, sig prot'
          close(10)
       end if

       if (ZeroPion_analysis) then
          call WriteHeader10(.true., 'neutrino_0pions.dat')
          call WriteHeader10(includeQE,'neutrino_0pions_QE.dat')
          call WriteHeader10(includeDELTA,'neutrino_0pions_Delta.dat')
          call WriteHeader10(includeRES,'neutrino_0pions_highRES.dat')
          call WriteHeader10(include1pi,'neutrino_0pions_1piBG.dat')
          call WriteHeader10(include2pi,'neutrino_0pions_2piBG.dat')
          call WriteHeader10(includeDIS,'neutrino_0pions_DIS.dat')
          call WriteHeader10(include2p2hQE .or. include2p2hDelta,&
          & 'neutrino_0pions_2p2h.dat')
       end if

       if (reconstruct_neutrino_energy .and. specificEvent_Analysis) then
          if (binsizeQ2.lt.0.001) binsizeQ2=0.02
          if (binsizeEnu.lt.0.001) binsizeEnu=0.01

          do m=1, max_speEvent
             if (.not.includeSpeEvent(m)) cycle
             do iHist=0, max_Hist
                if (.not.includeHist(iHist)) cycle

                ! Q2 reconstruction
                call CreateHist(H_Q2_real(m,iHist),  &
                   &'true Q2 for a specific event',0.,maxQ2, binsizeQ2)
                call CreateHist(H_Q2_rec(m,iHist),  'reconstructed Q2  &
                   & for a specific event',0., maxQ2,binsizeQ2)
                call CreateHist2D(H_Q2_rec_versus_real(m,iHist), &
                     'reconstructed Q2 versus real Q2 for a specific event', &
                     (/0.,0./),(/maxQ2,maxQ2/),  (/binsizeQ2,binsizeQ2/))

                ! neutrino energy reconstruction
                if (nuExp>0) call CreateHist(H_enu_real(m,iHist), &
                   & 'true Enu for a specific event',0.,maxEnu,binsizeEnu)
                call CreateHist(H_enu_rec(m,iHist),  'reconstructed Enu&
                   & for a specific event', 0.,maxEnu,binsizeEnu)
                if (nuExp>0) call CreateHist2D(H_enu_rec_versus_real(m,iHist), &
                     'reconstructed Enu versus real Enu for a specific event', &
                     (/0.,0./),(/maxEnu,maxEnu/),  (/binsizeEnu,binsizeEnu/))

                ! oscillations:  nu_mu survival, nu_e appearence
                if (isOSC()) then
                   call CreateHist(Oscmumu_enu_real(m,iHist), &
                        'nu_mu survival versus true energy',0.,maxEnu,binsizeEnu)
                   call CreateHist(Oscmumu_enu_rec(m,iHist) , &
                        'nu_mu survival versus reconstructed energy',0.,maxEnu, &
                        & binsizeEnu)
                   call CreateHist(Oscmuemax_enu_real(m,iHist), &
                        & 'nu_e appearence versus true energy for&
                        & delta_CP=pi/2',0., maxEnu,binsizeEnu)
                   call CreateHist(Oscmuemax_enu_rec(m,iHist),  &
                        &  'nu_e appearence versus reconstructed energy&
                        & for delta_CP=pi/2', 0.,maxEnu,binsizeEnu)
                   call CreateHist(Oscmue_enu_real(m,iHist), &
                        'nu_e appearence versus true energy',0.,maxEnu,binsizeEnu)
                   call CreateHist(Oscmue_enu_rec(m,iHist), &
                        'nu_e appearence versus reconstructed energy',0., &
                        & maxEnu,binsizeEnu)
                   call CreateHist(Oscmueantimax_enu_real(m,iHist), &
                        &  'nu_e appearence versus true energy for&
                        & delta_CP=-pi/2',0., maxEnu,binsizeEnu)
                   call CreateHist(Oscmueantimax_enu_rec(m,iHist),  &
                        &  'nu_e appearence versus reconstructed energy&
                        & for delta_CP=-pi/2', 0.,maxEnu,binsizeEnu)
                end if

             end do !iHist
          end do !m
       end if


       call CreateHist(hNucleonVacuumMass, 'mass of nucleons', 0.65,1.2,0.001)

       call Multiplicity_Reset

    end if

    call getNeutrinoInfo(raiseFlagVariable)

    numberOfCalls=numberOfCalls+1

    initHists = (numberOfCalls.eq.1) ! whether we initialize the histograms

    !==========================================================================
    !
    ! After all these preparations now the actual analysis of events starts
    !
    !==========================================================================

    write(*,*) '################### NEUTRINO ANALYSIS STARTS ##################'
    write(*,*) ' number of calls: ', numberofCalls, &
         ' number of finals: ', numberofFinals


    write(*,'(a,F12.4)')'radialScale',radialScale

    !==========================================================================
    ! First calculate
    ! total (tsigma) and kinetic energy differential (tksigma) cross section
    ! for pions and nucleons
    !==========================================================================

    !nullify everything
    tsigmapion=0.
    tksigmapion=0.
    tsigmanucleon=0.
    tksigmanucleon=0.
    numberofneutrons_tot=0
    numberofprotons_tot=0
    numberofANTIprotons_tot=0
    numberofneutrons_out=0
    numberofprotons_out=0
    numberofANTIprotons_out=0
    numberofnucleons_out=0
    Emiss_Fissum=0.

    Q2real=0.
    Enureal=0.
    Q2rec=0.
    Enurec=0.

    !==========================================================================
    ! First calculate total cross sections and kinetic energy distributions for
    ! pions and then numbers of nucleons in the final state
    ! acceptance cuts for leptons as well as bound-state criteria are taken into
    ! account
    ! Printout is controlled by switch Xsection_analysis
    !==========================================================================


    do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
       do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
          if (Particles(i,j)%ID.le.0) cycle

          Part=>Particles(i,j)

          if (Part%ID.eq.1) then ! == nucleon
             call AddHist(hNucleonVacuumMass,Part%mass,Part%perweight)
             select case (Part%charge)
             case (-1)
                numberofANTIprotons_tot=numberofANTIprotons_tot+1
             case (0)
                numberofneutrons_tot=numberofneutrons_tot+1
             case (1)
                numberofprotons_tot=numberofprotons_tot+1
             end select
          end if

          if (iand(applyCuts,2)==2) then
             if (IsBound(Part)) cycle
          end if
          if (iand(applyCuts,4)==4) then
             if (IsBelowThreshold(Part)) cycle
          end if

          ! The following call to function neutrinoProdInfo_Get retrieves
          ! information on the production process for particle j in ensemble i

          if (.not.neutrinoProdInfo_Get(Part%firstEvent,prod_id,&
            &  dummy,lepIn_mom,lep_mom,boson_mom)) then
             call TRACEBACK('error in getting production info')
          end if

          if (iand(applyCuts,1)==1) then
             if (.not.lepton_acceptance(lep_mom)) cycle
          end if


          tkin=Part%mom(0)-Part%mass
          if (tkin.lt.0.) tkin = 0.
          ntk=min(int(tkin/dEkin),200)


          generation=history_getGeneration(Part%history)
          if (generation.gt.max_generation) generation=max_generation

          select case (Part%ID)
          case (pion)

             !total cross section
             tsigmapion(Part%charge)=tsigmapion(Part%charge) &
                  +Part%perweight

             !kinetic energy differential xsection
             tksigmapion(Part%charge,ntk)= &
                  & Part%perweight/dEkin   &
                  & +tksigmapion(Part%charge,ntk)

          case (nucleon)

             select case (Part%charge)
             case (-1)
                numberofANTIprotons_out=numberofANTIprotons_out+1
             case (0)
                numberofneutrons_out=numberofneutrons_out+1
             case (1)
                numberofprotons_out=numberofprotons_out+1
             end select

             numberofnucleons_out(Part%charge,generation)  &
             & =numberofnucleons_out(Part%charge,generation)+1

             !total cross section

             tsigmanucleon(Part%charge)  &
             & =tsigmanucleon(Part%charge)+Part%perweight

             !kinetic energy differential xsection
             tksigmanucleon(Part%charge,ntk)=&
                  &Part%perweight/dEkin  &
                  & +tksigmanucleon(Part%charge,ntk)


             !-----------------------------------------------------------------
             ! Fissum analysis for O16(e,e'p) data, Phys.Rev.C70:034606,2004

             if (Fissum_analysis) then
                if (Part%charge.ne.1) cycle
                call get_phi_Theta(Part%mom(1:3),theta,phi)

                if (.not.(abs(cos(theta)-cos(pi/180.*38.45)).lt.0.02/2.)) cycle

                Emiss_Fissum(ntk)=Part%perweight/dEkin/0.02 &
                                 & +Emiss_Fissum(ntk)
             end if
             !-----------------------------------------------------------------

          end select

       end do ! loop over particle vector
    end do ! loop over particle vector: ensemble loop



    ! summing over several runs at the same energy

    sum_numberofneutrons_tot=numberofneutrons_tot+sum_numberofneutrons_tot
    sum_numberofprotons_tot=numberofprotons_tot+sum_numberofprotons_tot
    sum_numberofANTIprotons_tot=numberofprotons_tot+sum_numberofANTIprotons_tot
    sum_numberofneutrons_out=numberofneutrons_out+sum_numberofneutrons_out
    sum_numberofprotons_out=numberofprotons_out+sum_numberofprotons_out
    sum_numberofANTIprotons_out=numberofprotons_out+sum_numberofANTIprotons_out
    sum_numberofnucleons_out=numberofnucleons_out+sum_numberofnucleons_out

    sum_tsigmapion=sum_tsigmapion+tsigmapion
    sum_tksigmapion=sum_tksigmapion+tksigmapion
    sum_tsigmanucleon=sum_tsigmanucleon+tsigmanucleon
    sum_tksigmanucleon=sum_tksigmanucleon+tksigmanucleon

    sum_Emiss_Fissum=sum_Emiss_Fissum+Emiss_Fissum

    ! Now printout of calculated distributions, controlled by
    ! switch Xsection_analysis

    if (Xsection_analysis) then

       open(10,File='Neutrino_total_Xsection.dat',position='append')
       if (numberofcalls.ne.1) backspace(10)
       write(10,'(10g13.5)') raiseFlagVariable, &
            sum_tsigmapion/real(numberofcalls), &
            sum_tsigmanucleon/real(numberofcalls)
       close(10)

       open(11,File='Neutrino_kinetic_energy_Xsection.dat',position='append')
       if (numberofcalls.ne.1) then
          do i=0,200
             backspace(11)
          end do
       end if
       do i=0,200
          write(11,'(10g13.5)') raiseFlagVariable, &
               (float(i)+0.5)*dEkin,sum_tksigmapion(-1,i)/real(numberofcalls), &
               sum_tksigmapion(0,i)/real(numberofcalls), &
               sum_tksigmapion(1,i)/real(numberofcalls), &
               sum_tksigmanucleon(0,i)/real(numberofcalls), &
               sum_tksigmanucleon(1,i)/real(numberofcalls)
       end do
       close(11)

       open(12,File='Neutrino_Xsection_numbers.dat',position='append')
       if (numberofcalls.ne.1) backspace(12)
       write(12,'(30g13.5)') raiseFlagVariable, &
            sum_numberofprotons_tot/numberofcalls, &
            sum_numberofprotons_out/numberofcalls, &
            sum_numberofnucleons_out(1,:)/numberofcalls, &
            sum_numberofneutrons_tot/numberofcalls, &
            sum_numberofneutrons_out/numberofcalls, &
            sum_numberofnucleons_out(0,:)/numberofcalls, &
            sum_numberofANTIprotons_tot/numberofcalls, &
            sum_numberofANTIprotons_out/numberofcalls, &
            sum_numberofnucleons_out(-1,:)/numberofcalls
       close(12)

    end if

    !--------------------------------------------------------------------------
    ! Fissum analysis for O16(e,e'p)
    if (Fissum_analysis) then
       open(10,File='Neutrino_Emiss_spectrum.dat',position='append')
       if (numberofcalls.ne.1) then
          do i=0,200
             backspace(10)
          end do
       end if
       do i=0,200
          write(10,'(10g13.5)')raiseFlagVariable, &
             & (float(i)+0.5)*dEkin,sum_Emiss_Fissum(i)/real(numberofcalls)
       end do
       close(10)
    end if
    !--------------------------------------------------------------------------


    call WriteHist(hNucleonVacuumMass,file="NucleonVacuumMass.dat",  &
                  & mul = 1.0/numberofcalls)

    if (finalFlag) then
       sum_numberofneutrons_tot=0
       sum_numberofprotons_tot=0
       sum_numberofANTIprotons_tot=0
       sum_numberofneutrons_out=0
       sum_numberofprotons_out=0
       sum_numberofANTIprotons_out=0
       sum_numberofnucleons_out=0
       sum_tsigmapion=0.
       sum_tksigmapion=0.
       sum_tsigmanucleon=0.
       sum_tksigmanucleon=0.
       sum_Emiss_Fissum=0.
    end if

    !==========================================================================
    !==========================================================================
    !==========================================================================

    !==========================================================================
    ! now event analyses with more detailed information on multiplicities,
    ! generations and production process
    !
    ! (1) Setting up particles into the separate first events QE, Delta, DIS, ..
    !
    ! This is done with the help of %firstEvent:
    ! Particles stemming from the same event get in the init the same
    ! %firstEvent entry.
    ! During the run %firstEvent stays constant and is inherited during
    ! collisions.
    !==========================================================================

    firstEvents=getFirstEventRange()

    allocate(lepton(firstEvents(1):firstEvents(2)))
    allocate(struckNuc(firstEvents(1):firstEvents(2)))

    call setToDefault(lepton)
    call setToDefault(struckNuc)

    allocate(events(firstEvents(1):firstEvents(2)))
    allocate(events_QE(firstEvents(1):firstEvents(2)))
    allocate(events_Delta(firstEvents(1):firstEvents(2)))
    allocate(events_highRES(firstEvents(1):firstEvents(2)))
    allocate(events_1piBG(firstEvents(1):firstEvents(2)))
    allocate(events_2piBG(firstEvents(1):firstEvents(2)))
    allocate(events_DIS(firstEvents(1):firstEvents(2)))
    allocate(events_2p2h(firstEvents(1):firstEvents(2)))
    allocate(events_gen0(firstEvents(1):firstEvents(2)))
    allocate(events_gen1(firstEvents(1):firstEvents(2)))
    allocate(events_gen2(firstEvents(1):firstEvents(2)))
    allocate(events_gen3ormore(firstEvents(1):firstEvents(2)))
    allocate(events0(firstEvents(1):firstEvents(2)))

    do i=firstEvents(1),firstEvents(2)
       call event_init(events(i))
       call event_init(events_QE(i))
       call event_init(events_Delta(i))
       call event_init(events_highRES(i))
       call event_init(events_1piBG(i))
       call event_init(events_2piBG(i))
       call event_init(events_DIS(i))
       call event_init(events_2p2h(i))
       call event_init(events_gen0(i))
       call event_init(events_gen1(i))
       call event_init(events_gen2(i))
       call event_init(events_gen3ormore(i))
       call event_init(events0(i))
    end do

    !==========================================================================
    !==========================================================================
    !==========================================================================



    do i=firstEvents(1),firstEvents(2)
       lepton(i)%firstEvent=i
       call get_init_namelist(outLepton_ID=lepton(i)%ID, &
            outLepton_charge=lepton(i)%charge)
       if (.not.neutrinoProdInfo_Get(i, prod_id,lepton(i)%perweight,lepIn_mom, &
            lepton(i)%mom,boson_mom, &
            struckNuc(i)%mom, struckNuc(i)%charge)) then
          call TRACEBACK('error in getting perweight')
       end if


       ! Now acceptance cuts in kinetic energy and angle for outgoing lepton
       ! Events are accepted only if ekin_lepton > threshold energy
       ! and lepton-angle < threshold angle
       if (iand(applyCuts,1)==1) then
          if(.not.lepton_acceptance(lepton(i)%mom)) cycle
       end if

       ! These cuts affect FinalEvents.dat; events which do not fulfill this cut
       ! are not contained in FinalEvents.dat


       !=== put leptons into specific event types (QE, Delta, ..):

       Part=>lepton(i)
       call event_add(events(i),Part)

       !=== put struck nucleon into specific event types:

       Part=>struckNuc(i)
       struckNuc(i)%ID=1

       call event_add(events(i),Part)

       select case (prod_id)
       case (1)
          call event_add(events_QE(i),Part)
       case (2)
          call event_add(events_Delta(i),Part)
       case (3:31)
          call event_add(events_highRES(i),Part)
       case (32:33)
          call event_add(events_1piBG(i),Part)
       case (34)
          call event_add(events_DIS(i),Part)
       case (35:36)
          call event_add(events_2p2h(i),Part)
       case (37)
          call event_add(events_2piBG(i),Part)
       case default
          write(*,*) 'prod_id =', prod_id
          call TRACEBACK('strange prod_id')
       end select

       call event_add(events_gen0(i),Part)
       call event_add(events_gen1(i),Part)
       call event_add(events_gen2(i),Part)
       call event_add(events_gen3ormore(i),Part)

    end do

    !==========================================================================
    ! Now put all other particles into specific event types
    !==========================================================================

    do i=lbound(Particles,dim=1),ubound(Particles,dim=1)
       do j=lbound(Particles,dim=2),ubound(Particles,dim=2)
          if (Particles(i,j)%ID.le.0) cycle

          Part=>Particles(i,j)
          first=Part%firstEvent

          !  These acceptance switches do affect the FinalEvents.dat file
          if (iand(applyCuts,1)==1) then
             if(.not.lepton_acceptance(lepton(first)%mom)) cycle
          end if
          if (iand(applyCuts,2)==2) then
             if (IsBound(Part)) cycle
          end if
          if (iand(applyCuts,4)==4) then
             if (IsBelowThreshold(Part)) cycle
          end if

          call event_add(events(first),Part)
          call event_add(events0(first),Part)

          generation=history_getGeneration(Part%history)
          if (.not.neutrinoProdInfo_Get(Part%firstEvent,prod_id,dummy, &
               lepIn_mom, lep_mom,boson_mom)) then
             call TRACEBACK('error in getting production info')
          end if

          if (generation.ge.max_generation) generation=max_generation

          ! to simplify output individual resonance contributions are all lumped
          ! together in the highRes component, this is controlled by array
          ! K2Hist in initneutrino.f90

          ! now generation is the real particle generation and prod_id contains
          ! the information, in which process the particle was produced
          ! (QE=1, Delta=2, highRES=3:31, BG=32:33)


          select case (prod_id)
          case (1)
             call event_add(events_QE(first),Part)
          case (2)
             call event_add(events_Delta(first),Part)
          case (3:31)
             call event_add(events_highRES(first),Part)
          case (32:33)
             call event_add(events_1piBG(first),Part)
          case (34)
             call event_add(events_DIS(first),Part)
          case (35:36)
             call event_add(events_2p2h(first),Part)
          case (37)
             call event_add(events_2piBG(first),Part)
          case default
             write(*,*) 'prod_id =', prod_id
             call TRACEBACK('strange prod_id')
          end select

          select case (generation)
          case (0)
             call event_add(events_gen0(first),Part)
          case (1)
             call event_add(events_gen1(first),Part)
          case (2)
             call event_add(events_gen2(first),Part)
          case (3)
             call event_add(events_gen3ormore(first),Part)
          case default
             write(*,*) 'generation =', generation
             call TRACEBACK('strange generation')
          end select

       end do
    end do

    !==========================================================================
    ! (2) Use the list "events" to evaluate total cross sections
    !==========================================================================

    call event_sigma(events,sigma,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_QE,sigma_QE,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_Delta,sigma_Delta,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_highRES,sigma_highRES,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_1piBG,sigma_1piBG,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_2piBG,sigma_2piBG,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_DIS,sigma_DIS,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_2p2h,sigma_2p2h,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen0,sigma_gen0,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen1,sigma_gen1,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen2,sigma_gen2,initHists,numberOfCalls, &
         identifier=raiseFlagVariable)
    call event_sigma(events_gen3ormore,sigma_gen3ormore,initHists, &
         numberOfCalls, &
         identifier=raiseFlagVariable)


    if (DoOutChannels) then
       do i=firstEvents(1),firstEvents(2)
          if (CreateSortedPreEvent(events0(i),PreEv%preE)) then
             PreEv%weight = events0(i)%particleList%first%V%perweight
             call PreEvList_INSERT(ListPreEv,PreEv)
          end if
       end do

       open(141,file='OutChannels.FINAL.'//inttochar4(numberOfFinals)//'.dat', status='unknown')
       rewind(141)
       write(141,"(A,f8.4)") "# raiseFlagVariable =",raiseFlagVariable
       write(141,"(A,i6)") "# numberOfCalls = ",numberOfCalls
       call PreEvList_Print(141,ListPreEv,1.0/numberOfCalls)
       close(141)

    end if

    if (Xsection_analysis) then
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities.dat', sigma)
       call PrintVals10(includeQE, raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_QE.dat', sigma_QE)
       call PrintVals10(includeDelta, raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_Delta.dat', sigma_Delta)
       call PrintVals10(includeRes, raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_highRES.dat', sigma_highRES)
       call PrintVals10(include1pi, raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_1piBG.dat', sigma_1piBG)
       call PrintVals10(includeDIS, raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_DIS.dat', sigma_DIS)
       call PrintVals10((include2p2hQE.or.include2p2hDelta), &
            raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_2p2h.dat', sigma_2p2h)

       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_gen0.dat', sigma_gen0)
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_gen1.dat', sigma_gen1)
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_gen2.dat', sigma_gen2)
       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'Neutrino_total_Xsection_multiplicities_gen3ormore.dat', &
            sigma_gen3ormore)
    end if



    ! extra channels for ZeroPion analysis

    if (ZeroPion_analysis) then
       call event_sigma_0pions(events,sigma_0pions,initHists,numberOfCalls, &
            identifier=raiseFlagVariable)
       if (includeQE) call event_sigma_0pions(events_QE,sigma_QE_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (includeDELTA) call event_sigma_0pions(events_Delta,sigma_Delta_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (includeRES) call event_sigma_0pions(events_highRES,sigma_highRES_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (includeDIS) call event_sigma_0pions(events_DIS,sigma_DIS_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (include1pi) call event_sigma_0pions(events_1piBG,sigma_1piBG_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (include2pi) call event_sigma_0pions(events_2piBG,sigma_2piBG_0pions, &
            initHists,numberOfCalls,identifier=raiseFlagVariable)
       if (include2p2hQE .or. include2p2hDelta) call event_sigma_0pions(events_2p2h, &
          & sigma_2p2h_0pions,initHists,numberOfCalls,identifier=raiseFlagVariable)

       call PrintVals10(.true., raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions.dat', sigma_0pions)
       call PrintVals10(includeQE, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_QE.dat', sigma_QE_0pions)
       call PrintVals10(includeDelta, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_Delta.dat', sigma_Delta_0pions)
       call PrintVals10(includeRes, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_highRES.dat', sigma_highRES_0pions)
       call PrintVals10(include1pi, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_1piBG.dat', sigma_1piBG_0pions)
       call PrintVals10(include2pi, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_2piBG.dat', sigma_2piBG_0pions)
       call PrintVals10(includeDIS, raiseFlagVariable, numberOfCalls, &
            'neutrino_0pions_DIS.dat', sigma_DIS_0pions)
       call PrintVals10((include2p2hQE .or. include2p2hDelta), raiseFlagVariable, numberOfCalls,&
            'neutrino_0pions_2p2h.dat', sigma_2p2h_0pions)
    end if


    if (calorimetric_analysis) then
       call event_hadronicEnergy(events,numin,numax,nubin,Ehad_versusNu_hist, &
            & dSigdNu_hist, dSigdEhad_hist, &
            Enumin, Enumax, Enubin, Enurestored_versusEnu, dSigdEnu, dSigdEnurestored)
    end if



    ! Now differential cross sections



    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat
    ! NAME
    ! file diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat
    !
    ! with:
    ! * ZZZ denotes the origin (the first interaction vertex) of the event
    !   (see description in neutrino.EprimeCostplaneSX.ZZZ.dat).
    !
    !   If ZZZ is missing the total cross section
    !   (sum over all primary events) is given, otherwise ZZZ=Delta, DIS, ....
    ! * XXX=000, 001, 002 --- the first, second, third and so on values of
    !   the "raised variable"
    !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
    !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
    ! * HADRON = pion, N, K, K~, Lambda, SigmaResonance, eta
    !   only those are shown, for which the switches in the namelist
    !   "detailed_diff" are set to .true.
    !   possible hadrons are contained in the field particleIDs defined in
    !   AnaEvent.f90
    ! * CHARGE = charge of the outgoing hadron
    !
    ! PURPOSE
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "detailed_diff_output" in
    !   namelist "neutrinoAnalysis" the is set to .true.
    ! * eventtype=3=LowPhoto if switch "dE_switch" in
    !   namelist "LowElePhoto_Analysis" is set to .true.
    !
    ! The file shows the cross sections for ___1-HADRON___  final state
    ! versus kinetic energy of the outgoing hadron
    ! (one HADRON of a given CHARGE and no HADRONs with same flavor, but
    ! different charges;
    !  there could be additional hadrons with different flavor)
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV
    ! * For eventtype=3: microbarns=10^{-30}cm^2/GeV
    ! * All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: kinetic energy of the outgoing HADRON of given CHARGE [GeV]
    ! * #2: dsi/dEkin  xsec
    ! * #3: number of events contributed  (only for internal use, you can
    !   safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE_1X.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE_1X.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! for ___1-HADRON-X___  final state  (one HADRON of a given CHARGE and
    ! any number of  HADRONs of different charges)
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE_2X.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE_2X.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! for ___2-HADRON-X___  final state  (two HADRONs of a given CHARGE and
    ! any number of  HADRONs of different charges)
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEkin_HADRON_charge_CHARGE_MULTI.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE_MULTI.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! for ___MULTI-HADRON___  final state  (at least one HADRON of a given
    ! CHARGE and  any number of  HADRONs of different charges)
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE.dat
    !
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE_1X.dat
    !
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE_2X.dat
    !
    ! file diff_XXX_dSigma_dEcostheta_HADRON_charge_CHARGE_MULTI.dat)
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! dsigma/d(E(1-cosTheta))
    !
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "detailed_diff_output"
    !   in namelist "neutrinoAnalysis" the is set to .true.
    !
    ! The file shows the cross sections dsigma/d(E(1-cosTheta)) versus
    ! E*(1-costheta), where E is the energy (full energy, not kinetic) of the
    ! outgoing hadron, costheta its polar scattering angle (recall here that
    ! in neutrino runs neutrinos are moving along z-direction)
    !
    ! Columns:
    ! * #1: E*(1-cosTheta) [GeV]
    ! * #2: dsi/d(E(1-costheta))  xsec
    ! * #3: number of events contributed  (only for internal use, you can
    !   safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_ZZZ_XXX_dSigma_dTheta_HADRON_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE.dat
    !
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE_1X.dat
    !
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE_2X.dat
    !
    ! file diff_XXX_dSigma_dTheta_HADRON_charge_CHARGE_MULTI.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat but
    ! dsigma/dTheta
    !
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "detailed_diff_output"
    !   in namelist "neutrinoAnalysis" the is set to .true.
    ! * eventtype=3=LowPhoto if switch "dTheta_switch" in
    !   namelist "LowElePhoto_Analysis" is set to .true.
    !
    ! The file shows the cross sections dsigma/dTheta versus Theta,
    ! where Theta is the polar scattering angle (in radians) of the outgoing
    ! hadron (recall here that in neutrino runs neutrinos are moving along
    ! z-direction)
    !
    ! Columns:
    ! * #1: Theta [radians]
    ! * #2: dsi/dTheta  xsec
    ! * #3: number of events contributed  (only for internal use, you can
    !   safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************



    if (detailed_diff_output) then

       ! Here we fill the histograms with the appropriate event types
       ! The called routines are in module AnaEvent

       ! first sum of all events
       call event_dSigma_dE(events,EkinMin,EkinMax,dEkin,&
            & 'diff_'//trim(intToChar(numberofFinals)),numberOfCalls, &
            & dE_hists,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_Multi,hists1X=dE_hists_1X,hists2X=dE_hists_2X)
       call event_dSigma_dOmega(events,dTheta,dPhi,&
            & 'diff_'//trim(intToChar(numberofFinals)),numberOfCalls,  &
            & dTheta_hists, dPhi_hists, dOmega_hists,initHists,sameFileNameIn=.true.)
       call event_dSigma_dEcostheta(events,EcosthetaMin,EcosthetaMax,dEcostheta,'diff_'// &
            & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists,dcostheta_hists,&
            & dEcostheta_hists_MULTI,dcostheta_hists_MULTI, initHists,sameFileNameIn=.true.)

       if (includeQE) then
          call event_dSigma_dE(events_QE,EkinMin,EkinMax,dEkin,&
               & 'diff_QE_'//trim(intToChar(numberofFinals)),numberOfCalls, &
               & dE_hists_QE,initHists,sameFileNameIn=.true., histsMulti=dE_hists_QE_Multi,&
               & hists1X=dE_hists_QE_1X,hists2X=dE_hists_QE_2X)
          call event_dSigma_dOmega(events_QE,dTheta,dPhi,&
               & 'diff_QE_'//trim(intToChar(numberofFinals)),numberOfCalls,  &
               & dTheta_hists_QE, dPhi_hists_QE, dOmega_hists_QE,initHists,  &
               & sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_QE,EcosthetaMin,EcosthetaMax,dEcostheta, &
               &'diff_QE_'//trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_QE, &
               & dcostheta_hists_QE,dEcostheta_hists_QE_MULTI,dcostheta_hists_QE_MULTI, &
               & initHists,sameFileNameIn=.true.)
       end if

       if (includeDELTA) then
          call event_dSigma_dE(events_Delta,EkinMin,EkinMax,dEkin,&
               & 'diff_Delta_'//trim(intToChar(numberofFinals)),numberOfCalls, &
               & dE_hists_Delta,initHists,sameFileNameIn=.true.,  &
               & histsMulti=dE_hists_Delta_Multi,&
               & hists1X=dE_hists_Delta_1X,hists2X=dE_hists_Delta_2X)
          call event_dSigma_dOmega(events_Delta,dTheta,dPhi,&
               & 'diff_Delta_'//trim(intToChar(numberofFinals)),numberOfCalls,  &
               & dTheta_hists_Delta, dPhi_hists_Delta, dOmega_hists_Delta,initHists, &
               & sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_Delta,EcosthetaMin,EcosthetaMax,dEcostheta, &
               & 'diff_Delta_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_Delta, &
               & dcostheta_hists_Delta,&
               & dEcostheta_hists_Delta_MULTI,dcostheta_hists_Delta_MULTI,initHists, &
               & sameFileNameIn=.true.)
       end if

       if (includeRES) then
          call event_dSigma_dE(events_highRES,EkinMin,EkinMax,dEkin,&
               & 'diff_highRES_'//trim(intToChar(numberofFinals)),numberOfCalls,&
               & dE_hists_highRES,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_highRES_Multi,&
               & hists1X=dE_hists_highRES_1X,hists2X=dE_hists_highRES_2X)
          call event_dSigma_dOmega(events_highRES,dTheta,dPhi,&
               & 'diff_highRES_'//trim(intToChar(numberofFinals)),numberOfCalls,&
               & dTheta_hists_highRES, dPhi_hists_highRES, dOmega_hists_highRES,initHists, &
               & sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_highRES,EcosthetaMin,EcosthetaMax,dEcostheta,&
               & 'diff_highRES_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_highRES, &
               & dcostheta_hists_highRES,&
               & dEcostheta_hists_highRES_MULTI,dcostheta_hists_highRES_MULTI,initHists, &
               & sameFileNameIn=.true.)
       end if

       if (include1pi) then
          call event_dSigma_dE(events_1piBG,EkinMin,EkinMax,dEkin,&
               & 'diff_1piBG_'//trim(intToChar(numberofFinals)),&
               & numberOfCalls,dE_hists_1piBG,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_1piBG_Multi,&
               & hists1X=dE_hists_1piBG_1X,hists2X=dE_hists_1piBG_2X)
          call event_dSigma_dOmega(events_1piBG,dTheta,dPhi,'diff_1piBG_'// &
               & trim(intToChar(numberofFinals)), &
               & numberOfCalls,dTheta_hists_1piBG, dPhi_hists_1piBG, dOmega_hists_1piBG, &
               & initHists,sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_1piBG,EcosthetaMin,EcosthetaMax,dEcostheta, &
               & 'diff_1piBG_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_1piBG, &
               & dcostheta_hists_1piBG,&
               & dEcostheta_hists_1piBG_MULTI,dcostheta_hists_1piBG_MULTI,initHists, &
               & sameFileNameIn=.true.)
       end if

       if (include2pi) then
          call event_dSigma_dE(events_2piBG,EkinMin,EkinMax,dEkin,&
               & 'diff_2piBG_'//trim(intToChar(numberofFinals)),&
               & numberOfCalls,dE_hists_2piBG,initHists,sameFileNameIn=.true., &
               &  histsMulti=dE_hists_2piBG_Multi,&
               & hists1X=dE_hists_2piBG_1X,hists2X=dE_hists_2piBG_2X)
          call event_dSigma_dOmega(events_2piBG,dTheta,dPhi,'diff_2piBG_'//&
               &trim(intToChar(numberofFinals)), &
               & numberOfCalls,dTheta_hists_2piBG, dPhi_hists_2piBG, dOmega_hists_2piBG,&
               & initHists,sameFileNameIn=.true.)
          call event_dSigma_dEcostheta(events_2piBG,EcosthetaMin,EcosthetaMax,dEcostheta,&
               & 'diff_2piBG_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_2piBG, &
               & dcostheta_hists_2piBG,&
               & dEcostheta_hists_2piBG_MULTI,dcostheta_hists_2piBG_MULTI,initHists,&
               & sameFileNameIn=.true.)
       end if

       if (includeDIS) then
          call event_dSigma_dE(events_DIS,EkinMin,EkinMax,dEkin,&
               & 'diff_DIS_'//trim(intToChar(numberofFinals)), &
               & numberOfCalls,dE_hists_DIS,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_DIS_Multi,&
               & hists1X=dE_hists_DIS_1X,hists2X=dE_hists_DIS_2X)
          call event_dSigma_dEcostheta(events_DIS,EcosthetaMin,EcosthetaMax,dEcostheta,&
               & 'diff_DIS_'// &
               & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_DIS, &
               & dcostheta_hists_DIS,&
               & dEcostheta_hists_DIS_MULTI,dcostheta_hists_DIS_MULTI,initHists,&
               & sameFileNameIn=.true.)
       end if


       if (include2p2hQE .or. include2p2hDelta) then
          call event_dSigma_dE(events_2p2h,EkinMin,EkinMax,dEkin,&
               & 'diff_2p2h_'//trim(intToChar(numberofFinals)), &
               & numberOfCalls,dE_hists_2p2h,initHists,sameFileNameIn=.true., &
               & histsMulti=dE_hists_2p2h_Multi,&
               & hists1X=dE_hists_2p2h_1X,hists2X=dE_hists_2p2h_2X)
       end if



       call event_dSigma_dE(events_gen0,EkinMin,EkinMax,dEkin,&
            & 'diff_gen0_'//trim(intToChar(numberofFinals)), &
            & numberOfCalls,dE_hists_gen0,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_gen0_Multi,&
            & hists1X=dE_hists_gen0_1X,hists2X=dE_hists_gen0_2X)
       call event_dSigma_dOmega(events_gen0,dTheta,dPhi,'diff_gen0_'//&
            & trim(intToChar(numberofFinals)), &
            & numberOfCalls,dTheta_hists_gen0, dPhi_hists_gen0, dOmega_hists_gen0,initHists, &
            & sameFileNameIn=.true.)
       call event_dSigma_dEcostheta(events_gen0,EcosthetaMin,EcosthetaMax,dEcostheta, &
            & 'diff_gen0_'// &
            & trim(intToChar(numberofFinals)),numberOfCalls,dEcostheta_hists_gen0, &
            & dcostheta_hists_gen0,&
            & dEcostheta_hists_gen0_MULTI,dcostheta_hists_gen0_MULTI,initHists, &
            & sameFileNameIn=.true.)
    end if




    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dEkin_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dEkin_lepton_PPP.dat
    !
    ! with:
    ! * XXX=000, 001, 002 --- the first, second, third and so on values of
    !   the "raised variable"
    !   (e.g. Q^2 for nuXsectionMode=3=dSigmadQs mode,
    !   Elepton for nuXsectionMode=2=dSigmadQsdElepton  and so on)
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    !
    ! PURPOSE
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino if switch "specificEvent_Analysis" in namelist
    !   "neutrinoAnalysis" the is set to .true.
    !   and specific final states in the namelist "nl_specificEvent" are
    !   set to .true.
    !
    ! The file shows the cross sections for a specific final state versus
    ! kinetic energy of the outgoing lepton
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV
    !
    ! Columns:
    ! * #1: kinetic energy of the outgoing lepton  [GeV]
    ! * #2: dsi/dEkin  xsec   per nucleon
    ! * #3: number of events contributed  (only for internal use, you can
    !   safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dQ2_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dQ2_lepton_PPP.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_lepton_PPP.dat but for dsigma/dQ2
    !
    ! The file shows the cross sections for a specific final state versus Q2
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV^2
    ! * All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: Q2 [GeV^2]  squared transfer momentum (Q2 = -q_mu \cdot q^\mu)
    ! * #2: dsi/dQ2  xsec
    ! * #3: number of events contributed  (only for internal use, you can
    !   safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dQ2p_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dQ2p_lepton_PPP.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_lepton_PPP.dat but for dsigma/dQ2
    !
    ! The file shows the cross sections for a final state with 1 mu, 0 pi ,
    ! and (at least) 1 p versus Q2
    ! here Q2 is calculated from the kinematics of the outgoing leading
    ! proton: Q2prot = (mProt - epsB)**2 - mProt**2  &
    !            & + 2*(mNeut - epsB)*(Tp + mProt - mNeut + epsB)
    ! This distribution can be converted into a kinetic energy distribution
    ! of protons in 0pion events
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2/GeV^2
    ! * All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: Q2 [GeV^2]  squared transfer momentum, from proton kinematics
    ! * #2: dsi/dQ2  xsec
    ! * #3: number of events contributed  (only for internal use,
    !   you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dcos_lepton_PPP.dat
    ! NAME
    ! file diff_XXX_dSigma_dcos_lepton_PPP.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dEkin_lepton_PPP.dat but for
    ! dsigma/dcos(theta_l)
    !
    ! The file shows the cross sections for a specific final state
    ! versus cos of the scattering angle of the outgoin lepton
    ! (with respect to neutrino direction)
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2
    !
    ! Columns:
    ! * #1: cos(theta_l)  cos of the scattering angle of the outgoing lepton
    ! * #2: dsi/dcos(theta_l)  xsec
    ! * #3: number of events contributed
    !       (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_d2Sigma_dEdcost_lepton_no_pi.dat
    ! NAME
    ! file diff_XXX_d2Sigma_dEdcost_lepton_no_pi.dat
    !
    ! PURPOSE
    ! Double differential cross section d2sigma/(dE dcostheta) for outgoing
    ! lepton for 0-pi events
    !
    ! The file contains the cross section for the outgoing lepton
    ! versus cos of the scattering angle of the outgoing lepton
    ! (with respect to neutrino direction) and its total energy
    ! for 0-pion events
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2
    ! * All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: total energy of the outgoing lepton (in GeV)
    ! * #2: cos(theta_l)  cos of the scattering angle of the outgoing lepton
    ! * #3: d2si/(dcos(theta_l) dE_l)  dd-xsec
    ! the following columns not implemented yet
    ! * #4: number of events contributed  (only for internal use, you can
    !   safely neglect it)
    ! * #5: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_d2Sigma_dpLdpT_lepton_no_pi.dat
    ! NAME
    ! file diff_XXX_d2Sigma_dpLdpT_lepton_no_pi.dat
    !
    ! PURPOSE
    ! Double differential cross section d2sigma/(dpL dpT) for outgoing
    ! lepton for 0-pi events
    !
    ! The file contains the cross section for the outgoing lepton
    ! vs longitudinal and transverse momentum of the outgoing lepton
    ! for 0-pion events
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2
    ! * For eventtype=5 and process_ID=EM: nanobarns=10^{-33}cm^2
    ! * All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: p_T longitudinal momentum of the outgoing lepton (in GeV)
    ! * #2: p_L longitudinal momentum of the outgoing lepton (in GeV)
    ! * #3: d2si//dpL dpT) dd-xsec
    ! the following columns not implemented yet
    ! * #4: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #5: xsec-statistical-error
    !**************************************************************************


    ! Calculate differential cross sections dsigma/dx for lepton with
    ! specific events
    ! Now loop over specific events, such as 0 pion, 1p + Xn, ...



    if (specificEvent_Analysis) then
       do m=1, max_speEvent           ! loop over specific events
          if (.not.includeSpeEvent(m)) cycle

          do iHist=0, 0         ! loop over first event types, iHist=0: total
		                ! other event types here not implemented

             if (.not.includeHist(iHist)) cycle

             select case(iHist)

             case(0)
                string = 'diff_000'
             case(001)
                string = 'diff_QE'
             case(002)
                string = 'diff_Delta'
             case(003)
                string = 'diff_highRES'
             case(004)
                string = 'diff_1piBG'
             case(005)
                string = 'diff_DIS'
             case(006)
                string = 'diff_2p2h'
             case(008)
                string = 'diff_2piBG'
             case default
                write(*,*) '2p2h_Delta not implemented'
                stop

             end select
   !
   ! The following routine delivers cross sections as a function of lepton
   ! kinematics for the special events m
   ! These cross sections respect the angle and energy cuts
   ! for the outgoing lepton
   ! The called routine is located in module AnaEvent
   !
             call event_dSigma_dLeptonVariables(events,    &
                  maxQ2,binsizeQ2,      &
                  AngleUpperDetectionThresholdDegrees_lepton, &
                  kineticEnergyDetectionThreshold_lepton, &
                  trim(string),numberOfCalls, &
                  dEnu_hist(m,iHist), dElepton_hist(m,iHist),   &
                  dcoslepton_hist(m,iHist), &
                  dQ2lepton_hist(m,iHist),dQ2plepton_hist(m,iHist), &
                  dSigmaMC_EprimeCost_0pi(m,iHist),&
                  dsigmaMC_pLpT_0pi(m,iHist),&
                  initHists, sameFileNameIn=.true., specificEvent=m)

          end do   !iHist
       end do  !m
    end if





    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat
    !
    ! PURPOSE
    ! Similar to diff_XXX_dSigma_dEkin_HADRON_charge_CHARGE.dat
    ! (see notations XXX and CHARGE there) but
    ! gives pion-nucleon invariant mass (W) distribution
    ! (true invariant mass with the sum of the 4-momenta of the outgoing
    ! particles, W^2=(p_1out+p_2out)
    ! as opposed to widely used W2=mN2+2*mN*nu-Q2 from lepton kinematics )
    ! for events with ___1pion___  AND  ___1nucleon___ in the final state
    !
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dW_muon_nucleon_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dW_muon_nucleon_charge_CHARGE.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat but
    ! for the  muon-nucleon invariant mass distribution
    !
    ! NOTES
    ! Note that CHARGE is still a pion charge
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/diff_XXX_dSigma_dW_muon_pion_charge_CHARGE.dat
    ! NAME
    ! file diff_XXX_dSigma_dW_muon_pion_charge_CHARGE.dat
    !
    ! PURPOSE
    ! The same as diff_XXX_dSigma_dW_nucleon_pion_charge_CHARGE.dat but
    ! for the  muon-pion invariant mass distribution
    !
    ! NOTES
    ! Note  that CHARGE is still a pion charge
    !**************************************************************************



    if (include_W_dist) then
       call event_dSigma_dInvMass(events,'diff_'//trim(intToChar(numberofFinals)), &
            & numberOfCalls,&
            & dW_Npi,Wmin_Npi,Wmax_Npi,dW_Npi_hists, dW_mupi,Wmin_mupi,Wmax_mupi,  &
            & dW_mupi_hists, &
            & dW_muN,Wmin_muN,Wmax_muN,dW_muN_hists, initHists ,sameFileNameIn=.true.)
    end if


    ! kinetic energy distributions of nucleons in events with 0 pions
    if (detailed_diff_output .and. ZeroPion_analysis) then

       ! Here we initialize the histograms
       call event_dSigma_dE_0pions(events,EkinMin,EkinMax,dEkin,&
            & 'diff_'//trim(intToChar(numberofFinals)),numberOfCalls,&
            & dE_hists_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_Multi_0pions, &
            & hists1X=dE_hists_1X_0pions, hists2X=dE_hists_2X_0pions)

       if (includeQE) &
            & call event_dSigma_dE_0pions(events_QE,EkinMin,EkinMax,dEkin,&
            & 'diff_QE_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_QE_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_QE_Multi_0pions, &
            & hists1X=dE_hists_QE_1X_0pions, hists2X=dE_hists_QE_2X_0pions)

       if (includeDELTA) &
            & call event_dSigma_dE_0pions(events_Delta,EkinMin,EkinMax,dEkin,&
            & 'diff_Delta_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls ,dE_hists_Delta_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_Delta_Multi_0pions, &
            & hists1X=dE_hists_Delta_1X_0pions, hists2X=dE_hists_Delta_2X_0pions)

       if (includeRES) &
            & call event_dSigma_dE_0pions(events_highRES,EkinMin,EkinMax,dEkin,&
            & 'diff_highRES_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_highRES_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_highRES_Multi_0pions, &
            & hists1X=dE_hists_highRES_1X_0pions, hists2X=dE_hists_highRES_2X_0pions)

       if (includeDIS) &
            & call event_dSigma_dE_0pions(events_DIS,EkinMin,EkinMax,dEkin,&
            & 'diff_DIS_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_DIS_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_DIS_Multi_0pions, &
            & hists1X=dE_hists_DIS_1X_0pions, hists2X=dE_hists_DIS_2X_0pions)

       if (include1pi) &
            & call event_dSigma_dE_0pions(events_1piBG,EkinMin,EkinMax,dEkin,&
            & 'diff_1piBG_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_1piBG_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_1piBG_Multi_0pions, &
            & hists1X=dE_hists_1piBG_1X_0pions, hists2X=dE_hists_1piBG_2X_0pions)

       if (include2pi) &
            & call event_dSigma_dE_0pions(events_2piBG,EkinMin,EkinMax,dEkin,&
            & 'diff_2piBG_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_2piBG_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_2piBG_Multi_0pions, &
            & hists1X=dE_hists_2piBG_1X_0pions, hists2X=dE_hists_2piBG_2X_0pions)

       if (include2p2hQE .or. include2p2hDelta) &
            & call event_dSigma_dE_0pions(events_2p2h,EkinMin,EkinMax,dEkin,&
            & 'diff_2p2h_'//trim(intToChar(numberofFinals)),&
            & numberOfCalls, dE_hists_2p2h_0pions,initHists,sameFileNameIn=.true., &
            & histsMulti=dE_hists_2p2h_Multi_0pions, &
            & hists1X=dE_hists_2p2h_1X_0pions, hists2X=dE_hists_2p2h_2X_0pions)

    end if

    !multiplicity output
    write(Prefix_MultAna,'(f8.4)') raiseFlagVariable
    do i=firstEvents(1),firstEvents(2)
       call Multiplicity_AddEvent(events(i))
    end do


    call Multiplicity_Write(Prefix_MultAna)

    if (outputEvents) then
       write(*,*) 'Writing events to file'

       !***********************************************************************
       !****o* neutrinoAnalysis/FinalEvents.dat
       ! NAME
       ! file FinalEvents.dat
       ! PURPOSE
       ! The file contains positions and four-vectors of all outgoing
       ! particles from an event.
       !
       ! For a cross-section construction all of the events have to be
       ! weighted with the 'perweight' values in col. 5. The final
       ! differential cross section dsigma/dE, for example, then is obtained
       ! by first binning all the events with respect to the energy. Then
       ! the perweights have to be summed within each bit. To obtain the
       ! differential cross section the sum of weights in each bin has to be
       ! divided by the bin width. In addition, the number of runs at same
       ! energy (>1 for better statistics) has to be divided out.
       ! The resulting cross section is per nucleon (1/A). Its units are
       ! 10^{-38} cm^2 for neutrinos and 10^{-33} cm^2 for electrons.
       !
       ! The output of this file is switched on by the switch 'outputEvents'.
       !
       ! Columns:
       ! * #1: run number (from 1 to num_runs_SameEnergy)
       ! * #2: event number (from 1 to ... less then target_A*numEnsembles)
       ! * #3: ID of the outgoing particle (for antiparticles times factor -1)
       ! * #4: charge of the outgoing particle
       ! * #5: perweight of the event
       ! * #6-#8: position of the outgoing particle (=0 for outgoing lepton)
       ! * #9-#12: 4-momentum of the outgoing particle in GeV
       ! * #13: history of the outgoing particle (see history.f90 for def)
       ! * #14: production_ID (type of the first event: 1=QE, 2=Delta, 34=DIS)
       ! * #15: incoming neutrino energy
       !
       ! NOTES
       ! There is always in each event a particle with weight 0. This is the
       ! nucleon on which the initial interaction happened.
       !
       ! In the case of an initial 2p2h process the second initial-state nucleon
       ! is not written out. It is chosen to be at the same place as the first
       ! initial nucleon (the one with weight=0), with a randomly chosen
       ! momentum in the Fermi-sea.
       !
       ! For large target_A, numEnsembles and num_runs_SameEnergy, this file
       ! may become very large, e.g.:
       ! * target_A=12, 1000ens x 20runs, QE and Delta results in 41 Mb
       ! * target_A=56, 2000ens x 5runs, QE,Delta,highRES,1piBG,DIS results in
       !   500 Mb
       !***********************************************************************

       !***********************************************************************
       !****g* neutrinoAnalysis/production_ID
       ! SOURCE
       !
       ! PURPOSE
       ! contains info on the very first neutrino-interaction with the nucleus:
       ! * 1: nucleon (QE)
       ! * 2-31: non-strange baryon resonance (as in IdTable)
       ! * 32: pi neutron-background  (e.g. nu + n -> mu + pi+ + n)
       ! * 33: pi proton-background   (e.g. nu + n -> mu + pi0 + p)
       ! * 34: DIS
       ! * 35: 2p2h QE
       ! * 36: 2p2h Delta
       ! * 37: two pion background
       !***********************************************************************

       if (doPipe) then
          filename=fileNamePipe
          call event_dump(numberOfCalls,events,filename,writeNeutrinoProdID=.true.,doPipe=.TRUE.)
       else
          filename='FinalEvents.dat'
          call event_dump(numberOfCalls,events,filename,writeNeutrinoProdID=.true.)
       end if

       ! for Sigma_MC run huge files, and they are not necessary
       ! because the same info can be obtained from 'FinalEvents.dat' using column 14

       !filename='FinalEvents_QE.dat'
       !call event_dump(numberOfCalls,events_QE,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_Delta.dat'
       !call event_dump(numberOfCalls,events_Delta,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_highRES.dat'
       !call event_dump(numberOfCalls,events_highRES,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_1piBG.dat'
       !call event_dump(numberOfCalls,events_1piBG,filename,writeNeutrinoProdID=.true.)

       !filename='FinalEvents_DIS.dat'
       !call event_dump(numberOfCalls,events_DIS,filename,writeNeutrinoProdID=.true.)

	   !filename='FinalEvents_2p2h.dat'
       !call event_dump(numberOfCalls,events_2p2h,filename,writeNeutrinoProdID=.true.)
    end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! RECONSTRUCTION OF ENERGY AND Q^2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (reconstruct_neutrino_energy .and. specificEvent_Analysis) then

       do j=lbound(events,dim=1),ubound(events,dim=1)
          if (.not.neutrinoProdInfo_Get(j,prod_id,perweight,lepIn_mom,lep_mom,  &
            & boson_mom)) then
             write(*,*) j,prod_id,perweight
             call TRACEBACK('error in getting perweight, stop')
          end if
!
! All the files produced are not affected by cuts on outgoing lepton angle and energy
!

          do m=1, max_speEvent
             if (.not.includeSpeEvent(m)) cycle
             if (.not.IfPass_SpecificEvent(m,events(j))) cycle
             Q2real=-abs4Sq(boson_mom)
             Enureal=lep_mom(0)+boson_mom(0)

             select case (m)
             case (1) ! events with 0 pions,
                Q2rec=CCQE_recQs(lep_mom)           !CCQE QE-like reconstruction
                Enurec=CCQE_recEnergy(lep_mom)

             case (2,9:18) ! events with 0 pions and 1 proton and X neutrons,
                Q2rec=K2K_recQs(lep_mom)                 ! K2K QE-like reconstruction
                Enurec=K2K_recEnergy(lep_mom)

             case (3,4) ! events with 1 pi+,  the same for 1pi0
                Q2rec=CCQE_recQs_Delta(lep_mom)         ! CCQE energy reconstruction assuming
                Enurec=CCQE_recEnergy_Delta(lep_mom)    ! Delta mass

             case (5,6) ! events with at least 1 pi0 (as K2K),   the same for at least 1pi+
                Q2rec=K2K_recQs(lep_mom,1.483)               ! K2K  energy reconstruction
                Enurec=K2K_recEnergy(lep_mom,1.483)          ! assuming DIS with W=1.483

             case default
             end select


             call AddHist(H_Q2_real(m,0),H_Q2_real(m,K2Hist(prod_id)), Q2real,perweight/ &
                & float(num_runs_sameEnergy))
             call AddHist(H_Q2_rec(m,0), H_Q2_rec(m,K2Hist(prod_id)),  Q2rec, perweight/ &
                & float(num_runs_sameEnergy))

             if (nuEXP>0) call AddHist(H_enu_real(m,0),H_enu_real(m,K2Hist(prod_id)), &
                  & Enureal,perweight/float(num_runs_sameEnergy))
             call AddHist(H_enu_rec(m,0), H_enu_rec(m,K2Hist(prod_id)),  Enurec, perweight/ &
                & float(num_runs_sameEnergy))

             call AddHist2D(H_Q2_rec_versus_real(m,0),H_Q2_rec_versus_real(m,K2Hist(prod_id)), &
                  & (/Q2real,Q2rec/),  perweight/float(num_runs_sameEnergy))
             if (nuEXP>0) call AddHist2D(H_enu_rec_versus_real(m,0), &
                  & H_enu_rec_versus_real(m,K2Hist(prod_id)), &
                  & (/Enureal,Enurec/),perweight/float(num_runs_sameEnergy))

             !! add muon survival in Oscmumu
             !! add oscillated electron appearance for various CP violation phases
             if (isOSC()) then
                L=OSCLENGTH(nuEXP)
                call oscillationProbability(Enureal,L,0.,Posc_mumu,Posc_mue,Posc_mue_max, &
                   & Posc_mue_antimax)

                call AddHist(Oscmumu_enu_real(m,0),Oscmumu_enu_real(m,K2Hist(prod_id)), &
                     &  Enureal,Posc_mumu*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmumu_enu_rec(m,0), Oscmumu_enu_rec(m,K2Hist(prod_id)),  &
                     & Enurec, Posc_mumu*perweight/float(num_runs_sameEnergy))

                call AddHist(Oscmue_enu_real(m,0),Oscmue_enu_real(m,K2Hist(prod_id)), &
                     & Enureal,Posc_mue*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmue_enu_rec(m,0), Oscmue_enu_rec(m,K2Hist(prod_id)),  &
                     & Enurec, Posc_mue*perweight/float(num_runs_sameEnergy))

                call AddHist(Oscmuemax_enu_real(m,0),Oscmuemax_enu_real(m,K2Hist(prod_id)), &
                     & Enureal,Posc_mue_max*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmuemax_enu_rec(m,0), Oscmuemax_enu_rec(m,K2Hist(prod_id)),  &
                     & Enurec, Posc_mue_max*perweight/float(num_runs_sameEnergy))

                call AddHist(Oscmueantimax_enu_real(m,0), &
                     & Oscmueantimax_enu_real(m,K2Hist(prod_id)),&
                     & Enureal,Posc_mue_antimax*perweight/float(num_runs_sameEnergy))
                call AddHist(Oscmueantimax_enu_rec(m,0), &
                     & Oscmueantimax_enu_rec(m,K2Hist(prod_id)),&
                     & Enurec, Posc_mue_antimax*perweight/float(num_runs_sameEnergy))
             end if

          end do ! m

       end do ! j


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! energy reconstruction output in files

       do m=1, max_speEvent
          if (.not.includeSpeEvent(m)) cycle
          call SpecificEvent_Name(m,filename1)
          do iHist=0, max_Hist
             if (.not.includeHist(iHist)) cycle

             call writeHist(H_Q2_real(m,iHist), file='reconstruction_Q2real_'&
                & //trim(filename1)//"."//trim(intToChar(iHist))//'.dat')
             call writeHist(H_Q2_rec(m,iHist),  file='reconstruction_Q2rec_'// &
                     &  trim(filename1)//&
                & "."//trim(intToChar(iHist))//'.dat')
             if (nuEXP>0) call writeHist(H_enu_real(m,iHist), &
                  & file='reconstruction_Enureal_'//trim(filename1)//"."  &
                  & //trim(intToChar(iHist))//'.dat')
             call writeHist(H_enu_rec(m,iHist), file='reconstruction_Enurec_'  &
                  & //trim(filename1) &
                  &//"."//trim(intToChar(iHist))//'.dat')
             call writeHist2D_Gnuplot(H_Q2_rec_versus_real(m,iHist), &
                  & file='reconstruction_Q2_rec_versus_real_'    &
                  &     //trim(filename1)//"."// &
                  & trim(intToChar(iHist))//'.dat')
             if (nuEXP>0)  call writeHist2D_Gnuplot&
                  &(H_enu_rec_versus_real(m,iHist), &
                  & file='reconstruction_Enu_rec_versus_real_'&
                  & //trim(filename1)//"."// &
                  & trim(intToChar(iHist))//'.dat')
          end do ! iHist
       end do ! m


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! oscillation versus true and reconstructed energy  output in files

       do m=1, max_speEvent
          if (.not.includeSpeEvent(m) .or. (m>6 .AND. m /= 22)) cycle
          call SpecificEvent_Name(m,filename1)
          do iHist=0, max_Hist
             if (.not.includeHist(iHist)) cycle

             if (isOSC()) then

              call writeHist(Oscmumu_enu_real(m,iHist), &
                   & file='oscillations_mumu_real_'//trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmumu_enu_rec(m,iHist), &
                   & file='oscillations_mumu_rec_'//trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmuemax_enu_real(m,iHist), &
                   & file='oscillations_mue_max_real_'//trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmuemax_enu_rec(m,iHist), &
                   & file='oscillations_mue_max_rec_'//trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmue_enu_real(m,iHist), &
                   & file='oscillations_mue_real_'//trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmue_enu_rec(m,iHist), &
                   & file='oscillations_mue_rec_'//trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmueantimax_enu_real(m,iHist), &
                   & file='oscillations_mue_antimax_real_'//  &
                   & trim(filename1)//"."// &
                   & trim(intToChar(iHist))//'.dat')
              call writeHist(Oscmueantimax_enu_rec(m,iHist), &
                   & file='oscillations_mue_antimax_rec_'//trim(filename1)// &
                   & "."//trim(intToChar(iHist))//'.dat')

             end if
          end do ! iHist
       end do ! m

    end if ! reconstruct_neutrino_energy .and. specificEvent_Analysis






    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Enureal_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Enureal_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for true neutrino energy
    !
    ! The file shows the event distribution (flux times xsec) versus true
    ! neutrino energy for a specific final state
    ! The file is produced only for the runs with a neutrino flux
    ! (in the namelist "neutrino_induced" nuXsectionMode > 10, nuExp>0)
    !
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! Columns:
    ! * #1: true neutrino energy  [GeV]
    ! * #2: event distribution: normalized flux times xsec
    ! * #3: number of events contributed
    !      (only for internal use, you can safely neglect it)
    ! * #4: statistical-error of #2
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Enurec_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Enurec_PPP.ZZZ.dat
    !
    ! with:
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    ! * ZZZ=000 - 008 is the origin (the first interaction vertex) of the event
    !   (see description in  neutrino.EprimeCostplaneXS.ZZZ.dat)
    !
    ! PURPOSE
    ! The file shows the event distribution (normalized flux times xsec)
    ! versus reconstructed neutrino energy (see arXiv:1208.3678 [nucl-th])
    ! for a specific final state
    !
    ! Reconstruction method depends on a specific process;
    ! * generally for no_pion events it is based on QE-like kinematics
    ! * and for pion events it is based on on-shell-Delta-creation assumption
    !
    ! The file is produced in runs with:
    ! * eventtype=5=neutrino
    !   if switch "reconstruct_neutrino_energy" is set to .true.
    !   and switch "specificEvent_Analysis" in namelist "neutrinoAnalysis"
    !   is set to .true.
    !   and specific final states in the namelist "nl_specificEvent" are set
    !   to .true.
    !
    !   The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! Units:
    ! * For eventtype=5  and process_ID=CC and NC: 10^{-38} cm^2/GeV
    ! * For eventtype=5  and process_ID=EM: nanobarns=10^{-33}cm^2/GeV
    ! * All x-sec per particle (1/A)
    !
    ! Columns:
    ! * #1: reconstructed neutrino energy  [GeV]
    ! * #2: event distribution: flux times xsec)
    ! * #3: number of events contributed
    !      (only for internal use, you can safely neglect it)
    ! * #4: statistical-error of #2
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Enu_rec_versus_real_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Enu_rec_versus_real_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but double-differential
    !
    ! The file shows the 2-D density of the flux times xsec
    ! (see arXiv:1208.3678 [nucl-th])
    ! versus true and reconstructed neutrino energies for a specific final state
    ! The file is produced only for the runs with a neutrino flux
    ! (in the namelist "neutrino_induced" nuXsectionMode > 10, nuExp>0)
    !
    ! Units:
    ! * For event_type=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * All xsec per particle (1/A)
    !
    ! Columns:
    ! * #1: true neutrino energy  [GeV]
    ! * #2: reconstructed neutrino energy  [GeV]
    ! * #3: flux-folded xsec
    ! * #4: number of events contributed
    !      (only for internal use, you can safely neglect it)
    ! * #5: xsec-statistical-error
    !
    ! NOTES
    ! up to normalization this is "migration matrix" between true and
    ! reconstructed energies
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Q2real_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Q2real_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for the true Q2
    !
    ! The file contains the flux averaged cross section dsigma/dQ2,
    ! i.e. integral  \int \Phi(E) dsigma/dQ2 (E) dE,
    ! versus true Q2 for a specific final state
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC:  10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM (one can run it, but makes no
    !   physical sense): nanobarns=10^{-33}cm^2/GeV^2
    ! * All xsec per particle (1/A)
    !
    ! Columns:
    ! * #1: true Q2  [GeV^2] (squared momentum transfer Q2 = - q_mu \cdot q^mu)
    ! * #2: flux-averaged dsigma/dQ2true
    ! * #3: number of events contributed  (only for internal use,
    !       you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Q2rec_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Q2rec_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for the reconstructed Q2
    ! (see arXiv:1208.3678 [nucl-th])
    !
    ! The file contains the flux averaged cross section dsigma/dQ2,
    ! i.e. integral  \int \Phi(E) dsigma/dQ2 (E) dE,
    ! versus reconstructed Q2 for a specific final state
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^2
    ! * For eventtype=5 and process_ID=EM (one can run it, but makes no
    !   physical sense): nanobarns=10^{-33}cm^2/GeV^2
    ! * All xsec per particle (1/A)
    !
    ! Columns:
    ! * #1: reconstructed Q2  [GeV^2]  = - q_mu \cdot q^mu)
    ! * #2: flux-averaged dsigma/dQ2rec
    ! * #3: number of events contributed
    !       (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !**************************************************************************
    !****o* neutrinoAnalysis/reconstruction_Q2_rec_versus_real_PPP.ZZZ.dat
    ! NAME
    ! file reconstruction_Q2_rec_versus_real_PPP.ZZZ.dat
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but double-differential
    ! in Q2
    !
    ! The file shows the dsigma/dQ2true dQ2rec distribution
    ! versus true and reconstructed Q^2 for a specific final state
    ! The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! Units:
    ! * For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV^4
    ! * For eventtype=5 and process_ID=EM (one can run it, but makes no
    !   physical sense): nanobarns=10^{-33}cm^2/GeV^4
    ! * All cross sections are given per nucleon (1/A)
    !
    ! Columns:
    ! * #1: true Q2  [GeV^2]
    ! * #2: reconstructed Q2  [GeV^2]
    ! * #3: dsigma/dQ2true dQ2rec
    ! * #4: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #5: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/oscillations_UUU_real_PPP.ZZZ.dat
    ! NAME
    ! file oscillations_UUU_real_PPP.ZZZ.dat
    !
    ! with:
    ! * UUU =  "mumu", "mue", "mue_max", "mue_antimax"
    !   for oscillation of a given experimental flux
    !   for muon_neutrino survival ("mumu"),
    !   mu_e appearance with delta_CP=0 ("mue"),
    !   delta_CP=pi ("mue_max"),
    !   delta_CP=-pi ("mue_antimax")
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    ! * ZZZ=000 - 008 is the origin (the first interaction vertex) of the event
    !   (see description in neutrino.EprimeCostplaneXS.ZZZ.dat)
    !
    !   The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! PURPOSE
    ! Similar to reconstruction_Enureal_PPP.ZZZ.dat, but for oscillated
    ! neutrino flux
    !
    ! Columns:
    ! * #1: true neutrino energy  [GeV]
    ! * #2: flux-folded xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************


    !**************************************************************************
    !****o* neutrinoAnalysis/oscillations_UUU_rec_PPP.ZZZ.dat
    ! NAME
    ! file oscillations_UUU_rec_PPP.ZZZ.dat
    !
    ! with:
    ! * UUU =  "mumu", "mue", "mue_max", "mue_antimax"
    !   for oscillation of a given experimental flux
    !   for muon_neutrino survival ("mumu"),
    !   mu_e appearance with delta_CP=0 ("mue"),
    !   delta_CP=pi ("mue_max"),
    !   delta_CP=-pi ("mue_antimax")
    ! * PPP = no_pi, p_Xn_no_pi, piplus, pi0, ....
    !   (see namelist nl_specificEvent for the full list)
    !   standing for the specific final state under consideration
    ! * ZZZ=000 - 008 is the origin (the first interaction vertex) of the event
    !   (see description in neutrino.EprimeCostplaneXS.ZZZ.dat)
    !   The file is not affected by cuts on outgoing lepton's angle and energy
    !
    ! PURPOSE
    ! Similar to reconstruction_Enurec_PPP.ZZZ.dat, but for oscillated
    ! neutrino flux
    !
    ! Columns:
    ! * #1: reconstructed neutrino energy  [GeV]
    ! * #2: event rate: flux times xsec
    ! * #3: number of events contributed  (only for internal use, you can safely neglect it)
    ! * #4: xsec-statistical-error
    !**************************************************************************

    !##########################################################################
    ! energy-reconstruction/oscillation  ANALYSIS END
    !##########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! (3) Clear the list "events" to clear the memory
    do i=firstEvents(1),firstEvents(2)
       call event_clear(events(i))
       call event_clear(events_QE(i))
       call event_clear(events_Delta(i))
       call event_clear(events_highRES(i))
       call event_clear(events_1piBG(i))
       call event_clear(events_2piBG(i))
       call event_clear(events_DIS(i))
       call event_clear(events_2p2h(i))
       call event_clear(events_gen0(i))
       call event_clear(events_gen1(i))
       call event_clear(events_gen2(i))
       call event_clear(events_gen3ormore(i))
    end do
    deallocate(events)
    deallocate(events_QE)
    deallocate(events_Delta)
    deallocate(events_highRES)
    deallocate(events_1piBG)
    deallocate(events_2piBG)
    deallocate(events_DIS)
    deallocate(events_2p2h)
    deallocate(events_gen0)
    deallocate(events_gen1)
    deallocate(events_gen2)
    deallocate(events_gen3ormore)
    deallocate(events0)

    if (finalflag) then
       call Multiplicity_Reset
       call PreEvList_CLEAR(ListPreEv)
    end if

    !    call neutrinoProdInfo_clear
    call neutrinoInfoStorage_clear

    if (finalFlag) then
       numberOfCalls=0
       numberOfFinals=numberOfFinals+1
    end if


    write(*,*) '################### NEUTRINO ANALYSIS FINISHED #######################'

  end subroutine neutrino_Analyze






  !****************************************************************************
  !****f* neutrinoAnalysis/IsBound
  ! NAME
  ! logical function IsBound(part)
  ! PURPOSE
  ! return whether particle is bound (.true.) or not (.false.)
  ! This cut affects FinalEvents.dat
  !****************************************************************************
  logical function IsBound(part)
    use particleDefinition
    use potentialMain, only: potential_LRF
    use nucleusDefinition
    use nucleus, only: getTarget
    use vector, only: absVec

    type(particle), intent(in) :: part
    type(tnucleus),pointer :: TargetNuc
    real :: NucRadius

    IsBound=.false.

    targetNuc => getTarget()
    NucRadius=targetNuc%radius

    if (kineticEnergy(part)+potential_LRF(part).lt.0) IsBound=.true.

!   if (absVec(part%pos).lt.radialScale*NucRadius) IsBound=.true.

  end function IsBound



  !****************************************************************************
  !****f* neutrinoAnalysis/IsBelowThreshold
  ! NAME
  ! logical function IsBelowThreshold(part)
  ! PURPOSE
  ! Returns .true. when a particle is below a given detection threshold
  ! This routine can be used to remove particles with kinetic energies
  ! and angles below specified detection thresholds.
  ! This cut affects FinalEvents.dat
  !****************************************************************************
  logical function IsBelowThreshold(part)
    use particleDefinition
    use vector, only: absVec
    use degRad_conversion, only: radian

    type(particle), intent(in) :: part
    real :: cos_theta

    IsBelowThreshold=.false.

    cos_theta = part%mom(3)/absVec(part%mom(1:3))

     !kinetic energy and angle threshold for lepton
    if (part%id.eq.902 .and. ((part%mom(0)-part%mass) .lt.                &
         & kineticEnergyDetectionThreshold_lepton .or.                         &
         & cos_theta .lt.                                                      &
         & cos(radian(AngleUpperDetectionThresholdDegrees_lepton))))           &
         & IsBelowThreshold=.true.

    !kinetic energy and angle threshold for nucleons
    if (part%id.eq.1 .and. ((part%mom(0)-part%mass) .lt.                  &
         & kineticEnergyDetectionThreshold_nucleon .or.                        &
         & cos_theta .lt.                                                      &
         & cos(radian(AngleUpperDetectionThresholdDegrees_nucleon))))          &
         & IsBelowThreshold=.true.

    !kinetic energy and angle threshold for charged pions
    if (part%id.eq.101 .and. part%charge.ne.0.and.((part%mom(0)-part%mass) &
         & .lt. kineticEnergyDetectionThreshold_chargedpion .or.               &
         & cos_theta .lt.                                                      &
         & cos(radian(AngleUpperDetectionThresholdDegrees_chargedpion))))      &
         & IsBelowThreshold=.true.

    !kinetic energy and angle threshold for neutral pions
    if (part%id.eq.101.and.part%charge.eq.0.and.((part%mom(0)-part%mass).lt.&
         & kineticEnergyDetectionThreshold_neutralpion .or.                    &
         & cos_theta .lt.                                                      &
         & cos(radian(AngleUpperDetectionThresholdDegrees_neutralpion))))      &
         & IsBelowThreshold=.true.

  end function IsBelowThreshold

  !****************************************************************************
  !****f* neutrinoAnalysis/lepton_acceptance
  ! NAME
  ! logical function lepton_acceptance
  ! PURPOSE
  ! Returns .true. when the lepton is above a given kinetic energy and
  ! below a given angle, i.e. it is accepted
  !
  ! This routine can be used to remove all events with lepton kinetic energies
  ! below and angles below specified detection thresholds. If an event is not
  ! accepted it is not added to the full event in subroutine neutrino_Analyze.
  ! This cut affects FinalEvents.dat
  !****************************************************************************
  logical function lepton_acceptance(lep_mom)
    use vector, only: absVec
    use degRad_conversion, only: radian
    use minkowski, only: abs4Sq

    intent(in) :: lep_mom
    real :: cos_theta,ekin_lepton
    real, dimension (0:3) :: lep_mom

    lepton_acceptance = .true.

    cos_theta = lep_mom(3)/absVec(lep_mom(1:3))
    ekin_lepton = lep_mom(0) - sqrt(max(0.,abs4Sq(lep_mom)))

    if (ekin_lepton < kineticEnergyDetectionThreshold_lepton) &
         & lepton_acceptance = .false.
    if ( cos_theta    &
         & < cos(radian(AngleUpperDetectionThresholdDegrees_lepton))) &
         & lepton_acceptance = .false.

  end function lepton_acceptance




  !****************************************************************************
  !****s* neutrinoAnalysis/oscillationProbability
  ! NAME
  ! subroutine oscillationProbability(Enu,L,deltaCP,
  ! Posc_mumu,Posc_mue,Posc_mue_max,Posc_mue_antimax)
  ! PURPOSE
  ! Prepare Calculation of Oscillation Probability
  !****************************************************************************
  subroutine oscillationProbability(Enu,L,deltaCP, &
       Posc_mumu,Posc_mue,Posc_mue_max,Posc_mue_antimax)


    implicit none

    real, intent(in) :: enu ! neutrino energy  (GeV)
    real, intent(in)  :: L ! distance (in km) from the near to the far detector !
                           ! defined in initNeutrino
    real, optional, intent(in) :: deltaCP ! CP violating phase
    real, intent(out) :: Posc_mumu, Posc_mue, Posc_mue_max,Posc_mue_antimax

    real ::  sin_deltaCP, Posc_mutau, pid


    if (present(deltaCP)) then
       sin_deltaCP=sin(deltaCP)
    else
       sin_deltaCP=0.0
    end if

    !! In the following oscillation probabilities are calculated for
    !! deltaCP = 0,+pi/2 and -pi/2

    !! first switch for neutrino-antineutrino

    if (process_ID < -1) then
       pid = -1
    else if (process_ID > +1) then
       pid = +1
    else if (process_ID == +1 .or. process_ID == -1) then
       call TRACEBACK('oscillation for electrons makes no sense')
    end if

    call Posc(pid,L,Enu,sin_deltaCP,Posc_mue,Posc_mutau,Posc_mumu)

    sin_deltaCP=1.0
    call Posc(pid,L,Enu,sin_deltaCP,Posc_mue)
    Posc_mue_max = Posc_mue

    sin_deltaCP=-1.0
    call Posc(pid,L,Enu,sin_deltaCP,Posc_mue)
    Posc_mue_antimax = Posc_mue

  end subroutine oscillationProbability


  !****************************************************************************
  !****s* neutrinoAnalysis/PrintVals10
  ! NAME
  ! PURPOSE
  ! abbreviation used for many outputfiles: printing the vals stored in sigma,
  ! including its error
  !****************************************************************************
  subroutine PrintVals10(flag, raiseFlagVariable, numberOfCalls, fName, sigma)

    logical, intent(in) :: flag
    real, intent(in) :: raiseFlagVariable
    integer, intent(in) :: numberOfCalls
    character*(*), intent(in) :: fName
    real, dimension(1:dimSigma,1:2), intent(in) :: sigma

    if (.not.flag) return

    open(10,File=fName,position='append')
    if (numberOfCalls.ne.1) backspace(10)
    if (numberOfCalls.gt.1) then
       write(10,'(300E13.5)') raiseFlagVariable,sigma(2:,1)/float(numberOfCalls), &
            sqrt(max(0.,((sigma(2:,2)-sigma(2:,1)**2/float(numberOfCalls)) &
            / (float(numberOfCalls-1)*float(numberOfCalls)))))
    else
       write(10,'(300E13.5)') raiseFlagVariable,sigma(2:,1)
    end if
    close(10)

  end subroutine PrintVals10


  !****************************************************************************
  !****s* neutrinoAnalysis/PrintHeader10
  ! NAME
  ! PURPOSE
  ! abbreviation used for many outputfiles: printing a header line
  !****************************************************************************
  subroutine WriteHeader10(flag, fName)

    logical, intent(in) :: flag
    character*(*), intent(in) :: fName

    character*(*), parameter :: header = &
       "(' # total Xsections for defined finalstates (see sigma_0pions.dat)',/, &
       & ' # next 120 columns: error',/, &
       & ' # order of Xsections as in sigma_0pions.dat')"

    if (.not.flag) return

    open(10,file=fName)
    write(10,header)
    close(10)

  end subroutine WriteHeader10


end module neutrinoAnalysis
