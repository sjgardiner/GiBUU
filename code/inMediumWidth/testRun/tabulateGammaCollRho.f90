program tabulateGammaCollRho
  use inputGeneral, only: readinputGeneral
  use particleProperties, only: initParticleProperties
  use version, only: PrintVersion
  use inMediumWidth_rho
  
  call PrintVersion

  call readinputGeneral
  call initParticleProperties

  call tabulate_GammaColl_rho
  
end program tabulateGammaCollRho

