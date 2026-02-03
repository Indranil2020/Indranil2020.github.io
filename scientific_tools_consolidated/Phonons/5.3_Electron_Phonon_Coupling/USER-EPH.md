# USER-EPH

## Official Resources
- Homepage: https://github.com/LLNL/USER-EPH
- Source Repository: https://github.com/LLNL/USER-EPH
- License: MIT License

## Overview
USER-EPH is a LAMMPS extension package for incorporating electron-phonon coupling effects in molecular dynamics simulations. Developed at Lawrence Livermore National Laboratory, it enables simulation of energy exchange between electronic and ionic subsystems.

**Scientific domain**: Electron-phonon coupling, two-temperature model, radiation damage  
**Target user community**: LAMMPS users studying electron-phonon energy transfer

## Theoretical Methods
- Two-temperature model (TTM)
- Electronic heat diffusion
- Electron-phonon coupling parameter
- Langevin thermostat for electrons
- Energy exchange dynamics
- Non-equilibrium electron-ion dynamics

## Capabilities (CRITICAL)
- Electron-phonon energy transfer
- Two-temperature MD simulations
- Electronic heat diffusion
- Radiation damage simulations
- Laser-matter interactions
- Non-equilibrium dynamics
- LAMMPS integration

## Key Strengths

### LAMMPS Integration:
- Standard LAMMPS package
- Easy installation
- Familiar interface
- Large user base

### TTM Implementation:
- Robust two-temperature model
- Electronic subsystem
- Energy coupling
- Non-equilibrium physics

## Inputs & Outputs
- **Input formats**:
  - LAMMPS input scripts
  - Electronic parameters
  - Coupling coefficients
  
- **Output data types**:
  - Electronic temperature
  - Ionic temperature
  - Energy transfer rates
  - Standard LAMMPS output

## Interfaces & Ecosystem
- **LAMMPS**: Primary MD code
- **fix eph**: Main LAMMPS fix
- Compatible with LAMMPS potentials


## Advanced Features
- **Two-temperature model**: Electronic and ionic subsystems
- **Electronic heat diffusion**: Spatial temperature evolution
- **Langevin thermostat**: Stochastic electron-ion coupling
- **Energy exchange**: Non-equilibrium dynamics
- **LAMMPS fix**: Standard LAMMPS integration
- **Parallel support**: MPI-compatible

## Performance Characteristics
- LAMMPS-based: Efficient MD
- Scales with system size
- MPI parallelization available

## Computational Cost
- Standard LAMMPS MD cost
- Additional overhead for TTM: Moderate
- Electronic grid: Scales with spatial resolution
- Overall: Comparable to standard LAMMPS simulations

## Best Practices
- Validate electron-phonon coupling parameter
- Check electronic temperature evolution
- Use appropriate time steps for energy transfer
- Compare with experimental relaxation times

## Limitations & Known Constraints
- Requires coupling parameter input
- Simplified electronic model
- Classical ion dynamics
- Parameter fitting needed

## Application Areas
- Radiation damage
- Laser ablation
- Ion irradiation
- Ultrafast dynamics
- Materials under extreme conditions

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/LLNL/USER-EPH
2. A. Tamm et al., Phys. Rev. B 94, 024305 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- LLNL development
