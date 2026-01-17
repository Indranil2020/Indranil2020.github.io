# CE-TDDFT

## Official Resources
- Homepage: https://github.com/dceresoli/ce-tddft
- Documentation: https://github.com/dceresoli/ce-tddft/blob/master/README.md
- Source Repository: https://github.com/dceresoli/ce-tddft
- License: GNU General Public License (same as Quantum ESPRESSO)

## Overview
CE-TDDFT is a Real-Time Time-Dependent Density Functional Theory (RT-TDDFT) extension for Quantum ESPRESSO. It enables explicit time-propagation of Kohn-Sham orbitals to simulate non-equilibrium electron dynamics, strong-field phenomena, and Ehrenfest molecular dynamics.

**Scientific domain**: Ultrafast electron dynamics, strong-field physics, nonlinear optics, photochemistry  
**Target user community**: Researchers using Quantum ESPRESSO who need real-time propagation capabilities beyond linear-response TDDFT

## Theoretical Methods
- Real-Time Time-Dependent DFT (RT-TDDFT)
- Explicit time propagation of KS orbitals
- Ehrenfest molecular dynamics (coupled electron-ion dynamics)
- Electric field perturbations (delta-kick, continuous fields)
- Time-dependent charge density analysis
- Electron wavepacket dynamics (experimental)
- Plane-wave basis set (from QE)
- Norm-conserving pseudopotentials

## Capabilities
- **Optical absorption spectra** via Fourier transform of induced dipole
- **Strong-field dynamics** (beyond linear response)
- **Ehrenfest dynamics** for coupled electron-ion motion
- **Time-dependent observables** (dipole moment, charge density)
- **Electron wavepacket propagation** (experimental feature)
- **Time-dependent charge density** saving for post-processing
- **Compatible with QE-6.3+** (tested up to recent versions)

## Key Strengths

### Tight QE Integration:
- Direct extension of Quantum ESPRESSO
- Uses QE's plane-wave infrastructure
- Familiar input format for QE users
- Leverages QE pseudopotentials and XC functionals

### Ehrenfest Dynamics:
- Coupled electron-ion propagation
- Non-adiabatic dynamics
- Photochemical reaction pathways
- Energy transfer mechanisms

### Active Development:
- Regular updates (2017-2019 active development)
- Six official releases
- Responsive maintainer (D. Ceresoli)

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO input files
  - Pseudopotential files (UPF format)
  - Additional TDDFT control parameters
  
- **Output data types**:
  - Time-dependent dipole moments
  - Absorption spectra (via post-processing)
  - Time-resolved charge densities
  - Ehrenfest trajectories

## Interfaces & Ecosystem
- **Quantum ESPRESSO integration**:
  - Built on top of QE's PWscf module
  - Uses QE's parallelization (MPI)
  - Compatible with QE's post-processing tools

## Performance Characteristics
- **Speed**: Dependent on QE's plane-wave efficiency
- **Scalability**: Inherits QE's MPI parallelization
- **System size**: Limited by plane-wave scaling (typically <1000 atoms)
- **Time step**: Typically 0.01-0.1 atomic units

## Limitations & Known Constraints
- **QE version dependency**: Must match compatible QE version
- **Norm-conserving only**: No ultrasoft/PAW (RT propagation constraint)
- **System size**: Plane-wave scaling limits large systems
- **Memory**: Dense time-dependent data storage
- **Periodic systems**: Primarily designed for periodic/slab geometries

## Comparison with Other Codes
- **vs Octopus**: CE-TDDFT plane-wave, Octopus real-space; Octopus standalone, CE-TDDFT QE extension
- **vs SALMON**: Both plane-wave RT-TDDFT; SALMON standalone with more laser options
- **vs turboTDDFT**: turboTDDFT linear-response, CE-TDDFT real-time propagation
- **Unique strength**: Direct QE integration for existing QE users, Ehrenfest dynamics

## Application Areas
- Ultrafast spectroscopy simulations
- Strong-field ionization
- High-harmonic generation
- Photochemical dynamics
- Charge transfer processes
- Plasmon dynamics in nanostructures

## Best Practices
- Start with converged ground state from PWscf
- Use small time steps (0.01-0.05 a.u.) for stability
- Apply delta-kick perturbation for linear response spectra
- Monitor total energy conservation
- Use sufficient simulation time for frequency resolution

## Community and Support
- Open-source (GPL, same as QE)
- GitHub repository with issue tracker
- 6 official releases
- Author: Davide Ceresoli (CNR-ISTM, Milan)
- Integration with QE community

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/dceresoli/ce-tddft
2. Quantum ESPRESSO: https://www.quantum-espresso.org/

**Confidence**: VERIFIED
- Repository: ACCESSIBLE (GitHub)
- Documentation: Available in README
- Active releases: 6 releases (2017-2019)
- Integration: Works with QE 6.1-6.3+

**Verification status**: ✅ VERIFIED
