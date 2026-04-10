# UppASD

## Official Resources
- Source Repository: https://github.com/UppASD/UppASD
- GitLab: https://gitlab.com/UppASD/UppASD
- Documentation: https://uppasd.github.io/UppASD/
- License: GNU General Public License v3

## Overview
**UppASD** (Uppsala Atomistic Spin Dynamics) is a simulation tool for atomistic spin dynamics and Monte Carlo simulations of Heisenberg spin systems. It studies magnetization dynamics using the atomistic Landau-Lifshitz-Gilbert (LLG) equation and can compute magnon dispersion via linear spin-wave theory.

**Scientific domain**: Atomistic spin dynamics, Monte Carlo magnetism, magnon spectroscopy  
**Target user community**: Researchers studying magnetic materials at the atomistic level, magnon dispersion, and phase transitions

## Theoretical Methods
- Atomistic Landau-Lifshitz-Gilbert (LLG) equation
- Monte Carlo simulations (Metropolis, heat bath)
- Heisenberg Hamiltonian (isotropic, anisotropic, DM)
- Linear spin-wave theory for magnon dispersion
- Spin-spin correlation functions
- Dzyaloshinskii-Moriya interaction
- Uniaxial and biaxial anisotropy
- External magnetic fields

## Capabilities (CRITICAL)
- Atomistic spin dynamics simulation
- Monte Carlo simulation (equilibrium properties)
- Magnon dispersion calculation
- Spin-spin correlation functions
- Critical temperature (Tc) determination
- Magnetic ground state identification
- Temperature-dependent magnetization
- Hysteresis loops
- Non-collinear magnetism
- Spin-lattice coupling (external)

**Sources**: GitHub repository, O. Eriksson et al., "Atomistic Spin Dynamics" (Oxford University Press, 2017)

## Key Strengths

### Comprehensive Spin Dynamics:
- LLG equation integration
- Multiple integration schemes
- Stochastic thermal fluctuations
- Time-dependent magnetization

### Monte Carlo:
- Equilibrium properties
- Phase transitions
- Critical temperatures
- Ground state determination

### Magnon Spectroscopy:
- Linear spin-wave theory
- Magnon dispersion
- Dynamical structure factor
- Comparison with INS/RIXS

### DFT Integration:
- Exchange parameters from DFT
- DM vectors from DFT
- Anisotropy from DFT
- Combined DFT+ASD workflow

## Inputs & Outputs
- **Input formats**:
  - UppASD input files (jfile, kfile, etc.)
  - Exchange coupling parameters
  - DM interaction vectors
  - Anisotropy constants
  
- **Output data types**:
  - Magnetization vs temperature
  - Magnon dispersion
  - Spin-spin correlations
  - Hysteresis loops
  - Energy vs time

## Interfaces & Ecosystem
- **DFT codes**: Exchange parameter input
- **AiiDA**: aiida-uppasd plugin available
- **Python**: Post-processing scripts
- **VASP/QE**: Parameter extraction workflows

## Performance Characteristics
- **Speed**: Fast (LLG integration)
- **Accuracy**: Depends on exchange parameters
- **System size**: Millions of spins
- **Parallelization**: MPI + OpenMP

## Computational Cost
- **Equilibrium MC**: Hours
- **Spin dynamics**: Hours
- **Magnon dispersion**: Minutes
- **Typical**: Moderate

## Limitations & Known Constraints
- **Classical spins**: No quantum effects
- **Heisenberg model**: Limited Hamiltonian forms
- **Exchange parameters**: Need external DFT calculation
- **No orbital moments**: Spin-only dynamics
- **3D only**: No 2D-specific optimizations

## Comparison with Other Codes
- **vs Spirit**: UppASD is Fortran, Spirit is C++
- **vs VAMPIRE**: UppASD has magnon dispersion, VAMPIRE is more general
- **vs SpinW**: UppASD is dynamics, SpinW is spin-wave theory
- **Unique strength**: Atomistic spin dynamics with magnon dispersion, Monte Carlo, and DFT integration

## Application Areas

### Magnetic Materials:
- Transition metal ferromagnets
- Rare-earth magnets
- Antiferromagnets
- Ferrimagnets

### Magnon Spectroscopy:
- Magnon dispersion comparison
- Inelastic neutron scattering
- RIXS magnon spectra
- Spin waves in nanostructures

### Phase Transitions:
- Critical temperature calculation
- Order-disorder transitions
- Spin reorientation transitions
- Multiferroic transitions

### Spintronics:
- Domain wall dynamics
- Skyrmion dynamics
- Spin torque effects
- Ultrafast demagnetization

## Best Practices

### Exchange Parameters:
- Use well-converged DFT calculations
- Include sufficient neighbor shells
- Validate against experimental Tc
- Consider DM interaction for non-centrosymmetric

### Monte Carlo:
- Use sufficient thermalization steps
- Average over multiple runs
- Test finite-size effects
- Use heat bath for faster convergence

### Spin Dynamics:
- Choose appropriate time step
- Include thermal fluctuations
- Monitor energy conservation
- Use sufficient averaging

## Community and Support
- Open source (GPL v3)
- Developed at Uppsala University
- Published textbook: "Atomistic Spin Dynamics" (Oxford, 2017)
- Active development
- AiiDA plugin available

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/UppASD/UppASD
2. O. Eriksson et al., "Atomistic Spin Dynamics" (Oxford University Press, 2017)
3. A. Bergman et al., Phys. Rev. B 81, 144416 (2010)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub/GitLab)
- Documentation: ACCESSIBLE
- Published methodology: Oxford University Press
- Active development: Ongoing
- Specialized strength: Atomistic spin dynamics with magnon dispersion, Monte Carlo, DFT integration
