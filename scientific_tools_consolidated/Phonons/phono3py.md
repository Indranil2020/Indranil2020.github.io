# phono3py

## Official Resources
- Homepage: https://phonopy.github.io/phono3py/
- Documentation: https://phonopy.github.io/phono3py/
- Source Repository: https://github.com/phonopy/phono3py
- License: BSD 3-Clause License

## Overview
phono3py is a code for computing lattice thermal conductivity and related properties from first principles using three-phonon interactions. It extends Phonopy to include anharmonic effects via third-order force constants, enabling calculations of phonon lifetimes, thermal conductivity, and other anharmonic properties essential for thermoelectric materials and thermal management applications.

**Scientific domain**: Anharmonic lattice dynamics, thermal conductivity, phonon-phonon interactions  
**Target user community**: Researchers studying thermal transport, thermoelectric materials, and anharmonic phonon properties

## Theoretical Methods
- Third-order force constants (3rd IFCs)
- Phonon-phonon interaction strength
- Phonon Boltzmann transport equation (BTE)
- Relaxation time approximation (RTA)
- Iterative solution of BTE
- Phonon linewidth and lifetime
- Lattice thermal conductivity tensor
- Mode-dependent thermal conductivity
- Cumulative thermal conductivity
- Group velocities from phonon dispersion

## Capabilities (CRITICAL)
- Lattice thermal conductivity calculation
- Phonon lifetimes and linewidths
- Phonon-phonon scattering rates
- Temperature-dependent thermal conductivity
- Directional thermal conductivity (tensor)
- Mode-resolved contributions
- Cumulative thermal conductivity
- Spectral thermal conductivity
- Phonon-phonon interaction strengths
- Phonon linewidths and lifetimes
- Lattice thermal conductivity tensor (RTA and full BTE solution)
- Cumulative thermal conductivity analysis
- Mode-resolved thermal conductivity contributions
- Thermal conductivity vs temperature
- Isotope scattering effects
- Boundary scattering models
- Group velocity and mean free path
- Phase space volume for three-phonon processes
- Mode Gruneisen parameters from 3rd order IFCs
- Interface to multiple DFT codes (same as phonopy)
- HDF5 output for large datasets
- Parallelization support

**Sources**: Official phono3py documentation, GitHub, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - POSCAR (VASP structure format)
  - FORCES_FC3 (forces on displaced atoms for 3rd order)
  - FORCES_FC2 (forces for 2nd order force constants)
  - fc3.hdf5 (precomputed third-order force constants)
  - fc2.hdf5 (precomputed second-order force constants)
  - DFT code outputs via interfaces
  
- **Output data types**:
  - kappa-*.hdf5 (thermal conductivity data)
  - fc3.hdf5 (third-order force constants)
  - gamma-*.hdf5 (phonon linewidths)
  - Thermal conductivity vs temperature
  - Mode-resolved thermal conductivity
  - Cumulative thermal conductivity plots

## Interfaces & Ecosystem
- **DFT code interfaces** (inherited from phonopy):
  - VASP - most common
  - Quantum ESPRESSO
  - ABINIT
  - CRYSTAL
  - TURBOMOLE
  - CASTEP
  - CP2K
  - DFTB+
  - Elk
  - SIESTA
  - FHI-aims
  
- **Framework integrations**:
  - Phonopy - requires phonopy for 2nd order IFCs
  - ASE - structure manipulation
  - pymatgen - structure I/O
  - ShengBTE - alternative BTE solver (can compare results)
  
- **Workflow integration**:
  - Can be integrated into high-throughput workflows
  - Python API for custom analysis

## Limitations & Known Constraints
- **Computational cost**: Extremely expensive; requires forces for numerous displaced supercells (scales as N³)
- **Supercell size**: Large supercells needed for convergence; 2×2×2 or 3×3×3 minimum for many systems
- **Cutoff distance**: Third-order cutoff must be carefully converged; long-range interactions may be important
- **Three-phonon processes only**: Neglects four-phonon and higher-order processes (important at high T)
- **Harmonic phonon requirement**: Requires stable harmonic phonons; unstable modes cause failures
- **Classical treatment**: Uses classical Bose-Einstein statistics; quantum corrections not included
- **Isotope scattering**: Simplified model; detailed isotope configurations not considered
- **Boundary scattering**: Phenomenological models; not fully first-principles
- **Memory**: HDF5 files can become very large for fine q-point meshes
- **Convergence**: Requires extensive convergence testing (supercell size, cutoffs, q-mesh)

## Verification & Sources
**Primary sources**:
1. Official documentation: https://phonopy.github.io/phono3py/
2. GitHub repository: https://github.com/phonopy/phono3py
3. A. Togo et al., Phys. Rev. B 91, 094306 (2015) - phono3py methodology
4. L. Chaput, Phys. Rev. Lett. 110, 265506 (2013) - Direct BTE solution method

**Secondary sources**:
1. phono3py examples and tutorials
2. Comparison studies with ShengBTE and ALAMODE
3. High-throughput thermal conductivity databases using phono3py
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues, shared with phonopy)
- Academic citations: >400 (Google Scholar)
- DFT interfaces: Verified for 10+ codes
