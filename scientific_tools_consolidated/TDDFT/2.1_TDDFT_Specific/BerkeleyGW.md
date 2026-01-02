# BerkeleyGW

## Official Resources
- Homepage: https://berkeleygw.org/
- Documentation: https://berkeleygw.org/documentation/
- Source Repository: Available to users (registration required)
- License: BSD-like license (free for academic use)

## Overview
BerkeleyGW is a massively parallel code for computing the quasiparticle and optical properties of materials using many-body perturbation theory within the GW approximation and the Bethe-Salpeter equation (BSE). It is designed for large-scale calculations on leadership-class supercomputers and provides highly accurate band structures, band gaps, and optical spectra beyond DFT.

**Scientific domain**: Many-body perturbation theory, GW approximation, optical properties, excited states  
**Target user community**: Researchers studying electronic excitations, band structures, and optical properties of materials

## Theoretical Methods
- GW approximation (G₀W₀, eigenvalue self-consistent GW)
- Generalized Plasmon Pole (GPP) model
- Contour-deformation method (full frequency)
- Bethe-Salpeter equation (BSE) for optical spectra
- Time-dependent density functional theory (TDDFT)
- Electron-hole interaction (excitonic effects)
- Static and dynamic screening
- Coulomb-hole screened-exchange (COHSEX)
- Scissors operator corrections
- Contour-deformation method
- Bethe-Salpeter Equation (BSE)
- Static and dynamical screening
- Vertex corrections (selected cases)
- Hybrid functional starting points

## Capabilities (CRITICAL)
- Quasiparticle band structures via GW
- Quasiparticle band gaps and corrections to DFT
- Optical absorption spectra including excitonic effects (BSE)
- Electron-hole interaction analysis
- Exciton wavefunctions and binding energies
- Static and frequency-dependent dielectric functions
- Interface to multiple DFT codes (mean-field starting point)
- Massively parallel execution (MPI + OpenMP)
- GPU acceleration (experimental)
- Real-space grids and subspace methods
- Static subspace approximation for reduced cost
- Spin-orbit coupling support
- Non-collinear magnetism support
- Interpolation schemes for band structures

**Sources**: Official BerkeleyGW website, documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - WFN files (wavefunctions from DFT codes)
  - RHO files (charge density)
  - VXC files (exchange-correlation potential)
  - Input files (epsilon.inp, sigma.inp, kernel.inp, absorption.inp)
  - k-point and q-point grids
  
- **Output data types**:
  - eps0mat, epsmat (dielectric matrices)
  - eqp.dat (quasiparticle energies)
  - sigma.log (self-energy calculations)
  - absorption_eh.dat (optical absorption)
  - eigenvectors (exciton wavefunctions)

## Interfaces & Ecosystem
- **DFT code interfaces** (verified):
  - Quantum ESPRESSO - primary interface via pw2bgw.x
  - PARATEC - native support
  - PARSEC - supported
  - SIESTA - interface available
  - Octopus - interface available
  - ABINIT - can be interfaced
  
- **Pre-processing tools**:
  - kgrid.x - k-point grid generation
  - wfn2hdf.x - wavefunction format conversion
  - mf_convert_wrapper.sh - mean-field conversion scripts
  
- **Post-processing tools**:
  - plotxct.x - exciton wavefunction plotting
  - absorption.x - optical spectra
  - inteqp.x - band structure interpolation
  - offdiag.x - off-diagonal matrix elements
  
- **Analysis utilities**:
  - Python scripts for data analysis
  - Plotting utilities
  - Convergence checking tools

## Limitations & Known Constraints
- **Computational cost**: Extremely expensive; GW scales as O(N⁴), BSE even worse
- **Memory requirements**: Very large; dielectric matrices can be tens to hundreds of GB
- **Convergence**: Many parameters to converge (cutoffs, bands, k-points, q-points)
- **DFT dependency**: Quality depends on mean-field starting point
- **Registration required**: Free but requires registration for download
- **Learning curve**: Steep; requires understanding of GW theory and convergence procedures
- **System size**: Practical limit ~100-200 atoms for GW; smaller for BSE
- **Static approximation**: Plasmon-pole models introduce approximations
- **Parallelization**: Requires careful setup for optimal performance
- **GPU support**: Experimental; not all features GPU-accelerated

## Verification & Sources
**Primary sources**:
1. Official website: https://berkeleygw.org/
2. Documentation: https://berkeleygw.org/documentation/
3. J. Deslippe et al., Comput. Phys. Commun. 183, 1269 (2012) - BerkeleyGW code paper
4. M. S. Hybertsen and S. G. Louie, Phys. Rev. B 34, 5390 (1986) - GW method
5. M. Rohlfing and S. G. Louie, Phys. Rev. B 62, 4927 (2000) - BSE method

**Secondary sources**:
1. BerkeleyGW tutorials and workshops
2. Quantum ESPRESSO interface documentation
3. Published benchmarks and applications
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (requires registration)
- Community support: Active (mailing list, workshops)
- Academic citations: >600 (main code paper)
- HPC optimization: Extensively benchmarked on supercomputers
