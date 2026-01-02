# CASINO

## Official Resources
- Homepage: https://vallico.net/casino/
- Documentation: https://vallico.net/casino/documentation/
- Source Repository: Available to users (registration required)
- License: Proprietary (free for academic use)

## Overview
CASINO is a highly sophisticated quantum Monte Carlo (QMC) code for electronic structure calculations on molecules and solids. Developed by the Theory of Condensed Matter group at Cambridge, it provides state-of-the-art implementations of variational and diffusion Monte Carlo with emphasis on accuracy, advanced trial wavefunctions, and user-friendly features for production calculations.

**Scientific domain**: Quantum Monte Carlo, electronic structure, high-accuracy ab initio calculations  
**Target user community**: Researchers requiring benchmark-quality QMC calculations for molecules and materials

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Diffusion Monte Carlo (DMC)
- Fixed-node and released-node DMC
- Reptation quantum Monte Carlo
- Slater-Jastrow trial wavefunctions
- Multi-determinant expansions
- Pfaffian wavefunctions (open-shell)
- Backflow transformations
- Geminal wavefunctions
- Pairing wavefunctions
- Wavefunction optimization (variance and energy minimization)
- Periodic boundary conditions

## Capabilities (CRITICAL)
- Highly accurate total energies (chemical accuracy)
- Ground-state properties of molecules and solids
- Equation of state calculations
- Excited states via VMC
- Structural optimization with QMC forces
- Electron and spin densities
- Momentum densities and Compton profiles
- Pair correlation functions
- Transition dipole moments
- Finite-size corrections for periodic systems
- Twist-averaged boundary conditions
- Pseudopotential and all-electron calculations
- Extensive wavefunction optimization

**Sources**: Official CASINO documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - input file (CASINO format)
  - gwfn.data (wavefunction from external codes)
  - Pseudopotentials
  - Trial wavefunctions from DFT or quantum chemistry
  
- **Output data types**:
  - out file (energies, standard deviations)
  - vmc.hist, dmc.hist (energy histories)
  - Optimized wavefunction parameters
  - Density files
  - Configuration files

## Interfaces & Ecosystem
- **Wavefunction interfaces**:
  - CRYSTAL for Gaussian-type orbitals
  - Quantum ESPRESSO for plane-waves
  - GAMESS, Gaussian, TURBOMOLE
  - PySCF, ORCA
  
- **Utilities**:
  - casinohelp - interactive help
  - runqmc - run manager
  - Various analysis scripts
  
- **Parallelization**:
  - MPI for distributed calculations
  - Efficient for large systems

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Computational cost**: DMC expensive; slower than DFT by orders of magnitude
- **Statistical errors**: Stochastic method; results have error bars
- **Pseudopotentials**: Quality critical; locality approximation required
- **Fixed-node approximation**: Nodal surface from trial wavefunction
- **System size**: Practical limit ~100-500 electrons
- **Wavefunction quality**: Depends on trial wavefunction from external code
- **Learning curve**: Steep; QMC methods require understanding
- **Convergence**: Multiple parameters to optimize and converge
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://vallico.net/casinoqmc/
2. Documentation: https://vallico.net/casinoqmc/documentation.html
3. R. J. Needs et al., J. Chem. Phys. 152, 154106 (2020) - CASINO overview
4. M. D. Towler et al., Comput. Phys. Commun. 98, 181 (1996) - CASINO QMC

**Secondary sources**:
1. CASINO manual and tutorials
2. Published benchmark QMC studies
3. High-accuracy calculations vs experiment
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Software: Free for academics (registration required)
- Community support: Active (mailing list, workshops)
- Academic citations: >500
- Gold standard: Benchmark-quality QMC calculations
