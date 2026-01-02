# SIESTA

## Official Resources
- Homepage: https://siesta-project.org/siesta/
- Documentation: https://docs.siesta-project.org/
- Source Repository: https://gitlab.com/siesta-project/siesta
- License: GNU General Public License v3.0

## Overview
SIESTA (Spanish Initiative for Electronic Simulations with Thousands of Atoms) is an efficient DFT code using strictly localized numerical atomic orbital basis sets. It excels at large-scale calculations with linear-scaling capabilities and is particularly strong for low-dimensional systems and molecules.

**Scientific domain**: Large-scale materials, nanostructures, surfaces, molecules  
**Target user community**: Researchers needing efficient DFT for large systems (1000+ atoms)

## Theoretical Methods
- Density Functional Theory (DFT)
- Numerical atomic orbital (NAO) basis sets
- Strictly localized basis functions
- Norm-conserving pseudopotentials
- LDA, GGA functionals
- van der Waals corrections (DFT-D, VDW-DF)
- DFT+U for correlated systems
- Spin-orbit coupling
- Non-collinear magnetism
- Time-Dependent DFT (via TranSIESTA)
- Linear-scaling O(N) DFT

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and MD (NVE, NVT, NPT)
- Large-scale systems (1000+ atoms)
- Linear-scaling DFT for very large systems
- Band structure and DOS
- Forces and stress tensors
- Phonon calculations via finite differences
- Molecular dynamics
- Quantum transport (TranSIESTA)
- Non-equilibrium Green's function (NEGF) for transport
- STM image simulation
- Optical properties
- Electric polarization
- Wannier functions
- Constrained DFT
- Thermostats and barostats
- Variable cell dynamics

**Sources**: Official SIESTA documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - fdf input files (Flexible Data Format)
  - XV, XYZ coordinate files
  - Pseudopotential files (.psf, .vps)
  
- **Output data types**:
  - Standard output with energies, forces
  - XV files (structures)
  - Density matrices
  - DOS and band structure files
  - LDOS, PDOS files
  - Molecular dynamics trajectories

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface
  - Phonopy - phonon calculations
  - pymatgen - structure I/O
  - Deneb - GUI and workflow
  
- **Transport calculations**:
  - TranSIESTA - quantum transport
  - Smeagol - transport properties
  - inelastica - inelastic transport
  
- **Utilities**:
  - Util/ directory with analysis tools
  - sisl - Python interface to SIESTA
  
- **Post-processing**:
  - Denchar - charge density plotting
  - grid2cube - grid file conversion

## Limitations & Known Constraints
- **Basis sets**: NAO basis sets require careful convergence testing
- **Basis completeness**: Strictly localized basis less complete than plane-waves
- **Pseudopotentials**: Limited to norm-conserving; quality varies
- **Accuracy**: Generally less accurate than plane-wave codes for same computational cost
- **Overlap matrix**: Can become ill-conditioned for small basis cutoffs
- **Documentation**: Comprehensive but can be overwhelming

## Verification & Sources
**Primary sources**:
1. Official website: https://siesta-project.org/siesta/
2. Documentation: https://docs.siesta-project.org/
3. GitLab repository: https://gitlab.com/siesta-project/siesta
4. J. M. Soler et al., J. Phys. Condens. Matter 14, 2745 (2002) - SIESTA method
5. E. Artacho et al., Phys. Status Solidi B 215, 809 (1999) - Linear-scaling
6. A. García et al., J. Chem. Phys. 152, 204108 (2020) - Recent developments

**Secondary sources**:
1. SIESTA manual and tutorials
2. Published DFT studies using SIESTA (>8,000 citations)
3. Workshop materials
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in ALL 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab, GPL v3)
- Community support: Active mailing list, workshops
- Academic citations: >10,000
- Active development: Regular releases
- Benchmark validation: Extensively validated
- Specialized strength: O(N) methods, large systems, localized basis, quantum transport
- Source code: OPEN (GitLab)
- Community support: Very active (mailing list, GitLab)
- Academic citations: >10,000 (main paper)
- Active development: Regular releases, large community
