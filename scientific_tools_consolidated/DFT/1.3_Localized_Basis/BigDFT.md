# BigDFT

## Official Resources
- Homepage: https://bigdft.org/
- Documentation: https://bigdft.org/Wiki/index.php?title=BigDFT_website
- Source Repository: https://gitlab.com/l_sim/bigdft-suite
- License: GNU General Public License v2.0

## Overview
BigDFT is a DFT code using Daubechies wavelets as a basis set, providing systematic convergence and efficient treatment of systems in vacuum, surfaces, and periodic systems. It features linear-scaling capabilities and excellent support for massively parallel calculations with particular strength in isolated and low-dimensional systems.

**Scientific domain**: Molecules, nanostructures, surfaces, linear-scaling calculations  
**Target user community**: Researchers studying isolated systems, surfaces, and needing systematic basis convergence

## Theoretical Methods
- Density Functional Theory (DFT)
- Daubechies wavelet basis sets
- Systematic basis set convergence
- Norm-conserving and PAW pseudopotentials
- LDA, GGA functionals
- Hybrid functionals (experimental)
- van der Waals corrections
- DFT+U for correlated systems
- Linear-scaling DFT (O(N) method)
- Poisson solver for arbitrary boundary conditions

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Systematic convergence with single parameter
- Isolated molecules (no spurious interactions)
- Surfaces and low-dimensional systems
- Periodic systems (1D, 2D, 3D)
- Linear-scaling DFT for large systems
- Geometry optimization and transition states
- Molecular dynamics (NVE, NVT, NPT)
- Band structure and DOS
- Forces and stress tensors
- Massively parallel calculations
- Adaptive grid refinement
- Fragment-based calculations
- Excited states via TDDFT (in development)
- Poisson solver for mixed boundary conditions

**Sources**: Official BigDFT documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Input files (YAML format)
  - XYZ coordinate files
  - Pseudopotential files
  - Python API for input generation
  
- **Output data types**:
  - YAML output files
  - Energies and forces
  - Optimized structures
  - Wavefunction data
  - Density files
  - DOS outputs

## Interfaces & Ecosystem
- **Python integration**:
  - PyBigDFT - Python interface
  - Jupyter notebook support
  - High-level Python API
  
- **Framework integrations**:
  - Can interface with ASE
  - Compatible with workflow tools
  
- **Visualization**:
  - v_sim - visualization tool
  - Compatible with standard tools
  
- **Linear-scaling**:
  - Fragment approach
  - Excellent parallel scaling

## Limitations & Known Constraints
- **Wavelet basis**: Less familiar than plane-waves or orbitals
- **Pseudopotentials**: Requires specific formats
- **Hybrid functionals**: Limited implementation
- **Community**: Smaller than major codes
- **Documentation**: Good but evolving
- **Learning curve**: Wavelet methods require understanding
- **k-point sampling**: Best for systems where Γ-point sufficient
- **Installation**: Requires compilation and libraries
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://bigdft.org/
2. Documentation: https://bigdft.org/Wiki/
3. GitLab repository: https://gitlab.com/l_sim/bigdft-suite
4. L. Genovese et al., J. Chem. Phys. 129, 014109 (2008) - Daubechies wavelets
5. S. Mohr et al., J. Chem. Phys. 140, 204110 (2014) - BigDFT linear-scaling

**Secondary sources**:
1. BigDFT tutorials and workshops
2. PyBigDFT documentation
3. Published applications
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitLab)
- Community support: Active (developers, mailing list)
- Academic citations: >300 (main papers)
- Active development: Regular releases
