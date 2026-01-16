# GPAW

## Official Resources
- Homepage: https://wiki.fysik.dtu.dk/gpaw/
- Documentation: https://wiki.fysik.dtu.dk/gpaw/documentation/documentation.html
- Source Repository: https://gitlab.com/gpaw/gpaw
- License: GNU General Public License v3.0 or later

## Overview
GPAW is a density-functional theory Python code based on the projector-augmented wave (PAW) method. It combines the efficiency of real-space grids with the accuracy of plane-wave methods and integrates seamlessly with the Atomic Simulation Environment (ASE). GPAW is particularly strong for large-scale calculations, time-dependent DFT, and as a platform for method development.

**Scientific domain**: Real-space DFT, PAW method, Python-based calculations  
**Target user community**: Researchers needing flexible, programmable DFT calculations, method developers

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Projector augmented wave (PAW) method
- Real-space grid representation
- Plane-wave mode available
- Finite-difference mode
- LCAO (linear combination of atomic orbitals) mode
- Plane-wave mode
- Linear combination of atomic orbitals (LCAO)
- Time-Dependent DFT (TDDFT)
- GW approximation
- Bethe-Salpeter Equation (BSE)
- Random Phase Approximation (RPA)
- Exact exchange and hybrid functionals
- Non-collinear magnetism and spin-orbit coupling

## Capabilities (CRITICAL)
- Ground-state electronic structure calculations
- Geometry optimization and structure relaxation
- Molecular dynamics (NVE, NVT, NPT)
- Band structure and density of states
- Total energy calculations
- Forces and stress tensors
- Optical absorption spectra via TDDFT
- Quasiparticle energies via GW
- Optical spectra via BSE
- Linear response calculations (phonons, dielectric)
- Magnetic properties (moments, anisotropies)
- STM image simulation
- Electric field calculations
- Bader charge analysis
- Python scripting for workflows
- Real-time TDDFT propagation
- Delta-SCF for excited states
- Implicit solvation models

**Sources**: Official GPAW documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Python scripts (native interface)
  - ASE Atoms objects
  - Standard structure files via ASE (CIF, XYZ, POSCAR, etc.)
  - GPAW restart files (.gpw)
  
- **Output data types**:
  - GPAW binary files (.gpw) with wavefunctions
  - Text output with energies, forces
  - Cube files for densities and orbitals
  - Eigenvalues and DOS
  - Trajectories in ASE format

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - native integration (GPAW built on ASE)
  - pymatgen - structure conversion
  - AiiDA - workflow automation possible
  - Fireworks - workflow integration
  
- **Built-in tools**:
  - TDDFT module for optical spectra
  - GW module for quasiparticles
  - BSE module for excitonic effects
  - Linear response for phonons
  - Response function calculations
  
- **Post-processing**:
  - gpaw-tools for common analysis
  - Python-based custom analysis
  - Integration with matplotlib, numpy, scipy

## Limitations & Known Constraints
- **Python overhead**: Slower than pure Fortran/C++ codes for routine calculations
- **Real-space grids**: Require careful convergence testing for grid spacing
- **Memory**: Can be memory-intensive for large systems with fine grids
- **Pseudopotentials**: Limited to PAW datasets included with GPAW
- **Parallelization**: MPI + OpenMP supported but efficiency varies
- **LCAO mode**: Basis set quality depends on available atomic orbital sets
- **GW/BSE**: Computationally expensive; limited to smaller systems
- **Documentation**: Good but assumes Python/ASE familiarity
- **Installation**: Dependencies can be complex (libxc, BLAS/LAPACK, FFTW)

## Computational Cost
- **LCAO Mode**: Very fast ($O(N)$), comparable to SIESTA.
- **FD/PW Mode**: Slower than VASP/QE for small systems due to Python overhead, but scales well for large grids.
- **Memory**: Real-space grids can consume significant RAM if grid spacing ($h$) is very small (< 0.15 Å).

## Comparison with Other Codes
- **vs VASP/QE**: GPAW offers greater flexibility (3 modes: PW, FD, LCAO) and full Python scripting, but pure compute speed in PW mode is often lower than optimized Fortran codes.
- **vs ASE**: GPAW is the "reference" calculator for ASE and has the tightest integration.

## Best Practices
- **Mode Selection**: Use LCAO for quick screening, then Finite Difference (FD) or Plane Wave (PW) for final high-accuracy numbers.
- **Grid Spacing**: Standard is $h=0.18-0.20$ Å. Refine to $0.15$ Å for hard potentials.
- **Parallelization**: Parallelize over k-points first, then domains.

## Community and Support
- **Mailing List**: `gpaw-users` list is very active.
- **Development**: Hosted on GitLab; easy to contribute Python code.

## Verification & Sources
**Primary sources**:
1. Official website: https://wiki.fysik.dtu.dk/gpaw/
2. Documentation: https://wiki.fysik.dtu.dk/gpaw/documentation/
3. GitLab repository: https://gitlab.com/gpaw/gpaw
4. J. Phys.: Condens. Matter 22, 253202 (2010) - GPAW code paper
5. J. Chem. Phys. 152, 124101 (2020) - GPAW LCAO developments

**Secondary sources**:
1. ASE documentation (GPAW calculator)
2. GPAW tutorials and workshops
3. Published applications using GPAW
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab)
- Community support: Active (mailing list, GitLab issues)
- Academic citations: >1,000 (main papers)
- Active development: Regular releases
