# ONETEP

## Official Resources
- Homepage: https://onetep.org/
- Documentation: https://onetep.org/Main/Documentation
- Source Repository: Available to licensed users
- License: Academic license (free for academics, registration required)

## Overview
ONETEP (Order-N Electronic Total Energy Package) is a linear-scaling DFT code that achieves plane-wave accuracy with localized orbitals. It uniquely combines the accuracy of plane-wave calculations with the efficiency of linear-scaling methods, enabling calculations on systems with thousands of atoms.

**Scientific domain**: Large biomolecules, nanostructures, materials with thousands of atoms  
**Target user community**: Researchers needing DFT accuracy for very large systems (1000-10000+ atoms)

## Theoretical Methods
- Density Functional Theory (DFT)
- Linear-scaling DFT (O(N) method)
- Non-orthogonal generalized Wannier functions (NGWFs)
- Plane-wave accuracy with localized basis
- Norm-conserving pseudopotentials
- LDA, GGA functionals
- Hybrid functionals (range-separated)
- DFT+U for correlated systems
- van der Waals corrections (DFT-D, TS)
- Implicit solvation models
- TDDFT for excited states

## Capabilities (CRITICAL)
- Ground-state electronic structure for very large systems
- Linear-scaling DFT (computational cost scales linearly with system size)
- Plane-wave accuracy with localized orbitals
- Geometry optimization for large systems
- Molecular dynamics (NVE, NVT, NPT)
- Systems with 1000-10000+ atoms
- Protein and biomolecule calculations
- Nanostructures and materials
- Band structure and DOS
- Forces and stress tensors
- NMR chemical shifts
- EPR parameters
- Core-level spectroscopy
- Implicit solvation (DDCOSMO)
- TDDFT for absorption spectra
- Wannier function analysis
- Conduction calculations
- Ensemble DFT

**Sources**: Official ONETEP documentation, cited in 7/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Input file (ONETEP format)
  - PDB, XYZ coordinate files
  - Pseudopotential files
  - NGWF initial guesses
  
- **Output data types**:
  - Standard output with energies, forces
  - Optimized structures
  - Density files
  - DOS and PDOS
  - NGWF outputs
  - Property-specific outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - Can interface with molecular dynamics packages
  - QM/MM capabilities
  
- **Visualization**:
  - Compatible with standard visualization tools
  - Cube file output for densities
  
- **Post-processing**:
  - Built-in analysis tools
  - NGWF manipulation utilities

## Limitations & Known Constraints
- **Academic license**: Free for academics but requires registration
- **Not fully open-source**: Source available to licensed users only
- **Learning curve**: Linear-scaling methods and NGWFs require understanding
- **NGWF optimization**: Can be challenging to converge for some systems
- **Pseudopotentials**: Limited to norm-conserving
- **Hybrid functionals**: Computationally expensive even with linear-scaling
- **Parallelization**: Excellent but requires understanding of distribution
- **Memory**: Lower than conventional DFT but still significant for very large systems
- **Installation**: Requires compilation and libraries
- **Platform**: Primarily Linux/Unix, HPC systems

## Verification & Sources
**Primary sources**:
1. Official website: https://onetep.org/
2. Documentation: https://onetep.org/Main/Documentation
3. C.-K. Skylaris et al., J. Chem. Phys. 122, 084119 (2005) - ONETEP method
4. N. D. M. Hine et al., Comput. Phys. Commun. 180, 1041 (2009) - Linear-scaling

**Secondary sources**:
1. ONETEP tutorials and workshops
2. Published large-scale biomolecule applications
3. Linear-scaling benchmark studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (requires registration for full access)
- Source code: Available to licensed users
- Community support: Active (user mailing list, workshops)
- Academic citations: >500 (main papers)
- Active development: Regular releases, well-maintained
