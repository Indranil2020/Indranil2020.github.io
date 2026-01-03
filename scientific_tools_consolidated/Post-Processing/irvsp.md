# irvsp (Irreducible Representations for VASP)

## Official Resources
- Homepage: https://github.com/zjwang11/irvsp
- Documentation: https://github.com/zjwang11/irvsp/blob/master/README.md
- Source Repository: https://github.com/zjwang11/irvsp
- License: MIT License (implied/open)

## Overview
irvsp is a program to compute the irreducible representations (irrep) of electronic states in VASP calculations. It determines the symmetry of Bloch wavefunctions at high-symmetry points in the Brillouin zone, which is crucial for identifying topological phases of matter, enforcing selection rules, and analyzing band connectivity.

**Scientific domain**: Topological materials, symmetry analysis, electronic structure  
**Target user community**: Topological physics researchers, VASP users

## Theoretical Methods
- Group theory
- Irreducible representations of space groups
- Trace of symmetry operators
- Projection of wavefunctions
- Topological invariants (symmetry-based indicators)

## Capabilities (CRITICAL)
- Calculation of irreps for all 230 space groups (with and without SOC)
- Identification of topological invariants (Z2, Z4, Chern numbers) based on symmetry indicators
- Analysis of band connectivity
- Support for spin-orbit coupling (SOC)
- Interfaces with VASP WAVECAR

**Sources**: irvsp GitHub, Comp. Phys. Comm. 261, 107760 (2021)

## Inputs & Outputs
- **Input formats**: WAVECAR, POSCAR, OUTCAR (from VASP)
- **Output data types**: Irreps labels (e.g., GM1+, X2-), trace data, topological indices

## Interfaces & Ecosystem
- **VASP**: Specialized for VASP output
- **Phonopy**: Can use Phonopy for symmetry analysis
- **ir2tb**: Workflow integration for tight-binding models

## Workflow and Usage
1. Perform self-consistent VASP calculation (with symmetry on).
2. Prepare `tband.dat` (optional) or rely on `OUTCAR`.
3. Run `irvsp`: `irvsp -sg <SpaceGroupNumber>`
4. Analyze output for band symmetries and topological classification.

## Performance Characteristics
- Fast analysis of wavefunction symmetry
- Requires compatibility with VASP compiled version (macros)

## Application Areas
- Topological insulators (TI)
- Topological crystalline insulators (TCI)
- Weyl and Dirac semimetals
- High-throughput topological screening

## Community and Support
- Open-source
- Developed by Zhijun Wang group (IOP CAS)
- GitHub issues for support

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/zjwang11/irvsp
2. Publication: J. Gao, Q. Wu, C. Persson, Z. Wang, Comp. Phys. Comm. 261, 107760 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE (README)
- Source: OPEN (GitHub)
- Development: ACTIVE (Wang Group)
- Applications: Irreducible representations, topological materials, VASP
