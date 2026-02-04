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

## Key Strengths

### VASP Specialized:
- Optimized for VASP output
- Direct WAVECAR parsing
- Efficient implementation
- Production-ready

### Complete Coverage:
- All 230 space groups
- SOC and non-SOC
- Symmorphic and non-symmorphic
- Magnetic groups supported

### High-Throughput Ready:
- Automated workflow
- Batch processing
- CheckTopologicalMat integration
- Database screening

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

## Limitations & Known Constraints
- **VASP-only**: Specialized for VASP output
- **Compilation**: Requires proper VASP version compatibility
- **Space group input**: User must specify correct space group
- **Memory**: Large WAVECAR files require adequate RAM

## Comparison with Other Tools
- **vs IrRep**: irvsp VASP-specific, IrRep multi-code
- **vs qeirreps**: irvsp for VASP, qeirreps for QE
- **vs SpaceGroupIrep**: irvsp DFT-based, SpaceGroupIrep database
- **Unique strength**: Optimized VASP integration, high-throughput ready

## Application Areas
- Topological insulators (TI)
- Topological crystalline insulators (TCI)
- Weyl and Dirac semimetals
- High-throughput topological screening

## Best Practices
- Run Phonopy first to verify space group
- Use pos2aBR to standardize POSCAR
- Set MAGMOM explicitly for SOC calculations
- Verify against known topological materials

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
