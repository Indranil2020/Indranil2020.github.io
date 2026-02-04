# Hirshfeld

## Official Resources
- Homepage: Method - Implemented in various codes (Tonto, Q-Chem, ORCA, VASP via scripts)
- Documentation: https://github.com/theochem/tonto/wiki (Tonto implementation)
- Source Repository: https://github.com/theochem/tonto
- License: LGPL (Tonto) / Varies by implementation

## Overview
"Hirshfeld" refers to Hirshfeld Charge Analysis (and its variants like Hirshfeld-I, Iterative Hirshfeld), a method for partitioning the electron density of a molecule or crystal into atomic contributions. It assigns partial charges to atoms based on the ratio of the free-atom density to the total molecular density ("stockholder" partitioning). The term here likely refers to specific implementations or scripts (e.g., Tonto, or VASP scripts) rather than a single code named "Hirshfeld".

**Scientific domain**: Charge analysis, population analysis, electron density partitioning  
**Target user community**: Computational chemists, crystallographers

## Theoretical Methods
- Hirshfeld Partitioning (Stockholder)
- Iterative Hirshfeld (Hirshfeld-I)
- Extended Hirshfeld (Hirshfeld-E)
- Atomic Dipole Moments
- Electrostatic Potential Fitting

## Capabilities (CRITICAL)
- Calculation of net atomic charges
- Basis-set independent partial charges
- Atomic volumes and moments
- Analysis of electron density topology
- Implemented in: Tonto, Q-Chem, ORCA, ADF, Multiwfn, VASP (via scripts), Gaussian

**Sources**: F. L. Hirshfeld, Theor. Chim. Acta 44, 129 (1977)

## Key Strengths

### Stockholder Partitioning:
- Intuitive physical basis
- Basis-set independent
- Smooth charge distribution
- Widely accepted

### Multiple Variants:
- Standard Hirshfeld
- Hirshfeld-I (iterative)
- Hirshfeld-E (extended)
- Flexibility in choice

### Broad Implementation:
- Tonto (crystallography)
- Multiwfn (molecules)
- Q-Chem, ORCA, ADF
- VASP via scripts

## Inputs & Outputs
- **Input formats**: Electron density (cube, CHGCAR), wavefunction files
- **Output data types**: Partial charges, atomic populations

## Interfaces & Ecosystem
- **Tonto**: Powerful crystallographic toolbox implementing Hirshfeld-I
- **Multiwfn**: Popular analysis tool supporting Hirshfeld
- **VASP**: Scripts available for charge analysis

## Performance Characteristics
- Fast for standard Hirshfeld
- Iterative Hirshfeld (Hirshfeld-I) requires convergence loops

## Limitations & Known Constraints
- **Reference atoms**: Depends on free-atom reference choice
- **Charged systems**: Standard Hirshfeld less reliable
- **Implementation**: No single standard code
- **Convergence**: Hirshfeld-I requires iterations

## Comparison with Other Tools
- **vs Bader**: Hirshfeld smoother, Bader topological
- **vs DDEC**: DDEC more robust for charged systems
- **vs Mulliken**: Hirshfeld basis-set independent
- **Unique strength**: Intuitive stockholder partitioning

## Application Areas
- Force field parameterization
- Analysis of ionic vs covalent character
- Electrostatic potential modeling
- Reactivity indices

## Best Practices
- Use Hirshfeld-I for better convergence
- Compare with other charge methods
- Verify against chemical expectations
- Consider system charge state

## Community and Support
- Method is standard in quantum chemistry
- Tonto and Multiwfn communities

## Verification & Sources
**Primary sources**:
1. Tonto Wiki: https://github.com/theochem/tonto/wiki
2. Publication: F. L. Hirshfeld, Theor. Chim. Acta 44, 129 (1977)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Method: STANDARD
- Implementation: Common in major codes (Tonto, Multiwfn, ORCA)
- Applications: Charge analysis, density partitioning
