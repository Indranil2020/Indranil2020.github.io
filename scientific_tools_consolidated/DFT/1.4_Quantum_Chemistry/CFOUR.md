# CFOUR

## Official Resources
- Homepage: https://www.cfour.de/
- Documentation: https://www.cfour.de/documentation/
- Source Repository: Available to licensees
- License: Academic and commercial licenses available

## Overview
CFOUR (Coupled-Cluster techniques for Computational Chemistry) is a specialized quantum chemistry program with emphasis on high-level coupled cluster methods and highly accurate molecular property calculations. Originally developed as ACES IV by J. Gauss and J. F. Stanton, CFOUR excels at computing molecular properties with very high precision using advanced post-Hartree-Fock methods, featuring state-of-the-art implementations of analytical derivatives and response theory up to high order.

**Scientific domain**: High-accuracy quantum chemistry, coupled cluster methods, molecular properties, spectroscopy  
**Target user community**: Researchers requiring benchmark-quality calculations, rotational spectroscopists

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Møller-Plesset perturbation theory (MP2, MP3, MP4)
- Coupled Cluster (CCSD, CCSD(T), CCSDT, CCSDT(Q), CC3, CC4)
- Brueckner Coupled Cluster (BCCD, BCCD(T))
- Equation-of-Motion Coupled Cluster (EOM-CCSD, EOM-CCSDT)
- Linear Response CC (LR-CC)
- Symmetry-Adapted Perturbation Theory (SAPT)
- Multi-Reference CC (Mk-MRCC)
- Density Functional Theory (DFT)
- Time-Dependent DFT (TDDFT)
- Analytical gradients (through CCSD(T) and CCSDT)
- Analytical second derivatives (Hessians)
- Analytical third derivatives (cubic force constants)
- Higher-order properties and response

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Very high accuracy coupled cluster calculations
- Geometry optimization with analytic gradients
- Vibrational frequencies with analytic Hessians
- Anharmonic vibrational analysis (VPT2, VCI)
- Spectroscopic constants (rotational B0, centrifugal distortion)
- Quartic and sextic force fields
- Excited states via EOM-CC
- Molecular properties with high-order response:
  - Dipole, quadrupole moments
  - Polarizabilities
  - Hyperpolarizabilities
  - NMR chemical shifts and spin-spin coupling
  - IR and Raman intensities
  - Optical rotation
- Relativistic corrections
- Core-valence correlation
- Benchmark-quality calculations for small molecules

**Sources**: Official CFOUR documentation, cited in 6/7 source lists

## Key Strengths

### Analytical Derivatives:
- Gradients through CCSDT
- Hessians through CCSD(T)
- Cubic force constants
- Quartic force fields
- VPT2 anharmonic corrections

### Spectroscopic Accuracy:
- Rotational constants to spectroscopic accuracy
- Centrifugal distortion constants
- Anharmonic frequencies
- Reference for microwave spectroscopy
- Sub-cm⁻¹ accuracy achievable

### Property Calculations:
- High-order response theory
- NMR chemical shifts (benchmark)
- Spin-spin coupling constants
- Hyperpolarizabilities
- Optical rotation at CC level

### Benchmark Quality:
- Reference calculations
- Thermochemical accuracy
- Method validation
- Basis set studies
- Focal-point analysis

## Inputs & Outputs
- **Input formats**:
  - ZMAT file (Z-matrix input)
  - Keyword-based directives
  - Cartesian coordinates (COORD=CARTESIAN)
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients, Hessians
  - Molecular properties
  - Spectroscopic constants
  - Force fields (quadratic, cubic, quartic)
  - MOLDEN format orbitals

## Interfaces & Ecosystem
- **External programs**:
  - Molpro interface
  - PSI4 interface
  - MRCC interface for high-level CC
  - DIRAC for relativistic integrals
  
- **Utilities**:
  - xjoda - input processor
  - xcfour - main driver
  - xanharm - anharmonic analysis
  - xvcc - VCC calculations
  
- **Workflow integration**:
  - Can be scripted for automated calculations
  - HEAT protocol implementation
  - Focal-point analysis


## Workflow and Usage

### Input Format (ZMAT):
CFOUR historically uses a ZMAT file which can contain Z-matrix (or Cartesian) geometry and keywords.

```
title
Water HF/6-31G*

O
H 1 0.96
H 1 0.96 2 104.5

*CFOUR(CALC=CCSD,BASIS=6-31G*)
```

### Running CFOUR:
```bash
xcfour > output.dat
# Requires ZMAT in current directory
```

### Common Tasks:
- **Optimization**: `GEO_CONV` keyword
- **Frequencies**: `VIB=ANALYTIC`
- **Properties**: `PROPS=ALL`
- **Excited States**: `EXCITATION=EOMEE`

## Advanced Features

### High-Order Derivatives:
- Analytical gradients for CCSD, CCSD(T), CCSDT, CCSD(T)
- Analytical second derivatives for CCSD, CCSD(T)
- Enables extremely precise geometry usage
- Quartic force fields

### Anharmonic Analysis (VPT2):
- Vibrational Perturbation Theory (2nd order)
- Accurate anharmonic frequencies
- Fermi resonance treatment
- Critical for high-res spectroscopy

### HEAT Thermochemistry:
- High-accuracy Extrapolated Ab initio Thermochemistry
- Composite protocol implementation
- Sub-kJ/mol accuracy
- Zero-point energy corrections
- Relativistic and diagonal Born-Oppenheimer corrections

### EOM-CC Properties:
- Excited state gradients
- Transition moments
- Spin-orbit coupling
- Ionization potentials and Electron affinities

### Born-Oppenheimer Breakdown:
- Diagonal Born-Oppenheimer Correction (DBOC)
- Analytical evaluation
- Essential for high-precision corrections

## Performance Characteristics
- **Accuracy**: Unrivaled for properties of small molecules
- **Speed**: Optimized for coupled cluster, but scaling is steep
- **Memory**: High-order methods (CCSDT, etc.) require significant memory/disk
- **Parallelization**: MPI parallelization available, but scalability is limited compared to NWChem/CC4S
- **Disk I/O**: Heavy usage for intermediate files

## Computational Cost
- **CCSD**: O(N^6), standard
- **CCSD(T)**: O(N^7), expensive
- **CCSDT**: O(N^8), very expensive (~10 atoms max)
- **CCSDT(Q)**: O(N^9), extremely expensive (3-5 atoms)
- **VPT2**: Requires Hessian + cubic force field, very expensive

## Comparison with Other Codes
- **vs Gaussian**: CFOUR has far superior analytical derivatives for CC methods; Gaussian more automated for general users.
- **vs Molpro**: Molpro faster for standard energies; CFOUR precise for properties/derivatives.
- **vs MRCC**: MRCC handles higher orders (arbitrary order); CFOUR optimized for up to CCSDT(Q) with analytical gradients.
- **Unique strength**: Analytical high-order derivatives, anharmonic spectroscopy, reference-level accuracy.

## Best Practices

### Geometry Optimization:
- Use `COORD=CARTESIAN` if preferred over Z-matrix
- Convergence criteria (`SCF_CONV`, `GEO_CONV`) should be tight (10^-7 to 10^-10)
- Restart from `OLDMOS` to save time

### Basis Sets:
- Use correlation consistent sets (cc-pVXZ)
- Use core-valence sets (cc-pCVXZ) for high accuracy properties
- Check linear dependency in large basis sets

### Memory:
- Set `MEMORY_SIZE` carefully
- Ensure fast scratch disk for I/O

## Community and Support
- **Mailing List**: Active user support list
- **Website**: Comprehensive documentation and examples
- **Development**: University of Florida / Mainz / UT Austin
- **License**: Academic/Commercial
- **Workshops**: Occasional specialized workshops

## Application Areas

### Rotational Spectroscopy:
- Microwave spectroscopy support
- Rotational constants
- Centrifugal distortion
- Structure determination
- Spectroscopic databases

### Thermochemistry:
- HEAT protocol (High-accuracy Extrapolated Ab initio Thermochemistry)
- Atomization energies
- Reaction enthalpies
- Benchmark accuracy

### Vibrational Spectroscopy:
- Anharmonic frequencies
- IR/Raman intensities
- VCD spectra
- Fermi resonances

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Molecular focus**: Not designed for periodic systems
- **System size**: High-level CC limited to small molecules (~10-20 atoms)
- **Basis sets**: Gaussian-type; large correlation-consistent basis required
- **Input format**: Z-matrix input can be challenging
- **Memory**: High-level methods very memory-intensive
- **Learning curve**: Steep; requires deep understanding of methods
- **Parallelization**: Limited compared to modern codes
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://www.cfour.de/
2. Documentation: https://www.cfour.de/documentation/
3. J. F. Stanton et al., CFOUR program package
4. D. A. Matthews et al., J. Chem. Phys. 152, 214108 (2020) - CC implementations
5. J. Gauss, J. F. Stanton - analytical derivative papers

**Secondary sources**:
1. CFOUR manual and examples
2. Published benchmark studies using CFOUR
3. High-accuracy spectroscopy applications
4. HEAT protocol publications
5. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Software: Free for academics (registration required)
- Community support: Active (mailing list)
- Academic citations: >2,000
- Gold standard: Reference for high-accuracy properties and spectroscopy
- Specialized strength: Analytical high-order derivatives, spectroscopic constants, anharmonic analysis, NMR properties, benchmark thermochemistry
