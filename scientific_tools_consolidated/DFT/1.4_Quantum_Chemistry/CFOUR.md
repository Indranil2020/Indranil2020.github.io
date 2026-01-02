# CFOUR

## Official Resources
- Homepage: https://www.cfour.de/
- Documentation: https://www.cfour.de/documentation/
- Source Repository: Available to licensees
- License: Academic and commercial licenses available

## Overview
CFOUR (Coupled-Cluster techniques for Computational Chemistry) is a specialized quantum chemistry program with emphasis on high-level coupled cluster methods and highly accurate calculations. It excels at computing molecular properties with very high precision using advanced post-Hartree-Fock methods, featuring state-of-the-art implementations of analytical derivatives and response theory.

**Scientific domain**: High-accuracy quantum chemistry, coupled cluster methods, molecular properties  
**Target user community**: Researchers requiring benchmark-quality calculations and precise molecular properties

## Theoretical Methods
- Hartree-Fock (HF)
- Møller-Plesset perturbation theory (MP2-MP4)
- Coupled Cluster (CCSD, CCSD(T), CCSDT, CC3)
- Equation-of-Motion Coupled Cluster (EOM-CCSD)
- Linear Response CC (LR-CC)
- Symmetry-Adapted Perturbation Theory (SAPT)
- Multi-Reference CC (MRCC)
- Density Functional Theory (DFT)
- Time-Dependent DFT (TDDFT)
- Analytical gradients (up to CCSD(T))
- Analytical second derivatives (Hessians)
- Analytical third derivatives (cubic force constants)
- Higher-order properties
- Relativistic corrections

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Very high accuracy coupled cluster calculations
- Geometry optimization with analytic gradients
- Vibrational frequencies with analytic Hessians
- Anharmonic vibrational analysis (VPT2)
- Spectroscopic constants (rotational, vibrational)
- Excited states via EOM-CC
- Molecular properties (dipole, quadrupole, polarizability)
- NMR chemical shifts and spin-spin coupling
- IR and Raman intensities
- Optical rotation
- Hyperpolarizabilities
- Response properties
- Thermochemistry
- Benchmark-quality calculations for small molecules

**Sources**: Official CFOUR documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - ZMAT file (Z-matrix input)
  - Keyword-based directives
  - Cartesian coordinates
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients, Hessians
  - Molecular properties
  - Spectroscopic data
  - Wavefunction analysis

## Interfaces & Ecosystem
- **External programs**:
  - Molpro interface
  - PSI4 interface
  - MRCC interface for high-level CC
  
- **Utilities**:
  - xjoda - input processor
  - xcfour - main driver
  - Analysis tools for properties
  
- **Workflow integration**:
  - Can be scripted for automated calculations
  - Compatible with benchmark databases

## Limitations & Known Constraints
- **Registration required**: Free for academics but requires registration
- **Molecular focus**: Not designed for periodic systems
- **System size**: High-level CC limited to small molecules (~10-20 atoms)
- **Basis sets**: Gaussian-type; extensive basis set requirements
- **Input format**: Z-matrix can be challenging for complex molecules
- **Memory**: High-level methods very memory-intensive
- **Learning curve**: Steep; requires deep understanding of methods
- **Parallelization**: Limited compared to modern codes
- **Platform**: Primarily Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: https://www.cfour.de/
2. Documentation: http://www.cfour.de/documentation/
3. J. F. Stanton et al., CFOUR program package
4. D. A. Matthews et al., J. Chem. Phys. 152, 214108 (2020) - Coupled cluster implementations

**Secondary sources**:
1. CFOUR manual and examples
2. Published benchmark studies using CFOUR
3. High-accuracy spectroscopy applications
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Software: Free for academics (registration required)
- Community support: Active (mailing list)
- Academic citations: >1,000
- Gold standard: Reference for high-accuracy calculations
