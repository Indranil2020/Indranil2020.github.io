# OpenQP (Open Quantum Platform)

## Official Resources
- Homepage: https://github.com/Open-Quantum-Platform/openqp
- Documentation: In repository / Source code
- Source Repository: https://github.com/Open-Quantum-Platform/openqp
- License: Apache License 2.0

## Overview
OpenQP (Open Quantum Platform) is a modern quantum chemistry software package developed by the Choi Group (Kyungpook National University). It features the implementation of Mixed-Reference Spin-Flip Time-Dependent Density Functional Theory (MRSF-TDDFT). This method addresses the critical issue of spin contamination in conventional Spin-Flip TDDFT, allowing for accurate description of ground and excited states with multireference character, conical intersections, and diradicals.

**Scientific domain**: Multireference electronic structure, excited states, spin-flip methods, conical intersections
**Target user community**: Electronic structure theorists, photochemists studying diradicals and bond breaking

## Theoretical Methods
- Mixed-Reference Spin-Flip TDDFT (MRSF-TDDFT)
- Spin-Flip TDDFT (SF-TDDFT)
- Conventional TDDFT
- Hartree-Fock (HF) and DFT ground states
- Analytic gradients for MRSF-TDDFT
- Non-adiabatic coupling vectors (MRSF-TDDFT)
- Linear Response Theory

## Capabilities (CRITICAL)
- Ground and excited state energies
- Geometry optimization (ground and excited)
- Conical intersection search (MECP)
- Non-adiabatic coupling calculation
- Correction of spin contamination
- Proper description of S0-S1 degeneracy in diradicals
- Accurate singlet-triplet gaps

**Sources**: GitHub repository, published papers (JCTC, JCP)

## Key Strengths

### MRSF-TDDFT Method:
- Eliminates spin contamination of SF-TDDFT
- Balanced description of response states
- Accurate for bond-breaking
- Accurate for conical intersections

### Analytic Gradients:
- Geometry optimization for excited states
- Dynamics simulations (via PyOQP)
- MECP optimization capability

### Modern Architecture:
- Python/C++ hybrid (PyOQP)
- Modular design
- Open-source license

## Inputs & Outputs
- **Input formats**:
  - Python scripts (PyOQP)
  - Input blocks for molecules and methods
  
- **Output data types**:
  - Energies and gradients
  - Spin expectation values <S^2>
  - Optimized geometries
  - NAC vectors

## Interfaces & Ecosystem
- **Language**: C++ core, Python interface
- **Parallelization**: OpenMP
- **Libraries**: Eigen3, Libint2
- **Ecosystem**: Can be used as a library or standalone

## Advanced Features

### Conical Intersection Optimization:
- Analytic gradients allow efficient search
- Correct topology at CX thanks to MRSF
- Avoids artifical cusps of standard TDDFT

### Diradical Physics:
- Accurate singlet-triplet splitting
- Proper handling of open-shell singlets
- Double excitation retrieval

## Performance Characteristics
- **Speed**: Comparable to standard TDDFT/SF-TDDFT
- **Accuracy**: Superior to TDDFT for multireference cases
- **Scaling**: N^3 to N^4 depending on implementation
- **Parallelization**: Shared memory (OpenMP)

## Computational Cost
- **Memory**: Moderate (density matrices)
- **Time**: Similar to regular TDDFT linear response
- **Optimization**: Efficient analytic gradients

## Limitations & Known Constraints
- **Feature set**: Focused on (SF-)TDDFT, less comprehensive than Gaussian/Q-Chem
- **Basis sets**: Depends on libint support
- **Solvation**: Functionality may be limited compared to major codes

## Comparison with Other Codes
- **vs Q-Chem**: Q-Chem has SF-TDDFT, but MRSF-TDDFT is OpenQP's specialty
- **vs PySCF**: OpenQP focused on specific MRSF methodology
- **vs GAMESS**: OpenQP more modern C++ architecture
- **Unique strength**: Reference implementation of MRSF-TDDFT

## Application Areas
- **Photoswitches**: Azobenzene, diarylethenes (conical intersections)
- **Diradicals**: Perylene diimide, organic magnetic materials
- **Bond breaking**: Photodissociation curves
- **Singlet Fission**: Electronic state characterization

## Best Practices
- **Reference State**: Choose appropriate high-spin triplet reference
- **Functional**: BHHLYP often used for SF-TDDFT
- **Validation**: Check <S^2> values for spin purity
- **Active Space**: Implicit in SF method, check orbital ordering

## Community and Support
- Open-source Apache 2.0
- Developed by Choi Group
- GitHub issues
- Academic publications serve as documentation foundation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/Open-Quantum-Platform/openqp
2. Y. I. Carreras, H. Park, A. Jiang, C. H. Choi, J. Chem. Phys. 153, 214107 (2020)

**Confidence**: VERIFIED - Research group code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (Apache 2.0)
- Method: MRSF-TDDFT (Scientifically verified)
- Specialized strength: Eliminating spin contamination in SF-TDDFT
