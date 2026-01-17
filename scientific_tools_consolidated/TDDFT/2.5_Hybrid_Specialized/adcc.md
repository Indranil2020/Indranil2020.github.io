# adcc (ADC-connect)

## Official Resources
- Homepage: https://adc-connect.org/
- Documentation: https://adcc.readthedocs.io/
- Source Repository: https://github.com/adc-connect/adcc
- License: BSD 3-Clause License

## Overview
adcc (Algebraic Diagrammatic Construction-connect) is a Python-based hybrid library designed to perform excited-state calculations using the Algebraic Diagrammatic Construction (ADC) scheme for the polarization propagator. It serves as a flexible backend that can interface with various SCF codes (like PySCF and Psi4) to obtain the Hartree-Fock reference, enabling correlated excited-state calculations up to the third order of perturbation theory.

**Scientific domain**: Correlated electronic structure, excited states, propagator methods, ADC
**Target user community**: Developers and researchers needing modular ADC implementations or reliable excited states beyond TDDFT

## Theoretical Methods
- ADC(2): Second-order ADC
- ADC(2)-x: Extended second-order ADC
- ADC(3): Third-order ADC
- CVS-ADC: Core-Valence Separation for core states
- Intermediate State Representation (ISR)
- Dyson equation solvers
- Davidson diagonalization

## Capabilities (CRITICAL)
- Excitation energies (Singlet and Triplet)
- Transition moments and oscillator strengths
- One-particle density matrices
- Excited state dipole moments
- Core-excited states (XAS/XES) via CVS
- Interactive usage via Python/Jupyter
- Verification against Q-Chem/MOLCAS results

**Sources**: Official website, GitHub, J. Chem. Phys. 2020

## Key Strengths

### Modular Design:
- Host-code agnostic (needs matrix/tensor interface)
- Easy Python integration
- Tensor-based implementation (using opt_einsum)

### Methodological Range:
- From fast ADC(2) to accurate ADC(3)
- Rigorous handling of double excitations
- Size-consistent excited states

### Core Spectroscopy:
- CVS implementation essential for X-ray absorption
- Target specific edges
- Accurate relaxation effects

## Inputs & Outputs
- **Input formats**:
  - Python objects (SCF results from PySCF/Psi4)
  - Python scripts
  
- **Output data types**:
  - Excitation tables (Pandas dataframes)
  - Spectral data
  - Analyzing state character (doubles %, etc.)
  - Densities

## Interfaces & Ecosystem
- **Host interactions**: PySCF (native support), Psi4 (native support), VeloxChem
- **Language**: Python (frontend), C++ (backend kernels)
- **Dependencies**: NumPy, H5Py, opt_einsum, pybind11

## Advanced Features

### Interactive Analysis:
- Analyze state composition interactively
- Plot spectra directly in notebooks
- Export densities for visualization

### High-Order Methods:
- ADC(3) implementation available
- Useful for benchmarking TDDFT
- Captures double excitation character

## Performance Characteristics
- **Speed**: Optimized tensor contractions
- **Accuracy**: ADC(2) similar to CC2, ADC(3) higher
- **System size**: Moderate (tens to small hundreds of atoms)
- **Parallelization**: Shared memory (OpenMP/NumPy threads)

## Computational Cost
- **ADC(2)**: N^5 scaling (feasible for molecules)
- **ADC(3)**: N^6 scaling (smaller systems)
- **Memory**: Tensor storage can be significant
- **Intermediate**: Requires HF ground state

## Limitations & Known Constraints
- **Reference**: Restricted to HF (Canonical orbitals usually)
- **Open-shell**: UHF support exists but limited
- **Solvation**: No implicit solvent models in library itself
- **Gradients**: Analytical gradients not fully exposed/implemented for all methods

## Comparison with Other Codes
- **vs Q-Chem**: adcc is free/open-source, Q-Chem commercial
- **vs TURBOMOLE**: adcc offers Python flexibility, ADC(3)
- **vs EOM-CCSD**: ADC(2) cheaper than EOM-CCSD, often sufficient
- **Unique strength**: Open-source, hackable, Pythonic ADC(3) implementation

## Application Areas
- **Benchmarking**: Validating TDDFT results
- **Charge Transfer**: Handling CT states correctly
- **Core States**: X-ray spectroscopy simulation
- **Double Excitations**: States with significant non-single character

## Best Practices
- **Basis Set**: Requires decent basis (aug-cc-pVDZ/TZ)
- **Memory**: Ensure sufficient RAM for tensor intermediates
- **CVS**: Use large core-valence separation energy gap
- **Comparison**: Compare ADC(2) vs ADC(3) for convergence

## Community and Support
- Open-source BSD
- Developed by Heidelberg/Tubingen groups
- Documentation on ReadTheDocs
- Active GitHub

## Verification & Sources
**Primary sources**:
1. Website: https://adc-connect.org/
2. GitHub: https://github.com/adc-connect/adcc
3. M. F. Herbst et al., J. Chem. Phys. 152, 244119 (2020)

**Confidence**: VERIFIED - Active academic project

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Source code: OPEN (BSD)
- Method: ADC(n) (Scientifically verified)
- Specialized strength: Python library for high-level excited states
