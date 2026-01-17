# PyMEX

## Official Resources
- **Repository**: https://github.com/imaitygit/PyMEX
- **License**: MIT (Check repository)
- **Primary Citation**: Indrajit Maity et al., "Atomistic theory of twist-angle dependent intralayer and interlayer exciton properties in twisted bilayer materials"

## Overview
PyMEX (Python Moiré Exciton) is a specialized Python package designed to calculate exciton properties in moiré superlattices, such as twisted bilayer transition metal dichalcogenides (TMDs). It solves the Bethe-Salpeter Equation (BSE) using a Wannier function basis, enabling the efficient study of intralayer and interlayer excitons in large-unit-cell systems.

**Scientific domain**: 2D materials, Moiré superlattices, Excitonics
**Target user community**: Researchers in twistronics and 2D material optics

## Theoretical Methods
- **Bethe-Salpeter Equation (BSE)**: Solves the eigenproblem for neutral excitations
- **Wannier Basis**: Uses maximally localized Wannier functions to construct the Hamiltonian
- **Tight-Binding**: Input via tight-binding models
- **Zero-Momentum Excitons**: Current implementation focus (extensible to finite momentum)

## Capabilities
- **BSE Solver**: Computes exciton eigenvalues and eigenvectors
- **Optical Conductivity**: Calculates optical selection rules and conductivity spectra
- **Large Systems**: optimized for the large unit cells characteristic of Moiré patterns
- **Hybrid Implementation**: Performance-critical loops optimized with Cython

## Inputs & Outputs
- **Input formats**:
  - Tight-binding Hamiltonian (Wannier90 format or internal)
  - Coulomb interaction definitions
  - Parameter files for twist angles and lattice constants
- **Output data types**:
  - Exciton energies
  - Excitonic wavefunctions (real-space visualization)
  - Optical absorption/conductivity spectra

## Performance Characteristics
- **Parallelization**: MPI and OpenMP support via Cython and libraries
- **Efficiency**: Hybrid Python/Cython approach balances usability and speed
- **Memory**: Optimized with h5py-parallel for large data handling

## Comparison with Other Codes
- **vs Yambo/BerkeleyGW**: PyMEX is specialized for Moiré systems and Wannier basis, whereas Yambo/BGW are general-purpose plane-wave codes.
- **vs NanoGW**: Both treat confined/specific systems, but PyMEX focuses on 2D twistronics.

## Usage & Best Practices
- **Prerequisites**: Python 3.x, NumPy, SciPy, MPI4Py.
- **Workflow**: Generate Wannier Hamiltonian -> Define Moiré geometry -> Run PyMEX BSE solver -> Analyze spectra.

## Limitations & Known Constraints
- **Zero-Momentum**: Primary release focus on Q=0 excitons.
- **System**: Specialized for Moiré/Twisted systems.
