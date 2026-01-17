# pypwdft

## Official Resources
- Source Repository: https://github.com/ifilot/pypwdft
- License: Open Source (GPL-3.0)

## Overview
pypwdft is a pure Python implementation of a Plane-Wave Density Functional Theory (PW-DFT) code. It provides a complete, self-contained educational platform for understanding how plane-wave basis sets, pseudopotentials, and reciprocal space operations work in standard codes like VASP or Quantum ESPRESSO, but in a readable, non-HPC Python environment.

**Scientific domain**: Educational theory, Plane-Wave DFT
**Target user community**: Students, developers learning PW algorithms

## Theoretical Methods
- Kohn-Sham DFT
- Plane-Wave basis set
- Periodic Boundary Conditions (PBC)
- Pseudopotentials (Local and Non-local parts)
- Ewald Summation for electrostatics
- Conjugate Gradient minimization for electrons
- LDA/GGA Functionals

## Capabilities
- Periodic crystal calculations (e.g., Silicon, Aluminum)
- Electronic band structure generation
- Total energy minimization
- Charge density calculation
- Reciprocal space handling (k-points)

## Key Strengths

### Pedagogical Value:
- Demonstrates FFT grid operations explicitely
- Shows construction of the Hamiltonian in G-space
- Clear implementation of Ewald summation
- Handles pseudopotential integration grids

### Code Readability:
- Python classes map clearly to physical concepts (Cell, Atom, Wavefunction)
- Avoids complex parallelization logic

## Inputs & Outputs
- **Inputs**:
  - Lattice parameters
  - Atom positions
  - Plane-wave cutoff energy
  - Pseudopotential files
- **Outputs**:
  - Total energy
  - Forces (if implemented)
  - Stress (if implemented)
  - Band structure data

## Interfaces & Ecosystem
- **Python**: Core implementation.
- **NumPy/SciPy**: Heavy usage for FFTs and linear algebra.

## Advanced Features
- **CG Optimization**: Uses conjugate gradient for wavefunction optimization (common in large codes).
- **Pseudopotential parsing**: Simple readers for format handling.

## Performance Characteristics
- **Slow**: Pure Python + NumPy FFTs are significantly slower than C/Fortran.
- **Memory**: Stores dense grids, limited to small unit cells.

## Computational Cost
- **Moderate**: Can handle small unit cells (2-8 atoms) on a desktop, but scales poorly compared to compiled codes.

## Limitations & Known Constraints
- **Performance**: Not a production code.
- **Features**: Lacks advanced features like hybrid functionals, spin-orbit coupling, or complex constraints.
- **Pseudopotentials**: Limited library support compared to Quantum ESPRESSO.

## Comparison with Other Codes
- **vs SimpleDFT/eminus**: Similar scope, but pypwdft often focuses more explicitly on the "traditional" VASP-like algorithms (CG minimization).
- **vs Quantum ESPRESSO**: pypwdft is the "textbook model" of what QE does at high scale.
- **Unique strength**: Accessibility for understanding reciprocal space logic without Fortran.

## Application Areas
- **Teaching**: Solid State Physics courses.
- **Algorithm Development**: Testing new minimization schemes or functional forms in Python.

## Best Practices
- **Small Cutoffs**: Keep kinetic energy cutoffs low for testing.
- **Simple Systems**: Start with Silicon or Aluminum bulk.
- **Profile**: Use Python profilers to seeing where the time goes (usually FFTs).

## Community and Support
- **GitHub**: Repository by active maintainer (Ivo Filot).

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/ifilot/pypwdft

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GPL-3.0)
- Functionality: Functional PW-DFT prototype.
