# Pomerol

## Official Resources
- Homepage: https://github.com/pomerol-ed/pomerol
- Documentation: https://github.com/pomerol-ed/pomerol/wiki
- Source Repository: https://github.com/pomerol-ed/pomerol
- License: GNU General Public License v2.0

## Overview
Pomerol is an exact diagonalization (full ED) code written in C++ for solving condensed matter second-quantized models of interacting fermions and bosons on finite size lattices at finite temperatures. It is designed to produce single and two-particle Green's functions and can be used as an impurity solver in DMFT calculations. Pomerol uses Lehmann representation for efficient computation.

**Scientific domain**: Exact diagonalization, quantum impurity solvers, DMFT  
**Target user community**: Researchers needing exact solutions for small quantum systems

## Theoretical Methods
- Exact diagonalization (full ED)
- Lehmann representation
- Finite size lattices
- Finite temperature formalism
- Fermions and bosons
- Green's function calculations
- Two-particle correlation functions

## Capabilities (CRITICAL)
- Exact diagonalization for quantum models
- Single-particle Green's functions
- Two-particle Green's functions
- Finite temperature calculations
- Fermionic and bosonic systems
- Impurity solver for DMFT
- Real-frequency spectral functions (no analytical continuation)
- Lehmann representation for efficiency
- C++ implementation
- Python interface (pomerol2triqs)
- Integration with TRIQS

**Sources**: Official Pomerol repository (https://github.com/pomerol-ed/pomerol), confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- Hamiltonian definition (C++ or Python)
- System parameters
- Temperature
- Frequency grids

**Output data types**:
- Green's functions (real and Matsubara frequencies)
- Two-particle correlation functions
- Spectral functions
- Observables
- Energy eigenvalues

## Interfaces & Ecosystem
- **TRIQS**: pomerol2triqs provides TRIQS interface
- **DCore**: Can use Pomerol as impurity solver
- **C++**: Native C++ implementation
- **Python**: Python bindings available
- **libcommute**: Uses libcommute for operator algebra

## Limitations & Known Constraints
- Limited to small systems (exponential scaling)
- Hilbert space must be manageable
- Memory intensive
- Not suitable for large multi-orbital systems
- Bath discretization affects accuracy
- Exact solver: computationally expensive
- Best for small clusters or impurity problems

## Comparison with Other Codes
| Feature | Pomerol | HPhi |
| :--- | :--- | :--- |
| **Method** | Exact Diagonalization (Full ED) | Exact Diagonalization (Lanczos, TPQ) |
| **Primary Application** | DMFT Impurity Solver, Green's Functions | Lattice Models, Ground State/Thermal Properties |
| **Performance** | Optimized for GF calculation, automatic symmetries | Massively parallel, TPQ for finite T |
| **Key Strength** | Versatility for generic fermionic interactions (libcommute) | Scalability to larger lattice sizes |

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/pomerol-ed/pomerol
2. Documentation: https://github.com/pomerol-ed/pomerol/wiki
3. pomerol2triqs: https://github.com/pomerol-ed/pomerol2triqs

**Secondary sources**:
1. TRIQS benchmarks (includes Pomerol)
2. DCore documentation (lists Pomerol as solver)
3. Exact diagonalization literature
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official repository: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (Wiki)
- Source code: OPEN (GitHub, GPL v2)
- TRIQS integration: CONFIRMED (pomerol2triqs)
- Active development: Regular updates
- C++ implementation: Efficient
- Exact solver: No approximations
