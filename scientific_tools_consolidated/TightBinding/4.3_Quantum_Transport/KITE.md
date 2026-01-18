# KITE

## Official Resources
- **Homepage**: https://quantum-kite.com/
- **Repository**: https://github.com/QUANTUM-KITE/kite
- **License**: GPL-3.0

## Overview
**KITE** is a high-performance open-source software suite for accurate quantum transport and electronic structure calculations in **large-scale real-space tight-binding models**. It employs the **Kernel Polynomial Method (KPM)**, using Chebyshev polynomial expansions of the Green's function and density operator to avoid the cubic sealing of exact diagonalization. This enables KITE to simulate systems with **billions of atomic orbitals**, making it ideal for disordered materials, twisted bilayers (twistronics), and topological systems.

**Scientific domain**: Condensed Matter Physics, Twistronics, Disordered Systems
**Target user community**: Researchers studying spectral and transport properties of large systems

## Theoretical Methods
- **Kernel Polynomial Method (KPM)**: Expands spectral functions (DOS, Conductivity) in Chebyshev polynomials.
- **Real-Space Tight-Binding**: Handles arbitrary lattice geometries and multi-orbital Hamiltonians.
- **Disorder on the Fly**: Adds disorder (vacancies, Anderson) to the Hamiltonian application step without storing the full disordered matrix.
- **Linear Response**: Calculates Kubo-Bastin or Kubo-Greenwood formulas via polynomial moments.

## Capabilities
- **Observables**:
  - Density of States (DOS).
  - DC Conductivity ($\sigma_{xx}, \sigma_{xy}$).
  - Optical Conductivity $\sigma(\omega)$.
  - Quantum Hall and Anomalous Hall conductivity.
  - Time evolution of wavepackets.
- **Systems**:
  - Graphene and 2D materials (TMDs, Phosphorene).
  - Twisted Bilayer Graphene (Moiré lattices).
  - 3D Topological Insulators.
  - Disordered alloys/quasicrystals.

## Key Strengths
- **Scalability**: $O(N)$ scaling allows simulations of experimental-size samples ($10^9$ atoms).
- **Efficiency**: "Disorder on the fly" minimizes memory usage, enabling massive systems on single nodes or clusters.
- **Versatility**: Handles both spectral properties (DOS) and transport (Conductivity) in the same framework.

## Inputs & Outputs
- **Inputs**: Python scripts defining the lattice, tight-binding model (`PyBinding` style), and KPM parameters (number of moments).
- **Outputs**:
  - HDF5 files containing Chebyshev moments.
  - Post-processed spectral data (CSV/TXT).

## Interfaces & Ecosystem
- **PyBinding**: Compatible with PyBinding for easy model construction.
- **Python**: Core logic wrapped in Python for ease of use.

## Performance Characteristics
- **Speed**: Orders of magnitude faster than diagonalization for large $N$.
- **Parallelism**: Comparison-free parallelization strategies; MPI/OpenMP enabled.
- **Compute**: CPU optimized; GPU version often discussed/available.

## Comparison with Other Codes
- **vs. Kwant**: Kwant uses sparse matrix solvers (MUMPS) which scale worse ($> O(N)$) for 2D/3D bulk; KITE's KPM is stricter $O(N)$ and better for very large bulk systems (though Kwant is better for open boundaries/leads).
- **vs. TBPLaS**: Similar KPM capabilities; KITE is known for its "disorder on the fly" efficiency.

## Application Areas
- **Twistronics**: Simulating moiré patterns with thousands of atoms per unit cell.
- **Anderson Localization**: Studying localization transitions in very large disordered lattices.
- **Topological Invariants**: Computing Chern numbers via real-space markers (`KITE` has modules for this).

## Community and Support
- **Development**: University of York (Aires Ferreira) and University of Minho.
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://quantum-kite.com/](https://quantum-kite.com/)
- **Primary Publication**: S. M. Joao et al., R. Soc. Open Sci. 7, 191809 (2020).
- **Verification status**: ✅ VERIFIED
  - Active and modern code (released ~2020).
