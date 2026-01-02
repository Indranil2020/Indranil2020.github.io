# CONQUEST

## Official Resources
- Homepage: https://www.order-n.org/
- Documentation: https://www.order-n.org/tutorials.html
- Source Repository: https://github.com/OrderN/CONQUEST-release
- License: MIT License (open-source)

## Overview
CONQUEST is a linear-scaling (O(N)) DFT code designed for massively parallel calculations on extremely large systems (millions of atoms). It uses local orbital basis sets and achieves excellent parallel scalability, making it ideal for biomolecules, nanostructures, and large-scale materials simulations.

**Scientific domain**: Very large biomolecules, nanostructures, materials (10,000+ atoms)  
**Target user community**: Researchers needing DFT for extremely large systems on supercomputers

## Theoretical Methods
- Density Functional Theory (DFT)
- Linear-scaling DFT (O(N) method)
- Blip functions (B-spline functions on a grid)
- Pseudo-atomic orbital (PAO) basis sets
- Norm-conserving pseudopotentials
- LDA, GGA functionals
- van der Waals corrections
- DFT+U for correlated systems
- Multiple time-step integration for MD

## Capabilities (CRITICAL)
- Ground-state electronic structure for very large systems
- True linear-scaling DFT (O(N) cost)
- Systems with millions of atoms demonstrated
- Massively parallel (1000s of processors)
- Geometry optimization for large systems
- Molecular dynamics (NVE, NVT)
- Large biomolecules (proteins, DNA)
- Nanostructures and materials
- Band structure and DOS
- Forces and stress tensors
- Mixed basis (blips + PAOs)
- Multisite support functions for accuracy
- Excellent weak and strong scaling
- Adaptive grid refinement

**Sources**: Official CONQUEST documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Coord.dat (coordinates)
  - Input files (CONQUEST format)
  - Pseudopotential files
  - PAO basis definitions
  
- **Output data types**:
  - Standard output with energies, forces
  - Optimized structures
  - MD trajectories
  - Density files
  - DOS outputs

## Interfaces & Ecosystem
- **Parallelization**:
  - MPI for distributed memory
  - Excellent scaling to 1000s of processors
  - Tested on major supercomputers
  
- **Basis sets**:
  - Blip functions (smooth B-splines)
  - PAO basis sets
  - Mixed blip/PAO approach
  
- **Post-processing**:
  - Analysis utilities included
  - Compatible with standard visualization tools

## Limitations & Known Constraints
- **Specialized for large systems**: Not optimal for small systems (<1000 atoms)
- **Linear-scaling overhead**: O(N) methods have overhead; crossover point at ~1000 atoms
- **Basis sets**: Blip/PAO basis requires understanding
- **Pseudopotentials**: Limited to norm-conserving
- **Documentation**: Good but smaller community than major codes
- **Learning curve**: Linear-scaling methods require expertise
- **Installation**: Requires MPI and libraries
- **k-point sampling**: Limited; best for large supercells (Gamma-point)
- **Properties**: Fewer property calculations than general-purpose codes
- **Platform**: Primarily HPC systems

## Verification & Sources
**Primary sources**:
1. Official website: https://www.order-n.org/
2. Documentation: https://www.order-n.org/tutorials.html
3. GitHub repository: https://github.com/OrderN/CONQUEST-release
4. D. R. Bowler and T. Miyazaki, Rep. Prog. Phys. 75, 036503 (2012) - Linear-scaling review
5. A. Nakata et al., J. Chem. Phys. 152, 164112 (2020) - CONQUEST large-scale

**Secondary sources**:
1. CONQUEST tutorials and examples
2. Published applications on very large systems
3. HPC scaling benchmarks
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (developers, mailing list)
- Academic citations: >300
- Active development: Regular releases, HPC optimization
