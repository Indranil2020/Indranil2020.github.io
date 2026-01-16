# CONQUEST

## Official Resources
- Homepage: https://www.order-n.org/
- Documentation: https://conquest.readthedocs.io/
- Source Repository: https://github.com/OrderN/CONQUEST-release
- License: MIT License (open-source)

## Overview
CONQUEST is a linear-scaling (O(N)) DFT code designed for massively parallel calculations on extremely large systems (up to millions of atoms). It uses local orbital basis sets and achieves excellent parallel scalability, making it ideal for biomolecules, nanostructures, and large-scale materials simulations on supercomputers.

**Scientific domain**: Very large biomolecules, nanostructures, materials (10,000+ atoms)  
**Target user community**: Researchers needing DFT for extremely large systems on supercomputers

## Theoretical Methods
- Density Functional Theory (DFT)
- Linear-scaling DFT (O(N) method)
- Exact diagonalization (for smaller systems)
- Blip functions (B-spline functions on a grid)
- Pseudo-atomic orbital (PAO) basis sets
- Multi-site support functions (MSSF)
- Norm-conserving pseudopotentials
- LDA, GGA functionals (via LibXC)
- van der Waals corrections (DFT-D2/D3, TS, vdW-DF)
- DFT+U for correlated systems
- Multiple time-step integration for MD

## Capabilities (CRITICAL)
- Ground-state electronic structure for very large systems
- True linear-scaling DFT (O(N) cost)
- Systems with millions of atoms demonstrated
- Massively parallel (scaling to ~200,000 cores)
- Geometry optimization for large systems
- Molecular dynamics (NVE, NVT, NPT)
- Large biomolecules (proteins, DNA)
- Nanostructures and materials
- Band structure and DOS
- Forces and stress tensors
- Mixed basis (blips + PAOs)
- Multisite support functions for accuracy
- Excellent weak and strong scaling
- Adaptive grid refinement
- Tersoff-Hamann STM images

**Sources**: Official CONQUEST documentation, cited in 6/7 source lists

## Key Strengths

### Linear-Scaling (O(N)) Methodology:
- Solves for density matrix directly
- Avoids Kohn-Sham eigenstate diagonalization
- Computational cost scales linearly
- Demonstrated for 1,000,000+ atoms
- Essential for truly large-scale DFT

### Flexible Basis Sets:
- Pseudo-atomic orbitals (PAO)
- Blip functions (B-splines on grid)
- Multi-site support functions (MSSF)
- Systematically improvable
- Up to 10,000 atoms with exact diagonalization

### Massive Parallelism:
- Scales to ~200,000 cores
- Excellent weak scaling
- Distributed memory architecture
- Tested on major supercomputers
- Production-ready HPC code

### Sparse Matrix Operations:
- Locality exploited throughout
- Efficient sparse storage
- Sparse matrix algebra
- Memory efficient for large systems

### Open-Source MIT License:
- Fully open-source
- MIT license (permissive)
- Active GitHub development
- 11 releases available
- Community contributions welcome

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
  - Charge density files
  - STM images

## Interfaces & Ecosystem
- **Parallelization**:
  - MPI for distributed memory
  - Excellent scaling to 1000s of processors
  - Tested on major supercomputers
  - Production HPC optimization
  
- **Basis sets**:
  - Blip functions (smooth B-splines)
  - PAO basis sets
  - Mixed blip/PAO approach
  - MSSF for higher accuracy
  
- **LibXC Integration**:
  - LDA functionals
  - GGA functionals
  - metaGGA (in development)
  - Hybrid functionals (in development)
  
- **Dispersion Corrections**:
  - DFT-D2 (Grimme)
  - DFT-D3 (Grimme)
  - Tkatchenko-Scheffler (TS)
  - vdW-DF functionals
  
- **Post-processing**:
  - Analysis utilities included
  - Compatible with standard visualization tools

## Advanced Features

### O(N) vs Exact Diagonalization:
- O(N): >1,000 atoms typical crossover
- Exact: Up to 1,000 atoms (PAO), 10,000 (MSSF)
- Seamless switching between methods
- Automatic selection possible

### Multi-Site Support Functions:
- Enhanced accuracy over simple PAOs
- Systematic improvability
- Larger systems with exact diagonalization
- Bridging localized and delocalized

### Calculation Types:
- Static electronic structure
- Density of states
- Band structure
- Structural optimization (ions and cell)
- Molecular dynamics (NVE, NVT, NPT)
- STM image simulation

### Output Capabilities:
- Total energies and forces
- Stress tensors
- Charge density
- Orbital density
- Tersoff-Hamann STM images

## Performance Characteristics
- **Speed**: Excellent for O(N) large systems
- **Accuracy**: Comparable to conventional DFT for equivalent basis
- **System size**: Up to millions of atoms
- **Memory**: Efficient sparse storage
- **Parallelization**: Excellent scaling to ~200,000 cores

## Computational Cost
- **O(N) DFT**: Linear scaling with system size
- **Diagonalization**: Standard cubic scaling
- **Memory**: Scales linearly for O(N)
- **I/O**: Optimized for HPC
- **Crossover**: ~1,000 atoms for O(N) advantage

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

## Comparison with Other Codes
- **vs SIESTA**: Both O(N), CONQUEST more HPC-focused
- **vs ONETEP**: Both O(N), different basis (blips vs NGWFs)
- **vs VASP/QE**: CONQUEST for much larger systems
- **vs BigDFT**: Different O(N) approaches (blips vs wavelets)
- **Unique strength**: Massive parallelism, up to millions of atoms, MIT license, HPC focus

## Application Areas

### Biomolecules:
- Large proteins
- DNA/RNA structures
- Protein-ligand complexes
- Enzyme mechanisms
- Solvated biomolecules

### Nanostructures:
- Nanoparticles
- Nanowires
- Nanotubes
- Complex interfaces
- Core-shell structures

### Materials Defects:
- Point defects in supercells
- Dislocations
- Grain boundaries
- Surfaces with reconstructions
- Amorphous materials

### Device Simulation:
- Semiconductor nanostructures
- Heterostructures
- Large device models
- Interface engineering

## Best Practices

### Choosing O(N) vs Exact:
- O(N) for >1,000 atoms
- Exact for high accuracy on smaller systems
- Test convergence with both if possible
- Consider MSSF for intermediate sizes

### Basis Set Selection:
- Start with standard PAO sets
- Test blip-only for comparison
- MSSF for enhanced accuracy
- Check cutoff radius convergence

### Parallelization:
- Match atoms to processors
- Aim for good load balance
- Use appropriate block sizes
- Benchmark scaling

### Convergence:
- Grid spacing for blips
- Support function radii
- SCF convergence criteria
- L tolerance for O(N)

## Community and Support
- Open-source MIT license
- GitHub development
- Developer mailing list
- ReadTheDocs documentation
- Regular releases

## Verification & Sources
**Primary sources**:
1. Official website: https://www.order-n.org/
2. Documentation: https://conquest.readthedocs.io/
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
- Documentation: ACCESSIBLE (ReadTheDocs)
- Source code: OPEN (GitHub, MIT license)
- Community support: Active (developers, mailing list)
- Academic citations: >300
- Active development: Regular releases (11 on GitHub), HPC optimization
- Specialized strength: Linear-scaling O(N), millions of atoms, massive parallelism, MIT open-source, HPC focus
