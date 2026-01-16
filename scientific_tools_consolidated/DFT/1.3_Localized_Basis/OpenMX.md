# OpenMX

## Official Resources
- Homepage: http://www.openmx-square.org/
- Documentation: http://www.openmx-square.org/openmx_man3.9/
- Source Repository: http://www.openmx-square.org/ (download page)
- License: GNU General Public License v3.0

## Overview
OpenMX (Open source package for Material eXplorer) is an efficient DFT code using localized pseudo-atomic orbitals (PAO) with particular strength in large-scale calculations, non-collinear magnetism, and spin-orbit coupling. It provides excellent performance for complex magnetic systems, topological materials, and spintronics applications.

**Scientific domain**: Magnetism, spintronics, topological materials, large systems  
**Target user community**: Researchers studying magnetic materials, spin-orbit physics, topological properties, large-scale systems

## Theoretical Methods
- Density Functional Theory (DFT)
- Pseudo-atomic orbital (PAO) basis sets
- Norm-conserving pseudopotentials
- LDA, GGA functionals
- DFT+U for correlated systems
- van der Waals corrections
- Spin-orbit coupling (fully relativistic, self-consistent)
- Non-collinear magnetism (unconstrained)
- Constrained DFT
- Effective screening medium (ESM) method
- O(N) Krylov subspace method for large systems

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and MD
- Large-scale systems (1000+ atoms with O(N))
- Non-collinear magnetism with spin-orbit coupling
- Magnetic anisotropy energy
- Rashba and Dresselhaus spin splitting
- Topological properties (Z2 invariants, Chern numbers)
- Band structure including spin texture
- Wannier functions and maximally localized Wannier functions
- Quantum transport (NEGF method)
- STM image simulation
- Optical conductivity
- Berry phase calculations
- Electric polarization
- Orbital magnetization
- ESM method for slab calculations
- Linear-scaling DFT
- Anomalous Hall conductivity

**Sources**: Official OpenMX documentation, cited in 7/7 source lists

## Key Strengths

### Spin-Orbit Coupling:
- Full self-consistent SOC implementation
- Unconstrained non-collinear DFT
- Explicit spin-orbit in all calculations
- Accurate for heavy elements
- Essential for topological studies

### Topological Materials Analysis:
- Z2 topological invariant (Fukui-Hatsugai method)
- Chern number calculations
- Berry phase and curvature
- Wilson loop calculations
- Parity calculations
- Anomalous Hall conductivity

### Spin Texture Analysis:
- kSpin post-processing code
- K-space spin texture resolution
- Rashba/Dresselhaus spin splitting
- Atom-resolved spin contributions
- Orbital-resolved spin analysis

### Non-Collinear Magnetism:
- Fully unconstrained spins
- Complex magnetic structures
- Spin spirals
- Magnetic anisotropy
- Exchange interactions

### Boundary State Calculations:
- Slab models for surfaces
- Green's function methods
- Topological surface states
- Edge states in 2D materials

## Inputs & Outputs
- **Input formats**:
  - Input file (OpenMX format)
  - Coordinate files (XYZ, PDB)
  - PAO basis definitions
  - Pseudopotential files
  
- **Output data types**:
  - Standard output with energies, forces
  - Band structure files with spin information
  - DOS and PDOS files
  - Density and spin density files
  - Wannier function outputs
  - Transmission coefficients for transport
  - Topological invariant data

## Interfaces & Ecosystem
- **Post-processing tools**:
  - OpenMX Viewer - visualization
  - Z2FH - Z2 invariant calculation
  - kSpin - spin texture analysis
  - Band unfolding tools
  - Various analysis utilities included
  
- **Transport calculations**:
  - Built-in NEGF for quantum transport
  - Electrode-device-electrode setup
  - Multi-terminal configurations
  
- **Topological analysis**:
  - Berry phase module
  - Wilson loop calculations
  - Chern number computation
  - Anomalous Hall conductivity
  
- **Workflow integration**:
  - Can be interfaced with ASE
  - Compatible with standard workflow tools

## Advanced Features

### ADPACK:
- Pseudopotential and PAO generator
- Fully relativistic pseudopotentials
- Optimized for OpenMX
- User-customizable basis sets

### VPS/PAO Databases:
- Ver. 2019 standard database
- Core excitation specialized database
- Ready-to-use basis sets
- Validated for many elements

### Technical Notes:
- In-depth methodology documentation
- Application examples
- Best practices guides
- Algorithm explanations

### Video Lectures:
- Educational materials
- Tutorial walkthroughs
- Research presentations
- Workshop recordings

## Performance Characteristics
- **Speed**: Efficient for magnetic systems
- **Accuracy**: Good for PAO basis calculations
- **System size**: Up to thousands of atoms with O(N)
- **Memory**: Generally efficient
- **Parallelization**: MPI parallelization; good scaling

## Computational Cost
- **DFT**: Efficient PAO implementation
- **SOC**: Moderate additional cost
- **Non-collinear**: ~2x spin-polarized cost
- **Transport**: Moderate NEGF overhead
- **Topological**: Post-processing mostly

## Limitations & Known Constraints
- **Basis sets**: PAO basis requires convergence testing
- **Pseudopotentials**: Limited to norm-conserving
- **Documentation**: Comprehensive but English translations vary in quality
- **Community**: Smaller than VASP/QE; primarily Japan-based
- **Installation**: Requires compilation; dependencies (BLAS, LAPACK, FFT)
- **Parallelization**: MPI parallelization; efficiency varies
- **Memory**: Generally efficient but depends on basis size
- **k-point sampling**: Required for periodic systems
- **Platform**: Primarily Linux/Unix

## Comparison with Other Codes
- **vs VASP**: OpenMX localized basis, VASP plane-wave; OpenMX better for SOC details
- **vs SIESTA**: Similar approach, OpenMX stronger for magnetism/topology
- **vs FHI-aims**: OpenMX pseudopotential, FHI-aims all-electron
- **vs Quantum ESPRESSO**: OpenMX better for non-collinear SOC
- **Unique strength**: Comprehensive topological invariant tools, spin texture analysis, non-collinear SOC, open-source

## Application Areas

### Topological Materials:
- Topological insulators (TIs)
- Weyl semimetals
- Node-line semimetals
- Topological crystalline insulators
- Higher-order TIs

### Spintronics:
- Spin Hall effect
- Rashba/Dresselhaus systems
- Magnetic tunnel junctions
- Spin-orbit torque
- Magnetization dynamics

### Magnetic Materials:
- Complex magnets
- Frustrated systems
- Magnetic anisotropy
- Exchange coupling
- Spin spirals

### 2D Materials:
- Graphene spintronics
- TMD magnetism
- Van der Waals magnets
- Heterostructure topology

## Best Practices

### Basis Set Selection:
- Standard vs. precise PAOs
- Convergence with cutoff radius
- Semi-core state inclusion
- Reference to database recommendations

### SOC Calculations:
- Full relativistic pseudopotentials
- Converge without SOC first
- Check spin texture convergence
- Compare collinear vs non-collinear

### Topological Analysis:
- Sufficient k-point mesh
- Check gauge consistency
- Validate with multiple methods
- Surface/edge state verification

### Convergence:
- Energy cutoff for grid
- k-point sampling
- SCF convergence criteria
- Spin convergence for magnets

## Community and Support
- Open-source GPL v3
- OpenMX Forum for support
- Developer meetings (annual)
- Video lecture resources
- Japan-based core team

## Verification & Sources
**Primary sources**:
1. Official website: http://www.openmx-square.org/
2. Manual: http://www.openmx-square.org/openmx_man3.9/
3. T. Ozaki, Phys. Rev. B 67, 155108 (2003) - OpenMX method
4. T. Ozaki and H. Kino, Phys. Rev. B 69, 195113 (2004) - O(N) method
5. T. Ozaki et al., Phys. Rev. B 81, 035116 (2010) - Krylov subspace

**Secondary sources**:
1. OpenMX tutorials and examples
2. Published applications in spintronics and topology
3. Benchmark studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (download from website)
- Community support: Active (forum, Japanese community, meetings)
- Academic citations: >1,000 (main papers)
- Active development: Regular updates, Patch 3.9.9 (Oct 2021)
- Specialized strength: Spin-orbit coupling, topological invariants, non-collinear magnetism, spin texture analysis, open-source
