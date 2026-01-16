# ONETEP

## Official Resources
- Homepage: https://onetep.org/
- Documentation: https://onetep.org/documentation/
- Source Repository: Available to licensed users
- License: Academic license (free for academics, registration required)

## Overview
ONETEP (Order-N Electronic Total Energy Package) is a linear-scaling DFT code that achieves plane-wave accuracy with localized orbitals. It uniquely combines the accuracy of plane-wave calculations with the efficiency of linear-scaling methods through Non-orthogonal Generalized Wannier Functions (NGWFs), enabling calculations on systems with thousands to tens of thousands of atoms.

**Scientific domain**: Large biomolecules, nanostructures, materials with thousands of atoms  
**Target user community**: Researchers needing DFT accuracy for very large systems (1000-10000+ atoms)

## Theoretical Methods
- Density Functional Theory (DFT)
- Linear-scaling DFT (O(N) method)
- Non-orthogonal generalized Wannier functions (NGWFs)
- Periodic cardinal sine (psinc) basis equivalent to plane-waves
- Plane-wave accuracy with localized basis
- Norm-conserving pseudopotentials
- LDA, GGA functionals
- Hybrid functionals (range-separated)
- DFT+U for correlated systems
- van der Waals corrections (DFT-D, TS)
- Implicit solvation models (DDCOSMO)
- TDDFT for excited states
- Finite electronic temperature

## Capabilities (CRITICAL)
- Ground-state electronic structure for very large systems
- Linear-scaling DFT (computational cost scales linearly with system size)
- Plane-wave accuracy with localized orbitals
- Geometry optimization for large systems
- Molecular dynamics (NVE, NVT, NPT)
- Systems with 1000-10000+ atoms (14,000+ demonstrated)
- Protein and biomolecule calculations
- Nanostructures and materials
- Band structure and DOS
- Forces and stress tensors
- NMR chemical shifts
- EPR parameters
- Core-level spectroscopy
- Protein-ligand binding free energies
- Implicit solvation (DDCOSMO)
- TDDFT for absorption spectra
- Wannier function analysis
- Conduction calculations
- Ensemble DFT
- QM/MM capabilities

**Sources**: Official ONETEP documentation, cited in 7/7 source lists

## Key Strengths

### Non-orthogonal Generalized Wannier Functions (NGWFs):
- Spatially localized orbitals
- Optimized in situ during calculation
- Psinc basis (plane-wave equivalent)
- Systematic accuracy improvement
- Species-dependent localization radii

### Plane-Wave Accuracy:
- Equivalent to large plane-wave basis sets
- Systematic convergence
- No basis set superposition error
- Transferable across systems
- Benchmark-quality results

### Linear Scaling (O(N)):
- Density matrix optimization
- Avoid eigenstate diagonalization
- Computational cost linear with size
- Previously unattainable system sizes
- Thousands of cores parallelization

### Biomolecular Applications:
- Full QM treatment of proteins
- Protein-ligand binding energies
- Charge transfer and polarization
- DNA electronic structure calculations
- Enzyme catalysis studies

### Comprehensive Properties:
- NMR chemical shifts
- EPR g-tensors
- Core-level spectroscopy
- Optical absorption
- Electronic transport

## Inputs & Outputs
- **Input formats**:
  - Input file (ONETEP format)
  - PDB, XYZ coordinate files
  - Pseudopotential files
  - NGWF initial guesses
  
- **Output data types**:
  - Standard output with energies, forces
  - Optimized structures
  - Density files (cube format)
  - DOS and PDOS
  - NGWF outputs
  - Distributed multipole analysis
  - Property-specific outputs

## Interfaces & Ecosystem
- **Framework integrations**:
  - QM/MM capabilities (embedding)
  - Molecular dynamics interfaces
  - i-PI compatibility
  
- **Visualization**:
  - Compatible with standard visualization tools
  - Cube file output for densities
  - NGWF visualization
  
- **Analysis Tools**:
  - Fragment charge analysis
  - Distributed multipole analysis
  - Interaction energy decomposition
  - Population analysis
  
- **Solvation Models**:
  - DDCOSMO implicit solvation
  - Environment effects
  - Solvation free energies

## Advanced Features

### Protein-Ligand Binding:
- Full QM binding free energies
- Gas-phase and solvation contributions
- Charge transfer effects
- Polarization captured correctly
- Drug discovery applications

### Electronic Transport:
- Conductance calculations
- Non-equilibrium systems
- Nanoscale junctions
- Material interfaces

### TDDFT:
- Optical absorption spectra
- Excited state properties
- Large-system optical properties
- Material characterization

### ONETEP 7.2 (2024):
- Latest academic release
- Performance improvements
- New features added
- Active development

### Finite Electronic Temperature:
- Metallic systems
- Fractional occupations
- Convergence improvement
- Extended applications

## Performance Characteristics
- **Speed**: Efficient O(N) implementation
- **Accuracy**: Plane-wave equivalent
- **System size**: Up to 14,000+ atoms demonstrated
- **Memory**: Lower than conventional for large systems
- **Parallelization**: Excellent scaling to thousands of cores

## Computational Cost
- **O(N) DFT**: Linear scaling achieved
- **NGWF optimization**: Moderate overhead
- **Hybrid functionals**: More expensive but feasible
- **Properties**: Variable depending on type
- **Crossover**: Benefits at ~500+ atoms

## Limitations & Known Constraints
- **Academic license**: Free for academics but requires registration
- **Not fully open-source**: Source available to licensed users only
- **Learning curve**: Linear-scaling methods and NGWFs require understanding
- **NGWF optimization**: Can be challenging to converge for some systems
- **Pseudopotentials**: Limited to norm-conserving
- **Hybrid functionals**: Computationally expensive even with linear-scaling
- **Parallelization**: Excellent but requires understanding of distribution
- **Memory**: Lower than conventional DFT but still significant for very large systems
- **Installation**: Requires compilation and libraries
- **Platform**: Primarily Linux/Unix, HPC systems

## Comparison with Other Codes
- **vs CONQUEST**: Both O(N), ONETEP plane-wave accuracy, CONQUEST blip basis
- **vs SIESTA**: ONETEP more accurate per atom, SIESTA faster
- **vs VASP/QE**: ONETEP for larger systems with similar accuracy
- **vs FHI-aims**: Different O(N) approaches, both high accuracy
- **Unique strength**: Plane-wave accuracy at O(N) cost, NGWF technology, biomolecular applications

## Application Areas

### Drug Discovery:
- Protein-ligand binding
- Structure-based drug design
- Binding affinity prediction
- QM/MM studies
- Fragment-based analysis

### Biomolecular Science:
- Amyloid fibril structures
- DNA electronic properties
- Enzyme mechanisms
- Protein folding effects
- Large biomolecular complexes

### Nanomaterials:
- Nanoparticles
- Functionalized surfaces
- Graphene and 2D materials
- Carbon nanotubes
- Nanoscale devices

### Materials Science:
- Defects in semiconductors
- Interface properties
- Molecular crystals
- Charged adsorbates
- Battery materials

## Best Practices

### NGWF Optimization:
- Choose appropriate radii per species
- Start from good initial guess
- Monitor convergence carefully
- Adjust localization if needed

### System Setup:
- Use quality pseudopotentials
- Check NGWF coverage
- Appropriate psinc grid spacing
- Localization consistent with system

### Large Calculations:
- Parallel scaling tests
- Memory estimation
- I/O optimization
- Checkpoint strategies

### Convergence:
- NGWF convergence threshold
- Density kernel tolerance
- Energy convergence
- Grid cutoff energy

## Community and Support
- Academic license model
- User workshops (Masterclass events)
- Mailing list support
- Active development team
- Regular releases

## Verification & Sources
**Primary sources**:
1. Official website: https://onetep.org/
2. Documentation: https://onetep.org/documentation/
3. C.-K. Skylaris et al., J. Chem. Phys. 122, 084119 (2005) - ONETEP method
4. N. D. M. Hine et al., Comput. Phys. Commun. 180, 1041 (2009) - Linear-scaling

**Secondary sources**:
1. ONETEP tutorials and workshops (Masterclass 2024)
2. Published large-scale biomolecule applications
3. Linear-scaling benchmark studies
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE (requires registration for full access)
- Source code: Available to licensed users
- Community support: Active (user mailing list, workshops, Masterclass)
- Academic citations: >500 (main papers)
- Active development: Regular releases (v7.2 in 2024), well-maintained
- Specialized strength: Plane-wave accuracy at O(N) cost, NGWF technology, biomolecular applications, protein-ligand binding
