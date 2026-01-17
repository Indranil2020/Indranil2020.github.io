# FreeON (formerly MondoSCF)

## Official Resources
- Homepage: https://github.com/FreeON/freeon
- Documentation: https://github.com/FreeON/freeon/wiki
- Source Repository: https://github.com/FreeON/freeon
- License: GNU General Public License v3.0

## Overview
FreeON is an experimental, open-source suite of programs for linear-scaling quantum chemistry calculations. Formerly known as MondoSCF, it performs Hartree-Fock, pure DFT, and hybrid DFT calculations using a Cartesian-Gaussian LCAO basis. FreeON emphasizes O(N) and O(N log N) algorithms for non-metallic systems.

**Scientific domain**: Molecules, clusters, large organic systems, biomolecules  
**Target user community**: Researchers needing linear-scaling molecular quantum chemistry for large non-metallic systems

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- Hybrid HF/DFT (B3LYP and others)
- Cartesian-Gaussian LCAO basis
- Linear-scaling O(N) algorithms
- O(N log N) Coulomb evaluation
- DIIS convergence acceleration
- Norm-conserving pseudopotentials (optional)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Hartree-Fock calculations
- Pure DFT (LDA, GGA)
- Hybrid DFT (B3LYP, etc.)
- Linear-scaling O(N) Fock matrix construction
- Geometry optimization
- Molecular dynamics
- Forces and gradients
- Frequency calculations (Hessian)
- Large molecular systems

**Sources**: Wikipedia, GitHub repository, academic publications

## Key Strengths

### Linear-Scaling Focus:
- O(N) Fock matrix construction
- O(N log N) Coulomb operator
- Targeted at non-metallic systems
- Scales to large molecules
- Efficient for insulators

### Open Source GPL v3:
- Fully open source
- Community contributions welcome
- Transparent implementation
- Academic and commercial use

### Hybrid Functionals:
- B3LYP support
- PBE0 support
- Exact exchange with O(N)
- Accurate thermochemistry

### Robust Implementation:
- DIIS convergence
- Multiple SCF algorithms
- Fortran95 and C code
- Well-structured codebase

## Inputs & Outputs
- **Input formats**:
  - Native input file format
  - XYZ coordinates
  - Basis set specifications
  - Calculation parameters
  
- **Output data types**:
  - Total energies
  - Orbital energies
  - Forces and gradients
  - Optimized geometries
  - Molecular dynamics trajectories
  - Frequency data

## Interfaces & Ecosystem
- **Standalone operation**:
  - Self-contained package
  - Built-in basis sets
  - Integrated geometry optimizer
  
- **Build system**:
  - Autoconf/Automake
  - GNU toolchain
  - Linux/Unix platforms
  - FreeBSD support

## Advanced Features

### Linear-Scaling Methods:
- Fast Multipole Method (FMM)
- Hierarchical matrix algebra
- Sparse matrix techniques
- Locality exploitation

### Geometry Optimization:
- Quasi-Newton methods
- Conjugate gradient
- Constraints support
- Transition state search

### Molecular Dynamics:
- Born-Oppenheimer MD
- Velocity Verlet integrator
- Temperature control
- Trajectory analysis

### Frequency Calculations:
- Hessian evaluation
- Vibrational analysis
- IR intensities
- Thermodynamic properties

## Performance Characteristics
- **Speed**: O(N) for large non-metallic systems
- **Accuracy**: Standard Gaussian basis precision
- **System size**: Hundreds to thousands of atoms
- **Memory**: Efficient sparse matrix storage
- **Parallelization**: Some parallel capabilities

## Computational Cost
- **Linear scaling**: Achieved for insulators/semiconductors
- **Hybrid DFT**: Efficient O(N) exact exchange
- **Break-even**: ~100-500 atoms vs cubic codes
- **Typical**: Competitive for large molecules

## Limitations & Known Constraints
- **Metallic systems**: Linear scaling breaks down
- **Development status**: Experimental, less active recently
- **Community**: Smaller than major codes
- **Documentation**: Academic-level
- **Periodic systems**: Not primary focus
- **Platform**: Linux/Unix only
- **Basis sets**: Limited built-in selection

## Comparison with Other Codes
- **vs Gaussian**: FreeON O(N) focus vs Gaussian robustness
- **vs NWChem**: FreeON specialized linear-scaling
- **vs ONETEP**: FreeON molecular, ONETEP materials
- **vs BigDFT**: FreeON Gaussian, BigDFT wavelets
- **Unique strength**: Open-source O(N) molecular quantum chemistry

## Application Areas

### Large Organic Molecules:
- Natural products
- Pharmaceuticals
- Polymers
- Dendrimers

### Biomolecular Systems:
- Protein fragments
- Nucleic acid segments
- Enzyme models
- Drug-receptor complexes

### Molecular Clusters:
- Water clusters
- Hydrogen-bonded systems
- Molecular crystals (fragments)
- Host-guest complexes

### Materials Fragments:
- Molecular solids
- Organic semiconductors
- Self-assembled monolayers
- Nanoparticle ligands

## Best Practices

### System Selection:
- Best for non-metallic systems
- Check band gap for linear scaling
- Use for insulators/semiconductors
- Avoid highly metallic systems

### Basis Set Choice:
- Start with standard basis (6-31G*)
- Test convergence with larger basis
- Document basis for reproducibility

### SCF Convergence:
- Use DIIS for difficult cases
- Monitor energy convergence
- Adjust mixing if needed

### Linear-Scaling Settings:
- Tune cutoff parameters
- Balance accuracy and speed
- Verify against standard calculation

## Community and Support
- Open-source GPL v3
- GitHub repository
- Academic publications
- Legacy from LANL development
- Contributor community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/FreeON/freeon
2. Wikipedia: FreeON article
3. M. Challacombe et al., J. Chem. Phys. publications

**Secondary sources**:
1. Linear-scaling quantum chemistry literature
2. Open-source chemistry software surveys

**Confidence**: VERIFIED - Open source, published methodology

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Academic use: Published applications
- Documentation: Wiki and papers
- Development status: Experimental/Mature
- Specialty: O(N) linear-scaling molecular quantum chemistry
