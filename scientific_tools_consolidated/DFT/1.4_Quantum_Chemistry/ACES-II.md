# ACES II

## Official Resources
- Homepage: https://www.qtp.ufl.edu/ACES/
- Documentation: UFL Quantum Theory Project
- Source: Academic license via QTP
- License: Academic use

## Overview
ACES II (Advanced Concepts in Electronic Structure II) is an ab initio quantum chemistry program developed at the University of Florida's Quantum Theory Project (QTP) by Rodney Bartlett and collaborators. It pioneered many high-level coupled-cluster implementations and is the predecessor to both ACES III and CFOUR, representing a seminal contribution to coupled-cluster methodology.

**Scientific domain**: High-accuracy coupled-cluster quantum chemistry  
**Target user community**: Historic; methods and algorithms now in CFOUR and ACES III

## Theoretical Methods
- Coupled Cluster (CCSD, CCSDT, CCSDTQ, full CC)
- Many-body perturbation theory (MBPT2-4)
- Equation-of-Motion CC (EOM-CC) for excited states
- Analytical first and second derivatives
- Property calculations
- IP/EA-EOM for ionized/attached states
- Analytic gradient-based optimizations

## Capabilities (CRITICAL)
- High-accuracy coupled cluster energies
- Full CCSDT and CCSDTQ implementations
- Analytical gradients for geometry optimization
- EOM-CC for excited states
- Property calculations via response theory
- Bartlett group methods
- Benchmark-quality results
- Foundation for CFOUR development
- Parallel capabilities in later versions

## Key Strengths

### Coupled Cluster Theory:
- High-level truncations (CCSDT, CCSDTQ)
- Rigorous implementations
- Analytical derivatives
- Benchmark accuracy
- Size-extensivity

### Excited States:
- EOM-CCSD
- EOM-CCSDT
- IP/EA-EOM variants
- Transition properties
- State-of-the-art methods

### Properties:
- Response theory
- Dipole moments
- Polarizabilities
- NMR parameters
- Various molecular properties

### Research Foundation:
- Bartlett group development
- Led to CFOUR
- Led to ACES III
- Training ground for developers

## Inputs & Outputs
- **Input formats**:
  - ACES II input files
  - ZMAT coordinates
  - Basis set specification
  
- **Output data types**:
  - Correlation energies
  - Optimized geometries
  - Properties
  - Excited state data

## Interfaces & Ecosystem
- **Standalone**: Complete program
- **Successors**: CFOUR, ACES III
- **Integration**: Academic research use

## Advanced Features

### High-Level CC:
- Full CCSDT without approximations
- CCSDTQ for benchmarks
- Arbitrary truncations
- Orbital optimization

### EOM Methods:
- Excitation energies
- Ionization potentials
- Electron affinities
- Transition moments

### Analytical Derivatives:
- First derivatives (gradients)
- Second derivatives (Hessians)
- Response properties
- Geometry optimization

### Method Development:
- Research platform
- New CC variants
- Testing ground
- Publication vehicle

## Performance Characteristics
- **Speed**: Standard CC scaling
- **Accuracy**: Spectroscopic accuracy possible
- **System size**: Limited by CC scaling
- **Memory**: Large for high-level CC
- **Era**: Competitive for 1990s-2000s

## Computational Cost
- **CCSD**: O(N^6)
- **CCSDT**: O(N^8)
- **CCSDTQ**: O(N^10)
- **Typical**: Small molecules for high accuracy

## Limitations & Known Constraints
- **Succeeded**: By CFOUR and ACES III
- **Support**: Limited for ACES II itself
- **Modern features**: In successor codes
- **Community**: Use successors now
- **Historic**: Important legacy

## Comparison with Other Codes
- **vs CFOUR**: CFOUR is direct successor/evolution
- **vs ACES III**: ACES III is parallel extension
- **vs Gaussian**: ACES II more CC-focused
- **vs MOLPRO**: Different origins, overlapping methods
- **Legacy**: Foundational for modern CC codes

## Application Areas

### Spectroscopic Accuracy:
- Bond energies
- Reaction barriers
- Equilibrium geometries
- Vibrational frequencies

### Excited States:
- Electronic spectra
- Transition moments
- State orderings
- Photochemistry

### Benchmarking:
- Reference calculations
- Method validation
- Basis set limits
- Accurate thermochemistry

## Historical Context

### Development Timeline:
- Early 1990s: Initial development
- Mid 1990s: EOM-CC implementations
- Late 1990s: ACES II-MAB branch (→CFOUR)
- 2000s: ACES III for parallel
- Present: Use CFOUR or ACES III

### Key Contributors:
- Rodney Bartlett (PI)
- John Stanton (→CFOUR)
- Jürgen Gauss (→CFOUR)
- Many postdocs and students

### Seminal Papers:
- Coupled cluster theory papers
- EOM-CC methodology
- Analytical derivative theory
- Benchmark applications

## Community and Support
- Historic QTP development
- Academic licensing
- CFOUR as successor
- Extensive publication record
- Bartlett group legacy

## Verification & Sources
**Primary sources**:
1. UFL QTP: https://www.qtp.ufl.edu/
2. Stanton, Bartlett et al., ACES II papers
3. Review of CC theory (Bartlett, Musial)
4. CFOUR website references ACES II history

**Confidence**: VERIFIED (Historic)
- Status: Historic, succeeded by CFOUR/ACES III
- Significance: Pioneer in high-level CC
- Impact: CFOUR, ACES III descended from it
- Methods: State-of-the-art CC implementations
