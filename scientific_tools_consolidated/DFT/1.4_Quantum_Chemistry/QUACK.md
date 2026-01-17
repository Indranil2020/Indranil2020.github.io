# QUACK

## Official Resources
- Homepage: https://github.com/pfloos/QuACK
- Documentation: https://github.com/pfloos/QuACK/wiki
- Source Repository: https://github.com/pfloos/QuACK
- License: GNU General Public License v3.0

## Overview
QUACK (Quantum Chemistry in ACKnowledgement) is an open-source software for emerging quantum electronic structure methods. It specializes in Green's function methods including GW approximation and Bethe-Salpeter equation (BSE) for molecular systems, along with coupled cluster and other advanced correlation methods.

**Scientific domain**: Many-body perturbation theory, GW, BSE, coupled cluster for molecules  
**Target user community**: Researchers developing and applying GW/BSE methods to molecular systems

## Theoretical Methods
- Hartree-Fock (RHF, UHF, GHF)
- Random Phase Approximation (RPA)
- GW approximation (G0W0, evGW, qsGW)
- Bethe-Salpeter Equation (BSE)
- Coupled Cluster (CCD, CCSD)
- Equation-of-Motion Coupled Cluster
- ADC methods
- Configuration Interaction
- Møller-Plesset Perturbation Theory

## Capabilities (CRITICAL)
- One-shot GW (G0W0)
- Eigenvalue self-consistent GW (evGW)
- Quasiparticle self-consistent GW (qsGW)
- BSE for optical excitations
- Molecular GW calculations
- Ionization potentials and electron affinities
- Excitation energies
- Coupled cluster correlation
- Educational/development platform

## Key Strengths

### GW Methods:
- Multiple self-consistency schemes
- Molecular focus
- Quasiparticle energies
- Spectral functions
- Benchmark quality

### BSE Implementation:
- Optical excitations
- Electron-hole interactions
- Absorption spectra
- Exciton binding energies

### Research Platform:
- Method development
- Testing new approximations
- Benchmark calculations
- Educational purposes

### Clean Implementation:
- Readable Fortran code
- Well-structured
- Easy to modify
- Good documentation

## Inputs & Outputs
- **Input formats**:
  - QUACK input files
  - Molecular coordinates
  - Basis specifications
  
- **Output data types**:
  - Quasiparticle energies
  - Excitation energies
  - Correlation energies
  - Spectral data

## Interfaces & Ecosystem
- **Integral packages**: External integral generation
- **Basis sets**: Standard basis sets
- **Post-processing**: Text output parsing

## Advanced Features

### Self-Consistent GW:
- Partial self-consistency
- Full self-consistency
- Vertex corrections
- Beyond G0W0

### BSE Optical Properties:
- Singlet/triplet excitations
- Oscillator strengths
- Exciton analysis
- Optical gaps

### Coupled Cluster:
- Ground state correlation
- EOM for excitations
- Comparison with GW/BSE

## Performance Characteristics
- **Speed**: Efficient for molecules
- **Accuracy**: High-level many-body methods
- **System size**: Small to medium molecules
- **Memory**: Standard requirements
- **Parallelization**: OpenMP threading

## Computational Cost
- **HF**: Fast baseline
- **G0W0**: O(N^4) with approximations
- **evGW/qsGW**: Multiple iterations
- **BSE**: Depends on excitation number
- **Typical**: Suitable for benchmarks

## Limitations & Known Constraints
- **System size**: Best for small-medium molecules
- **Periodicity**: Molecular focus only
- **Features**: Specialized scope
- **Community**: Research group centered
- **Production**: More for development/benchmarks
- **Gradients**: Not available

## Comparison with Other Codes
- **vs molgw**: Similar GW focus, different implementations
- **vs BerkeleyGW**: QUACK molecular, BerkeleyGW periodic
- **vs Turbomole**: QUACK more GW variants
- **vs ORCA**: QUACK specialized GW, ORCA general purpose
- **Unique strength**: GW method variety, educational platform

## Application Areas

### Spectroscopy:
- Photoelectron spectra
- Optical absorption
- Ionization energies
- Electron affinities

### Method Development:
- Testing GW variants
- BSE approximations
- Vertex corrections
- Benchmark sets

### Molecular Properties:
- HOMO-LUMO gaps
- Fundamental gaps
- Optical gaps
- Quasiparticle renormalization

## Best Practices

### GW Calculations:
- Appropriate starting point
- Basis set convergence
- Self-consistency level
- Frequency treatment

### BSE:
- Number of states
- Kernel approximations
- TDA vs full BSE
- Convergence checks

## Community and Support
- Open-source GPL v3
- Active GitHub development
- Academic research group
- Publications and benchmarks
- Growing community

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/pfloos/QuACK
2. Loos et al. publications on GW/BSE
3. Benchmark studies (GW100, etc.)
4. Active development

**Confidence**: VERIFIED
- Source code: OPEN (GitHub, GPL v3)
- Documentation: Wiki available
- Active development: Yes
- Academic citations: Yes
