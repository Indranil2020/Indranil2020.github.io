# molgw (Molecular GW)

## Official Resources
- Homepage: https://www.molgw.org/
- Documentation: https://www.molgw.org/doc.html
- Source Repository: https://github.com/molgw/molgw
- License: GNU GPL v3

## Overview
molgw is an open-source quantum chemistry code implementing many-body perturbation theory (GW approximation and Bethe-Salpeter equation) for molecules and clusters using Gaussian basis sets. Developed by Fabien Bruneval (CEA, France), molgw focuses on efficient GW calculations for finite systems with emphasis on accurate ionization potentials, electron affinities, and optical excitations. It provides a user-friendly, lightweight alternative to larger quantum chemistry packages for MBPT calculations.

**Scientific domain**: GW approximation, BSE, molecular systems, Gaussian basis  
**Target user community**: Molecular physicists, quantum chemists, GW method developers

## Theoretical Methods
- GW approximation (G₀W₀, evGW, qsGW)
- Bethe-Salpeter Equation (BSE)
- Hartree-Fock (HF)
- Hybrid functionals (PBE0, B3LYP)
- Range-separated hybrids (RSH)
- MP2 and RPA
- Gaussian basis sets
- Resolution-of-identity (RI)
- Auxiliary basis sets

## Capabilities (CRITICAL)
- GW quasiparticle energies
- Ionization potentials (IP)
- Electron affinities (EA)
- Fundamental gaps
- Optical excitations (BSE)
- Absorption spectra
- Self-consistent GW variants
- Small to medium molecules
- Clusters and finite systems
- Open-source and free
- Efficient algorithms
- Educational tool

**Sources**: molgw website (https://www.molgw.org/)

## Key Strengths

### Molecular Focus:
- Optimized for molecules
- Finite system algorithms
- No periodicity
- Cluster calculations
- Size-appropriate methods

### Gaussian Basis:
- Standard quantum chemistry basis
- Extensive basis libraries
- Systematic convergence
- Efficient for molecules
- Well-established

### GW Variants:
- G₀W₀ (perturbative)
- evGW (eigenvalue self-consistent)
- qsGW (quasiparticle self-consistent)
- User choice
- Method comparison

### Open Source:
- GNU GPL license
- Free software
- Transparent code
- Community contributions
- Educational value

### Efficiency:
- Resolution-of-identity
- Auxiliary basis optimization
- Efficient algorithms
- Reasonable cost
- Production calculations

## Inputs & Outputs
- **Input formats**:
  - Simple text input
  - Molecular geometries (xyz, etc.)
  - Basis set specifications
  - Calculation parameters
  
- **Output data types**:
  - Quasiparticle energies
  - Ionization potentials
  - Electron affinities
  - Excitation energies
  - Absorption spectra
  - Detailed output files

## Interfaces & Ecosystem
- **Standard Formats**:
  - xyz geometries
  - Gaussian basis libraries
  - Standard output
  
- **Visualization**:
  - Plotting tools
  - Spectral data
  - Standard formats

## Workflow and Usage

### Typical Workflow:
1. Prepare molecular geometry
2. Choose basis set
3. Run HF or DFT calculation
4. Perform GW calculation
5. Analyze quasiparticle energies
6. Optional: BSE for excitations

### Input Example:
```fortran
# molgw input
molecule benzene
basis cc-pVTZ
scf PBE0
postscf G0W0
```

### Running molgw:
```bash
molgw input_file
```

## Advanced Features

### Self-Consistent GW:
- evGW: Eigenvalue SC
- qsGW: Quasiparticle SC
- Improved accuracy
- Systematic approach
- Research capability

### BSE Implementation:
- Optical excitations
- Singlet and triplet
- Absorption spectra
- Oscillator strengths
- Molecular spectroscopy

### Resolution-of-Identity:
- RI-GW
- Auxiliary basis sets
- Computational efficiency
- Controlled approximation
- Standard technique

### Range-Separated Hybrids:
- RSH starting point
- Optimal for GW
- Reduced error
- Systematic improvement
- Modern approach

## Performance Characteristics
- **Speed**: Good for molecules (RI acceleration)
- **Accuracy**: Excellent for IPs/EAs
- **System size**: Small to medium molecules
- **Scaling**: Standard GW scaling
- **Typical**: Research and benchmarks

## Computational Cost
- **GW**: Moderate for molecules
- **RI**: Significantly faster
- **Basis**: Controlled by choice
- **Production**: Feasible up to ~100 atoms
- **Benchmarks**: Standard molecules

## Limitations & Known Constraints
- **System size**: Limited to medium molecules
- **Periodic systems**: Not designed for solids
- **Features**: Focused on GW/BSE
- **GUI**: Command-line only
- **Platform**: Linux primarily

## Comparison with Other Codes
- **vs FHI-aims**: molgw more specialized GW
- **vs Turbomole**: molgw focused on MBPT
- **vs BerkeleyGW**: molgw for molecules, BerkeleyGW for solids
- **Unique strength**: Open-source molecular GW, Gaussian basis, self-consistent variants, educational clarity

## Application Areas

### Ionization Potentials:
- Molecular IPs
- Photoemission
- Benchmark accuracy
- Experimental comparison
- Systematic studies

### Electron Affinities:
- Molecular EAs
- Reduction potentials
- Electron attachment
- Accurate predictions

### Fundamental Gaps:
- HOMO-LUMO gaps
- Optical vs fundamental
- Gap renormalization
- Theoretical predictions

### Optical Excitations:
- BSE excitations
- Absorption spectra
- Exciton binding
- Molecular spectroscopy

### Method Development:
- GW algorithm testing
- Benchmark calculations
- Method comparison
- Educational use

## Best Practices

### Basis Set Selection:
- cc-pVTZ or larger
- Augmented for anions
- Systematic convergence
- Auxiliary basis

### Starting Point:
- HF or hybrid DFT
- PBE0 recommended
- Range-separated hybrids
- Consistent approach

### GW Convergence:
- Basis set convergence
- Auxiliary basis
- Self-consistency effects
- Numerical parameters

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Documentation available
- Developer support (Fabien Bruneval)
- User community
- Active development

## Educational Resources
- molgw documentation
- Tutorial examples
- GW method papers
- Benchmark datasets
- Source code (educational)

## Development
- Fabien Bruneval (CEA, France)
- Open development on GitHub
- Regular updates
- Community contributions
- Research-driven

## Research Applications
- Molecular ionization
- Optical properties
- GW benchmarks
- Method validation
- Spectroscopy theory

## Open-Source Advantage
- Free software
- Transparent algorithms
- Modifiable code
- Educational value
- Community development

## Molecular GW Expertise
- Optimized for molecules
- Gaussian basis efficiency
- RI acceleration
- Self-consistent methods
- Production quality

## Verification & Sources
**Primary sources**:
1. molgw website: https://www.molgw.org/
2. GitHub: https://github.com/molgw/molgw
3. F. Bruneval et al., Comp. Phys. Comm. (2016)
4. Documentation and tutorials

**Secondary sources**:
1. GW approximation literature
2. Molecular MBPT papers
3. Benchmark studies
4. Application publications

**Confidence**: VERIFIED - Established open-source code

**Verification status**: ✅ VERIFIED
- Website: ACCESSIBLE
- GitHub: Active repository
- Documentation: Comprehensive
- Community support: Developer and users
- Active development: Regular commits
- Specialized strength: Open-source molecular GW, Gaussian basis, self-consistent GW variants, ionization potentials, electron affinities, BSE excitations, RI acceleration, educational clarity
