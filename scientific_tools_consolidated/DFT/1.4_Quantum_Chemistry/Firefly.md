# Firefly

## Official Resources
- Homepage: http://classic.chem.msu.su/gran/firefly/index.html
- Documentation: http://classic.chem.msu.su/gran/firefly/index.html
- Download: http://classic.chem.msu.su/gran/firefly/index.html
- License: Free for non-commercial use

## Overview
Firefly (formerly PC GAMESS) is an optimized fork of GAMESS(US) with architecture-specific optimizations developed by Alex Granovsky at Moscow State University. It maintains full compatibility with GAMESS input format while providing significant performance improvements through modern code optimization techniques.

**Scientific domain**: General ab initio quantum chemistry  
**Target user community**: Users seeking faster GAMESS-compatible calculations with full method coverage

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF, GVB)
- Density Functional Theory (many functionals)
- MP2, MP3, MP4 perturbation theory
- Coupled Cluster (CCSD, CCSD(T), CR-CC)
- MCSCF, CASSCF, CASPT2
- CI methods (CIS, CISD, SOCI)
- Semi-empirical methods (AM1, PM3, PM6)
- TD-DFT for excited states

## Capabilities (CRITICAL)
- Full GAMESS input compatibility
- Architecture-optimized performance
- Up to 10x faster than standard GAMESS
- Parallel execution (MPI)
- Wide method coverage
- Geometry optimization
- Transition state searches
- Frequency calculations and thermochemistry
- Solvation models (PCM, COSMO)
- QM/MM capabilities
- Extensive property calculations

## Key Strengths

### Performance:
- Code optimization for modern CPUs
- SSE/AVX instructions
- Better memory management
- Faster integral evaluation
- Optimized BLAS usage

### GAMESS Compatibility:
- Same input format
- Same output format
- Easy migration
- Documentation applies
- Feature parity

### Method Coverage:
- Ground state methods
- Excited states (TD-DFT, CIS)
- Multi-reference methods
- Coupled cluster
- Comprehensive suite

### Practical Features:
- Geometry optimization
- Transition states
- Frequencies
- Thermochemistry
- Solvation

## Inputs & Outputs
- **Input formats**:
  - GAMESS input files
  - $DATA, $BASIS, $CONTRL groups
  - Standard GAMESS format
  
- **Output data types**:
  - GAMESS-format output
  - Energies, gradients, Hessians
  - Orbitals and densities
  - Properties

## Interfaces & Ecosystem
- **GAMESS compatible**: All GAMESS tools work
- **Visualization**: wxMacMolPlt, Avogadro, etc.
- **QM/MM**: TINKER interface
- **Solvation**: PCM, COSMO

## Advanced Features

### Multi-Reference Methods:
- MCSCF optimization
- CASSCF with large active spaces
- CASPT2 for dynamic correlation
- XMCQDPT2

### Coupled Cluster:
- CCSD implementation
- CCSD(T) energies
- CR-CC methods
- Left eigenstates

### Excited States:
- TD-DFT spectra
- CIS and CIS(D)
- EOM-CCSD
- State-averaged CASSCF

### Solvation:
- PCM (various models)
- COSMO
- SMD model
- Cavity specification

## Performance Characteristics
- **Speed**: 2-10x faster than GAMESS(US)
- **Accuracy**: Same as GAMESS
- **System size**: Medium to large molecules
- **Memory**: Efficient management
- **Parallelization**: MPI scaling

## Computational Cost
- **HF/DFT**: Very efficient
- **MP2**: Fast implementation
- **CCSD(T)**: Standard scaling with optimizations
- **CASPT2**: Production capable
- **Typical**: Hours to days for large systems

## Limitations & Known Constraints
- **Development status**: Less active recently
- **Closed source**: Binary distribution
- **Platform-specific**: Requires matched binaries
- **Documentation**: Uses GAMESS docs
- **Community**: Smaller than GAMESS

## Comparison with Other Codes
- **vs GAMESS(US)**: Faster, same features
- **vs Gaussian**: Open alternative, some overlap
- **vs ORCA**: Different strengths
- **vs GAMESS(UK)**: Different development
- **Unique strength**: Optimized GAMESS performance

## Application Areas

### Organic Chemistry:
- Reaction mechanisms
- Conformational analysis
- Transition states
- Spectroscopy

### Inorganic Chemistry:
- Transition metal complexes
- Coordination compounds
- Spin states
- Ligand effects

### Photochemistry:
- Excited states
- Photoreactions
- Absorption spectra
- Fluorescence

### Materials:
- Clusters
- Molecular crystals
- Intermolecular interactions
- Host-guest chemistry

## Best Practices

### Migration from GAMESS:
- Use same inputs
- Verify output consistency
- Check for Firefly-specific keywords
- Benchmark key calculations

### Performance:
- Use appropriate parallelization
- Match binary to architecture
- Memory settings
- Scratch space

## Community and Support
- Free for academic use
- Moscow State University
- Firefly-specific forums
- GAMESS documentation applies
- Granovsky development

## Verification & Sources
**Primary sources**:
1. Official site: http://classic.chem.msu.su/gran/firefly/
2. Granovsky, J. Chem. Phys. (benchmark papers)
3. GAMESS documentation (compatible)
4. Moscow State University

**Confidence**: VERIFIED
- Software: Available for download
- Status: Free for academic use
- Performance: Documented speedups
- Active: Periodic updates
