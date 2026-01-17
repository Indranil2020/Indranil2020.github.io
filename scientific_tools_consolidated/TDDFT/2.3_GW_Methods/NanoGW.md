# NanoGW

## Official Resources
- Homepage: https://codebase.helmholtz.cloud/nanogw/nanogw (or LBL hosted)
- Documentation: Bundled with code
- Source Repository: Available through Berkeley Lab / NERSC
- License: Open Source (BSD-like)

## Overview
NanoGW is an open-source software package for linear-response TDDFT, GW, and Bethe-Salpeter equation (BSE) calculations using a real-space grid. Designed specifically for confined systems such as molecules and nanoclusters, it performs full-frequency GW calculations with optional LDA vertex corrections.

**Scientific domain**: Molecules, nanoclusters, quantum dots, excited states  
**Target user community**: Researchers studying finite systems requiring accurate quasiparticle and optical properties

## Theoretical Methods
- Full-frequency GW approximation
- Bethe-Salpeter equation (BSE)
- Linear-response TDDFT
- LDA vertex function corrections
- Real-space grid representation
- Quasiparticle self-energy
- Optical excitations

## Capabilities (CRITICAL)
- Full-frequency GW calculations
- G0W0 quasiparticle energies
- Bethe-Salpeter equation for optical spectra
- Linear-response TDDFT
- LDA vertex corrections
- Real-space grid discretization
- Molecules and nanoclusters (<30 atoms optimized)
- Ionization potentials and electron affinities
- Optical absorption spectra
- Convergence studies

**Sources**: LBL NanoGW pages, PARSEC integration documentation

## Key Strengths

### Real-Space Grid:
- Systematic convergence
- No basis set artifacts
- Natural for confined systems
- Flexible boundary conditions

### Full-Frequency GW:
- No plasmon-pole approximation
- Accurate self-energy
- Frequency-dependent screening
- Precise spectral functions

### Finite System Focus:
- Optimized for molecules
- Nanoclusters up to ~30 atoms
- Quantum dots
- No periodic image artifacts

### BSE Capability:
- Optical excitations
- Excitonic effects
- Neutral excitations
- Absorption spectra

## Inputs & Outputs
- **Input formats**:
  - PARSEC wavefunctions and energies
  - PARATEC plane-wave (converted)
  - Real-space grid data
  
- **Output data types**:
  - Quasiparticle energies
  - Self-energy matrices
  - Optical absorption spectra
  - BSE excitation energies
  - Oscillator strengths

## Interfaces & Ecosystem
- **DFT Integration**:
  - PARSEC (real-space DFT code)
  - PARATEC (plane-wave, with conversion)
  - Kohn-Sham wavefunctions required
  
- **Post-processing**:
  - Spectrum analysis tools
  - Convergence utilities

## Advanced Features

### LDA Vertex Corrections:
- Beyond standard GW
- Improved electron-electron description
- Optional enhancement

### Crystalline Systems:
- Can handle crystals (less tested)
- Periodic boundary support
- Primary focus remains finite systems

## Performance Characteristics
- **Speed**: Grid-based efficiency
- **Accuracy**: Full-frequency precision
- **System size**: Optimized for <30 atoms
- **Memory**: Grid point dependent

## Computational Cost
- **GW**: Full frequency integration
- **BSE**: Two-particle calculations
- **Scaling**: Depends on grid density
- **Typical**: Hours for small molecules

## Limitations & Known Constraints
- **System size**: Best for small molecules/clusters
- **Spin-orbit**: Not supported
- **Crystalline**: Less thoroughly tested
- **DFT input**: Requires PARSEC or PARATEC

## Comparison with Other Codes
- **vs BerkeleyGW**: NanoGW real-space, BerkeleyGW plane-wave
- **vs molgw**: Both molecular, different basis approaches
- **vs Yambo**: NanoGW finite systems, Yambo periodic
- **Unique strength**: Real-space grid for molecules, full-frequency GW

## Application Areas

### Molecular Spectroscopy:
- IP/EA calculations
- Optical absorption
- Electronic excitations
- Photoemission

### Nanoclusters:
- Quantum dot excitations
- Cluster electronic structure
- Size-dependent properties

### Quantum Chemistry:
- Beyond-DFT corrections
- Accurate HOMO-LUMO gaps
- Excitonic binding

## Best Practices

### Grid Convergence:
- Systematic grid refinement
- Check energy convergence
- Balance accuracy vs cost

### System Selection:
- Best for <30 atoms
- Finite systems preferred
- Use PARSEC for input generation

## Community and Support
- Berkeley Lab development
- Academic user community
- PARSEC integration documented
- Open-source availability

## Verification & Sources
**Primary sources**:
1. LBL NanoGW pages
2. PARSEC code documentation
3. Published applications

**Confidence**: VERIFIED
- Code availability: Through LBL/NERSC
- Documentation: Available
- Active use: Academic community
