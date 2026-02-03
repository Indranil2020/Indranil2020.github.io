# phonon-sed (pSED)

## Official Resources
- Homepage: https://github.com/tyst3273/phonon-sed
- Source Repository: https://github.com/tyst3273/phonon-sed
- Documentation: Included manual.pdf
- License: Open Source

## Overview
phonon-sed (pSED) is a Python code to calculate the phonon spectral energy density (SED) from molecular dynamics simulations. It uses the normal-mode-decomposition technique to extract phonon frequencies, lifetimes, and linewidths from MD trajectories.

**Scientific domain**: Phonon spectral analysis, molecular dynamics post-processing  
**Target user community**: Researchers analyzing phonon properties from MD simulations

## Theoretical Methods
- Spectral energy density (SED) method
- Normal mode decomposition
- Fourier transform of velocities
- Phonon frequency extraction
- Linewidth and lifetime analysis
- Temperature-dependent phonon properties

## Capabilities (CRITICAL)
- Phonon spectral energy density calculation
- Frequency extraction from MD
- Phonon lifetime estimation
- Linewidth analysis
- Temperature-dependent properties
- Compatible with LAMMPS output
- Command-line interface (pSED)

## Key Strengths

### SED Method:
- Direct from MD trajectories
- Includes anharmonic effects
- Temperature-dependent
- No perturbation theory needed

### MD Compatibility:
- Works with LAMMPS
- General trajectory format
- Post-processing tool
- Flexible input

## Inputs & Outputs
- **Input formats**:
  - MD trajectory files
  - LAMMPS dump files
  - Structure information
  
- **Output data types**:
  - Spectral energy density
  - Phonon frequencies
  - Linewidths
  - Lifetimes

## Interfaces & Ecosystem
- **LAMMPS**: Primary MD code support
- **Python**: Pure Python implementation
- **Command-line**: pSED executable


## Advanced Features
- **Normal mode decomposition**: Phonon-resolved analysis
- **Fourier transform**: Velocity autocorrelation processing
- **Temperature dependence**: Finite-temperature phonon properties
- **Linewidth extraction**: Phonon lifetime analysis
- **LAMMPS integration**: Direct trajectory processing

## Performance Characteristics
- Python-based: Moderate speed
- FFT-limited: Scales with trajectory length
- Memory: Depends on system size and trajectory

## Computational Cost
- MD simulation: Dominant cost (external)
- SED calculation: Minutes to hours
- Depends on trajectory length and system size
- Long trajectories needed for frequency resolution

## Best Practices
- Use sufficiently long MD trajectories (>100 ps)
- Equilibrate system before production run
- Check frequency resolution vs trajectory length
- Validate against harmonic phonon calculations
- Compare with experimental spectroscopy when available

## Limitations & Known Constraints
- Requires long MD trajectories
- Computational cost for large systems
- Resolution depends on simulation time
- Classical MD limitations

## Application Areas
- Anharmonic phonon analysis
- Temperature-dependent phonon properties
- Thermal transport studies
- Phonon lifetime extraction
- MD validation against experiment

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/tyst3273/phonon-sed
2. Based on: T. Sun et al., J. Appl. Phys. 117, 135104 (2015)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Documentation: Manual included
