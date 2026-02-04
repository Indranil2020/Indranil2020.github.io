# ezSpectra

## Official Resources
- Homepage (Krylov): https://iopenshell.usc.edu/downloads/
- Homepage (Mosey): https://github.com/mosey-group/ezSpectra
- Publication: WIREs Comput. Mol. Sci. 12, e1546 (2022)
- License: Varies by package

## Overview
**Note**: Two distinct software packages share this name:
1. **ezSpectra Suite (Krylov Group)**: Toolkit for electronic spectroscopy (photoelectron, absorption) with vibronic structure via Franck-Condon factors.
2. **ezSpectra (Mosey Group)**: Python package for vibrational spectra from MD trajectories.

**Scientific domain**: Electronic/vibrational spectroscopy simulation
**Target user community**: Spectroscopists, photochemists

## Theoretical Methods
### Krylov Suite:
- Franck-Condon factors (ezFCF)
- Dyson orbitals for photoionization (ezDyson)
- Vibronic coupling

### Mosey Package:
- MD trajectory analysis
- Dipole/polarizability correlation
- IR, Raman, VCD spectra

## Capabilities (CRITICAL)
- **ezFCF**: Franck-Condon factors for vibronic spectra
- **ezDyson**: Photoionization cross-sections
- **Vibrational**: IR, Raman, VCD from MD
- **Electronic**: Photoelectron spectra
- **Absorption**: UV-Vis with vibronic structure

**Sources**: Krylov group website, GitHub repositories

## Key Strengths

### Krylov Suite:
- Vibronic structure
- Photoionization
- Q-Chem integration
- Well-documented

### Mosey Package:
- MD-based spectra
- Python native
- CP2K compatible
- Open source

## Inputs & Outputs
- **Input formats**: Q-Chem output (Krylov), MD trajectories (Mosey)
- **Output data types**: Simulated spectra, Franck-Condon factors

## Limitations & Known Constraints
- **Two packages**: Name ambiguity
- **Code-specific**: Different input requirements
- **Documentation**: Varies by package
- **Scope**: Each specialized for its domain

## Comparison with Other Tools
- **vs Gaussian FC**: ezFCF more specialized
- **vs Phonopy-Spectroscopy**: Different approaches (MD vs phonons)
- **Unique strength**: Vibronic structure (Krylov), MD spectra (Mosey)

## Application Areas
- Photoelectron spectroscopy
- UV-Vis absorption with vibronic structure
- IR/Raman from dynamics
- Photochemistry

## Best Practices
- Identify correct package for your needs
- Validate against experimental spectra
- Use appropriate theory level
- Check convergence parameters

## Community and Support
- Krylov group (USC)
- Mosey group (Queen's)
- Active development
- Published methodologies

## Verification & Sources
**Primary sources**:
1. Krylov: https://iopenshell.usc.edu/downloads/
2. Mosey: https://github.com/mosey-group/ezSpectra
3. WIREs Comput. Mol. Sci. 12, e1546 (2022)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Websites: ACCESSIBLE
- Source: AVAILABLE
- Method: Vibronic/MD spectroscopy
