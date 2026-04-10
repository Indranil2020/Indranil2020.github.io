# VASP-Raman

## Official Resources
- Source Repository: https://github.com/raman-sc/VASP
- License: Open source

## Overview
**VASP-Raman** is a Python program for evaluating off-resonance Raman activity using VASP as the backend. It computes Raman spectra by calculating the derivative of the macroscopic dielectric tensor (polarizability) with respect to phonon normal mode coordinates using finite displacements.

**Scientific domain**: Off-resonance Raman spectroscopy  
**Target user community**: VASP users needing Raman spectra from first-principles calculations

## Theoretical Methods
- Off-resonance Raman theory (Placzek approximation)
- Finite displacement method for dielectric tensor derivatives
- Density Functional Perturbation Theory (via VASP)
- Phonon normal modes from finite differences
- Polarizability derivative method (dP/dQ)

## Capabilities (CRITICAL)
- Off-resonance Raman spectra
- Raman activities per mode
- Polarization-resolved Raman
- Temperature-dependent Raman (Bose-Einstein factors)
- Automated VASP job management
- Support for both DFPT and finite-difference phonons

**Sources**: GitHub repository, Comput. Phys. Commun. 184, 1863 (2013)

## Key Strengths

### VASP Integration:
- Direct use of VASP dielectric calculations
- No VASP modifications needed
- Uses standard VASP capabilities (LEPSILON, LCALCEOL)
- Compatible with VASP parallelization

### Automated Workflow:
- Generates displaced structures
- Submits VASP calculations
- Collects and processes results
- Produces Raman spectra

### Well-Established:
- Widely used in VASP community
- Published methodology
- Multiple validation studies
- Active maintenance

## Inputs & Outputs
- **Input formats**:
  - VASP POSCAR, INCAR, KPOINTS, POTCAR
  - Phonon mode data (OUTCAR from IBRION=5/6 or DFPT)
  
- **Output data types**:
  - Raman activities per mode
  - Raman spectra (convoluted)
  - Polarization-resolved spectra
  - Mode assignments

## Interfaces & Ecosystem
- **VASP**: Required backend
- **Python**: Scripting
- **Phonopy**: Can use Phonopy for mode generation

## Performance Characteristics
- **Speed**: Limited by VASP dielectric calculations
- **Accuracy**: DFT-level
- **System size**: Limited by VASP
- **Memory**: VASP-dependent

## Computational Cost
- **Per mode**: One VASP dielectric calculation per displacement
- **3N modes**: 6N+1 VASP calculations needed
- **Typical**: Hours to days for moderate systems
- **vs ML methods**: Much slower than ramannoodle

## Limitations & Known Constraints
- **Off-resonance only**: No resonance Raman
- **VASP only**: No other DFT code support
- **Expensive**: Requires many VASP calculations
- **No ML acceleration**: Pure finite-difference
- **No excitonic effects**: DFT-level dielectric

## Comparison with Other Codes
- **vs ramannoodle**: VASP-Raman is exact finite-difference, ramannoodle is ML-accelerated
- **vs Phonopy-Spectroscopy**: VASP-Raman uses VASP, Phonopy-Spectroscopy uses Phonopy
- **vs QERaman**: VASP-Raman is off-resonance, QERaman is resonance
- **Unique strength**: Direct off-resonance Raman from VASP dielectric tensor, well-established

## Application Areas

### Crystalline Raman:
- Phonon mode identification
- Raman active mode characterization
- Polarization analysis
- Temperature-dependent Raman

### 2D Materials:
- Graphene and TMDs Raman
- Layer number dependence
- Strain characterization
- Defect signatures

### Phase Transitions:
- Raman across phase boundaries
- Soft mode identification
- Order parameter tracking
- Pressure-dependent Raman

## Best Practices

### Displacement Size:
- Test convergence with displacement amplitude
- Typical: 0.01-0.03 Å
- Too large: anharmonic contamination
- Too small: numerical noise

### VASP Settings:
- Use LEPSILON = .TRUE. for dielectric
- Adequate k-point density
- Well-converged ENCUT
- Consistent POTCAR

## Community and Support
- Open source on GitHub
- Widely used in VASP community
- Published methodology
- Active maintenance

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/raman-sc/VASP
2. A. Fonari and S. Stauffer, Comput. Phys. Commun. 184, 1863 (2013)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Maintained
- Specialized strength: Off-resonance Raman spectra from VASP dielectric tensor
