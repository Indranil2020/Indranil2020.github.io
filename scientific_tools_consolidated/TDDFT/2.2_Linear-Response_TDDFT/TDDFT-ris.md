# TDDFT-ris

## Official Resources
- Homepage: https://github.com/John-zzh/pyscf_TDDFT_ris
- Source Repository: https://github.com/John-zzh/pyscf_TDDFT_ris
- Documentation: README with theory and usage
- License: Open Source

## Overview
TDDFT-ris is a high-performance Python implementation of a semiempirical Linear-Response TDDFT method that achieves ~300x speedup over traditional ab initio TDDFT while maintaining accuracy within 0.06 eV for excitation energies of organic molecules. The method uses Resolution-of-the-Identity (RI) approximation with minimal auxiliary basis (single s-type orbital per atom) and disables the XC kernel, providing an excellent balance of speed and accuracy for UV-Vis spectroscopy.

**Scientific domain**: Molecular excited states, UV-Vis absorption, large organic molecules  
**Target user community**: Researchers needing fast, accurate TDDFT for large molecular systems

## Theoretical Methods
- Linear-Response Time-Dependent DFT (LR-TDDFT)
- Resolution-of-the-Identity (RI) approximation
- Minimal auxiliary basis (ris/risp/rispd)
- Semi-empirical atomic radii parameterization
- Hybrid XC functional compatibility
- MOKIT integral interface
- PySCF backend

## Capabilities
- Extremely fast TDDFT calculations (~300x speedup)
- Excitation energy calculations
- UV-Vis absorption spectra
- Excited state analysis
- Multiple auxiliary basis levels (s, sp, spd)
- Hybrid functional support
- Large molecule capability
- Command-line interface
- PySCF object interface

## Key Strengths

### Exceptional Performance:
- ~300x faster than ab initio TDDFT
- 0.06 eV average deviation
- 4x more accurate than sTDDFT
- Scales to large molecules

### Flexible Auxiliary Basis:
- ris: single s-function per atom (fastest)
- risp: s+p functions (default, more accurate)
- rispd: s+p+d functions (most accurate)
- Atomic radii-based exponents

### Production Integration:
- Built into TURBOMOLE 7.7dev
- Built into AMESP v2.1dev
- ORCA 6.0 compound script support
- PySCF native integration

### Easy Workflow:
- .fch file input support
- Command-line interface
- Python API
- Automatic spectrum plotting

## Inputs & Outputs
- **Input formats**:
  - Formatted checkpoint files (.fch)
  - PySCF mean-field objects
  - MOKIT-compatible files
  
- **Output data types**:
  - Excitation energies
  - Oscillator strengths
  - UV-Vis spectra (plotted)
  - Excited state analysis
  - Natural transition orbitals

## Interfaces & Ecosystem
- **Dependencies**:
  - PySCF (required)
  - MOKIT (required)
  - NumPy, SciPy

- **Compatible software**:
  - TURBOMOLE
  - AMESP
  - ORCA 6.0
  - Gaussian (via .fch files)

## Command-Line Usage
```bash
# Basic calculation
python -m TDDFT_ris input.fch --nstates 10

# Plot spectrum
python -m TDDFT_ris input.fch --nstates 20 --plot
```

## Performance Characteristics
- **Speed**: ~300x faster than ab initio TDDFT
- **Accuracy**: 0.06 eV average deviation
- **System size**: Large organic molecules
- **Memory**: Reduced compared to full TDDFT
- **Basis sets**: All PySCF-supported basis sets

## Limitations & Known Constraints
- **Accuracy**: Slightly lower than full ab initio TDDFT
- **Core excitations**: Not optimized for core-level
- **Charge-transfer**: May need careful validation
- **Metals**: Parameterized for organic molecules
- **Dependencies**: Requires MOKIT installation

## Comparison with Other Codes
- **vs ab initio TDDFT**: 300x faster, 0.06 eV deviation
- **vs sTDDFT (Grimme)**: 4x more accurate (0.06 vs 0.24 eV)
- **vs sTDA**: Similar speed, better accuracy
- **Unique strength**: Best speed/accuracy trade-off for organics

## Application Areas

### UV-Vis Spectroscopy:
- Absorption spectra of large molecules
- Chromophore design
- Photophysical properties
- Dye molecules

### High-Throughput Screening:
- Virtual screening of chromophores
- Materials discovery
- Large dataset calculations

### Photochemistry:
- Initial excited state characterization
- Photosensitizer evaluation
- OLED material screening

## Best Practices
- Use risp (default) for balance of speed/accuracy
- Validate against full TDDFT for new molecule classes
- Use .fch files for Gaussian interoperability
- Enable spectrum plotting for visualization

## Community and Support
- Open-source on GitHub
- 2 contributors
- 1 release (v2.0)
- Active development
- Python 100%

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/John-zzh/pyscf_TDDFT_ris
2. Built into TURBOMOLE 7.7dev
3. Built into AMESP v2.1dev
4. ORCA 6.0 compound script support

**Confidence**: VERIFIED - Production software integration confirms validity

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Documentation: Comprehensive README
- Production use: TURBOMOLE, AMESP, ORCA
- Language: Python 100%
