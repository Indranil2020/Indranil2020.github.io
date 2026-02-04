# fold2Bloch-VASP

## Official Resources
- **GitHub**: https://github.com/rubel75/fold2Bloch-VASP
- **Publication**: Phys. Rev. B 90, 115202 (2014)
- **License**: MIT License

## Overview
fold2Bloch-VASP is a Fortran utility for unfolding supercell band structures from VASP calculations into primitive cell representation, computing effective band structures with Bloch spectral weights.

**Scientific domain**: Band unfolding, supercell calculations
**Target user community**: VASP users studying defects, alloys, interfaces

## Capabilities (CRITICAL)
- **Band Unfolding**: Supercell to primitive cell transformation
- **Effective Band Structure**: Primitive BZ representation
- **Spectral Weights**: Bloch character calculation
- **Large Supercells**: Efficient handling of complex systems

## Key Strengths
- Fortran implementation (fast)
- VASP WAVECAR support
- Command-line interface
- Well-documented

## Inputs & Outputs
- **Input formats**: VASP WAVECAR, POSCAR
- **Output data types**: Unfolded band structure data

## Application Areas
- Alloy band structures
- Defect calculations
- Interface/heterostructure studies

## Installation
```bash
git clone https://github.com/rubel75/fold2Bloch-VASP.git
cd fold2Bloch-VASP
make
```

## Limitations & Known Constraints
- **VASP-specific**: Only processes VASP WAVECAR
- **Compilation**: Requires Fortran compiler
- **Memory**: WAVECAR loading needs sufficient RAM

## Comparison with Other Tools
- **vs easyunfold**: fold2Bloch-VASP Fortran (fast), easyunfold Python (easier)
- **vs BandUP**: Different unfolding approaches
- **vs VaspBandUnfolding**: fold2Bloch-VASP compiled, VaspBandUnfolding Python
- **Unique strength**: Fast Fortran implementation, established methodology

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: Oleg Rubel (McMaster)
