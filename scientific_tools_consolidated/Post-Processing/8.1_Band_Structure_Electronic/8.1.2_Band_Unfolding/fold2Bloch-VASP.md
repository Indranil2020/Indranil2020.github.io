# fold2Bloch-VASP

## Overview
**fold2Bloch-VASP** is a Fortran utility designed to unfold the band structure of a supercell obtained with VASP and compute an effective band structure in a primitive representation.

## Official Resources
- **GitHub**: https://github.com/rubel75/fold2Bloch-VASP
- **Publication**: Phys. Rev. B 90, 115202 (2014)

## Capabilities
- **Band Unfolding**: Supercell to primitive cell
- **Effective Band Structure**: Primitive BZ representation
- **Large-scale Calculations**: Handle complex supercells
- **Alloy/Defect Studies**: Interpret disordered systems

## Key Features
- Fortran implementation (fast)
- VASP WAVECAR support
- Bloch spectral weights
- Command-line interface

## Applications
- Alloy band structures
- Defect calculations
- Interface/heterostructure studies
- Disordered systems

## Inputs & Outputs
- **Inputs**: VASP WAVECAR, POSCAR
- **Outputs**: Unfolded band structure data

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
- **Developer**: Oleg Rubel
