# vaspy

## Official Resources
- **GitHub**: https://github.com/arafune/vaspy
- **PyPI**: https://pypi.org/project/vaspy/
- **License**: MIT License

## Overview
vaspy is a VASP pre- and post-processing system written in Python providing utilities for manipulating VASP input/output files and analyzing calculation results including WAVECAR parsing.

**Scientific domain**: VASP file manipulation, wavefunction analysis
**Target user community**: VASP users needing Python-based file handling

## Capabilities (CRITICAL)
- **File Parsing**: Read/write VASP files (POSCAR, CONTCAR, OUTCAR)
- **WAVECAR Analysis**: Wavefunction manipulation and extraction
- **Band Structure**: Electronic structure analysis
- **DOS**: Density of states processing
- **Structure Manipulation**: Atomic structure operations

## Key Strengths
- Comprehensive VASP file support
- WAVECAR parsing (including gamma-only)
- Python scripting interface
- Structure manipulation utilities

## Inputs & Outputs
- **Input formats**: POSCAR, WAVECAR, OUTCAR, EIGENVAL, DOSCAR
- **Output data types**: Processed data, modified input files

## Installation
```bash
pip install vaspy
```

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Note: Different from vasppy (separate package)
