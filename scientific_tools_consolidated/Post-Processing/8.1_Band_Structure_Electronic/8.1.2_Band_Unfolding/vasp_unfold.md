# vasp_unfold

## Official Resources
- **GitHub**: https://github.com/tomkeus/vasp_unfold
- **License**: MIT License

## Overview
vasp_unfold is a Python tool for unfolding VASP band structures using PROCAR files with phase information (LORBIT=12), providing spectral weight calculation for supercell systems.

**Scientific domain**: Band unfolding, VASP post-processing
**Target user community**: VASP users with PROCAR-based workflows

## Capabilities (CRITICAL)
- **PROCAR-based Unfolding**: Uses projected character data
- **Phase Information**: Requires LORBIT=12
- **Fractional Translations**: Handles supercell translations
- **Spectral Weights**: Calculate unfolding weights

## Key Strengths
- PROCAR file parsing
- Phase-aware unfolding
- Command-line interface

## Requirements
- VASP calculation with LORBIT=12
- PROCAR file with phase information

## Usage
```bash
vasp_unfold PROCAR
```

## Installation
```bash
git clone https://github.com/tomkeus/vasp_unfold.git
cd vasp_unfold
pip install -e .
```

## Limitations & Known Constraints
- **LORBIT=12 required**: Needs VASP phase information in PROCAR
- **PROCAR-based**: Different approach from WAVECAR-based tools
- **VASP-specific**: Only works with VASP output

## Comparison with Other Tools
- **vs easyunfold**: vasp_unfold PROCAR-based, easyunfold WAVECAR-based
- **vs fold2Bloch**: Different input requirements (PROCAR vs WAVECAR)
- **Unique strength**: PROCAR-based unfolding with phase information

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Target Code: VASP
