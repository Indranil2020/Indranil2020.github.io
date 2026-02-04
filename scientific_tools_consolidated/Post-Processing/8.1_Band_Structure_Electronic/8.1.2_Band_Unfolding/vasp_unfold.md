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

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Target Code: VASP
