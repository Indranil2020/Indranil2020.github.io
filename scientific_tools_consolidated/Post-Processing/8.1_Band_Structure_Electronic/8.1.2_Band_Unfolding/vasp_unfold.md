# vasp_unfold

## Overview
**vasp_unfold** is a Python tool for unfolding VASP band structures using PROCAR files. It requires LORBIT=12 to include phase information.

## Official Resources
- **GitHub**: https://github.com/tomkeus/vasp_unfold

## Capabilities
- **PROCAR-based Unfolding**: Uses projected character data
- **Phase Information**: Requires LORBIT=12
- **Fractional Translations**: Handles supercell translations
- **Spectral Weights**: Calculate unfolding weights

## Key Features
- PROCAR file parsing
- Phase-aware unfolding
- Supercell translation handling
- Command-line interface

## Requirements
- VASP calculation with LORBIT=12
- PROCAR file with phase information

## Usage
After VASP bandstructure calculation:
```bash
vasp_unfold PROCAR
```

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
- **Target Code**: VASP
