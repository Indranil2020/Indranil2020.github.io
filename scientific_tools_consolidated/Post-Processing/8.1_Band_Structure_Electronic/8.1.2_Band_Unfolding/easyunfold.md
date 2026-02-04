# easyunfold

## Overview
**easyunfold** is a Python package for band structure unfolding, making it easy to obtain effective band structures from supercell calculations. It properly accounts for symmetry-breaking by sampling necessary additional k-points.

## Official Resources
- **GitHub**: https://github.com/SMTG-Bham/easyunfold
- **Documentation**: https://smtg-bham.github.io/easyunfold/
- **JOSS Paper**: Published in Journal of Open Source Software

## Capabilities
- **K-point Generation**: Generate k-points for supercell sampling
- **Band Unfolding**: Extract effective band structure from supercells
- **Symmetry Handling**: Proper treatment of symmetry-breaking
- **Multi-code Support**: VASP, CASTEP, Quantum ESPRESSO

## Key Features
- Command-line interface
- Symmetry-aware k-point sampling
- Spectral weight calculation
- Publication-ready plots
- Based on PyVaspwfc

## Supported Codes
- VASP
- CASTEP
- Quantum ESPRESSO

## Installation
```bash
pip install easyunfold
```

## Usage
```bash
easyunfold generate  # Generate k-points
easyunfold unfold    # Perform unfolding
```

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
- **Developer**: SMTG Birmingham
