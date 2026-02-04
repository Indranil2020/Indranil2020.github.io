# VaspBandUnfolding

## Official Resources
- **GitHub**: https://github.com/QijingZheng/VaspBandUnfolding
- **Tutorial**: http://staff.ustc.edu.cn/~zqj/posts/Band-unfolding-tutorial/
- **License**: MIT License

## Overview
VaspBandUnfolding is a Python toolkit for VASP band unfolding from supercell calculations, providing WAVECAR parsing and spectral weight calculation with well-documented tutorials.

**Scientific domain**: Band unfolding, VASP post-processing
**Target user community**: VASP users studying supercell systems

## Capabilities (CRITICAL)
- **Band Unfolding**: Supercell band structure unfolding
- **WAVECAR Parsing**: Read VASP wavefunction files
- **Spectral Functions**: Calculate spectral weights
- **K-path Generation**: Create unfolding k-paths

## Key Strengths
- WAVECAR reading utilities
- Spectral weight calculation
- Well-documented tutorials
- Based on Popescu & Zunger methodology

## Inputs & Outputs
- **Input formats**: VASP WAVECAR, POSCAR
- **Output data types**: Unfolded band structures, spectral functions

## Installation
```bash
git clone https://github.com/QijingZheng/VaspBandUnfolding.git
cd VaspBandUnfolding
pip install -e .
```

## Limitations & Known Constraints
- **VASP-specific**: Only processes VASP WAVECAR files
- **Memory**: Large WAVECAR files need substantial RAM
- **Documentation**: Tutorial-based, less formal documentation

## Comparison with Other Tools
- **vs easyunfold**: VaspBandUnfolding more manual, easyunfold more automated
- **vs BandUP**: Different implementation approaches
- **vs fold2Bloch**: Both VASP unfolding, different interfaces
- **Unique strength**: Well-documented tutorials, WAVECAR parsing utilities

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Developer: Qijing Zheng (USTC)
