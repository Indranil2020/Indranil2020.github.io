# QEPlotter

## Official Resources
- **GitHub**: https://github.com/romerogroup/QEPlotter
- **Compatible with**: Quantum ESPRESSO 5.4.0+
- **License**: MIT License

## Overview
QEPlotter is a Quantum ESPRESSO band structure, DOS/PDOS, and fatband plotting toolkit with automated band gap detection. It provides both command-line and GUI interfaces for easy visualization of QE calculation results.

**Scientific domain**: Quantum ESPRESSO post-processing, electronic structure visualization
**Target user community**: Quantum ESPRESSO users needing quick visualization

## Capabilities (CRITICAL)
- **Band Structure**: Electronic band plotting with high-symmetry paths
- **DOS/PDOS**: Density of states visualization
- **Fatbands**: Orbital-projected band structures
- **Band Gap Detection**: Automated direct/indirect gap identification
- **Spin-Polarized**: Support for magnetic calculations

## Key Strengths
- GUI interface for easy use
- Automated band gap detection
- Publication-ready figures
- Quantum ESPRESSO native support

## Inputs & Outputs
- **Input formats**: QE output files (bands.dat, dos.dat, pdos files)
- **Output data types**: Matplotlib figures, data files

## Installation
```bash
git clone https://github.com/romerogroup/QEPlotter.git
cd QEPlotter
pip install -e .
```

## Limitations & Known Constraints
- **QE-specific**: Only works with Quantum ESPRESSO outputs
- **Version compatibility**: Best with QE 5.4.0+
- **Dependencies**: Requires matplotlib, numpy

## Comparison with Other Tools
- **vs postqe**: QEPlotter GUI-focused, postqe Python API
- **vs abipy**: QEPlotter for QE, abipy for ABINIT
- **vs sumo**: QEPlotter QE-native, sumo multi-code via pymatgen
- **Unique strength**: GUI interface, automated band gap detection

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Target Code: Quantum ESPRESSO
