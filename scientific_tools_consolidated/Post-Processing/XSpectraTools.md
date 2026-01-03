# XSpectraTools

## Official Resources
- Homepage: https://github.com/calandra-group/XSpectraTools
- Documentation: README in repository
- Source Repository: https://github.com/calandra-group/XSpectraTools
- License: GNU General Public License v3.0

## Overview
XSpectraTools is a set of Python scripts designed to facilitate the calculation and analysis of X-ray Absorption Spectra (XAS) using the `xspectra.x` code of Quantum ESPRESSO. It automates the generation of input files, especially for core-hole supercell calculations, and provides utilities for post-processing and plotting the results.

**Scientific domain**: X-ray spectroscopy, XANES, workflow automation  
**Target user community**: Quantum ESPRESSO users, XAS researchers

## Theoretical Methods
- Core-hole supercell generation
- Spectral broadening (convolution)
- XANES calculation setup
- Configuration averaging (for disordered systems)

## Capabilities (CRITICAL)
- **Supercell Generation**: Creates supercells with a core-hole on a specific atom
- **Input Setup**: Generates `xspectra.in` and `pw.x` input files
- **Broadening**: Applies Lorentzian/Gaussian broadening to raw spectra
- **Alignment**: Aligns spectra based on core-level shifts
- **Plotting**: Simple visualization of spectra

**Sources**: XSpectraTools GitHub repository

## Inputs & Outputs
- **Input formats**: Structure file (CIF/POSCAR/QE input), parameters
- **Output data types**: QE input files, processed spectral data (xy format)

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Dedicated frontend for `xspectra.x`
- **Python**: Native implementation
- **ASE**: Uses ASE for structure manipulation

## Workflow and Usage
1. Define primitive structure.
2. Run script to generate supercell with core-hole: `python make_supercell.py ...`
3. Run SCF calculation (`pw.x`).
4. Run XSpectra calculation (`xspectra.x`).
5. Post-process: `python broaden.py spectrum.dat`

## Performance Characteristics
- Lightweight scripting tool
- Reduces manual error in complex input setup

## Application Areas
- XAS of defects and dopants
- High-throughput XANES calculations
- Teaching and tutorials for XSpectra

## Community and Support
- Open-source (GPL v3)
- Maintained by Calandra Group (Sorbonne Université)
- GitHub issues

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/calandra-group/XSpectraTools

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Documentation: AVAILABLE (README)
- Source: OPEN (GPL)
- Development: ACTIVE (Research group)
- Applications: XAS workflow, QE integration
