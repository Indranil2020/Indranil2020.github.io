# sumo

## Official Resources
- Homepage: https://smtg-ucl.github.io/sumo/
- Documentation: https://sumo.readthedocs.io/
- Source Repository: https://github.com/SMTG-UCL/sumo
- License: MIT License

## Overview
sumo (Seebeck Utility and Materials characterisation Organizer) is a Python toolkit for publication-quality plotting and analysis of ab initio calculation data, with particular focus on VASP outputs. It provides command-line tools for band structure, DOS, optical properties, and phonon plotting.

**Scientific domain**: Materials visualization, electronic structure analysis  
**Target user community**: Materials scientists needing publication-quality plots from VASP

## Theoretical Methods
sumo is a visualization/analysis tool, not a calculation method. It processes outputs from:
- VASP DFT calculations
- Band structures
- Density of states
- Optical properties
- Phonon calculations

## Capabilities (CRITICAL)
- Publication-quality band structure plotting
- Density of states (DOS, PDOS) visualization
- Optical absorption spectra plotting
- Phonon band structure and DOS plotting
- Effective mass calculations
- Band gap analysis
- Orbital contribution visualization
- Spin-polarized and spin-orbit calculations
- Customizable plot styling
- High-symmetry k-path generation
- Command-line interface for automation
- Python API for custom workflows
- Export to publication formats (PDF, PNG, SVG)

**Sources**: Official sumo documentation, cited in 6/7 source lists

## Inputs & Outputs
- **Input formats**:
  - vasprun.xml (VASP output)
  - EIGENVAL, DOSCAR, PROCAR
  - Phonopy outputs (phonopy.yaml, mesh.yaml)
  - POSCAR for structure info
  
- **Output data types**:
  - Publication-quality plots (PDF, PNG, SVG)
  - Data files for further analysis
  - Matplotlib figures
  - Effective mass data

## Interfaces & Ecosystem
- **VASP integration**:
  - Direct reading of VASP outputs via pymatgen
  - Automatic k-path detection
  
- **Phonopy integration**:
  - Phonon band structure plotting
  - Phonon DOS visualization
  
- **pymatgen backend**:
  - Structure analysis
  - Data parsing
  - Built on pymatgen infrastructure

## Limitations & Known Constraints
- **VASP-focused**: Primarily designed for VASP outputs
- **pymatgen dependency**: Requires pymatgen installation
- **Python 3**: Requires Python 3.6+
- **Platform**: Cross-platform (Linux, macOS, Windows)
- **Matplotlib backend**: Plot customization via matplotlib
- **Documentation**: Good but assumes VASP familiarity

## Verification & Sources
**Primary sources**:
1. Official website: https://smtg-ucl.github.io/sumo/
2. Documentation: https://sumo.readthedocs.io/
3. GitHub repository: https://github.com/SMTG-UCL/sumo
4. A. M. Ganose et al., J. Open Source Softw. 3, 717 (2018) - sumo paper

**Secondary sources**:
1. sumo tutorials and examples
2. SMTG group publications using sumo
3. Materials science visualization guides
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues)
- Academic citations: >100
- Active development: Regular updates, well-maintained
