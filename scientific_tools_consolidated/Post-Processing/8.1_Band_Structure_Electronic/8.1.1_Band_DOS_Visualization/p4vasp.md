# p4vasp

## Official Resources
- **Homepage**: http://www.p4vasp.at/
- **SourceForge**: https://sourceforge.net/projects/p4vasp/
- **License**: GPL v2

## Overview
p4vasp is a classic graphical visualization tool for VASP calculations. It provides an intuitive GUI for plotting band structures, density of states, and visualizing crystal structures from VASP output files. **Note**: This is legacy software; py4vasp is the recommended modern replacement.

**Scientific domain**: VASP visualization, electronic structure
**Target user community**: VASP users preferring GUI-based analysis

## Capabilities (CRITICAL)
- **Band Structure**: Electronic band structure plotting
- **DOS/PDOS**: Total and projected density of states
- **Structure Visualization**: Crystal structure display
- **Charge Density**: Electron density visualization
- **Local DOS**: Site-projected electronic structure
- **Animation**: MD trajectory visualization

## Key Strengths
- Graphical user interface
- Direct VASP file parsing
- Interactive visualization
- No programming required

## Inputs & Outputs
- **Input formats**: vasprun.xml, EIGENVAL, DOSCAR, PROCAR, CHGCAR
- **Output data types**: Plots, images (PNG, PS)

## Limitations & Known Constraints
- **Legacy software**: No longer actively developed
- **Python 2**: Requires older Python version
- **Superseded**: py4vasp is the modern replacement

## Verification & Sources
**Confidence**: VERIFIED - Legacy but functional

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Status: Legacy (use py4vasp for new projects)
