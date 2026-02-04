# OrbVis

## Official Resources
- **GitHub**: https://github.com/staradutt/OrbVis
- **PyPI**: https://pypi.org/project/orbvis/
- **License**: MIT License

## Overview
OrbVis is a lightweight Python package for plotting orbital-projected band structures and density of states from VASP output files (PROCAR & DOSCAR).

**Scientific domain**: VASP visualization, orbital projections
**Target user community**: VASP users needing quick orbital analysis

## Capabilities (CRITICAL)
- **Orbital-Projected Bands**: s, p, d, f orbital contributions
- **Atom-Projected Bands**: Element-specific band character
- **DOS/PDOS**: Orbital-resolved density of states
- **Fatbands**: Band width proportional to orbital character

## Key Strengths
- Lightweight and fast
- Easy-to-use API
- Publication-ready figures
- Customizable color schemes

## Inputs & Outputs
- **Input formats**: PROCAR, DOSCAR, EIGENVAL
- **Output data types**: Matplotlib figures

## Installation
```bash
pip install orbvis
```

## Limitations & Known Constraints
- **VASP-specific**: Only works with VASP output files
- **PROCAR required**: Needs PROCAR file for orbital projections
- **Limited features**: Focused tool, less comprehensive than pyprocar

## Comparison with Other Tools
- **vs pyprocar**: OrbVis simpler, pyprocar more comprehensive
- **vs vaspvis**: Similar capabilities, different API styles
- **vs sumo**: OrbVis VASP-only, sumo more general
- **Unique strength**: Lightweight, easy-to-use orbital projection plotting

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, MIT)
- Target Code: VASP
