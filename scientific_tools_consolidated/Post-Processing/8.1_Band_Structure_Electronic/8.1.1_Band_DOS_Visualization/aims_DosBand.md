# aims_DosBand

## Official Resources
- **FHI-aims**: https://fhi-aims.org/
- **Target Code**: FHI-aims

## Overview
aims_DosBand is a plotting tool for band structures and density of states from FHI-aims calculations, providing publication-quality figures with orbital projections.

**Scientific domain**: FHI-aims post-processing, electronic structure visualization
**Target user community**: FHI-aims users

## Capabilities (CRITICAL)
- **Band Structure**: Electronic band plotting along k-paths
- **DOS/PDOS**: Total and projected density of states
- **Fatbands**: Orbital-projected band structures
- **Spin-Polarized**: Support for magnetic calculations

## Key Strengths
- Native FHI-aims support
- Publication-quality figures
- Orbital projections
- Spin-resolved plotting

## Inputs & Outputs
- **Input formats**: FHI-aims output files (band*.out, KS_DOS*.dat)
- **Output data types**: Matplotlib figures

## Limitations & Known Constraints
- **FHI-aims specific**: Only processes FHI-aims output
- **File formats**: Requires specific FHI-aims output files
- **FHI-aims license**: FHI-aims code requires academic/commercial license

## Comparison with Other Tools
- **vs aimstools**: Similar scope, different implementations
- **vs sumo**: aims_DosBand FHI-aims native, sumo multi-code
- **Unique strength**: Native FHI-aims band/DOS visualization

## Verification & Sources
**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Target Code: FHI-aims
