# effectivemass (effmass)

## Official Resources
- Homepage: https://github.com/lucydot/effmass
- Documentation: https://effmass.readthedocs.io/
- Source Repository: https://github.com/lucydot/effmass
- License: MIT License

## Overview
effmass is a Python package for calculating the effective mass of charge carriers (electrons and holes) from electronic band structures. It supports both parabolic and non-parabolic definitions of effective mass and provides tools for selecting band segments extrema. effmass is designed to work with outputs from major DFT codes like VASP, FHI-aims, CASTEP, and Quantum ESPRESSO.

**Scientific domain**: Semiconductor physics, transport properties, electronic structure  
**Target user community**: Materials scientists, device physicists, PV researchers

## Theoretical Methods
- Parabolic effective mass (curvature at band extrema)
- Non-parabolic kane dispersion fitting
- Finite difference method (3-point, 5-point)
- Least-squares fitting
- Density-of-states (DOS) effective mass
- Conductivity effective mass
- Optical effective mass

## Capabilities (CRITICAL)
- Automatic identification of CBM and VBM
- Calculation of effective mass tensors
- Support for multiple definitions of effective mass (transport, optical, curvature)
- Handling of multiple bands and degeneracy
- Integration with ASE and vasppy
- Command-line interface and Python API
- Visualization of band segments and fits

**Sources**: effmass documentation, J. Open Source Softw. 3, 797 (2018)

## Inputs & Outputs
- **Input formats**: VASP (vasprun.xml, OUTCAR/PROCAR), CASTEP, FHI-aims, QE, cube files
- **Output data types**: Effective mass values (tensor/scalar), plots of fits, JSON summary

## Interfaces & Ecosystem
- **VASP**: Primary target via `vasprun.xml`
- **ASE**: Compatible with ASE Atoms objects
- **Cli**: `effmass` command line tool
- **Python**: Library for custom analysis

## Workflow and Usage
1. Perform DFT band structure calculation (dense k-mesh near extrema).
2. Run effmass: `effmass-cli` (interactive) or script.
3. Select extrema and directions.
4. Choose effective mass definition (e.g., curvature).
5. Output effective masses and plots.

## Performance Characteristics
- Very fast post-processing
- Depends on file I/O speed for large vasprun.xml

## Application Areas
- Photovoltaics (carrier mobility)
- Thermoelectrics
- Transparent conducting oxides
- Semiconductor device modeling

## Community and Support
- Open-source (MIT)
- GitHub repository
- Developed by Lucy Whalley (Imperial College London / Northumbria University)

## Verification & Sources
**Primary sources**:
1. Homepage: https://github.com/lucydot/effmass
2. Documentation: https://effmass.readthedocs.io/
3. Publication: L. D. Whalley, J. Open Source Softw. 3, 797 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Whalley)
- Applications: Effective mass, carrier transport, band fitting
