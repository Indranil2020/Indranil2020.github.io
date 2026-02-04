# yambopy

## Official Resources
- Homepage: https://www.yambo-code.org/
- Documentation: https://www.yambo-code.org/wiki/index.php?title=Yambopy
- Source Repository: https://github.com/yambo-code/yambopy
- License: GNU General Public License v3.0

## Overview
yambopy is a python package to read, write, and analyze Yambo input and output files. It facilitates the workflow of GW and Bethe-Salpeter Equation (BSE) calculations using the Yambo code. yambopy provides tools for convergence testing, excitonic wavefunction analysis, band structure interpolation, and handling of large-scale GW-BSE data.

**Scientific domain**: Many-body perturbation theory (GW, BSE), workflow automation, post-processing  
**Target user community**: Yambo users, researchers in excited states and optics

## Capabilities (CRITICAL)
- **Input Generation**: Pythonic creation of Yambo input files (yambo.in)
- **Output Parsing**: Reading of netCDF and ASCII output (QP energies, spectra)
- **Convergence**: Automated convergence tests for GW and BSE parameters
- **Analysis**: Plotting of quasiparticle band structures, absorption spectra
- **Excitons**: Visualization of excitonic wavefunctions (weights, real-space plots)
- **Interpolation**: Fourier interpolation of GW bands
- **Database**: Tools to manage Yambo databases (ndb files)

**Sources**: Yambopy documentation, Yambo website

## Inputs & Outputs
- **Input formats**: Yambo input structures (Python dictionaries), Quantum ESPRESSO save folders
- **Output data types**: Matplotlib plots, JSON data, netCDF databases

## Interfaces & Ecosystem
- **Yambo**: The primary backend code
- **Quantum ESPRESSO**: Interface for generating input from QE
- **Pymatgen**: Integration for structure handling
- **Python**: Full scripting environment

## Workflow and Usage
1. Perform DFT calculation (Quantum ESPRESSO).
2. Initialize Yambo databases (`p2y`).
3. Use yambopy script:
   ```python
   from yambopy import *
   y = YamboIn.from_runlevel(' -o b -k sex -y d', folder='.')
   y['FFTGvecs'] = [10, 'Ry']
   y.write('yambo_run.in')
   ```
4. Run Yambo.
5. Analyze results:
   ```python
   yo = YamboOut('.')
   yo.pack.get_energies()
   ```

## Performance Characteristics
- Python wrapper; performance depends on Yambo execution
- Efficient parsing of netCDF files
- Automates tedious convergence studies

## Application Areas
- Photovoltaics (quasiparticle gaps)
- 2D materials (excitons)
- Optical properties of solids
- High-throughput GW calculations

## Community and Support
- Part of Yambo-code project
- Active development on GitHub
- Tutorials available on Yambo wiki

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.yambo-code.org/
2. GitHub: https://github.com/yambo-code/yambopy
3. Documentation: https://www.yambo-code.org/wiki/index.php?title=Yambopy

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Yambo Team)
- Applications: GW/BSE workflow, analysis, plotting
