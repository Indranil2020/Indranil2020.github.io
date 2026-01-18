# LinReTraCe

## Official Resources
- Homepage: https://github.com/linretracedev/linretrace
- Documentation: https://github.com/linretracedev/linretrace/wiki (or minimal docs on repo)
- Source Repository: https://github.com/linretracedev/linretrace
- License: GPLv3

## Overview
LinReTraCe (Linear Response Transport Centre) is a massively parallel code for calculating transport properties of solids. It is specifically designed to work with spectral functions from many-body calculations (like DMFT), capturing lifetime effects and renormalization beyond the constant relaxation time approximation.

**Scientific domain**: Quantum transport, Thermoelectrics, Correlated materials
**Target user community**: DMFT practitioners, Thermoelectricity researchers

## Theoretical Methods
- Linear Response Theory (Kubo formalism)
- Boltzmann Transport Equation (with energy/momentum dependent lifetimes)
- Optical conductivity
- Seebeck coefficient
- Thermal conductivity
- Hall effect (some versions)
- Integration over Brillouin Zone

## Capabilities (CRITICAL)
- **Transport Coefficients**: DC/Optical conductivity, Seebeck, Thermal conductivity.
- **Many-Body Inputs**: Accepts self-energies or spectral functions from DMFT.
- **Interfaces**: Works with Wien2k, VASP, and Wannier90 inputs.
- **Parallelization**: MPI parallelization for efficient k-space integration.
- **Verification**: Artifact-free integration schemes.

## Key Features

### DMFT-Transport Link:
- Bridges the gap between DMFT electronic structure and experimental transport observables.
- Handles frequency-dependent self-energies.

### Efficiency:
- Written in Fortran/C++.
- Optimized for large k-grids required for transport convergence.

### HDF5 Output:
- All results are stored in structured HDF5 files for efficient data management and analysis.

## Inputs & Outputs
- **Input formats**:
  - `LINRETRACE.in` (aka `config.lrtc`): Free-format configuration file defining calculation parameters, temperature range, and desired observables.
  - Electronic structure data (eigenvalues, velocities) from DFT/Wannier.
  - Self-energy $\Sigma(\omega)$ from DMFT.
- **Output data types**:
  - HDF5 files (`*.h5`) containing:
    - Transport tensors ($\sigma$, $S$, $\kappa$) vs Temperature or Chemical Potential.
    - Optical conductivity spectra $\sigma(\omega)$.
  - `lprint` tool provided to extract/plot data from HDF5 files.

## Interfaces & Ecosystem
- **Upstream**: Reads data from Wien2k, VASP, Wannier90.
- **Downstream**: Produces HDF5 data; analysis via `lprint` or `h5py` (Python).

## Workflow and Usage
1.  Perform DFT+DMFT calculation.
2.  Generate Wannier functions or use DFT velocities.
3.  Export Self-energy.
4.  Configure `LINRETRACE.in` with desired T-mesh and chemical potential range.
5.  Run LinReTraCe to integrate Kubo formulas over the BZ.
6.  Use `lprint` to visualize the transport coefficients from the HDF5 output.

## Performance Characteristics
- **Scaling**: Scales well with k-points via MPI.
- **Accuracy**: Sophisticated integration tetrahedrons/adaptive schemes to handle Fermi surface complexity.


## Comparison with Other Codes
- **vs BoltzTraP**: BoltzTraP uses constant relaxation time (semi-classical); LinReTraCe captures lifetime effects and renormalization from self-energies.
- **vs BoltzWann**: Similar Wannier-based approach, but LinReTraCe focuses on linear response with full self-energy inputs.
- **vs TRIQS/transport**: LinReTraCe is a dedicated standalone transport code optimized for large k-grids.
- **Unique strength**: Handling of frequency-dependent self-energies (lifetimes) in transport coefficients.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/linretracedev/linretrace
2. Publication: "LinReTraCe: The Linear Response Transport Centre", arXiv/PRB.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Integration: Standard link to DMFT and Wannier workflows
- Focus: Dedicated transport post-processing
