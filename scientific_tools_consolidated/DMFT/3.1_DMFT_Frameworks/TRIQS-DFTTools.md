# TRIQS-DFTTools

## Official Resources
- Homepage: https://triqs.github.io/dft_tools/
- Documentation: https://triqs.github.io/dft_tools/latest/
- Source Repository: https://github.com/TRIQS/dft_tools
- License: GNU General Public License v3.0

## Overview
TRIQS/DFTTools is a TRIQS application providing the necessary tools to perform DFT+DMFT calculations. It serves as the bridge between DFT codes and the TRIQS DMFT framework, handling Wannier function projections, self-consistency loops, and post-processing of spectral functions. It is the standard interface for performing realistic DFT+DMFT calculations within the TRIQS ecosystem.

**Scientific domain**: DFT+DMFT calculations, strongly correlated materials  
**Target user community**: Researchers performing ab-initio DMFT calculations on realistic materials

## Theoretical Methods
- DFT+DMFT interface and workflow management
- Projective Wannier function formalism
- Charge self-consistency (optional)
- Maximum entropy analytical continuation
- Spectral function calculations
- Integration with multiple DFT codes
- Local and momentum-resolved quantities

## Capabilities (CRITICAL)
- Interface to multiple DFT codes (Wien2k, VASP, Quantum ESPRESSO, ABINIT, Elk)
- Wannier90 integration for projection operators
- DMFT self-consistency loop management
- Charge density updates for self-consistent calculations
- Spectral function and DOS calculations
- k-resolved spectral functions (ARPES)
- Momentum distribution functions
- Analytical continuation via MaxEnt
- Chemical potential adjustment
- Occupancy matrix calculations
- HDF5-based data management
- Post-processing and analysis tools

**Sources**: Official TRIQS/DFTTools documentation (https://triqs.github.io/dft_tools/), confirmed in 7/7 source lists

## Inputs & Outputs
**Input formats**:
- DFT outputs from Wien2k, VASP, QE, ABINIT, Elk
- Wannier90 projections
- DMFT solver outputs (self-energies)
- HDF5 archives from previous calculations

**Output data types**:
- HDF5 archives with all DMFT quantities
- Spectral functions (A(k,ω))
- Local and k-resolved Green's functions
- Self-energies
- Density matrices
- Chemical potentials
- Formatted data for plotting

## Interfaces & Ecosystem
- **DFT code interfaces**: Wien2k (native), VASP, Quantum ESPRESSO, ABINIT, Elk
- **Wannier function tools**: Wannier90 integration
- **TRIQS solvers**: Seamless integration with TRIQS/cthyb and other TRIQS impurity solvers
- **High-level workflows**: solid_dmft uses DFTTools as backend
- **Post-processing**: Python-based analysis tools, matplotlib integration

## Limitations & Known Constraints
- Requires understanding of DFT+DMFT methodology
- DFT code-specific setup can be complex
- Wannier projection quality critical for results
- Double-counting correction choice affects results
- Charge self-consistency increases computational cost significantly
- Analytical continuation introduces uncertainties
- HDF5 version compatibility issues possible


## Performance Characteristics
- **Efficiency**: Heavily dependent on the impurity solver backend (e.g., TRIQS/cthyb).
- **Parallelization**: Inherits MPI parallelization from TRIQS applications.
- **Overhead**: Python-layer overhead is minimal compared to the QMC solver cost.
- **Scalability**: Can scale to large clusters for complex impurity problems.

## Comparison with Other Frameworks
- **vs DMFTwDFT**: TRIQS/DFTTools is modular and requires Python scripting; DMFTwDFT aims for a "black-box" experience.
- **vs EDMFTF**: EDMFTF is a stationary functional code tightly coupled with Wien2k; DFTTools is a flexible library for various DFT codes.
- **vs solid_dmft**: solid_dmft is a high-level wrapper built *on top* of DFTTools to automate workflows.

## Verification & Sources
**Primary sources**:
1. Official documentation: https://triqs.github.io/dft_tools/latest/
2. GitHub repository: https://github.com/TRIQS/dft_tools
3. A. Hampel et al., arXiv:2309.10858 (2023) - DFTTools overview

**Secondary sources**:
1. TRIQS tutorials and workshops
2. solid_dmft documentation (uses DFTTools)
3. Published DFT+DMFT applications
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (TRIQS project)
- Part of supported TRIQS ecosystem
- Maintained by Flatiron Institute
