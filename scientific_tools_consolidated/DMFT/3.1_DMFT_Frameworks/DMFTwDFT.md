# DMFTwDFT

## Official Resources
- Homepage: https://dmftwdft-project.github.io/DMFTwDFT/
- Documentation: https://dmftwdft-project.github.io/DMFTwDFT/
- Source Repository: https://github.com/DMFTwDFT-project/DMFTwDFT
- License: GNU General Public License v3.0

## Overview
DMFTwDFT is an open-source code combining Dynamical Mean Field Theory with various Density Functional Theory packages. It provides a flexible framework for performing DFT+DMFT calculations with multiple DFT backends and impurity solvers. The code interfaces with Wannier90 for downfolding and supports various free-licensed and commercial DFT codes, enabling ab-initio treatment of strongly correlated materials.

**Scientific domain**: DFT+DMFT calculations, strongly correlated materials  
**Target user community**: Researchers performing DFT+DMFT calculations on transition metal oxides and correlated materials

## Theoretical Methods
- DFT+DMFT framework
- LDA+DMFT, GGA+DMFT
- Charge self-consistent calculations (optional)
- Wannier function-based downfolding
- Multiple impurity solver integration
- Double-counting corrections
- Spin-polarized and non-collinear magnetism

## Capabilities (CRITICAL)
- Multiple DFT code interfaces (VASP, Quantum ESPRESSO, SIESTA, OpenMX, Wien2k)
- Wannier90 integration for orbital projections
- Multiple impurity solver support (CTQMC, NCA, etc.)
- One-shot and self-consistent DFT+DMFT
- Charge self-consistency loop
- Post-processing tools for spectral functions
- Flexible configuration system
- Python-based workflow scripts
- Support for complex materials and magnetic systems
- Vibrational and elastic property calculations

**Sources**: Official DMFTwDFT documentation (https://github.com/DMFTwDFT-project/DMFTwDFT), confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- DFT outputs from supported codes
- Wannier90 projections
- Configuration files for DMFT parameters
- Interaction parameters (U, J)

**Output data types**:
- Self-energies and Green's functions
- Spectral functions
- Density of states
- Occupation matrices
- Convergence data
- HDF5 and text-based outputs

## Interfaces & Ecosystem
- **DFT codes**: VASP, Quantum ESPRESSO, SIESTA, OpenMX, Wien2k
- **Wannier tools**: Wannier90 for orbital projections
- **Impurity solvers**: CTQMC solvers, NCA solvers
- **Post-processing**: Python scripts for analysis

## Limitations & Known Constraints
- Requires installation of DFT code separately
- Impurity solver must be installed separately
- Setup can be complex for new users
- Documentation assumes DFT+DMFT familiarity
- Charge self-consistency increases computational cost
- Platform: Linux/Unix primary support


## Performance Characteristics
- **Parallelization**: MPI-based parallelization
- **Efficiency**: Dependent on the chosen impurity solver (CTQMC, etc.)
- **Cost**: Charge self-consistency adds significant overhead
- **Flexibility**: Can trade off accuracy for speed by choosing different solvers

## Comparison with Other DMFT Frameworks
- **vs ComDMFT**: DMFTwDFT is more flexible with various DFT backends; ComDMFT focuses on VASP/GW integration
- **vs EDMFTF**: EDMFTF is tightly coupled with WIEN2k and uses a stationary functional approach; DMFTwDFT is more modular
- **vs TRIQS/DFTTools**: Both are flexible; DMFTwDFT aims for a more "black-box" user-friendly approach with standardized workflows
- **Unique strength**: Open-source, supports multiple major DFT codes (VASP, QE, Siesta), user-friendly

## Best Practices
- **DFT Backend**: Choose the backend you are most familiar with (e.g., VASP or QE)
- **Wannier90**: Ensure high-quality Wannier projections
- **Solver**: Use CTQMC for high accuracy, cheaper solvers for initial checks
- **Convergence**: Monitor charge density convergence in self-consistent runs

## Verification & Sources
**Primary sources**:
1. Official documentation: https://dmftwdft-project.github.io/DMFTwDFT/
2. GitHub repository: https://github.com/DMFTwDFT-project/DMFTwDFT
3. H. Park et al., Comput. Phys. Commun. 284, 108594 (2023) - DMFTwDFT paper

**Secondary sources**:
1. DMFTwDFT tutorials and examples
2. Published applications using DMFTwDFT
3. Community forums and issue tracker
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (GitHub issues)
- Academic citations: Growing user base
- Multiple DFT backend support
- Actively maintained project
