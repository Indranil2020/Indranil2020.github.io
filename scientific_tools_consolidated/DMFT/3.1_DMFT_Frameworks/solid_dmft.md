# solid_dmft

## Official Resources
- Homepage: https://triqs.github.io/solid_dmft/
- Documentation: https://triqs.github.io/solid_dmft/latest/
- Source Repository: https://github.com/TRIQS/solid_dmft
- License: GNU General Public License v3.0

## Overview
solid_dmft is a versatile Python wrapper to perform DFT+DMFT calculations utilizing the TRIQS software library. It provides a high-level, user-friendly interface for one-shot and charge self-consistent DFT+DMFT calculations, supporting multiple DFT codes (VASP, Quantum ESPRESSO) and impurity solvers. Designed for accessibility, it automates many workflow steps while maintaining flexibility for advanced users.

**Scientific domain**: DFT+DMFT calculations, strongly correlated materials  
**Target user community**: Researchers performing realistic DFT+DMFT calculations on correlated materials

## Theoretical Methods
- DFT+DMFT (one-shot and charge self-consistent)
- LDA+DMFT, GGA+DMFT
- Integration with TRIQS/DFTTools backend
- Multiple impurity solver support
- Wannier function downfolding
- Double-counting corrections (FLL, AMF, etc.)
- Spectral function calculations

## Capabilities (CRITICAL)
- High-level Python interface for DFT+DMFT workflows
- VASP interface (native support)
- Quantum ESPRESSO interface
- Wannier90 integration for projections
- Multiple impurity solvers (TRIQS/cthyb, TRIQS/Hubbard-I, w2dynamics)
- One-shot DFT+DMFT calculations
- Charge self-consistent calculations
- Automated convergence checking
- Flexible configuration via TOML files
- Post-processing and analysis tools
- Spectral function calculations
- k-resolved quantities
- Tutorial-driven documentation
- Example calculations included

**Sources**: Official solid_dmft documentation (https://triqs.github.io/solid_dmft/), A. Hampel et al., arXiv:2103.13522 (2021), confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- TOML configuration files
- DFT outputs from VASP or Quantum ESPRESSO
- Wannier90 projections
- HDF5 archives

**Output data types**:
- HDF5 archives with all DMFT quantities
- Spectral functions
- Self-energies and Green's functions
- Occupancies and observables
- Convergence histories
- Formatted output for plotting

## Interfaces & Ecosystem
- **DFT codes**: VASP (primary), Quantum ESPRESSO
- **TRIQS ecosystem**: Uses TRIQS/DFTTools and TRIQS/cthyb
- **Impurity solvers**: TRIQS/cthyb, Hubbard-I, w2dynamics
- **Wannier tools**: Wannier90 for orbital projections
- **Python ecosystem**: Modern Python-based workflow

## Limitations & Known Constraints
- Requires TRIQS installation and dependencies
- DFT code dependency (VASP or QE required)
- Impurity solver must be separately installed
- Learning curve for DFT+DMFT methodology
- Computational cost from DMFT iterations
- Charge self-consistency significantly increases cost
- Documentation assumes DMFT familiarity


## Performance Characteristics
- **Efficiency**: Inherits the performance of the underlying TRIQS/cthyb solver.
- **Overhead**: Python overhead is negligible; cost is dominated by the impurity solver and DFT runs.
- **Parallelization**: Fully supports MPI parallelization for both DFT (via VASP/QE support) and DMFT parts.

## Comparison with Other Frameworks
- **vs DMFTwDFT**: Both wrap DFT+DMFT workflows; solid_dmft relies on the TRIQS ecosystem, while DMFTwDFT is standalone/Wannier-based.
- **vs Zen**: solid_dmft is a Python wrapper for TRIQS; Zen is a native Julia/Fortran package.
- **vs Manual TRIQS**: solid_dmft automates the complex scripting usually required for TRIQS, acting as a user-friendly layer.

## Verification & Sources
**Primary sources**:
1. Official documentation: https://triqs.github.io/solid_dmft/
2. GitHub repository: https://github.com/TRIQS/solid_dmft
3. A. Hampel et al., J. Open Source Softw. 6(65), 3278 (2021) - solid_dmft paper

**Secondary sources**:
1. TRIQS documentation and tutorials
2. solid_dmft examples and tutorials
3. Published applications
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: Active (TRIQS project)
- Maintained by Flatiron Institute
- Growing user base
- Designed for accessibility
