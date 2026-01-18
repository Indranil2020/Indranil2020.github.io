# opendf

## Official Resources
- Homepage: https://github.com/CQMP/opendf
- Documentation: https://github.com/CQMP/opendf/blob/master/README.md
- Source Repository: https://github.com/CQMP/opendf
- License: GNU General Public License v2.0

## Overview
opendf is a condensed matter physics code that solves strongly correlated lattice problems (such as the Hubbard model) in finite dimensions using the dual fermion method. It extends DMFT by including non-local correlations through diagrammatic extensions, providing a systematic way to treat spatial correlations beyond local DMFT approximations.

**Scientific domain**: Strongly correlated systems, dual fermion method, diagrammatic extensions of DMFT
**Target user community**: Researchers studying non-local correlations in strongly correlated materials

## Theoretical Methods
- Dual fermion (DF) method
- Diagrammatic extensions of DMFT
- Non-local correlation effects
- Ladder approximation
- Hubbard model in finite dimensions
- Spatial correlations beyond DMFT
- Vertex function calculations

## Capabilities (CRITICAL)
- **Dual fermion calculations**: Systematic inclusion of non-local correlations.
- **Hubbard model solutions**: Supports 2D and 3D lattice models.
- **Diagrammatic extensions**: Includes ladder diagrams and higher-order correlations.
- **Data Standard**: Utilizes the Open Data Format (ODF) for data exchange.
- **Parallelization**: MPI parallelization for intensive vertex calculations.

## Inputs & Outputs
- **Input formats**:
  - **ODF Packages (`.odf.zip`)**: A zipped archive containing:
    - **Data (`.csv`)**: Primary numerical data.
    - **Metadata (`.xml`)**: Descriptive metadata using DDI-Codebook schema.
    - **Version (`.json`)**: Versioning information.
  - Configuration files for solver parameters.
- **Output data types**:
  - Self-energies (momentum-dependent)
  - Green's functions
  - Vertex functions
  - Observables
  - ODF-compliant output archives.

## Interfaces & Ecosystem
- **ALPSCore**: Built on ALPSCore libraries for Monte Carlo and grid utilities.
- **C++**: Modern C++ implementation.
- **HDF5**: Standard data format for internal storage.
- **MPI**: Parallel execution.

## Limitations & Known Constraints
- **Computational Cost**: Higher than standard DMFT due to vertex function calculations.
- **Dependency**: Requires ALPSCore libraries, which can be complex to install.
- **Scope**: Primarily focused on model Hamiltonians (Hubbard) rather than full ab-initio materials (unless interfaced).


## Performance Characteristics
- **Cost**: Significantly higher than standard DMFT due to vertex function calculation ($O(N^4)$ complexity).
- **Parallelization**: MPI parallelization crucial for diagrammatic summations.
- **Scaling**: Scales well on large clusters.

## Comparison with Other Methods
- **vs Standard DMFT**: opendf includes non-local spatial correlations via Dual Fermions; standard DMFT is purely local.
- **vs GW+DMFT**: Dual Fermion is a diagrammatic extension of DMFT; GW+DMFT combines GW results with DMFT.
- **Unique strength**: Systematic inclusion of non-local correlations using the Dual Fermion method.

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/CQMP/opendf
2. README and code documentation.
3. Open Data Format specifications.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Method: Dual Fermion implementation
- Data Standard: Adopts ODF
