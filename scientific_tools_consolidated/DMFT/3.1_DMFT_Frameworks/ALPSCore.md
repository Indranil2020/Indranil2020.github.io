# ALPSCore

## Official Resources
- Homepage: https://alpscore.org/
- Documentation: https://alpscore.org/documentation.html
- Source Repository: https://github.com/ALPSCore/ALPSCore
- License: GNU General Public License v2.0

## Overview
ALPSCore (ALPS Core Libraries) represents the modernized core libraries extracted from the ALPS project. It provides a set of maintained, well-documented, and reusable C++ libraries for condensed matter physics simulations, with a focus on strongly correlated electron systems. ALPSCore libraries are designed to be lightweight, easy to integrate, and provide essential functionality for physics applications.

**Scientific domain**: Condensed matter physics, strongly correlated systems  
**Target user community**: Developers and users of condensed matter physics simulation codes

## Theoretical Methods
- Core library infrastructure for physics simulations
- Support for Monte Carlo algorithms
- Lattice model utilities
- Statistical analysis tools
- Parallel computing infrastructure

## Capabilities (CRITICAL)
- Generic C++ physics libraries
- HDF5 I/O support
- MPI parallelization utilities
- Statistical accumulation and analysis
- Parameter parsing and management
- Random number generation
- Lattice utilities
- Used by multiple physics codes (CT-INT, opendf, etc.)
- CMake-based build system
- Well-documented API

**Sources**: Official ALPSCore website (https://alpscore.org/), GitHub repository, confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- C++ API for library integration
- HDF5 data files
- Parameter files

**Output data types**:
- HDF5 archives
- Statistical results
- Library functions return standard C++ types

## Interfaces & Ecosystem
- **Applications using ALPSCore**: CT-INT, opendf, various QMC codes
- **Programming language**: C++ with Python bindings (some components)
- **Parallel computing**: MPI support
- **Data format**: HDF5 for efficient I/O

## Limitations & Known Constraints
- Primarily a library, not a standalone application
- Requires C++ compilation knowledge
- Documentation assumes programming experience
- Some features still under development
- Smaller community than original ALPS project

## Verification & Sources
**Primary sources**:
1. Official website: https://alpscore.org/
2. GitHub repository: https://github.com/ALPSCore/ALPSCore
3. G. Guertler et al., Comput. Phys. Commun. 228, 216 (2018) - ALPSCore paper

**Secondary sources**:
1. ALPSCore tutorials
2. Applications using ALPSCore
3. ALPS project documentation
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: OPEN (GitHub, GPL v2)
- Community support: Active development
- Modern successor to ALPS libraries
- Used by multiple physics applications
