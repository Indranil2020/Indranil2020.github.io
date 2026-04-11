# spglib

## Official Resources
- Homepage: https://spglib.github.io/spglib/
- Documentation: https://spglib.github.io/spglib/
- Source Repository: https://github.com/spglib/spglib
- License: BSD 3-Clause License

## Overview
Spglib is a C library for finding and handling crystal symmetries. It provides algorithms for finding space groups, symmetry operations, and primitive unit cells. Spglib is widely used as the underlying symmetry engine for many major materials science codes (including Phonopy, Phono3py, Pymatgen, ASE, and Quantum ESPRESSO) due to its robustness and efficiency.

**Scientific domain**: Crystallography, symmetry analysis, space groups  
**Target user community**: Developers of materials science software, crystallographers

## Capabilities (CRITICAL)
- **Space Group Determination**: Finds International Table number, Hall symbol, and Hermann-Mauguin symbol
- **Symmetry Operations**: Identifies rotation and translation operations
- **Standardization**: Converts arbitrary unit cells to standardized primitive or conventional cells
- **Wyckoff Positions**: Identifies Wyckoff positions of atoms
- **K-point Reduction**: Reduces k-point grids based on symmetry (for Brillouin zone sampling)
- **Magnetic Symmetry**: Support for magnetic space groups (msg)
- **Tolerance**: Robust handling of numerical noise with adjustable tolerance

**Sources**: Spglib documentation, arXiv:1805.05143

## Inputs & Outputs
- **Input formats**: Crystal structure (lattice vectors, atomic positions, atomic numbers) passed via API
- **Output data types**: Symmetry dataset (space group number, operations, standardized cell)

## Interfaces & Ecosystem
- **C API**: Core library interface
- **Python**: `spglib` python package
- **Fortran**: Interface available
- **Ruby**: Interface available
- **Julia**: LibSymspg.jl wrapper

## Workflow and Usage
1. Define crystal structure (lattice, positions, numbers).
2. Call `get_symmetry_dataset(cell, symprec=1e-5)`.
3. Retrieve space group number or operations.
4. Call `standardize_cell()` to get primitive cell.

## Performance Characteristics
- Highly optimized C implementation
- Extremely fast (milliseconds for typical cells)
- Robust against small structural distortions

## Application Areas
- High-throughput screening (structure classification)
- Band structure calculations (k-path standardization)
- Phonon calculations (force constant symmetry)
- Structure database creation

## Community and Support
- Open-source (BSD)
- Developed by Atsushi Togo (NIMS/Kyoto University)
- Widely adopted standard

## Verification & Sources
**Primary sources**:
1. Homepage: https://spglib.github.io/spglib/
2. GitHub: https://github.com/spglib/spglib

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Togo)
- Applications: Symmetry analysis, space groups, standardization
