# sisl

## Official Resources
- **Homepage**: https://zerothi.github.io/sisl/
- **Documentation**: https://zerothi.github.io/sisl/docs/latest/
- **Source Repository**: https://github.com/zerothi/sisl
- **License**: Mozilla Public License 2.0 (MPL-2.0)

## Overview
**sisl** is a high-performance, modern Python framework for electronic structure calculations and large-scale tight-binding modeling. Developed as a successor to the utility scripts of TranSIESTA, it has evolved into a general-purpose API that interfaces with multiple DFT codes (SIESTA, VASP, OpenMX, Wannier90, BigDFT) to manipulate Hamiltonians, geometries, and real-space grids. It is particularly renowned for its ability to handle extremely large sparse matrices, enabling tight-binding calculations on systems with millions of orbitals which are inaccessible to standard dense-matrix tools.

**Scientific domain**: Quantum Transport, Large-scale Tight-Binding, DFT Post-processing
**Target user community**: Users of SIESTA/TranSIESTA, researchers modeling large nanostructures, and developers needing a robust localized-basis API

## Theoretical Methods
- **Tight-Binding (TB) / LCAO**: Core data structure is the sparse Hamiltonian matrix $H_{ij}$ in a localized basis.
- **Green's Functions**: Interfaces with **TBtrans** to compute transport properties using Non-Equilibrium Green's Functions (NEGF).
- **Real-Space Grid Algebra**: Efficient manipulation of charge densities and potentials on 3D grids.
- **Geometry Operations**: Algorithms for constructing supercells, adding strain, creating vacancies, and building complex heterostructures.

## Capabilities
- **Input/Output**:
  - Reads/Writes: **SIESTA** (`.fdf`, `.XV`, `.TSHS`), **VASP** (`POSCAR`, `KPOINTS`), **OpenMX**, **Wannier90** (`_hr.dat`), **BigDFT**, **XYZ**, **PDB**.
  - Binary file support for high-speed I/O (NetCDF).
- **Analysis**:
  - **Electronic**: Band structure, Density of States (DOS), Projected DOS (PDOS), COOP/COHP.
  - **Wavefunctions**: Real-space plotting of orbitals and wavefunctions.
  - **Transport**: Analysis of transmission spectra and currents (with TBtrans).
- **Model Construction**:
  - Create Graphene/TMD nanoribbons, nanotubes, and flakes.
  - Twist bilayers, apply strain, and create defects.
  - Construct tight-binding models from scratch with nearest-neighbor rules.

## Key Strengths
- **Scalability**: Optimized C/Cython backend allows handling of sparse matrices with millions of rows/columns, far exceeding the capacity of `numpy` or `scipy.sparse` for physics workflows.
- **Interoperability**: Acts as a "Rosetta Stone" converting between different DFT formats (e.g., VASP structure to SIESTA input).
- **Transport Workflow**: Tightly integrated with **TBtrans**, providing a seamless pythonic interface for setting up and analyzing complex transport calculations.

## Inputs & Outputs
- **Inputs**: Almost any standard structure or Hamiltonian file from supported DFT codes.
- **Outputs**:
  - Normalized files for other codes.
  - Data arrays (Bands, DOS) for plotting.
  - Visualization files (VMD, XSF, Cube).

## Interfaces & Ecosystem
- **Upstream**:
  - **SIESTA**: Deep integration (reads dense/sparse matrices, grids).
  - **Wannier90**: Can read Hamiltonians for large-scale transport setup.
- **Downstream**:
  - **TBtrans**: The primary transport solver associated with sisl.
  - **KITE**: Interface available for quantum transport.

## Performance Characteristics
- **Speed**: Critical sections (Hamiltonian construction, neighbor searching) are written in Cython/C.
- **Memory**: Sparse matrix storage (CSR/CSC) ensures minimal memory footprint for large empty systems.
- **Parallelism**: OpenMP threading for computationally intensive matrix operations.

## Limitations & Known Constraints
- **Basis Limitation**: Strictly localized basis sets (LCAO, Wannier, TB); not suitable for plane-wave data (except for grid operations).
- **Solver**: sisl itself is primarily a manipulator/analyzer; it diagonalizes small/medium matrices but relies on external codes (TBtrans) for heavy inversion/transport tasks on large systems.

## Comparison with Other Codes
- **vs. ASE (Atomic Simulation Environment)**: ASE focuses on geometry and calculators; sisl focuses on the *Hamiltonian* and *Electronic Structure* matrices themselves, offering much deeper access to the physics of the model.
- **vs. PythTB**: PythTB is excellent for small toy models; sisl is built for "production" large-scale DFT-derived tight-binding models.

## Application Areas
- **Nanodevices**: Modeling FETs, molecular junctions, and interconnects.
- **2D Materials**: Twisted bilayers (Moiré physics) and defect engineering.
- **Topology**: Constructing large hamiltonians for topological invariant analysis (e.g., with Kite).

## Community and Support
- **Development**: Maintained by Nick Papior (DTU/eScience).
- **Documentation**: Extensive documentation and tutorials.
- **Forum**: Active Discord channel and GitHub discussions.

## Verification & Sources
- **Official Website**: [https://zerothi.github.io/sisl/](https://zerothi.github.io/sisl/)
- **Primary Publication**: N. Papior et al., Comp. Phys. Comm. 212, 8 (2017) (TBtrans/sisl paper).
- **Verification status**: ✅ VERIFIED
  - widely used in the TranSIESTA community.
  - >400 citations.
