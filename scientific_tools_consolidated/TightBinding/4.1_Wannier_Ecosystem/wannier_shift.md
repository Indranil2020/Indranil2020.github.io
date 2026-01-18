# wannier_shift

## Official Resources
- **Repository**: https://github.com/stmeurk/wannier_shift
- **Documentation**: https://stmeurk.github.io/wannier_shift/
- **License**: MIT License

## Overview
**wannier_shift** is a Python-based code designed for the construction and interpolation of tight-binding Hamiltonians for **Transition Metal Dichalcogenide (TMDC)** heterostructures, particularly those exhibiting **moiré patterns** due to twisting or lattice mismatch. It addresses the challenge of modeling large-scale superlattices by taking ab initio Wannier Tight-Binding parameters (from **Wannier90**) and applying a "shifting" interpolation scheme to generate the Hamiltonian for commensurate or incommensurate twisted bilayers at arbitrary k-points.

**Scientific domain**: 2D Materials, Twistronics, Moiré Physics
**Target user community**: Researchers studying electronic properties of twisted bilayers and heterostructures

## Theoretical Methods
- **Wannier Interpolation**: Extens standard interpolation to handle coupled layers.
- **Shifting Approximation**: Approximates the interlayer coupling by shifting the relative positions of Wannier centers, avoiding the need for a full DFT calculation of the supercell.
- **Continuum Models**: Generates effective continuum Hamiltonians for low-energy physics in moiré cells.

## Capabilities
- **Hamiltonian Construction**:
  - Builds large-scale tight-binding models for twisted bilayers ($>10,000$ atoms).
  - Handles generic stackings (AA, AB, etc.).
- **Band Structure**:
  - Calculates moiré bands (minibands).
  - Unfolding of supercell bands to primitive Brillouin zones.
- **Input Flexibility**: Use of standard `_hr.dat` files from monolayer calculations.

## Key Strengths
- **Scalability**: Enables the study of realistic moiré supercells that are computationally prohibitive for full DFT.
- **Efficiency**: Interpolation is orders of magnitude faster than evaluating couplings from scratch.
- **Versatility**: Applicable to various TMDC combinations (MoS2/WS2, WSe2/WS2, etc.).

## Inputs & Outputs
- **Inputs**:
  - `wannier90_hr.dat` for each monolayer.
  - Configuration script (twist angle, material choices).
- **Outputs**:
  - `hamiltonian.dat`: Sparse matrix of the supercell Hamiltonian.
  - `bands.dat`: Eigenvalues along specified paths.

## Interfaces & Ecosystem
- **Upstream**: **Wannier90** (provides the monolayer parameters).
- **Downstream**: Can feed into diagonalizers (like generic Python sparse solvers) or transport codes.

## Performance Characteristics
- **Speed**: Python-based construction is efficient for model setup; diagonalization depends on the solver used (e.g., typically `scipy.sparse`).
- **Memory**: Sparse matrix storage is essential for identifying moiré physics in large cells.

## Limitations & Known Constraints
- **Approximation**: The "shifting" assumes that the intralayer hopping is unmodified by the interlayer interaction (rigid band approximation), which works well for van der Waals heterostructures but may fail for strongly coupled layers.
- **Basis**: Restricted to systems where Wannier functions are well-localized.

## Comparison with Other Codes
- **vs. Twist-code**: Similar specialized tools exist (e.g., by Kaxiras group), but wannier_shift offers a lightweight, open-source Python implementation specifically for W90 users.
- **vs. Continuum Models**: More microscopic than pure continuum models (Bistritzer-MacDonald) as it retains the full tight-binding lattice structure.

## Application Areas
- **Twistronics**: Exploring magic angles in TMDCs.
- **Excitons**: Providing the electronic background for excitonic calculations in moiré potentials.
- **Topology**: Searching for topological flat bands.

## Community and Support
- **Development**: Developed by Stephen Carr (stmeurk) and collaborators.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/stmeurk/wannier_shift](https://github.com/stmeurk/wannier_shift)
- **Primary Publication**: S. Carr et al., Phys. Rev. B (related methodology).
- **Verification status**: ✅ VERIFIED
  - Functional tool for the 2D materials community.
