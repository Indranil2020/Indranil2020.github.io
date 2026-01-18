# ccao-unfold

## Official Resources
- **Repository**: https://github.com/ccao/unfold
- **License**: GNU General Public License v3.0

## Overview
**ccao-unfold** is a lightweight but powerful utility for **band unfolding**. It allows researchers to recover the effective primitive-cell band structure from calculations performed in a supercell. This is essential for analyzing systems with broken translational symmetry, such as defects, alloys, or twisted bilayers, where the "folded" supercell bands are often too dense to interpret. Unlike plane-wave unfolding codes, ccao-unfold operates on the **Wannier90** tight-binding Hamiltonian, making it highly efficient and code-agnostic (as long as the code interfaces with Wannier90).

**Scientific domain**: Band Structure Analysis, Disordered Systems, Defects
**Target user community**: Researchers studying impurities, alloys, or surface reconstructions using supercells

## Theoretical Methods
- **Band Unfolding**: Projects eigenstates of the supercell Hamiltonian onto the Bloch states of the primitive cell.
- **Spectral Weight**: Calculates the spectral function $A(k_{prim}, E)$, representing the probability of finding a primitive Bloch state at a given energy.
- **Wannier Basis**: Performs the unfolding algebra using the localized Wannier basis, which avoids the need for large wavefunction files.

## Capabilities
- **Unfolding**:
  - From Supercell to Primitive Cell.
  - From large Interface cells to bulk-like projections.
- **Analysis**:
  - Resolves "shadow bands" (ghost bands) from folding.
  - Identifies the localized nature of defect states by their lack of spectral weight in the primitive BZ.
- **Input Flexibility**: 
  - Works with any Hamiltonian in the `_hr.dat` format (Wannier90 standard).

## Key Strengths
- **Efficiency**: Post-processing takes seconds to minutes, compared to the hours typical for wavefunction-based unfolding.
- **Portability**: Does not depend on the specific DFT code (VASP, QE, Siesta) used to generate the Wannier functions.
- **Simplicity**: Single-purpose tool that does one thing well.

## Inputs & Outputs
- **Inputs**:
  - `wannier90_hr.dat`: The supercell Hamiltonian.
  - `KPATH`: Specification of the primitive cell k-path.
  - `POSCAR` (or equivalent): Implementation dependent structural info for mapping.
- **Outputs**:
  - Spectral weight data (k-point, Energy, Weight) suitable for plotting "fat bands" or heatmaps.

## Interfaces & Ecosystem
- **Upstream**: **Wannier90** (and any DFT code that feeds it).
- **Visualization**: Outputs are typically plotted with Python (Matplotlib) or Gnuplot.

## Performance Characteristics
- **Speed**: Very fast; limited only by the size of the Hamiltonian matrix multiplication.
- **Memory**: Moderate usage, proportional to the number of orbitals in the supercell.

## Limitations & Known Constraints
- **Mapping**: Requires a clear geometric mapping between the supercell and the primitive cell (integer transformation matrix).
- **Basis Consistency**: Assumes the Wannier functions in the supercell can be essentially mapped to those in the primitive cell (valid for standard localized orbitals).

## Comparison with Other Codes
- **vs. BandUP**: BandUP is a more general plane-wave folder/unfolder; ccao-unfold is specific to Wannier models and thus faster for this specific workflow.
- **vs. VASP Unfolding**: VASP has built-in unfolding (via `LORBIT` extensions in some patches), but ccao-unfold works as a post-processor for any code.

## Application Areas
- **Doped Semiconductors**: Visualizing how dopant bands emerge within the gap.
- **High-Entropy Alloys**: Understanding the broadening of bands due to chemical disorder.
- **Surface States**: Differentiating surface resonances from bulk projected bands.

## Community and Support
- **Development**: Developed by Changmiao Cao.
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/ccao/unfold](https://github.com/ccao/unfold)
- **Verification status**: ✅ VERIFIED
  - Functional repository.
  - Methodologically standard approach to TB unfolding.
