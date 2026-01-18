# Z2Pack

## Official Resources
- **Homepage**: https://z2pack.ethz.ch/
- **Documentation**: https://z2pack.greschd.ch/
- **Repository**: https://github.com/Z2PackDev/Z2Pack
- **License**: GPL-3.0

## Overview
**Z2Pack** is the premier tool for the automated calculation of topological invariants and Berry phases. Developed at ETH Zurich, it implements the **Wilson Loop** and **Wannier Charge Center (WCC)** tracking methods to calculate $\mathbb{Z}_2$ invariants, Chern numbers, and Weyl point chiralities without the need for manual inspection or gauge fixing. It interfaces seamlessly with both tight-binding models (via **TBmodels**) and ab-initio codes (VASP, Quantum ESPRESSO, Wannier90).

**Scientific domain**: Topological Classification, Berry Phase
**Target user community**: Theorists and computational researchers characterizing topological matter

## Theoretical Methods
- **Wilson Loops**: Non-Abelian Berry phase over a closed k-path.
- **Wannier Charge Centers (WCC)**: Hybrid Wannier function centers $\bar{x}(k_y, k_z)$.
- **$\mathbb{Z}_2$ Invariant**: Tracking the flow of WCCs (time-reversal polarization pumping).
- **Chern Number**: Winding of WCCs over the BZ.
- **Surface Green's Functions**: Iterative renormalization (via surface module).

## Capabilities
- **Invariants**:
  - $\mathbb{Z}_2$ indices (3D TI, 2D QSH).
  - Chern Numbers (QAH, Weyl Chirality).
- **Methods**:
  - Automated convergence of Wilson loops (adaptive grid).
  - Surface States spectral function.
- **Integration**:
  - **FP (First Principles)**: VASP, QE, Abinit.
  - **TB (Tight Binding)**: TBmodels, PythTB.

## Key Strengths
- **Automation**: Z2Pack's "Surface" tracking algorithm adaptively adds k-points where the Berry curvature is high (e.g., near gap closings), ensuring integer convergence of topological numbers without wasting resources on trivial parts of the BZ.
- **Rigorousness**: Unlike parity-based indicators (which fail without inversion symmetry), Z2Pack calculates the invariant directly from the Berry phase evolution, making it valid for *any* symmetry class (including generic Weyl semimetals).
- **Visualization**: Built-in plots for WCC evolution lines are publication-ready and crucial for debugging topological calculations.

## Inputs & Outputs
- **Inputs**:
  - System definitions (Python scripts).
  - DFT output files or TB model files.
- **Outputs**:
  - `result.json` (numerical invariants).
  - Plots of WCC flow.

## Interfaces & Ecosystem
- **TBmodels**: Tight coupling; models generated in TBmodels flow directly into Z2Pack.
- **VASP**: Specialized interface (`vasp.py`) to drive VASP calculations automatically.

## Performance Characteristics
- **Efficiency**: The adaptive algorithm is orders of magnitude more efficient than uniform grid integration for singular quantities (like Berry flux near a Weyl point).
- **Cost**: For DFT, cost is proportional to the number of SCF/Non-SCF calls triggered. For TB, it is very fast.

## Comparison with Other Codes
- **vs. WannierTools**: WannierTools computes surface *states* (Green's functions) very well. Z2Pack computes the *invariant* number itself more rigorously (tracking WCC lines). They are complementary.
- **vs. IrRep**: IrRep is faster but restricted to symmetry-protected phases. Z2Pack solves the general case (no symmetry required).

## Application Areas
- **Weyl Semimetals**: Determining the charge (chirality) of a node by enclosing it in a sphere.
- **Topological Insulators**: Confirming the non-trivial nature of Bi2Se3, inversion-asymmetric TIs.

## Community and Support
- **Development**: Dominik Gresch (Microsoft Quantum / ETH Zurich Alumni).
- **Source**: GitHub.

## Verification & Sources
- **Primary Publication**: D. Gresch et al., Phys. Rev. B 95, 075146 (2017).
- **Verification status**: ✅ VERIFIED
  - The gold standard for numerical Z2 calculations.
