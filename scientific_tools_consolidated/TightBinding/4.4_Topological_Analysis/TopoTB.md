# TopoTB

## Official Resources
- **Homepage**: https://github.com/xlhuang-phy/TopoTB
- **Repository**: https://github.com/xlhuang-phy/TopoTB
- **License**: GPL-3.0

## Overview
**TopoTB** is a software package written in **Mathematica** for the interactive calculation and visualization of topological properties in tight-binding models. Its standout feature is its use of Mathematica's `Manipulate` functionality to create **real-time interactive phase diagrams**, allowing users to dynamically tune Hamiltonian parameters (like mass or hopping strength) and immediately see the effect on band structures, Berry curvature, and edge states.

**Scientific domain**: Topological Phases, Education, Model Building
**Target user community**: Theorists and Students exploring phase transitions

## Theoretical Methods
- **Tight-Binding**: Symbolic or numerical Hamiltonian definitions.
- **Topological Invariants**:
  - **Chern Number**: Integration of Berry curvature $\Omega_{xy}$.
  - **$\mathbb{Z}_2$ Invariant**: Via Wilson loops or Parity criteria.
- **Surface States**: Slab calculations to reveal bulk-boundary correspondence.

## Capabilities
- **Interactive UI**: Sliders to vary $t_1, t_2, \phi, M$, etc.
- **Vizualization**:
  - 3D Band structures.
  - Berry curvature density plots.
  - Edge state dispersion.
- **Invariants**: Automatic calculation of Chern/Z2 numbers for the current parameter set.
- **Phase Diagrams**: Can scan parameters to map out phase boundaries (gap closings).

## Key Strengths
- **Interactivity**: The ability to *watch* a gap close and the edge states change chirality in real-time is unmatched for pedagogical purposes and building intuition.
- **Symbolic Power**: Being based in Mathematica, it can handle analytic derivations of models before numerical evaluation.
- **Ease of Entry**: Requires no compilation, just dragging sliders.

## Inputs & Outputs
- **Inputs**: Mathematica Notebook cells defining the lattice and Hamiltonian.
- **Outputs**: Interactive plots and numerical invariant values.

## Interfaces & Ecosystem
- **Platform**: Requires Wolfram Mathematica.
- **Language**: Wolfram Language.

## Performance Characteristics
- **Speed**: Fast for small "toy models" (2-8 bands) typical in topological physics. Slower for large supercells compared to compiled codes like C++ or Julia.
- **Scalability**: Not intended for large-scale atomistic simulations, but for model Hamiltonians.

## Comparison with Other Codes
- **vs. PythTB**: PythTB is a Python library (script based). TopoTB is a Mathematica UI (interactive). TopoTB is better for "playing" with a model; PythTB is better for systematic scripting.
- **vs. Kwant**: Kwant is for transport. TopoTB is for band topology.

## Application Areas
- **Education**: Demonstrating the Haldane model or Kane-Mele model in class.
- **Research**: Rapidly checking the topology of a new effective Hamiltonian.

## Community and Support
- **Development**: Xielin Huang.
- **Source**: GitHub.

## Verification & Sources
- **Primary Publication**: X. Huang et al., arXiv:2403.08615 (2024).
- **Verification status**: ✅ VERIFIED
  - Active interactive tool.
