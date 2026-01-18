# Berry

## Official Resources
- **Homepage**: https://github.com/ricardoribeiro-2020/berry
- **Repository**: https://github.com/ricardoribeiro-2020/berry
- **License**: MIT License

## Overview
**Berry** is a Python package designed to calculate the Berry phase, Berry curvature, and related topological properties of crystalline materials directly from **Density Functional Theory (DFT)** wavefunctions, without relying on Maximally Localized Wannier Functions (MLWFs). It specifically interfaces with **Quantum ESPRESSO** to extract Bloch states and compute topological invariants and optical responses, such as the Anomalous Hall Conductivity (AHC) and Circular Dichroism.

**Scientific domain**: Topological Physics, Non-linear Optics
**Target user community**: DFT users studying topology and optical responses

## Theoretical Methods
- **Berry Connection & Curvature**: Direct calculation using finite differences of Bloch states in k-space.
- **Wavefunction Gauge**: Uses graph-based or AI-enhanced algorithms to ensure smooth gauge choices across the Brillouin Zone.
- **Optical Response**:
  - Anomalous Hall Conductivity (linear response).
  - Circular Dichroism (differential absorption of polarized light).
  - Second Harmonic Generation (SHG) ($\chi^{(2)}$ tensor).

## Capabilities
- **DFT Interface**:
  - Reads Quantum ESPRESSO `save` directories and XML files.
  - Handles spinor wavefunctions (SOC).
- **Calculations**:
  - Berry curvature vectors $\mathbf{\Omega}(\mathbf{k})$.
  - Chern numbers of 2D planes.
  - Optical conductivity spectra $\sigma_{\alpha\beta}(\omega)$.
  - SHG susceptibility tensor.
- **Tools**:
  - Band unfolding for supercells.
  - Visualization of curvature fields.

## Key Strengths
- **Direct Basis**: By working directly with DFT wavefunctions, it avoids potential artifacts or loss of information associated with Wannierization, especially for entangled bands where Wannier projection is difficult.
- **Optical Focus**: Specialized modules for Circular Dichroism and SHG make it unique among topological tools, which often focus only on static invariants.
- **Gauge Fixing**: Implements robust algorithms to handle the "random gauge" problem inherent in DFT codes, essential for numerical differentiation.

## Inputs & Outputs
- **Inputs**:
  - Quantum ESPRESSO output (`prefix.save`, `wfc` files).
  - Python configuration script.
- **Outputs**:
  - NumPy arrays of curvature and conductivity.
  - Plots of 2D/3D curvature distribution.

## Interfaces & Ecosystem
- **Upstream**: Quantum ESPRESSO.
- **Dependencies**: NumPy, SciPy.

## Performance Characteristics
- **Memory**: High. Processing full DFT wavefunctions requires significant RAM compared to tight-binding models.
- **Connectivty**: Requires dense k-grids in the DFT step for accurate derivatives.

## Comparison with Other Codes
- **vs. Wannier90**: Wannier90 interpolates bands to ultra-dense grids cheaply. Berry works on the original grid (or requires a dense DFT run), which is more expensive but stays "closer to the truth" of the DFT basis.
- **vs. Yambo**: Yambo is a full Many-Body/Optics code. Berry is a lighter, specialized tool for topological optics.

## Application Areas
- **Chiral Semimetals**: Calculating circular dichroism in enantiomorphic crystals.
- **Magnetic Topology**: AHC in ferromagnetic Weyl semimetals.

## Community and Support
- **Development**: Ricardo Mendes Ribeiro (University of Minho).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/ricardoribeiro-2020/berry](https://github.com/ricardoribeiro-2020/berry)
- **Verification status**: ✅ VERIFIED
  - Active research code.
