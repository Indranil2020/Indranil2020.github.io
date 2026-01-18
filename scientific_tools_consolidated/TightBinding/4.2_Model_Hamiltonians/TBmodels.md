# TBmodels

## Official Resources
- **Homepage**: https://tbmodels.greschd.ch/
- **Documentation**: https://tbmodels.greschd.ch/en/latest/
- **Repository**: https://github.com/Z2PackDev/TBmodels
- **License**: GPL-3.0

## Overview
**TBmodels** is a Python library developed as part of the **Z2Pack** ecosystem, designed for reading, creating, and manipulating tight-binding models. Its distinct feature is its comprehensive support for **symmetry operations** and seamless integration with **Wannier90**. It allows users to read Wannier Hamiltonians, symmetrize them to enforce crystal symmetries, and export them for topological invariant calculations, making it a critical tool in the workflow for identifying topological materials.

**Scientific domain**: Topological Materials, Symmetry Analysis
**Target user community**: Researchers using Wannier90 for topological characterization

## Theoretical Methods
- **Tight-Binding**: Evaluation of $H(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}} H(\mathbf{R})$.
- **Symmetry Analysis**: 
  - Application of space group operations to effective models.
  - Symmetrization: $\tilde{H} = \frac{1}{|G|} \sum_{g \in G} D(g) H(g^{-1} \mathbf{k}) D^\dagger(g)$.
- **Interpolation**: Fourier transform from real-space Wannier representation to k-space.

## Capabilities
- **Model Construction**:
  - Direct import from `wannier90_hr.dat`.
  - Manual construction from hopping terms.
  - HDF5 storage for efficient I/O.
- **Manipulation**:
  - Supercell creation.
  - Dimensionality reduction (slab cutting).
  - Basis transformation.
- **Observables**:
  - Band structures.
  - Eigenvalues/Eigenvectors along paths.
- **Symmetry**: Extensive tools to check and enforce symmetries on numerical models.

## Key Strengths
- **Wannier90 Integration**: The `from_wannier_files()` method is robust and handles the `wsvec.dat` (Wigner-Seitz vectors) correctly, solving phase ambiguity issues common in other parsers.
- **Topology Ready**: Directly feeds into **Z2Pack** for calculating Wilson loops and Chern numbers.
- **Sparse Storage**: Efficiently handles large models with many hopping terms using sparse matrix logic.

## Inputs & Outputs
- **Inputs**:
  - Wannier90 output files (`_hr.dat`, `_wsvec.dat`).
  - Python scripts.
- **Outputs**:
  - HDF5 model files.
  - Band structure arrays.

## Interfaces & Ecosystem
- **Z2Pack**: Primary consumer of TBmodels objects.
- **Wannier90**: Primary source of input data.
- **PythTB**: TBmodels can convert to/from PythTB formats.

## Performance Characteristics
- **Efficiency**: Optimized for evaluating Hamiltonians at many k-points (vectorized via NumPy).
- **Scalability**: Handles models with hundreds of orbitals (standard Wannier output).

## Comparison with Other Codes
- **vs. PythTB**: PythTB is excellent for pedagogy and simple models. TBmodels is built for the "read Wannier90 $\to$ symmetrize $\to$ calculate Z2" research pipeline.
- **vs. Pybinding**: Pybinding is better for massive disordered systems (KPM). TBmodels is better for accurate bulk topology of crystals.

## Application Areas
- **Material Discovery**: Screening databases for topological insulators.
- **Surface States**: Constructing slab models from bulk Wannier functions to see surface bands.

## Community and Support
- **Development**: Dominik Gresch (ETH Zurich alumni).
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://tbmodels.greschd.ch/](https://tbmodels.greschd.ch/)
- **Verification status**: ✅ VERIFIED
  - Mature library in the topological physics stack.
