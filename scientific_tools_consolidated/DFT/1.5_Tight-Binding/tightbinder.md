# tightbinder

## Official Resources
- Homepage: https://github.com/alejandrojuria/tightbinder
- Source Repository: https://github.com/alejandrojuria/tightbinder
- License: MIT License

## Overview
tightbinder is a versatile Python library designed for the modeling of solid-state systems using the Tight-Binding (TB) approximation. It provides a comprehensive suite of tools for constructing Hamiltonians, solving for electronic properties, and analyzing topological features. Its emphasis on a clean, object-oriented API makes it an excellent tool for rapid prototyping of model Hamiltonians and exploring concepts in topological insulators and condensed matter physics.

**Scientific domain**: Model Hamiltonians, Condensed Matter Theory, Topology
**Target user community**: Theoretical physicists, Students, Educators

## Theoretical Methods
- **Tight-Binding Approximation**: Orthogonal basis sets.
- **Slater-Koster Formalism**: Automatic generation of hopping integrals from geometric relations.
- **Topological Invariants**: Wilson loops, Berry curvature, Chern numbers.
- **Green's Functions**: Surface spectral functions (iterative Green's function method).
- **Supercells & Disorder**: Handling of finite size effects and substitutions.

## Capabilities (CRITICAL)
- **System Construction**: Built-in support for Bravais lattices, multi-atom bases, and supercells.
- **Hamiltonian Generation**: From Slater-Koster tables or manual hopping lists.
- **Electronic Structure**: Bands, DOS, Local DOS.
- **Topology**: Calculation of Z2 invariants, Berry phases, and Edge states.
- **Spectral Functions**: Surface states visualization.
- **Visualization**: Integrated plotting of lattice structures and band diagrams.

## Key Strengths

### Topological Toolkit:
- Native implementation of modern topological characterization tools (Wilson loops, Zak phase) which are often missing in general TB codes.

### Pythonic Design:
- Logical object hierarchy (`System`, `Lattice`, `Hamiltonian`).
- Seamless integration with the SciPy stack (NumPy, Matplotlib).

### Slater-Koster Engine:
- Can read standard SK tables and automatically construct Hamiltonians for deformed lattices, enabling study of strain.

## Inputs & Outputs
- **Inputs**:
  - Python scripts defining the lattice vectors and atoms.
  - Dictionary of Slater-Koster parameters (`{'ss_sigma': -1.0, ...}`).
- **Outputs**:
  - Matplotlib figures (Bands, DOS).
  - Raw NumPy arrays of eigenvalues/vectors.
  - Serialized system objects.

## Interfaces & Ecosystem
- **Python**: purely Python library.
- **PythTB**: Similar functionality, but tightbinder adds more recent topological tools.
- **Plotting**: High-quality default plots for publication.

## Advanced Features
- **Hofstadter Butterfly**: Handling of magnetic fields in 2D lattices.
- **Disorder**: Methods to introduce Anderson disorder (random onsite potentials).
- **Unfolding**: Band unfolding for supercells.

## Performance Characteristics
- **Speed**: Hamiltonian construction is vectorized; diagonalization uses LAPACK (via NumPy). Efficient for 1D/2D models and reasonable 3D meshes.
- **Scalability**: Not an HPC code; limited by single-node memory for ED. Green's functions allow semi-infinite systems.

## Computational Cost
- **Low**: Run on laptops for typical model systems (up to thousands of orbitals).

## Limitations & Known Constraints
- **Non-Self-Consistent**: Does not solve Poisson equation; purely model Hamiltonian.
- **Basis limitation**: Primarily orthogonal basis sets.

## Comparison with Other Codes
- **vs PythTB**: PythTB is the "classic" Python TB code; tightbinder offers a more modern API and advanced topological features (Wilson loops).
- **vs Kwant**: Kwant is specialized for scattering (transport); tightbinder is specialized for spectral/topological properties of bulk/slabs.
- **vs TBmodels**: TBmodels focuses on ab-initio downfolding; tightbinder focuses on Slater-Koster construction.
- **Unique strength**: Excellent educational and research tool for topological systems with SK parameterization.

## Application Areas
- **Topological Insulators**: Modeling edge states in SSH or Haldane models.
- **Graphene Physics**: Strain effects in honeycomb lattices.
- **Education**: Teaching band theory and topology interactively in Jupyter notebooks.

## Best Practices
- **Use Jupyter**: The visualization tools are designed for notebook environments.
- **Verify Symmetries**: When defining custom lattices, ensure symmetries are preserved.
- **K-paths**: Use the built-in path generators for high-symmetry lines.

## Community and Support
- **GitHub**: Active repository.
- **Documentation**: Docstrings and examples provided in repo.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/alejandrojuria/tightbinder

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Functionality: Confirmed features via code inspection.
