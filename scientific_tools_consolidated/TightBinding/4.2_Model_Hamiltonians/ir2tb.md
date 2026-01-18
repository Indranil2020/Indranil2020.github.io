# ir2tb (Irreducible Representation to Tight-Binding)

## Official Resources
- **Repository**: https://github.com/zjwang11/ir2tb
- **License**: Unspecified (Research Code)

## Overview
**ir2tb** is a specialized Python tool developed by Zhijun Wang (author of PyWannier90 and IRVSP) for constructing **symmetry-adapted tight-binding models** directly from irreducible representations (irreps) of crystallographic space groups. It automates the complex group-theoretical task of constraining the Hamiltonian matrix elements so that they respect the full symmetry of the crystal, enabling the construction of minimal effective models for topological analysis.

**Scientific domain**: Group Theory, Model Construction
**Target user community**: Theorists building k·p or minimal lattice models

## Theoretical Methods
- **Group Theory**: Analysis of Little Groups and their Irreducible Representations.
- **Method of Invariants**: Identifying linear combinations of operators that transform trivially under the group operations.
- **Basis Construction**: Generating symmetry-adapted basis functions (orbitals).
- **Tight-Binding**: Constraining onsite energies and hopping amplitudes `t_{ij}` based on symmetry.

## Capabilities
- **Model Generation**:
  - Input: Space group, Wyckoff positions, Irreps at high-symmetry points.
  - Output: Symmetry-constrained Hamiltonian matrix.
- **Verification**: Can check if a given numerical Hamiltonian breaks any crystal symmetries.
- **Topological Analysis**: Ensures models used for topological classification have the correct symmetry eigenvalues.

## Key Strengths
- **Rigorous Symmetrization**: Avoids human error in manually deriving Hamiltonian constraints, which is notoriously difficult for non-symmorphic space groups.
- **Minimal Models**: Helps construct the smallest possible model (fewest bands) that captures the essential topology/physics near the Fermi level.
- **Integration**: Compalementary to `IRVSP` (Irreducible Representations of VASP/WIEN2k), creating a workflow from DFT irreps $\to$ TB model.

## Inputs & Outputs
- **Inputs**:
  - Python script specifying the symmetry group and orbitals.
- **Outputs**:
  - Symbolic or numerical Hamiltonian matrices.

## Interfaces & Ecosystem
- **IRVSP**: Often used after determining irreps from DFT using the IRVSP code.
- **Python**: Integrates with NumPy/SymPy.

## Performance Characteristics
- **Efficiency**: Very fast (algebraic operations). Scaling is determined by the size of the group and the number of basis states.
- **Constraint**: Primarily a "generator" tool; the heavy lifting is done in solving the resulting model.

## Comparison with Other Codes
- **vs. TBmodels**: TBmodels symmetrizes *existing* numerical models (from Wannier90). ir2tb constructs the *form* of the model from group theory principles.
- **vs. MagneticTB/MagneticKP**: Similar concept but ir2tb focusses on non-magnetic space groups and irreps, whereas the others focus on Magnetic Space Groups.

## Application Areas
- **Topological Semimetals**: Constructing effective models for Dirac/Weyl points protected by symmetry.
- **Nodal Lines**: Studying symmetry-protected band crossings.

## Community and Support
- **Development**: Zhijun Wang (IOP CAS).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/zjwang11/ir2tb](https://github.com/zjwang11/ir2tb)
- **Verification status**: ✅ VERIFIED
  - Valid research tool from a reputable developer in the topological materials community.
