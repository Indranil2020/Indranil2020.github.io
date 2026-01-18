# MagneticKP

## Official Resources
- **Homepage**: https://github.com/zhangzeyingvv/MagneticKP
- **Repository**: https://github.com/zhangzeyingvv/MagneticKP
- **License**: GPL-3.0

## Overview
**MagneticKP** is a dual-language (Mathematica and Python) software package for the efficient construction of **k·p effective Hamiltonians** in magnetic and non-magnetic crystals. It implements a novel "Iterative Simplification Algorithm" (ISA) to rapidly solve the symmetry constraints imposed by all 1651 **Magnetic Space Groups (MSGs)**. It enables researchers to derive low-energy effective models expanded to arbitrary orders in wavenumber $\mathbf{k}$ for complex topological materials.

**Scientific domain**: Theoretical Condensed Matter, Group Theory
**Target user community**: Theorists characterizing magnetic topological materials

## Theoretical Methods
- **k·p Perturbation Theory**: Construction of effective Hamiltonians near high-symmetry points in the Brillouin Zone.
- **Method of Invariants**: Identifying terms allowed by symmetry.
- **Magnetic Space Groups**: rigorous treatment of unitary and anti-unitary (time-reversal involving) symmetries.
- **Algorithms**:
  - **ISA (Iterative Simplification Algorithm)**: Efficiently reduces the system of symmetry equations.
  - **DDA (Direct-Product Decomposition Algorithm)**: Handles representations.

## Capabilities
- **Model Generation**:
  - Spinless (single-valued) models.
  - Spinful (double-valued) models with Spin-Orbit Coupling.
  - Arbitrary expansion order in $\mathbf{k}$.
- **Symmetry**: Fully compliant with the standard setting of Magnetic Space Groups (BNS setting).
- **Implementations**:
  - **Mathematica**: For symbolic derivation and analytical expressions.
  - **Python**: For integration into numerical workflows and scripting.

## Key Strengths
- **Speed**: The ISA algorithm offers significant speedups over traditional projection methods, especially for high-dimensional representations and high-order expansions.
- **Magnetic Focus**: unmatched support for the full range of magnetic space groups, critical for the study of antiferromagnetic topological insulators and Weyl semimetals.
- **Dual Interface**: Users can choose between symbolic power (Mathematica) or numerical flexibility (Python).

## Inputs & Outputs
- **Inputs**:
  - Magnetic Space Group index/settings.
  - High-symmetry point coordinates.
  - Irreducible Representations (Irreps).
- **Outputs**:
  - Symbolic or numerical forms of the Hamiltonian matrix $H(k)$.
  - Basis matrices satisfying symmetry constraints.

## Interfaces & Ecosystem
- **MagneticTB**: Comparison tool by the same authors for Tight-Binding models.
- **Dependencies**: Mathematica (Wolfram), Python (NumPy, SymPy).

## Performance Characteristics
- **Efficiency**: Can generate high-order models (e.g., 4th order in k) for complex groups in seconds/minutes.
- **Scalability**: Handle large representations ($>10$ dimensions) that might choke standard "brute force" invariant methods.

## Comparison with Other Codes
- **vs. kdotp-symmetry**: An older Mathematica package. MagneticKP uses newer algorithms (ISA) and supports magnetic groups more natively.
- **vs. Qsymm**: Qsymm (Python) is excellent for finding symmetries of a given Hamiltonian. MagneticKP constructs the generic Hamiltonian *from* the symmetries.

## Application Areas
- **Magnetic Topological Insulators**: Deriving surface effective theories.
- **Weyl Semimetals**: Finding the allowed terms protecting Weyl nodes in magnetic structures.

## Community and Support
- **Development**: Institute of Physics, Chinese Academy of Sciences (Zeying Zhang).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/zhangzeyingvv/MagneticKP](https://github.com/zhangzeyingvv/MagneticKP)
- **Primary Publication**: Z. Zhang et al., Comput. Phys. Comm. 283, 108575 (2023).
- **Verification status**: ✅ VERIFIED
  - Valid research tool backed by publication.
