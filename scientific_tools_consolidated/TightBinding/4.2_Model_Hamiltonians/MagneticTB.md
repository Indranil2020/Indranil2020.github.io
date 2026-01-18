# MagneticTB

## Official Resources
- **Homepage**: https://github.com/zhangzeyingvv/MagneticTB
- **Repository**: https://github.com/zhangzeyingvv/MagneticTB
- **License**: GPL-3.0

## Overview
**MagneticTB** is a Mathematica package designed for the automated construction of tight-binding models for materials with any of the 1651 **Magnetic Space Groups (MSGs)**. It greatly simplifies the theoretical modeling of magnetic topological materials by automatically generating symmetry-allowed Hamiltonian matrices based on user-supplied Wyckoff positions and orbital characters.

**Scientific domain**: Topological Magnetism, Symmetry Analysis
**Target user community**: Theorists building effective models for magnetic systems

## Theoretical Methods
- **Magnetic Space Group Theory**: Full implementation of the representation theory for all 1651 magnetic space groups (Type I, II, III, IV).
- **Symmetry-Adapted Models**: Construction of invariant Hamiltonian terms.
- **Corepresentations**: Handling of anti-unitary symmetries (Time-Reversal * Geometric operation).
- **Wyckoff Positions**: Input atoms based on standard crystallographic sites.

## Capabilities
- **Model Generation**:
  - Spinless and Spinful models.
  - Inclusion of Spin-Orbit Coupling (SOC).
- **Input**:
  - MSG number (BNS or OG setting).
  - Orbitals (s, p, d, f) on specific sites.
- **Output**:
  - Parameterized Hamiltonian matrix $H(\mathbf{k})$.
  - Relationships between hopping parameters enforced by symmetry.
- **Interoperability**: Can export models for use in other codes or numerical analysis.

## Key Strengths
- **Automation**: Removes the tedious and error-prone process of manually deriving symmetry constraints for complex magnetic Unit cells.
- **Completeness**: Supports the full table of magnetic space groups, not just the non-magnetic ones.
- **Integration**: Works within Mathematica, allowing immediate symbolic manipulation of the result.

## Inputs & Outputs
- **Inputs**: Mathematica function calls specifying the group and orbitals.
- **Outputs**: Symbolic matrices and rules for parameters.

## Interfaces & Ecosystem
- **MagneticKP**: Companion package for k·p models.
- **WannierTools**: Models can be exported for topological surface state calculation.

## Performance Characteristics
- **Speed**: Symbolic generation is generally fast (seconds to minutes) for standard unit cells.
- **Scaling**: Capability depends on the number of orbitals; very large bases might slow down symbolic algebra.

## Comparison with Other Codes
- **vs. Pybinding**: Pybinding builds numerical models for known Hamiltonians. MagneticTB *derives* the Hamiltonian from symmetry.
- **vs. TB2J**: TB2J extracts magnetic parameters from DFT. MagneticTB constructs the model form from symmetry principles.

## Application Areas
- **Antiferromagnetic Spintronics**: Modeling collinear and non-collinear antiferromagnets.
- **Chiral Magnets**: Studying topological band crossings in magnetic systems.

## Community and Support
- **Development**: Institute of Physics, CAS (Zeying Zhang).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/zhangzeyingvv/MagneticTB](https://github.com/zhangzeyingvv/MagneticTB)
- **Primary Publication**: Z. Zhang et al., Comput. Phys. Comm. 270, 108153 (2022).
- **Verification status**: ✅ VERIFIED
  - Active research code.
