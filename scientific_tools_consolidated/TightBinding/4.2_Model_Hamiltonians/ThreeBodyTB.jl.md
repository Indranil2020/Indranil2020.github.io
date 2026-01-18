# ThreeBodyTB.jl

## Official Resources
- **Homepage**: https://pages.nist.gov/ThreeBodyTB.jl/
- **Repository**: https://github.com/usnistgov/ThreeBodyTB.jl
- **License**: MIT License

## Overview
**ThreeBodyTB.jl** is a high-accuracy tight-binding package developed by NIST. It distinguishes itself from standard Slater-Koster codes by including pre-fit **three-body interaction terms**, which dramatically improves transferability and accuracy for structures far from equilibrium (e.g., surfaces, defects, high pressure). Implemented in pure Julia, it provides a self-consistent field (SCF) solver that rivals DFT accuracy (specifically PBEsol) for a vast range of elemental and binary systems, but at a fraction of the computational cost.

**Scientific domain**: Materials Science, High-Throughput Screening
**Target user community**: Researchers constructing phase diagrams or running large-scale MD

## Theoretical Methods
- **Tight-Binding with 3-Body terms (TB3)**: $H = H_{2body} + H_{3body}$. The 3-body terms capture environmental dependence of bonding.
- **Self-Consistent Field (SCF)**: Solves for charge redistribution/transfer, essential for ionic systems and surfaces.
- **Parameterization**:
  - Uses a massive pre-computed database of parameters fitted to DFT data.
  - Covers >99% of ICSD prototypes for supported elements.
- **Magnetic Moments**: Collinear spin-polarized calculations.

## Capabilities
- **Simulations**:
  - SCF Ground state energy and forces.
  - Band structures and Density of States (DOS).
  - Molecular Dynamics (MD) trajectories.
  - Phonon spectra (via finite displacement).
- **System Support**:
  - Bulk crystals (Metals, Insulators, Semiconductors).
  - Low-dimensional systems (Surfaces, 2D materials).
  - Charged defects.
- **Ease of Use**: "Automatic" mode that guesses initial parameters and symmetry.

## Key Strengths
- **Accuracy**: Benchmarks show energy/force errors comparable to DFT-PBEsol, far superior to traditional non-SCF tight-binding.
- **Database**: Comes "batteries included" with parameters for most common elements, removing the need for users to perform their own fitting.
- **Speed**: $O(N)$ sparse matrix operations allow routine simulation of 1000+ atom supercells.

## Inputs & Outputs
- **Inputs**:
  - Crystal structure (POSCAR, CIF).
  - List of elements.
- **Outputs**:
  - Energies, Forces, Stress tensor.
  - Band plots.

## Interfaces & Ecosystem
- **Input Generation**: Uniquely integrated with `CrystalStructure.jl` logic (internal) for handling symmetry.
- **Visualization**: Built-in plotting recipes.

## Performance Characteristics
- **Efficiency**: Exploits Julia's Type system and SIMD; heavily optimized sparse solvers.
- **Parallelism**: Multi-threaded execution.

## Comparison with Other Codes
- **vs. DFTB+**: Both aim for "DFT quality" tight-binding. ThreeBodyTB.jl's unique selling point is the explicit 3-body term (better for structural relaxation) and its modern, hackable Julia codebase.
- **vs. xTB**: xTB is semi-empirical and great for molecules/organics. ThreeBodyTB.jl is parameterized specifically for solid-state crystals and materials science.

## Application Areas
- **Phase Stability**: Calculating formation energies of competing crystal polymorphs.
- **Defects**: Simulating large supercells to study vacancy/interstitial formation energies without finite-size errors.

## Community and Support
- **Development**: Kevin Garrity (NIST).
- **Source**: GitHub.

## Verification & Sources
- **Website**: [https://pages.nist.gov/ThreeBodyTB.jl/](https://pages.nist.gov/ThreeBodyTB.jl/)
- **Verification status**: ✅ VERIFIED
  - NIST-developed standard tool.
