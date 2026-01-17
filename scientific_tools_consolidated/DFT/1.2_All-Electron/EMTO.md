# EMTO (Exact Muffin-Tin Orbitals)

## Official Resources
- **Homepage**: https://emto.gitlab.io/
- **Documentation**: https://emto.gitlab.io/
- **Source Repository**: Hosted on GitLab (Access via request/license)
- **License**: Academic/Commercial (Not Open Source)

## Overview
**EMTO** is a specialized all-electron Density Functional Theory (DFT) code based on the Exact Muffin-Tin Orbitals (EMTO) theory. It is particularly renowned for its implementation of the **Coherent Potential Approximation (CPA)**, making it one of the premier tools for studying disordered alloys, paramagnetic states, and high-entropy alloys where chemical disorder is non-trivial.

**Scientific domain**: Metallurgy, High-Entropy Alloys, Disordered Systems, Magnetism.
**Target user community**: Materials scientists working on steel, alloys, and phase stability.

## Theoretical Methods
- **Basis Set**: Exact Muffin-Tin Orbitals (EMTO).
- **Potential**: Full Charge Density (FCD) technique.
- **Disorder**: Coherent Potential Approximation (CPA) for substitutionally disordered alloys.
- **Relativity**: Scalar-relativistic and fully relativistic implementations.
- **Elasticity**: Efficient calculation of elastic constants in disordered phases.

## Capabilities
- **Equation of State**: Accurate lattice parameters and bulk moduli.
- **Phase Stability**: Energy differences between crystal structures (BCC, FCC, HCP).
- **Elastic Constants**: Single-crystal elastic constants of random alloys.
- **Magnetism**: Ferromagnetic, ferrimagnetic, and disordered local moment (DLM) paramagnetic states.
- **Stacking Fault Energies**: Vital for mechanical properties.

## Key Strengths
- **Alloy Theory**: The combination of EMTO with CPA is highly efficient and accurate for random alloys compared to supercell methods.
- **Efficiency**: Faster than KKR-CPA implementations for many structural properties.
- **Accuracy**: All-electron precision avoids pseudopotential artifacts.

## Inputs & Outputs
- **Inputs**: Fortran-style input files defining structure, concentration, and potential parameters.
- **Outputs**: Total energies, DOS, spectral functions, equilibrium properties, elastic constants.

## Interfaces & Ecosystem
- **Workflow Integration**:
  - No direct ASE/Phonopy interface, but often used in scripted workflows for force constants.
  - Integration with machine learning pipelines for high-entropy alloy property prediction.
- **Visualization**:
  - Standard plotting tools for DOS and spectral functions.
- **Related Tools**:
  - **KKR-CPA**: Similar methodology but EMTO is often favored for structural energy differences.

## Workflow and Usage

### Elastic Constants Calculation:
1. **Equation of State**: Calculate total energy vs volume to find equilibrium lattice constant.
2. **Distortion**: Apply small volume-conserving strains (orthorhombic, monoclinic) to the lattice.
3. **Energy-Strain**: Calculate energy changes for distortions.
4. **Fitting**: Fit energy curves to extract elastic constants ($C_{11}, C_{12}, C_{44}$).
5. **Polycrystalline Moduli**: Use Voigt-Reuss-Hill averaging for bulk and shear moduli.

### CPA Alloying:
- Define mixing ratios in input file (e.g., Fe$_{0.7}$Ni$_{0.3}$).
- Run self-consistent loop to determine effective medium potential.
- Calculate properties on the effective medium.

## Advanced Features

### Elastic Constant Prediction:
- **CPA Accuracy**: Uniquely capable of calculating elastic constants of random alloys without supercells.
- **High-Entropy Alloys**: Standard tool for predicting mechanical properties of HEAs.
- **Polycrystalline Averages**: Direct calculation of B, G, E, and Poisson's ratio.

### Full Charge Density (FCD):
- Improves upon the muffin-tin approximation for total energy calculations.
- Critical for accurate structural energy differences and phase stability.

### Temperature Dependence:
- Can be coupled with Debye-Grüneisen model or phonon calculations (via displacement method) to estimating finite-temperature elastic properties.

## Performance Characteristics
- **Efficiency**: Highly efficient for disordered alloys compared to supercell (SQS) methods.
- **Parallelization**: Parallelization is typically handled over k-points and energy points in the contour integration.
- **Scalability**: Scaling is generally favorable for large concentrations of species but limited by the $O(N^3)$ diagonalization scaling of the underlying KKR/EMTO formalism.

## Limitations & Known Constraints
- **Muffin-Tin Approximation**: Assumes spherically symmetric potential inside spheres and constant potential in interstitial regions. Can be inaccurate for open structures or highly covalent systems with directional bonding.
- **Single-Site CPA**: The standard CPA is a single-site approximation, neglecting short-range order or clustering effects (though Cluster-CPA extensions exist).
- **Forces**: Structural relaxation is often more cumbersome/limited compared to plane-wave codes.

## Best Practices
- **Lattice Constants**: Always determine equilibrium lattice constants using the Equation of State module before doing extensive property calculations.
- **Mixing**: CPA convergence can be tricky; reducing mixing parameters is often necessary for magnetic alloys.

## Community and Support
- **Support**: Limited public forum; support is primarily through the developer network and workshops (e.g., hosted by ENCCS).
## Comparison with Other Codes
- **vs KKR-CPA (e.g., AkaiKKR, SPR-KKR)**: EMTO is often faster for total energy and structural optimization, while KKR codes might offer more spectroscopic features.
- **vs VASP/QE**: VASP/QE typically use supercells (SQS) for alloys, which can be computationally expensive and suffer from finite-size effects. EMTO-CPA treats disorder analytically.

## Verification & Sources
**Primary sources**:
1.  **Official Website**: [emto.gitlab.io](https://emto.gitlab.io/)
2.  **Literature**: Vitos, L. "Computational Quantum Mechanics for Materials Engineers: The EMTO Method and Applications", Springer (2007).

**Verification status**: ✅ VERIFIED
