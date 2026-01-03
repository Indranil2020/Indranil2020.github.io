# Magnon codes (SpinW, VAMPIRE, Spirit)

## Official Resources
- Homepage: Varies (SpinW, VAMPIRE, Spirit homepages)
- Documentation: Varies
- Source Repository: Varies
- License: Varies

## Overview
This entry serves as a category placeholder for magnon simulation and spin dynamics codes. The primary modern tools for magnon calculations are **Spirit**, **VAMPIRE**, **SpinW**, and **UppASD**. These codes typically perform atomistic spin dynamics (ASD) simulations or linear spin wave theory (LSWT) calculations to determine magnon dispersion relations, lifetimes, and thermodynamic magnetic properties.

**Scientific domain**: Magnetism, spin dynamics, magnons, spintronics  
**Target user community**: Magnetism researchers, spintronics engineers

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Atomistic Spin Dynamics (ASD)
- Linear Spin Wave Theory (LSWT)
- Monte Carlo simulations (for magnetic phase transitions)
- Heisenberg Hamiltonian with anisotropy and DMI
- Geodesic Nudged Elastic Band (GNEB) for magnetic transitions

## Capabilities (CRITICAL)
- **SpinW**: MATLAB/Python code for linear spin wave theory, fitting experimental neutron scattering data.
- **Spirit**: Framework for spin dynamics and transition state finding (GNEB), visualizations.
- **VAMPIRE**: High-performance atomistic spin dynamics, temperature dependence, recording media.
- **UppASD**: Atomistic spin dynamics, thermodynamics, magnon dispersion.
- **Calculation of**: Magnon dispersion relations, magnetic susceptibility, Curie temperature, hysteresis loops, skyrmion stability.

**Sources**: Respective tool documentations

## Inputs & Outputs
- **Input**: Crystal structure, magnetic exchange interactions (J_ij), anisotropy (K), DMI (D)
- **Output**: Spin configurations, magnon spectra (S(q,w)), magnetization vs temperature

## Interfaces & Ecosystem
- **TB2J**: Calculates parameters (J, D, K) from DFT for these codes
- **DFT**: VASP, QE, etc. provide magnetic moments and energies
- **Visualization**: ParaView, POV-Ray, built-in visualizers

## Application Areas
- Magnonics and spintronics
- Skyrmions and topological magnetism
- Magnetic storage media
- Frustrated magnetism
- Neutron scattering data analysis

## Limitations & Known Constraints
[TO BE COMPLETED - Requires official documentation review]

## Verification & Sources
**Primary sources**:
1. Spirit: https://spirit-code.github.io/
2. VAMPIRE: https://vampire.york.ac.uk/
3. SpinW: https://spinw.org/
4. UppASD: https://github.com/UppASD/UppASD

**Secondary sources**: [TO BE VERIFIED]

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Category: Placeholder for verified codes
- Recommendation: See individual entries for **Spirit**, **VAMPIRE**, **TB2J**. Use **SpinW** for experimental fitting.
