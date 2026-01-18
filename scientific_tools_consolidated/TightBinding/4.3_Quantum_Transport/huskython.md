# huskython

## Official Resources
- **Repository**: https://github.com/NicoRenaud/huskython
- **License**: Open Source

## Overview
**huskython** is a Python-based code for simulating quantum transport in molecular junctions. It is notable for implementing two complementary formalisms: the **Equivalent Scattering-State Quantum Conductance (ESQC)** method and the standard **Non-Equilibrium Green's Function (NEGF)** method. This dual approach allows users to calculate transport properties from both a scattering state perspective (wavefunction matching) and a Green's function perspective.

**Scientific domain**: Molecular Electronics, Chemical Physics
**Target user community**: Researchers studying electron transfer in organic molecules

## Theoretical Methods
- **ESQC (Equivalent Scattering-State Quantum Conductance)**:
  - Treats the molecule as a scattering center.
  - Solving the Schrödinger equation for scattering states $\Psi_S$ propagating from leads.
  - Direct calculation of the $S$-matrix.
- **NEGF**:
  - Calculation of Retarded/Advanced Green's functions.
  - Transmission via the Fisher-Lee relation $T = \text{Tr}[\Gamma_L G \Gamma_R G^\dagger]$.
- **Electronic Structure**:
  - Semi-empirical **Extended Hückel** method (EHT) for Hamiltonian construction.

## Capabilities
- **Observables**:
  - Zero-bias Conductance.
  - Transmission Spectra $T(E)$.
  - Scattering wavefunctions (visualization of transmission channels).
- **Systems**:
  - Single molecules (alkanes, benzene rings) bridging metal electrodes.
  - Constructive/Destructive interference mapping.

## Key Strengths
- **Scattering Insight**: The ESQC method provides direct access to the scattering wavefunctions, offering intuitive pictures of how electrons traverse the molecule (e.g., through $\sigma$ vs $\pi$ systems).
- **Comparison**: Unique ability to benchmark ESQC results directly against NEGF within the same code.
- **Pythonic**: Easy to script and integrate with other Python chemical tools.

## Inputs & Outputs
- **Inputs**:
  - Molecular geometry (XYZ).
  - Extended Hückel parameters.
- **Outputs**:
  - Transmission data.
  - Wavefunction coefficients.

## Interfaces & Ecosystem
- **Dependencies**: Standard SciPy stack.
- **Chemistry**: Can work with geometries from RDKit or OpenBabel.

## Performance Characteristics
- **Efficiency**: ESQC can be numerically efficient for zero-bias conductance as it avoids full matrix inversion at every energy point in the same way as RGF, solving instead a linear system for boundary conditions.
- **Parallelism**: Serial execution (typical for single-molecule model codes).

## Comparison with Other Codes
- **vs. Gollum**: Gollum interacts with DFT codes (SIESTA); huskython is self-contained with semi-empirical EHT.
- **vs. Kwant**: Kwant is a general tight-binding solver; huskython is specialized for molecular chemistry (orbitals, chemical species).

## Application Areas
- **Molecular Wires**: Length dependence of conductance (beta factor).
- **Interference**: Destructive quantum interference in meta-substituted benzene.

## Community and Support
- **Development**: Nico Renaud (Netherlands eScience Center / TU Delft).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/NicoRenaud/huskython](https://github.com/NicoRenaud/huskython)
- **Verification status**: ✅ VERIFIED
  - Research code for specific scattering methods.
