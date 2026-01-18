# Smeagol

## Official Resources
- **Homepage**: http://www.tcd.ie/Physics/Smeagol/
- **License**: Academic / Open Source (Requires SIESTA license/install)

## Overview
**Smeagol** (Spin and Molecular Electronics Algorithm on a Generalized atomic Orbital Landscape) is a leading computational package for simulating **quantum transport** in nanoscale devices using the **Non-Equilibrium Green's Function (NEGF)** formalism combined with **Density Functional Theory (DFT)**. It is specifically designed to handle spin-polarized transport, non-collinear magnetism, and strong electron correlations in two-terminal devices under finite bias voltage.

**Scientific domain**: Spintronics, Molecular Electronics, Nanotechnology
**Target user community**: Researchers simulating magnetic tunnel junctions, molecular spin valves, and atomic wires

## Theoretical Methods
- **DFT-NEGF**: Combines the Kohn-Sham Hamiltonian (from SIESTA) with the NEGF formalism to solve for the density matrix under non-equilibrium conditions.
- **Basis Set**: Uses Numerical Atomic Orbitals (NAOs) provided by SIESTA, which are efficient for linear scaling.
- **Spin**: Full treatment of collinear and non-collinear spin, including Spin-Orbit Coupling (SOC).
- **Correlations**: Implements LDA+U and Self-Interaction Correction (SIC) for localized d- and f-electrons.

## Capabilities
- **Transport Observables**:
  - Current-Voltage (I-V) characteristics.
  - Magnetoresistance (TMR/GMR ratios).
  - Spin-Transfer Torque (STT).
  - Transmission coefficients $T(E, V)$.
- **System Types**:
  - Magnetic Tunnel Junctions (Fe/MgO/Fe).
  - Single-Molecule Magnets (SMMs) in junctions.
  - Metallic atomic chains.
  - Graphene spintronic devices.

## Key Strengths
- **Spintronics**: Arguably the most advanced code available for ab initio spintronic transport (STT, non-collinear textures).
- **Self-Consistency**: Fully self-consistent charge and spin density under bias.
- **Accuracy**: Inherits the efficiency and accuracy of SIESTA's NAO basis.

## Inputs & Outputs
- **Inputs**:
  - SIESTA-style `.fdf` input files.
  - Pseudopotentials (`.vps`).
  - Lead Hamiltonians (`.bulk`).
- **Outputs**:
  - `Results.dat`: Current and Transmission data.
  - `Density.DM`: Non-equilibrium density matrix.
  - `Orbitalmom.dat`: Orbital magnetic moments.

## Interfaces & Ecosystem
- **SIESTA**: Tightly integrated; Smeagol effectively replaces the SIESTA executable for transport runs.
- **Wannier90**: Can be used via downfolding tools (though native SIESTA usage is standard).

## Performance Characteristics
- **Computational Cost**: More expensive than equilibrium codes (Gollum) due to the SCF cycle at each bias point. $O(N^3)$ scaling.
- **Parallelism**: Efficient MPI parallelization over energy contour integration (dozens to hundreds of cores).

## Comparison with Other Codes
- **vs. TranSIESTA**: Both use SIESTA + NEGF. Smeagol historically had stronger support for non-collinear spin and corrections (SIC/LDA+U), while TranSIESTA is the official SIESTA module.
- **vs. QuantumATK**: Comparable capabilities; Smeagol is open-source (academic), whereas ATK is commercial.

## Application Areas
- **MRAM**: Optimizing tunnel barriers for magnetic random access memory.
- **Molecular Spintronics**: Filtering spin through organic molecules.
- **STM**: Simulating scanning tunneling microscopy images with magnetic tips.

## Community and Support
- **Development**: Sanvito Group (Trinity College Dublin).
- **Source**: Available via request/download form.

## Verification & Sources
- **Website**: [http://www.tcd.ie/Physics/Smeagol/](http://www.tcd.ie/Physics/Smeagol/)
- **Primary Publication**: A. R. Rocha et al., Phys. Rev. B 73, 085414 (2006).
- **Verification status**: ✅ VERIFIED
  - A major code in the spintronics community for over 15 years.
