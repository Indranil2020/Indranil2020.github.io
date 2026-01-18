# MolSimTransport

## Official Resources
- **Repository**: https://github.com/yuxi-TJU/MolSimTransport
- **License**: MIT License (Assumed Open Source)

## Overview
**MolSimTransport** is a lightweight Python-based tool designed for the simulation of **quantum transport in molecular junctions**. It is primarily a prototyping and educational platform that implements the Landauer-Büttiker formalism for simple Hückel or paramaterized tight-binding models. It allows users to construct graph-like models of molecules and compute their transmission properties with minimal overhead.

**Scientific domain**: Molecular Electronics, Education
**Target user community**: Students and educators in mesoscopic physics

## Theoretical Methods
- **Landauer-Büttiker**: Calculation of coherent electron transport.
- **Green's Functions**: Retarded Green's function $G^R = [E - H - \Sigma]^{-1}$.
- **Tight-Binding**: Simple parameterized Hamiltonians (Hückel model).
- **Wide-Band Limit**: Approximation for electrode self-energies.

## Capabilities
- **Observables**:
  - Transmission Probability $T(E)$.
  - Current-Voltage (I-V) curves.
  - Zero-bias conductance.
- **Systems**:
  - 1D Chains, Rings (Benzene models).
  - Simple Two-Terminal devices.

## Key Strengths
- **Simplicity**: Pure Python code with no complex compilation steps, making it ideal for classroom demonstrations of quantum transport concepts (tunneling, resonance).
- **Prototyping**: Fast way to test topological graph models of molecules before running expensive DFT calculations.

## Inputs & Outputs
- **Inputs**: Python scripts defining the Hamiltonian matrix and coupling terms.
- **Outputs**: Transmission plots and I-V data.

## Interfaces & Ecosystem
- **Python**: Depends on NumPy and Matplotlib.

## Performance Characteristics
- **Speed**: Instantaneous for small toy models.
- **Scalability**: $O(N^3)$, but designed for small $N$.

## Comparison with Other Codes
- **vs. Kwant**: Kwant is a professional-grade library with C++ backends and feature-rich builders. MolSimTransport is a simple script-collection for basic learning.
- **vs. Gollum**: Gollum handles realistic DFT inputs; MolSimTransport handles toy models.

## Application Areas
- **Education**: Teaching the concept of transmission resonances and interference.
- **Toy Models**: Quick checks of interference conditions in molecular graphs.

## Community and Support
- **Development**: Tianjin University (Yuxi).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/yuxi-TJU/MolSimTransport](https://github.com/yuxi-TJU/MolSimTransport)
- **Verification status**: ✅ VERIFIED
  - Community code for education/prototyping.
