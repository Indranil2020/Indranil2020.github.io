# MagnetoTransport.jl

## Official Resources
- **Repository**: https://github.com/ChristopheBerthod/MagnetoTransport.jl
- **License**: MIT License

## Overview
**MagnetoTransport.jl** is a specialized Julia package for calculating the linear magneto-transport properties of **two-dimensional quantum systems**. It implements the **Kubo-Bastin formalism** to evaluate the longitudinal ($\sigma_{xx}$) and transverse/Hall ($\sigma_{xy}$) conductivity tensors. Unlike general transport codes, it specifically targets the accurate computation of orbital magnetic effects, including Hofstadter butterfly physics and quantum Hall plateaus, taking into account energy-dependent self-energies.

**Scientific domain**: Quantum Hall Effect, Topological Matter, 2D Materials
**Target user community**: Theorists investigating topological transport in lattice models

## Theoretical Methods
- **Kubo-Bastin Formula**: A formulation of the Kubo formula that is robust for dissipative systems and captures both Fermi sea and Fermi surface contributions.
- **Streda Formula**: Used for the thermodynamic part of the Hall conductivity.
- **Peierls Substitution**: Models the orbital magnetic field effects on lattice hoppings.
- **Self-Energy**: Incorporates disorder/scattering via a complex self-energy $\Sigma(E)$.

## Capabilities
- **Observables**:
  - Longitudinal Conductivity ($\sigma_{xx}$).
  - Hall Conductivity ($\sigma_{xy}$, quantized in units of $e^2/h$).
  - Density of States (DOS).
  - Integrated Berry curvature.
- **Models**:
  - Square, Honeycomb (Graphene), and Triangular lattices.
  - User-defined tight-binding Hamiltonians.

## Key Strengths
- **Julia Efficiency**: Leverages Julia's JIT compilation for fast numerical integration over the Brillouin zone.
- **Precision**: Capable of resolving fine spectral features (Landau levels) using high-resolution grids or adaptive integration (`QuadGK`).
- **Focus**: Dedicated solely to the conductivity tensor, ensuring correct implementation of the often-tricky Kubo terms.

## Inputs & Outputs
- **Inputs**:
  - Julia scripts defining the lattice parameters ($a, t$).
  - Magnetic flux per plaquette ($\phi$).
  - Scattering rate ($\eta$ or $\Sigma$).
- **Outputs**:
  - Conductivity tensors as a function of Fermi energy or Chemical potential.
  - Plots of $\sigma_{xy}$ showing integer quantization.

## Interfaces & Ecosystem
- **Julia**: Fully integrated with the Julia scientific stack (`LinearAlgebra`, `Plots`).
- **Quantica.jl**: Can potentially use Hamiltonians constructed in Quantica (with conversion).

## Performance Characteristics
- **Speed**: Vectorsized operations over k-space make it highly efficient for 2D bulk systems.
- **Scalability**: Limited to effective single-particle models (non-interacting or mean-field).

## Comparisons with Other Codes
- **vs. Kwant**: Kwant calculates Hall conductance for *finite* bars using leads. MagnetoTransport.jl calculates the *bulk* Hall conductivity tensor using the Kubo formula.
- **vs. LinReTraCe**: Both usage Kubo formulas; MagnetoTransport.jl is Julia-based and specialized for 2D magnetic response, while LinReTraCe is C++ and more general (3D, thermoelectrics).

## Community and Support
- **Development**: Christophe Berthod (University of Geneva).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/ChristopheBerthod/MagnetoTransport.jl](https://github.com/ChristopheBerthod/MagnetoTransport.jl)
- **Verification status**: ✅ VERIFIED
  - Functional academic package.
