# Gollum

## Official Resources
- **Homepage**: https://www.gollumcode.com/
- **License**: Academic License (Free registration required)

## Overview
**Gollum** is a highly efficient quantum transport simulation code designed to compute the charge, spin, and thermal transport properties of nanoscale devices. It primarily employs the **equilibrium transport theory** (Landauer-Büttiker formalism in the linear response limit) using Green's functions, which makes it significantly faster than fully self-consistent NEGF codes for many applications. It acts as a bridge between ab initio electronic structure codes (DFT) and transport observables.

**Scientific domain**: Molecular Electronics, Spintronics, Thermoelectrics
**Target user community**: Researchers bridging DFT and device physics

## Theoretical Methods
- **Equilibrium Green's Functions (EGF)**: Calculates the retarded Green's function of the scattering region coupled to semi-infinite leads.
- **Landauer-Büttiker Formalism**: Computes transmission $T(E)$ and integrates it to find conductance $G$, Seebeck $S$, and thermal conductance $\kappa$.
- **Hamiltonians**: Works with Tight-Binding models or DFT Hamiltonians (SIESTA, Wannier90).
- **Environment**: Can model gating, simple magnetic fields, and superconducting proximity effects (Bogoliubov-de Gennes).

## Capabilities
- **Transport Coefficients**:
  - Electrical Conductance ($G_0$).
  - Thermopower (Seebeck $S$).
  - Peltier coefficient ($\Pi$).
  - Electronic thermal conductance ($\kappa_e$).
- **Spintronics**:
  - Spin-polarized transmission ($T_{\uparrow}, T_{\downarrow}$).
  - Spin currents and spin-transfer torques.
- **Systems**:
  - Single molecules (break junctions).
  - Nanotubes and Nanoribbons.
  - Multiterminal devices.

## Key Strengths
- **Speed**: By avoiding the self-consistent loop (using the equilibrium Hamiltonian), it can screen thousands of molecular geometries rapidly.
- **Versatility**: Interfaces with almost all major DFT codes via Wannier90 or native inputs (SIESTA).
- **Ease of Use**: Simple input files and extensive examples/tutorials.

## Inputs & Outputs
- **Inputs**:
  - Hamiltonian ($H$) and Overlap ($S$) matrices (from DFT/TB).
  - Lead definitions.
  - `input` file controlling energy range and temperature.
- **Outputs**:
  - `transmission.dat`: Transmission probability vs Energy.
  - `conductance.dat`, `thermopower.dat`: Transport coefficients vs Chemical Potential/Temperature.

## Interfaces & Ecosystem
- **Upstream**:
  - **SIESTA**: Native support for `.HS` files.
  - **Wannier90**: Universal interface via `_hr.dat`.
  - **VASP/QE/ABINIT/OpenMX**: Via Wannier90.
- **Downstream**: Simple text outputs plotting in Gnuplot/Python.

## Performance Characteristics
- **Efficiency**: Calculations are typically limited by matrix inversion $O(N^3)$. Pre-calculation of leads allows fast sweeps over the scattering region.
- **Parallelism**: MPI parallelization over energy points.

## Limitations & Known Constraints
- **Bias**: Primarily for **zero-bias** (linear response) or low-bias; does not capture non-equilibrium Stark effects or charge rearrangement at high bias (unlike Smeagol/Transiesta).
- **Correlations**: Limited to the level of theory of the input DFT (LDA/GGA); no native DMFT transport (though can use DFT+U Hamiltonians).

## Comparison with Other Codes
- **vs. Smeagol/Transiesta**: These are non-equilibrium solvers (finite bias, I-V curves). Gollum is faster but limited to linear response.
- **vs. Artaios**: Both are post-processing tools; Gollum has broader Hamiltonian support (SIESTA native) and focuses more on thermoelectrics/spintronics.

## Application Areas
- **Molecular Junctions**: Conductance histograms of break junction experiments.
- **Thermoic Devices**: Designing quantum dot heat engines.
- **2D Material Contacts**: Calculating contact resistance in Graphene/TMD devices.

## Community and Support
- **Development**: Lancaster University (Colin Lambert group).
- **Source**: Available upon registration.

## Verification & Sources
- **Website**: [https://www.gollumcode.com/](https://www.gollumcode.com/)
- **Primary Publication**: J. Ferrer et al., New J. Phys. 16, 093029 (2014).
- **Verification status**: ✅ VERIFIED
  - Widely cited and used in the molecular electronics community.
