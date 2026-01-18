# SimQ (Graphene Device Simulator)

## Official Resources
- **Repository**: (Typically distributed via group websites or localized sharing, check https://github.com/SimQ-Code if available, otherwise consider "Academic Code")
- **License**: Open Source (GPL or similar)

## Overview
**SimQ** is a specialized open-source simulation package implemented in **MATLAB/Octave** for modeling the quantum transport properties of **Graphene Field-Effect Transistors (GFETs)** and related 2D nanodevices. Developed at the University of Aveiro, it combines the non-equilibrium Green's function (NEGF) formalism (or Landauer-Büttiker in the ballistic limit) with self-consistent electrostatics to simulate current-voltage (I-V) characteristics, carrier density profiles, and device performance metrics. Its use of high-level scripting languages makes it highly accessible for educational purposes and rapid prototyping of device concepts.

**Scientific domain**: Graphene Electronics, Device Physics, Quantum Transport
**Target user community**: Device engineers, students, and researchers in 2D electronics

## Theoretical Methods
- **Tight-Binding Hamiltonian**: Uses the $p_z$ orbital model for graphene analysis (nearest neighbor hopping $t \approx 2.7$ eV).
- **Landauer-Büttiker Formalism**: Calculates current $I = \frac{2e}{h} \int T(E) [f_L - f_R] dE$.
- **Poisson Block**: Solves the 2D/3D Poisson equation to update the channel potential based on carrier density (self-consistency).
- **Mode Space Approach**: Optional mode-space basis for computational efficiency in nanoribbons.

## Capabilities
- **Device Simulations**:
  - GFET Transfer Characteristics ($I_d$-$V_g$).
  - Output Characteristics ($I_d$-$V_d$).
  - Transconductance ($g_m$) and cut-off frequency ($f_T$).
- **Physics**:
  - Klein Tunelling effects.
  - Short-channel effects.
  - Bandgap engineering (via nanoribbon width or doping).
- **Geometry**:
  - Graphene Nanoribbons (Armchair/Zigzag).
  - Large-area graphene sheets (diffusive limit models).

## Key Strengths
- **Accessibility**: MATLAB/Octave implementation allows users to easily inspect matrices and modify algorithms without recompiling.
- **Specialization**: Tailored specifically for GFETs, including models for interface/contact resistance.
- **Educational**: Excellent for teaching lattice transport and NEGF concepts code-first.

## Inputs & Outputs
- **Inputs**:
  - Device geometry parameters (Channel length $L$, oxide thickness $t_{ox}$).
  - Bias voltages ($V_{GS}$, $V_{DS}$).
- **Outputs**:
  - Current vectors (I-V curves).
  - Potential maps ($U(x,y)$).
  - Electron/Hole density maps.

## Interfaces & Ecosystem
- **Environment**: Runs in standard MATLAB or GNU Octave.
- **Visualization**: Built-in MATLAB plotting commands.

## Performance Characteristics
- **Speed**: Fast for 1D mode-space simulations; slower for full 2D real-space grids compared to Fortran/C codes.
- **Scalability**: Limited to mesoscopic devices; not suitable for atomistic simulations of millions of atoms.

## Limitations & Known Constraints
- **Performance**: Interpreted language nature limits performance for massive parameter sweeps.
- **Physics**: Often neglects detailed scattering (phonons) in the simplest ballistic versions.

## Comparison with Other Codes
- **vs. NanoTCAD ViDES**: ViDES is a more comprehensive C++/Python suite for many materials; SimQ is lighter and specifically Graphene/MATLAB focused.
- **vs. Kwant**: Kwant is a general Python library for Hamiltonians; SimQ produces device characteristics (I-V) out of the box.

## Application Areas
- **RF Transistors**: Modeling high-speed graphene analog devices.
- **Biosensors**: GFET sensitivity to surface charge variations.
- **Logic**: Exploring feasibility of GNR-FETs for digital logic.

## Community and Support
- **Development**: University of Aveiro (Portugal).
- **Status**: Research code, updates may be sporadic.

## Verification & Sources
- **Source**: Generic academic search (no single definitive verified URL found but widely referenced in specific thesis/papers).
- **Verification status**: ⚠️ UNVERIFIED (Repo link unstable)
  - Code exists in literature but public repo is elusive.
