# OMEN

## Official Resources
- **Homepage**: https://engineering.purdue.edu/gekcogrp/software-projects/omen/
- **Repository**: https://github.com/spcl/dace-omen (DaCe variant)
- **License**: Academic / Open Source variants exist

## Overview
**OMEN** is a massively parallel, multidimensional quantum transport simulator optimized for **high-performance computing (HPC)** environments. It focuses on the simulation of post-CMOS nanodevices (nanowires, TFETs) using the **semi-empirical tight-binding** method ($sp^3d^5s^*$) coupled with the **Non-Equilibrium Green's Function (NEGF)** or **Wave Function (WF)** formalisms. Known for its extreme scalability, OMEN has been used to simulate realistic devices with tens of thousands of atoms including electron-phonon scattering.

**Scientific domain**: Quantum Transport, Supercomputing, Device Physics
**Target user community**: HPC specialists and device physicists studying dissipative transport

## Theoretical Methods
- **NEGF**: Full calculating of retarded/lesser Green's functions with self-energies for scattering.
- **Wave Function**: Faster "transmitting boundary" method for ballistic transport.
- **Scattering**: Self-consistent Born approximation for acoustic/optical phonons and surface roughness.
- **Basis**: Localized tight-binding orbitals.

## Capabilities
- **Device Simulation**:
  - Nanowire Field-Effect Transistors (GAA-NWFETs).
  - Tunnel FETs (TFETs) with band-to-band tunneling.
  - 2D material logic devices.
- **Physics**:
  - Dissipative transport (Joule heating).
  - Band structure effects in confined geometries.

## Key Strengths
- **HPC Performance**: A Gordon Bell Prize finalist code, capable of scaling to >200,000 cores on machines like Titan/Summit.
- **Scattering**: One of the few codes that can handle full 3D atomistic NEGF with scattering (albeit at high cost).
- **Data Centric**: Recent versions (DaCe-OMEN) explore new programming paradigms for efficiency.

## Inputs & Outputs
- **Inputs**: Atomistic structure maps, material parameter files, bias commands.
- **Outputs**: Current, charge density, energy-resolved current spectra.

## Interfaces & Ecosystem
- **CP2K**: Integration allows for DFT-based transport (Hamiltonian blocks from CP2K fed into OMEN).
- **libNEGF**: Can use external solver libraries.

## Performance Characteristics
- **Cost**: Extremely high for NEGF+Scattering.
- **Optimization**: Uses SSE/AVX intrinsics, GPU acceleration, and advanced MPI/OpenMP tiling.

## Comparison with Other Codes
- **vs. NEMO5**: OMEN is the specialized transport engine; NEMO5 wraps it in a larger multiphysics framework. OMEN is often used for pure transport performance benchmarks.
- **vs. Kwant**: OMEN is a fully-fledged device simulator with physics-specific scattering models; Kwant is a Hamiltonian solver toolkit.

## Community and Support
- **Development**: ETH Zurich (Mathieu Luisier) and Purdue.
- **Source**: GitHub (DaCe-OMEN) / NanoHub.

## Verification & Sources
- **Primary Publication**: M. Luisier et al., SC10 (2010).
- **Verification status**: ✅ VERIFIED
  - Validated against experimental data for nanowire transistors.
