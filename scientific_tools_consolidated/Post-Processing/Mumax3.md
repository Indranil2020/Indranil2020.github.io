# Mumax3

## Official Resources
- Homepage: https://mumax.github.io/
- Documentation: https://mumax.github.io/api.html
- Source Repository: https://github.com/mumax/3
- License: GNU General Public License v3.0

## Overview
Mumax3 is a GPU-accelerated micromagnetic simulation program. It solves the time-dependent Landau-Lifshitz-Gilbert (LLG) equation using finite difference discretization. It is designed to be highly efficient, running exclusively on NVIDIA GPUs, and provides a rich scripting interface (Go-like syntax) for defining complex geometries, time-dependent fields, and various magnetic interactions.

**Scientific domain**: Micromagnetics, spintronics, magnetic dynamics  
**Target user community**: Magnetism researchers, spintronics engineers, device physicists

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Finite Difference Method (FDM)
- Exchange interaction (Heisenberg)
- Dipole-Dipole interaction (Demagnetization)
- Magnetic Anisotropy (Uniaxial, Cubic)
- Dzyaloshinskii-Moriya Interaction (DMI) - Bulk and Interfacial
- Spin-Transfer Torque (Zhang-Li, Slonczewski)
- Thermal fluctuations (Langevin dynamics)

## Capabilities (CRITICAL)
- GPU-accelerated micromagnetic simulations (very fast)
- Simulation of domain walls, skyrmions, vortices
- Spin-transfer torque and spin-orbit torque
- Time-dependent excitation (RF fields, pulses)
- Complex geometries via constructive solid geometry (CSG) or image masks
- Voronoi tessellation for grains
- Periodic boundary conditions

**Sources**: Mumax3 website, AIP Advances 4, 107133 (2014)

## Inputs & Outputs
- **Input formats**: `.mx3` script files, image files (png) for masks
- **Output data types**: `.ovf` (vector fields), `.txt` (tables), `.png` (snapshots)

## Interfaces & Ecosystem
- **Go**: Scripting language based on Go
- **Python**: `mumax3-python` wrapper available
- **Visualization**: Output compatible with ParaView, Muview, OOMMF tools
- **Web Interface**: Mumax3-server

## Workflow and Usage
1. Write `.mx3` script: Define geometry, material parameters, and simulation protocol.
2. Run simulation: `mumax3 script.mx3`
3. Convert output: `mumax3-convert` to standard formats (OVF/VTK).
4. Analyze tables and visualize vector fields.

## Performance Characteristics
- Optimized for NVIDIA CUDA
- Speedup of 10x-100x compared to CPU codes (like OOMMF)
- Performance scales with number of CUDA cores

## Application Areas
- Magnetic storage devices (MRAM, HDD)
- Skyrmionics and topological magnetism
- Magnonics (spin waves)
- Nanomagnetism and logic devices

## Community and Support
- Open-source (GPL v3)
- Developed by DyNaMat group (Ghent University)
- Active Google Group for support
- Large user community in micromagnetics

## Verification & Sources
**Primary sources**:
1. Homepage: https://mumax.github.io/
2. GitHub: https://github.com/mumax/3
3. Publication: A. Vansteenkiste et al., AIP Advances 4, 107133 (2014)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE/MAINTAINED
- Applications: Micromagnetics, GPU acceleration, LLG
