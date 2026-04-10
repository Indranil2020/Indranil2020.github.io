# MicroMagnetic.jl

## Official Resources
- Source Repository: https://github.com/MagneticSimulation/MicroMagnetic.jl
- Documentation: https://magneticimulation.github.io/MicroMagnetic.jl/
- License: Open source

## Overview
**MicroMagnetic.jl** is a Julia package for classical spin dynamics and micromagnetic simulations with multi-platform GPU support (NVIDIA, AMD, Intel, Apple). It supports atomistic and continuum spin simulations, Monte Carlo, NEB energy barriers, and spin-transfer torque effects.

**Scientific domain**: Micromagnetic and atomistic spin simulation, GPU computing  
**Target user community**: Researchers needing fast GPU-accelerated magnetic simulations across different GPU platforms

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Finite-difference method
- Atomistic and continuum spin models
- Monte Carlo simulation
- Nudged elastic band (NEB)
- Spin-transfer torque (Zhang-Li, Slonczewski)
- Exchange, anisotropy, Zeeman, demagnetization
- Dzyaloshinskii-Moriya interaction
- Constructive solid geometry (CSG)

## Capabilities (CRITICAL)
- Micromagnetic simulation (continuum)
- Atomistic spin simulation
- Monte Carlo simulation
- NEB energy barrier calculation
- Multi-GPU support (CUDA, AMDGPU, oneAPI, Metal)
- Spin-transfer torque
- Periodic boundary conditions
- Constructive solid geometry
- Thermal fluctuations
- Julia scripting interface

**Sources**: GitHub repository

## Key Strengths

### Multi-Platform GPU:
- NVIDIA (CUDA)
- AMD (AMDGPU.jl)
- Intel (oneAPI.jl)
- Apple (Metal.jl)
- CPU fallback

### Julia Performance:
- Near-C performance
- JIT compilation
- Automatic differentiation
- Composable with Julia ecosystem

### Comprehensive Physics:
- Atomistic + continuum
- Monte Carlo + dynamics
- NEB energy barriers
- Spin-torque effects
- CSG geometry definition

## Inputs & Outputs
- **Input formats**:
  - Julia scripts
  - Material parameters
  - Mesh specifications (CSG)
  
- **Output data types**:
  - Magnetization fields
  - Energy vs time
  - NEB paths and barriers
  - VTK output

## Interfaces & Ecosystem
- **Julia**: Primary language
- **CUDA.jl/AMDGPU.jl**: GPU backends
- **Makie.jl**: Visualization
- **JLD2.jl**: Data storage

## Performance Characteristics
- **Speed**: Fast (GPU-accelerated)
- **Accuracy**: Good (validated)
- **System size**: Millions of cells (GPU)
- **Multi-GPU**: Supported

## Computational Cost
- **Small systems**: Seconds (GPU)
- **Large systems**: Minutes to hours
- **Typical**: Very efficient with GPU

## Limitations & Known Constraints
- **Julia language**: Less common than Python/C++
- **Newer code**: Less established than OOMMF
- **Documentation**: Growing but limited
- **Community**: Small but active

## Comparison with Other Codes
- **vs Mumax3**: MicroMagnetic.jl supports more GPU platforms, Mumax3 is NVIDIA-only
- **vs OOMMF**: MicroMagnetic.jl is GPU, OOMMF is CPU
- **vs Spirit**: MicroMagnetic.jl is Julia/GPU, Spirit is C++
- **Unique strength**: Multi-platform GPU support (NVIDIA/AMD/Intel/Apple), Julia performance, atomistic+continuum

## Application Areas

### GPU-Accelerated Micromagnetics:
- Large-scale domain dynamics
- Fast hysteresis loops
- Parameter sweeps
- Real-time simulation

### Skyrmion Dynamics:
- Skyrmion creation and motion
- Skyrmion lattice dynamics
- Current-driven skyrmions
- Temperature effects

### Energy Barriers:
- Switching field determination
- Thermal stability analysis
- NEB transition paths
- Coercivity estimation

### Novel Architectures:
- Apple Silicon GPU
- AMD GPU clusters
- Intel GPU nodes
- Cross-platform portability

## Best Practices

### GPU Selection:
- Use available GPU backend
- Test CPU vs GPU for small systems
- Monitor GPU memory usage
- Use appropriate precision

### Julia Setup:
- Install Julia with GPU support
- Pre-compile for first-run speed
- Use Makie for visualization
- Leverage Julia ecosystem

## Community and Support
- Open source on GitHub
- Developed by MicroMagnetic.jl team
- Active development
- Documentation growing

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/MagneticSimulation/MicroMagnetic.jl
2. Documentation: https://magneticimulation.github.io/MicroMagnetic.jl/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE
- Active development: Ongoing
- Specialized strength: Multi-platform GPU micromagnetic simulation (NVIDIA/AMD/Intel/Apple), Julia performance
