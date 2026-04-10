# magnum.af

## Official Resources
- Source Repository: https://github.com/magnum-af/magnum.af
- Documentation: https://magnum-af.github.io/
- License: Open source

## Overview
**magnum.af** is a finite-difference/finite-element micromagnetic simulation package that combines CPU and GPU solvers. It supports standard micromagnetic energy terms plus spin-transfer torque, spin-orbit torque, and true periodic boundary conditions for stray field calculation.

**Scientific domain**: Micromagnetic simulation, spintronics, domain dynamics  
**Target user community**: Researchers simulating magnetization dynamics with advanced boundary conditions and spin-torque effects

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Finite-difference and finite-element methods
- FFT-based demagnetization
- Spin-transfer torque (Zhang-Li, Slonczewski)
- Spin-orbit torque
- True periodic boundary conditions
- Thermal fluctuations

## Capabilities (CRITICAL)
- Micromagnetic simulation (LLG dynamics)
- Energy minimization
- Spin-transfer torque simulation
- Spin-orbit torque simulation
- True periodic boundary conditions
- GPU acceleration (CUDA)
- Finite-difference and finite-element solvers
- Domain wall dynamics
- Skyrmion simulation

**Sources**: GitHub repository, published in J. Magn. Magn. Mater.

## Key Strengths

### Advanced Boundary Conditions:
- True periodic BC for stray field
- Eliminates edge effects
- Accurate for infinite thin films
- Novel approach for periodic systems

### Spin-Torque Effects:
- Spin-transfer torque (STT)
- Spin-orbit torque (SOT)
- Zhang-Li and Slonczewski models
- Current-driven dynamics

### GPU Acceleration:
- CUDA implementation
- Fast demagnetization calculation
- Large system sizes feasible
- Efficient time integration

## Inputs & Outputs
- **Input formats**:
  - JSON configuration files
  - Mesh files (for FEM)
  - Material parameters
  
- **Output data types**:
  - Magnetization fields
  - Energy vs time
  - Hysteresis loops
  - VTK output for visualization

## Interfaces & Ecosystem
- **Python**: Scripting interface
- **C++**: Core computation
- **CUDA**: GPU acceleration
- **ParaView**: VTK visualization

## Performance Characteristics
- **Speed**: Fast with GPU
- **Accuracy**: High (validated)
- **System size**: Millions of cells (GPU)
- **Parallelization**: GPU (CUDA)

## Computational Cost
- **Small systems**: Minutes
- **Large systems**: Hours (GPU)
- **Typical**: Moderate with GPU

## Limitations & Known Constraints
- **CUDA required**: GPU acceleration needs NVIDIA
- **Limited documentation**: Could be more extensive
- **Community**: Smaller than OOMMF/Mumax3
- **No DFT integration**: Standalone micromagnetics

## Comparison with Other Codes
- **vs OOMMF**: magnum.af has GPU and true PBC, OOMMF is NIST standard
- **vs Mumax3**: magnum.af has FEM + true PBC, Mumax3 is pure FD/GPU
- **vs Spirit**: magnum.af is micromagnetic, Spirit is atomistic
- **Unique strength**: True periodic boundary conditions for stray field, FEM+FD solvers, spin-torque support

## Application Areas

### Spintronics:
- STT-MRAM switching
- SOT switching
- Domain wall motion by current
- Skyrmion Hall effect

### Thin Films:
- Periodic domain structures
- Stripe domains
- Skyrmion lattices
- Domain wall pinning

### Standard Problems:
- µMAG benchmarks
- Method comparison
- PBC-specific problems

## Best Practices

### Mesh Selection:
- Use cell size ≤ exchange length
- Test PBC effects
- Validate against analytical solutions
- Compare with non-PBC results

### GPU Usage:
- Ensure CUDA compatibility
- Monitor GPU memory usage
- Use appropriate block sizes
- Compare CPU/GPU results

## Community and Support
- Open source on GitHub
- Developed at TU Wien / Danube University Krems
- Published methodology
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/magnum-af/magnum.af
2. P. Heistracher et al., J. Magn. Magn. Mater. 548, 168875 (2022)
3. F. Bruckner et al., Scientific Reports 11, 9202 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE
- Published methodology: J. Magn. Magn. Mater.
- Active development: Ongoing
- Specialized strength: True periodic boundary conditions, FEM+FD solvers, spin-torque support
