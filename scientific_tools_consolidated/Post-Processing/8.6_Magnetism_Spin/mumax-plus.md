# mumax+

## Official Resources
- Source Repository: https://github.com/mumax/plus
- Documentation: Included in repository
- License: Open source

## Overview
**mumax+** (mumax plus) is a versatile and extensible GPU-accelerated micromagnetic simulator written in C++ and CUDA with a Python interface. It is the successor to mumax3, offering more extensibility, Python scripting, and additional physics capabilities beyond standard micromagnetics.

**Scientific domain**: GPU-accelerated micromagnetic simulation  
**Target user community**: Researchers needing fast, extensible GPU micromagnetic simulations with Python control

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Finite-difference method on GPU
- FFT-based demagnetization
- Exchange, anisotropy, Zeeman energies
- Dzyaloshinskii-Moriya interaction
- Spin-transfer torque
- Thermal fluctuations
- Reduced units system

## Capabilities (CRITICAL)
- GPU-accelerated micromagnetic simulation
- Python scripting interface
- Extensible C++ architecture
- All standard micromagnetic energy terms
- DMI support
- Spin-transfer torque
- Thermal fluctuations
- Custom energy terms (extensible)
- Batch simulation mode
- Real-time visualization

**Sources**: GitHub repository

## Key Strengths

### GPU Performance:
- CUDA-accelerated
- Fast FFT demagnetization
- Large system sizes
- Efficient memory usage

### Extensibility:
- C++ plugin architecture
- Custom energy terms
- Custom evolvers
- Python scripting
- More flexible than mumax3

### Python Interface:
- Full Python control
- Jupyter notebook compatible
- Easy post-processing
- Automation and scripting

### mumax3 Compatibility:
- Familiar input format
- Migration path from mumax3
- Same reduced units
- Enhanced capabilities

## Inputs & Outputs
- **Input formats**:
  - mumax+ input files
  - Python scripts
  - Material parameters
  
- **Output data types**:
  - Magnetization fields (OVF)
  - Energy vs time
  - Average magnetization
  - Table data
  - PNG snapshots

## Interfaces & Ecosystem
- **Python**: Scripting interface
- **C++/CUDA**: Core computation
- **OVF format**: Standard micromagnetic output
- **ParaView/MumaxView**: Visualization

## Performance Characteristics
- **Speed**: Very fast (GPU)
- **Accuracy**: High (validated)
- **System size**: Millions of cells
- **Memory**: GPU memory limited

## Computational Cost
- **Small systems**: Seconds
- **Large systems**: Minutes to hours
- **Typical**: Very efficient with GPU

## Limitations & Known Constraints
- **NVIDIA GPU only**: Requires CUDA
- **Finite differences**: No FEM
- **Regular grids**: No adaptive meshing
- **In development**: Not as mature as mumax3
- **Documentation**: Limited

## Comparison with Other Codes
- **vs mumax3**: mumax+ is extensible, has Python interface, still GPU
- **vs OOMMF**: mumax+ is GPU, OOMMF is CPU
- **vs MicroMagnetic.jl**: mumax+ is NVIDIA-only, MicroMagnetic.jl is multi-platform
- **Unique strength**: Extensible GPU micromagnetic simulator with Python interface, mumax3 successor

## Application Areas

### Fast Micromagnetics:
- Large-scale domain dynamics
- Parameter sweeps
- Hysteresis loops
- Real-time simulation

### Spintronics:
- STT-MRAM
- Domain wall motion
- Skyrmion dynamics
- Spin-orbit torque

### Custom Physics:
- Extended Hamiltonians
- Coupled systems
- Novel energy terms
- Research code development

## Best Practices

### GPU Setup:
- Use NVIDIA GPU with sufficient memory
- Monitor GPU utilization
- Choose appropriate cell size
- Use reduced units

### Python Scripting:
- Use Python for automation
- Post-process in Python
- Generate input files programmatically
- Use Jupyter for interactive analysis

## Community and Support
- Open source on GitHub
- Developed alongside mumax3
- Active development
- GitHub Discussions for support

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/mumax/plus
2. Related: A. Vansteenkiste et al., AIP Adv. 4, 107133 (2014) (mumax3 paper)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Ongoing
- Specialized strength: Extensible GPU micromagnetic simulator with Python interface, mumax3 successor
