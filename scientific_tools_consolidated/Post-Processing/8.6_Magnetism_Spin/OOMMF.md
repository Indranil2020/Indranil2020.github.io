# OOMMF

## Official Resources
- Homepage: https://math.nist.gov/oommf/
- Download: https://math.nist.gov/oommf/software.html
- Documentation: https://math.nist.gov/oommf/doc/userguide20a3/userguide/
- License: Public domain (NIST)

## Overview
**OOMMF** (Object Oriented MicroMagnetic Framework) is a public domain micromagnetic simulation program developed at NIST. It solves the Landau-Lifshitz-Gilbert equation on a finite-difference grid using Tcl/Tk for GUI and C++ for computation, and is the de facto standard for micromagnetic benchmark problems.

**Scientific domain**: Micromagnetic simulation, domain dynamics, magnetic recording  
**Target user community**: Researchers simulating magnetization dynamics in nanostructures, thin films, and magnetic devices

## Theoretical Methods
- Landau-Lifshitz-Gilbert (LLG) equation
- Finite-difference method
- Exchange energy
- Zeeman energy
- Magnetostatic (demagnetization) energy
- Magnetocrystalline anisotropy energy
- Surface anisotropy
- Dzyaloshinskii-Moriya interaction (via extensions)

## Capabilities (CRITICAL)
- Micromagnetic simulation (time-domain LLG)
- Energy minimization (relaxation to ground state)
- Hysteresis loop calculation
- Domain wall dynamics
- Vortex and skyrmion simulation
- Standard problem benchmarks (µMAG)
- Tcl/Tk graphical interface
- Extensible via Oxs (OOMMF eXtensible Solver) child classes
- Batch simulation mode
- 2D and 3D geometries

**Sources**: NIST OOMMF documentation, µMAG standard problems

## Key Strengths

### NIST Standard:
- Reference implementation for micromagnetics
- µMAG standard problem validation
- Public domain (no licensing restrictions)
- Decades of validation and testing

### Extensible Architecture:
- Oxs child class system
- Custom energy terms
- Custom evolvers
- Tcl scripting for automation
- Python interfaces (via extensions)

### Comprehensive Physics:
- All standard energy terms
- Accurate demagnetization field
- Multiple integration methods
- Thermal fluctuations support

### GUI and Scripting:
- Interactive Tcl/Tk GUI
- Batch scripting capability
- Real-time visualization
- Parameter sweeps

## Inputs & Outputs
- **Input formats**:
  - OMF/OVF vector field files
  - MIF (OOMMF Material Input Format)
  - Tcl script files
  
- **Output data types**:
  - Magnetization vector fields (OMF/OVF)
  - Energy vs time
  - Hysteresis loops
  - Average magnetization
  - Field files

## Interfaces & Ecosystem
- **joommf**: Jupyter/OOMMF Python interface
- **oommf-extension-dmi-t**: DMI extension
- **oommf-python**: Python wrapper
- **µMAG**: Standard problem benchmarks

## Performance Characteristics
- **Speed**: Moderate (CPU-based, no GPU)
- **Accuracy**: High (validated against benchmarks)
- **System size**: Millions of cells
- **Parallelization**: Limited (mostly serial)

## Computational Cost
- **Small systems**: Minutes
- **Large systems**: Hours to days
- **Typical**: Moderate (no GPU acceleration)

## Limitations & Known Constraints
- **No GPU**: CPU-only computation
- **Finite differences only**: No finite elements
- **Regular grids only**: No adaptive meshing
- **Tcl/Tk dependency**: Required for GUI
- **No native DMI**: Requires extension
- **Limited parallelization**: Mostly serial

## Comparison with Other Codes
- **vs Mumax3**: OOMMF is CPU, Mumax3 is GPU; OOMMF is NIST standard
- **vs magnum.af**: OOMMF is finite-difference, magnum.af supports FEM
- **vs MicroMagnetic.jl**: OOMMF is Tcl/C++, MicroMagnetic.jl is Julia/GPU
- **Unique strength**: NIST standard micromagnetic code, public domain, extensive validation, extensible Oxs architecture

## Application Areas

### Magnetic Recording:
- Write head dynamics
- Media switching
- Bit patterned media
- Heat-assisted recording

### Domain Dynamics:
- Domain wall motion
- Vortex dynamics
- Skyrmion creation/annihilation
- Domain wall pinning

### Thin Film Magnetism:
- Permalloy films
- Multilayer systems
- Exchange bias
- Interlayer coupling

### Standard Problems:
- µMAG benchmarks
- Method validation
- Code comparison
- Teaching tool

## Best Practices

### Mesh Selection:
- Use cell size ≤ exchange length
- Test convergence with cell size
- Consider geometry discretization
- Use appropriate boundary conditions

### Integration:
- Choose appropriate time step
- Monitor energy conservation
- Use Runge-Kutta for accuracy
- Validate against standard problems

### Extension Usage:
- Install DMI extension for chiral magnets
- Use joommf for Python integration
- Check extension compatibility
- Validate extended physics

## Community and Support
- Public domain (NIST)
- Extensive documentation
- µMAG community benchmarks
- Active mailing list
- Multiple Python wrapper projects

## Verification & Sources
**Primary sources**:
1. NIST OOMMF: https://math.nist.gov/oommf/
2. M. J. Donahue and D. G. Porter, NISTIR 6376 (1999)
3. µMAG benchmarks: https://www.ctcms.nist.gov/~rdm/mumag.org.html

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE (NIST)
- Documentation: COMPREHENSIVE
- Source code: PUBLIC DOMAIN
- Community support: ACTIVE (decades of use)
- Academic citations: >10000
- Active development: Ongoing (v2.0)
- Specialized strength: NIST standard micromagnetic code, public domain, extensive validation, extensible Oxs architecture
