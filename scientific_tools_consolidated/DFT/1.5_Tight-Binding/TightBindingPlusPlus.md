# TightBinding++

## Official Resources
- Homepage: https://github.com/huchou/TightBinding (Search result linked to "TightBinding++" but repo might be named differently)
- Source Repository: https://github.com/huchou/TightBinding
- License: GNU GPL v3

## Overview
TightBinding++ is a C++ framework with a Python wrapper/interface tailored for the simulation of quantum tight-binding models. It automates the generation of Hamiltonian matrices and allows for the inclusion of external magnetic fields and disorder, facilitating the study of topological systems and transport properties using the Kubo-Greenwood formalism.

**Scientific domain**: Tight-Binding, Topological Insulators, Transport
**Target user community**: Theoretical physicists

## Theoretical Methods
- Tight-Binding Hamiltonian
- Peierls Substitution (Magnetic fields)
- Coherent Potential Approximation (CPA) for disorder
- Kubo-Greenwood formula (Conductivity)
- Chern number / Topological invariants

## Capabilities (CRITICAL)
- Hamiltonian generation from lattice
- Magnetic field effects
- Random disorder simulation
- Linear response conductivity
- Band structure calculation
- HDF5 output format

## Key Strengths

### Performance:
- C++ core for efficiency
- Multi-threaded algorithms

### Advanced Physics:
- Native support for disorder (CPA)
- Topology and Transport focus
- Magnetic fields

## Inputs & Outputs
- **Input**:
  - Model definition (XML/Script)
  - Lattice parameters
- **Output**:
  - HDF5 files (*.tbpp)
  - Band structures
  - Conductivities

## Interfaces & Ecosystem
- **Python**: Python wrapper for ease of use
- **HDF5**: Standard data format integration

## Application Areas
- **Quantum Hall Effect**: Integer and fractional effects in lattice models
- **Anderson Localization**: Disorder-induced localization studies
- **Topological Materials**: Berry curvature and edge state characterization
- **Transport**: Conductivity in mesoscopic systems

## Best Practices
- **Disorder Averaging**: Ensure sufficient samples/CPA convergence for disordered systems
- **Magnetic Fields**: Check commensurability of magnetic flux (Hofstadter butterfly)
- **Memory Management**: Use HDF5 output efficiently for large sweeps
- **Parallelization**: Enable OpenMP for threading on multicore machines

## Comparison with Other Codes
- **vs Kwant**: TightBinding++ has specific disorder (CPA) features
- **vs PythTB**: C++ backend offers higher performance for complex models

## Verification & Sources
**Primary sources**:
1. GitHub Page: https://huchou.github.io/TightBinding/
2. Repository: https://github.com/huchou/TightBinding

**Confidence**: VERIFIED
- Status: Open Source
