# HPhi

## Official Resources
- Homepage: https://github.com/QLMS/HPhi
- Documentation: https://github.com/QLMS/HPhi/wiki
- Source Repository: https://github.com/QLMS/HPhi
- License: GNU General Public License v3.0

## Overview
HPhi is a software package for solving quantum lattice models using the exact diagonalization method. It supports a wide range of quantum lattice models, including Hubbard, Heisenberg, and Kondo lattice models. HPhi can calculate ground state and excited state properties, as well as thermal averages using the thermal pure quantum (TPQ) state method.

**Scientific domain**: Quantum lattice models, exact diagonalization, quantum many-body physics  
**Target user community**: Condensed matter physicists, quantum magnetism researchers

## Theoretical Methods
- Exact diagonalization (Lanczos, LOBPCG)
- Thermal Pure Quantum (TPQ) state method
- Quantum lattice models (Hubbard, Heisenberg, Kondo)
- Real-space parallelization

## Capabilities (CRITICAL)
- Ground state energy and wavefunction
- Excited states
- Finite temperature properties (specific heat, susceptibility)
- Green's functions
- Time evolution
- Spin dynamics
- Parallel computing support (MPI/OpenMP)

**Sources**: HPhi documentation, Comp. Phys. Comm. 217, 180 (2017)

## Inputs & Outputs
- **Input formats**: Standard input file format (text), model parameters
- **Output data types**: Energy spectra, correlation functions, thermodynamic quantities

## Interfaces & Ecosystem
- **Python**: HPhi inputs can be generated via Python scripts
- **C interface**: Library interface available

## Performance Characteristics
- Highly parallelized for distributed memory systems
- Efficient for intermediate system sizes

## Application Areas
- Quantum magnetism
- Strongly correlated electron systems
- Frustrated spin systems
- Quantum phase transitions

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Active development by ISSP, University of Tokyo

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/QLMS/HPhi
2. Documentation: https://github.com/QLMS/HPhi/wiki
3. Publication: Comp. Phys. Comm. 217, 180 (2017)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Applications: Exact diagonalization, TPQ method, quantum lattice models
