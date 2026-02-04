# KITE

## Official Resources
- Homepage: https://quantum-kite.com/
- Documentation: https://quantum-kite.com/documentation/
- Source Repository: https://github.com/quantum-kite/kite
- License: GNU General Public License v3.0

## Overview
KITE is an open-source quantum transport software optimized for large-scale tight-binding models (up to billions of orbitals). It uses Chebyshev polynomial expansion of the Green's function and the kernel polynomial method (KPM) to calculate electronic structure and transport properties of disordered systems, topological materials, and 2D heterostructures with high efficiency.

**Scientific domain**: Quantum transport, tight-binding, kernel polynomial method, large-scale simulations  
**Target user community**: Condensed matter physicists, 2D materials researchers

## Theoretical Methods
- Tight-Binding Model
- Kernel Polynomial Method (KPM)
- Chebyshev Expansion of Green's functions
- Kubo-Bastin formula (linear response)
- Stochastic trace evaluation (random vectors)
- Real-space recursion

## Capabilities (CRITICAL)
- **Scalability**: Simulations of systems with $10^9$ orbitals
- **Transport**: Longitudinal and transverse (Hall) conductivity
- **Optical**: Optical conductivity (AC transport)
- **Electronic**: Density of States (DOS)
- **Disorder**: Efficient handling of Anderson disorder, vacancies, and structural defects
- **Topology**: Chern numbers and topological invariants
- **Non-linear**: Non-linear optical response (second harmonic generation)

**Sources**: KITE website, Sci. Adv. 6, eaay0343 (2020)

## Inputs & Outputs
- **Input formats**: Python script defining lattice and Hamiltonian (`kite.Lattice`, `kite.System`)
- **Output data types**: HDF5 files containing spectral data, conductivities, and moments

## Interfaces & Ecosystem
- **Python**: Interface for setting up and running simulations
- **C++**: High-performance backend
- **PyBinding**: Compatible for model construction
- **Graphene/TMDs**: Pre-defined lattices available

## Workflow and Usage
1. Define lattice and hopping parameters in Python script.
2. Configure calculation: `configuration = kite.Configuration(...)`.
3. Run calculation: `kite.Calculation(configuration)`.
4. Post-process: `kite.Dos(...)` or `kite.Conductivity(...)` using Chebyshev moments.

## Performance Characteristics
- Linear scaling O(N) with system size
- Highly parallelized (OpenMP/MPI)
- Memory efficient (sparse matrices and stochastic trace)

## Application Areas
- Disordered graphene and 2D materials
- Quantum Hall effect in large samples
- Twisted bilayer graphene (Moiré patterns)
- Topological insulators
- Quasicrystals

## Community and Support
- Developed by Yuan Group (University of York / Wuhan University)
- Open-source (GPL v3)
- Active development and workshops

## Verification & Sources
**Primary sources**:
1. Homepage: https://quantum-kite.com/
2. GitHub: https://github.com/quantum-kite/kite
3. Publication: S. M. Joao et al., Royal Society Open Science 7, eaax5847 (2020)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Large-scale tight-binding, KPM, transport
