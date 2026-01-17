# RT-tddft (sheyua)

## Official Resources
- Homepage: https://github.com/sheyua/RT-tddft
- Source Repository: https://github.com/sheyua/RT-tddft
- License: Open Source (based on Quantum ESPRESSO)

## Overview
RT-tddft is a Real-Time Plane-Wave Time-Dependent Density Functional Theory code designed for simulating ultrafast electron dynamics, particularly focused on discharging nanostructures and non-equilibrium phenomena. It is built as an extension to Quantum ESPRESSO.

**Scientific domain**: Ultrafast electron dynamics, nanostructure discharge, non-equilibrium processes  
**Target user community**: Researchers studying time-dependent phenomena in nanoscale systems using plane-wave DFT

## Theoretical Methods
- Real-Time Time-Dependent DFT (RT-TDDFT)
- Plane-wave basis set implementation
- Explicit time propagation of Kohn-Sham equations
- Non-equilibrium electron dynamics
- Strong-field response simulations

## Capabilities
- **Ultrafast dynamics simulation** of nanostructures
- **Discharge dynamics** modeling
- **Real-time electron propagation** in extended systems
- **Non-equilibrium phenomena** simulation
- **Plane-wave accuracy** for periodic and slab systems

## Key Strengths

### Specialized Focus:
- Designed for nanostructure dynamics
- Discharge process simulation
- Non-equilibrium electron behavior

### QE Foundation:
- Built on Quantum ESPRESSO infrastructure
- Plane-wave basis accuracy
- Established pseudopotential library

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO-style input files
  - Structure files
  - Pseudopotentials (QE format)
  
- **Output data types**:
  - Time-dependent electronic properties
  - Dynamics trajectories
  - Charge evolution data

## Interfaces & Ecosystem
- **Quantum ESPRESSO** base code integration
- Standard QE pseudopotential compatibility
- QE post-processing tools applicable

## Performance Characteristics
- **Basis**: Plane-wave (systematic convergence)
- **Parallelization**: MPI (inherited from QE)
- **Accuracy**: Controlled by energy cutoff
- **Time step**: Requires small steps (attosecond scale) for stability

## Advanced Features
- **Discharge Dynamics**: Specialized algorithms for simulating charge loss/gain in nanostructures.
- **Strong Field Interaction**: Coupling with external laser fields.
- **Time-Dependent Charge Analysis**: Tools to monitor charge fluctuations in real-time.

## Computational Cost
- **High**: RT-TDDFT is computationally demanding, typically 10-100x more than ground state DFT.
- **Scaling**: Scales similarly to standard Plane-Wave DFT ($O(N^3)$) with system size, multiplied by thousands of time steps.
- **Memory**: Stores time-dependent wavefunctions, requiring significant RAM for large systems.

## Best Practices
- **Ground State**: Ensure a well-converged ground state before propagation.
- **Time Step**: Use conservative time steps ($< 0.1$ atomic units) to avoid instability.
- **Vacuum Padding**: For nanostructures, ensure sufficient vacuum to avoid periodic image interactions, especially for discharge simulations.

## Community and Support
- **Source**: GitHub repository (sheyua/RT-tddft).
- **Issues**: Use GitHub Issues for bug reports.
- **Relation to QE**: Leverages the broader Quantum ESPRESSO community for underlying DFT questions.


## Limitations & Known Constraints
- **Documentation**: Limited README/docs
- **Maintenance**: Research code status
- **Generality**: Focused on specific applications

## Comparison with Other Codes
- **vs CE-TDDFT**: Both QE-based RT-TDDFT; different focus areas
- **vs SALMON**: SALMON more mature with broader documentation
- **Unique aspect**: Nanostructure discharge specialization

## Application Areas
- Nanostructure electronics
- Discharge dynamics
- Non-equilibrium transport
- Ultrafast phenomena in nanomaterials

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/sheyua/RT-tddft

**Confidence**: VERIFIED
- Repository: ACCESSIBLE (GitHub)
- Code: Available and complete
- Status: Research-grade implementation

**Verification status**: ✅ VERIFIED
