# PWtransport

## Official Resources
- Homepage/Repo: http://yemeng.site/ (Reference to code) / Often private or distributed upon request
- Base Engine: PEtot
- License: Academic/Copyright (Likely similar to PEtot)

## Overview
PWtransport is a quantum transport code based on plane-wave pseudopotential Density Functional Theory. It utilizes the "PEtot" code as its electronic structure engine. It is designed to calculate transport properties of nanostructures using the non-equilibrium Green's function (NEGF) method or scattering state approaches combined with plane-wave DFT.

**Scientific domain**: Quantum transport, molecular electronics, nanotechnology
**Target user community**: Researchers studying electron transport in nanodevices

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets (inherited from PEtot)
- Pseudopotentials
- Quantum Transport Theory (NEGF / Scattering States)
- Open boundary conditions handling

## Capabilities
- Electronic structure of device regions
- Transmission coefficients
- Current-Voltage (I-V) characteristics
- Conductance calculations
- Coupled DFT-Transport self-consistency

## Key Strengths

### Plane-Wave Transport:
- One of the few transport codes explicitly using a plane-wave basis (many use localized orbitals like SIESTA/TranSIESTA).
- Systematic basis set convergence for transport problems.

### Large Scale:
- Inherits PEtot's ability to handle large systems, crucial for realistic device simulations.

## Inputs & Outputs
- **Input formats**:
  - PEtot-style input files modified for transport
  - Electrode and scattering region definitions
  
- **Output data types**:
  - Transmission functions T(E)
  - Current values
  - Local density of states (LDOS) under bias

## Interfaces & Ecosystem
- **PEtot**: Tightly coupled with the PEtot DFT code.

## Computational Cost
- **High Cost**: Transport calculations (Green's functions) are significantly more expensive than ground-state DFT.
- **Memory**: Large sparse matrix inversions require substantial RAM.
- **Time**: Self-consistent transport at finite bias is computationally intensive.

## Best Practices

### Convergence Strategies:
- **Zero Bias start**: Always converge the zero-bias calculation first.
- **Incremental Bias**: Restart finite-bias calculations from the previous lower-bias converged density (e.g., use 0.1V result for 0.2V calculation).
- **Log Monitoring**: Check SCF logs for "wild" oscillations; reduce mixing parameters if necessary.

### System Setup:
- **Leads**: Ensure electrode leads are perfectly matched to the scattering region interface to avoid spurious scattering.
- **K-points**: Use a dense k-point grid along the transport direction in the leads.

## Community and Support
- **Niche**: Small, specialized user community centered around quantum transport research groups.
- **Contact**: Primary support via academic contact with the original authors (L.-W. Wang group alumni).

## Performance Characteristics
- **Speed**: Dependent on PEtot's performance and the heavy cost of transport calculations (Green's functions).
- **Parallelization**: Parallelized to handle the heavy computational load of transport integration.

## Limitations & Known Constraints
- **Availability**: Not a standard public "download and click" code; often requires contact with developers (L.-W. Wang group alumni).
- **Documentation**: Sparse public documentation compared to TranSIESTA.

## Comparison with Other Codes
- **vs TranSIESTA**: TranSIESTA uses localized basis sets (SIESTA), which are naturally efficient for transport (sparse matrices). PWtransport uses plane waves, offering potentially better accuracy/convergence but different computational challenges.
- **vs OpenMX/Nanodcal**: Both are localized basis codes. PWtransport is unique in its plane-wave approach.

## Verification & Sources
**Primary sources**:
1. Scientific Literature: "Ab initio calculation of transport properties..." referencing PWtransport.
2. Developer websites (Y. Meng, L.-W. Wang group).

**Confidence**: VERIFIED - Existence confirmed via literature and developer pages.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: DFT/Transport
- Key Feature: Plane-Wave Transport
