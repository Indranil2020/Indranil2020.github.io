# ACE1pack

## Official Resources
- Source Repository: https://github.com/ACEsuit/ACE1pack
- Documentation: https://acesuit.github.io/ACE1pack.jl/
- License: Open source (MIT)

## Overview
**ACE1pack.jl** is a Julia package providing convenience functionality for fitting Atomic Cluster Expansion (ACE) interatomic potentials. It integrates ACE1.jl, ACEfit.jl, and JuLIP.jl for a complete fitting workflow.

**Scientific domain**: ACE potential fitting in Julia  
**Target user community**: Researchers fitting ACE potentials with Julia ecosystem

## Theoretical Methods
- Atomic Cluster Expansion (ACE)
- Linear and nonlinear ACE fitting
- Bayesian regression
- JuLIP atomistic simulation
- ACE basis construction

## Capabilities (CRITICAL)
- ACE basis construction
- Linear/nonlinear fitting
- Bayesian regression
- JuLIP integration
- LAMMPS export
- Multi-species support

**Sources**: GitHub repository

## Key Strengths

### ACE Framework:
- Complete ACE fitting workflow
- Linear and nonlinear models
- Systematic convergence
- Julia performance

### Julia Ecosystem:
- JuLIP atomistic simulation
- ACE1.jl core
- ACEfit.jl fitting
- Efficient computation

### Integration:
- LAMMPS potential export
- ASE data reading
- Multi-format I/O

## Inputs & Outputs
- **Input formats**: Training data (extxyz, JSON)
- **Output data types**: ACE potentials, LAMMPS files

## Interfaces & Ecosystem
- **JuLIP**: Atomistic simulation
- **Julia**: Core language
- **LAMMPS**: MD engine

## Performance Characteristics
- **Speed**: Fast (Julia)
- **Accuracy**: ACE-level (systematic)
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Fitting**: Minutes to hours
- **MD**: Fast (linear ACE)

## Limitations & Known Constraints
- **Julia required**: Not Python
- **Learning curve**: Julia ecosystem
- **LAMMPS export**: Format conversion needed
- **Documentation**: Could be more extensive

## Comparison with Other Codes
- **vs ACEpot (Python)**: ACE1pack is Julia, ACEpot is Python
- **vs PACE (LAMMPS)**: ACE1pack is fitting, PACE is evaluation
- **vs MACE**: ACE1pack is ACE only, MACE is equivariant NN
- **Unique strength**: Complete Julia-based ACE fitting workflow with systematic convergence

## Application Areas

### ACE Fitting:
- Metallic systems
- Molecular systems
- Multi-component alloys
- Systematic accuracy improvement

### Research:
- ACE basis development
- Fitting methodology
- Benchmark studies

## Best Practices
- Start with linear ACE
- Increase basis size systematically
- Validate with phonons and elastic constants
- Use Bayesian fitting for uncertainty

## Community and Support
- Open source (MIT)
- ACEsuit maintained
- Julia ecosystem
- Documentation available

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ACEsuit/ACE1pack

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Complete Julia-based ACE fitting workflow with systematic convergence
