# Jx_DMFT

## Official Resources
- Source Repository: https://github.com/KAIST-ELST/Jx_DMFT
- Documentation: Included in repository
- License: Open source

## Overview
**Jx_DMFT** is a software for calculating magnetic exchange parameters (Jx) from the magnetic force theorem, combined with both DFT and dynamical mean-field theory (DMFT). It can compute magnon dispersion and spectral functions from the exchange parameters, including correlation effects beyond DFT.

**Scientific domain**: Magnetic exchange parameters with DMFT, correlated magnetism  
**Target user community**: Researchers studying magnetic properties of strongly correlated materials where DFT alone is insufficient

## Theoretical Methods
- Magnetic force theorem for exchange parameters
- Density Functional Theory (DFT)
- Dynamical Mean-Field Theory (DMFT)
- Heisenberg exchange coupling
- Magnon dispersion calculation
- Momentum-dependent spectral functions
- Self-consistent DMFT

## Capabilities (CRITICAL)
- Exchange parameter calculation with DFT
- Exchange parameter calculation with DMFT
- Self-consistent DMFT calculation
- Magnon dispersion from exchange parameters
- Local and momentum-dependent spectral functions
- Post-processing scripts for measurable quantities
- Multiple DFT code interfaces

**Sources**: GitHub repository

## Key Strengths

### DMFT for Correlated Systems:
- Beyond-DFT exchange parameters
- Correlation effects on magnetism
- Temperature-dependent exchange
- Mott physics included

### Comprehensive Output:
- Magnon dispersion
- Spectral functions
- Exchange parameters
- Temperature dependence

### DFT+DMFT Integration:
- Combines both methods
- Self-consistent calculation
- Multiple DFT backends
- Systematic improvement

## Inputs & Outputs
- **Input formats**:
  - DFT Hamiltonian files
  - DMFT self-energy data
  - Exchange calculation parameters
  
- **Output data types**:
  - Exchange parameters (Jij)
  - Magnon dispersion
  - Spectral functions
  - Temperature-dependent quantities

## Interfaces & Ecosystem
- **DFT codes**: Multiple interfaces
- **DMFT codes**: Self-energy input
- **Python**: Post-processing scripts
- **Fortran**: Core computation

## Performance Characteristics
- **Speed**: Depends on DMFT convergence
- **Accuracy**: Beyond DFT for correlated systems
- **System size**: Limited by DMFT
- **Memory**: Moderate to high

## Computational Cost
- **DFT Jij**: Hours
- **DMFT Jij**: Hours to days
- **Typical**: Expensive for DMFT

## Limitations & Known Constraints
- **DMFT complexity**: Requires DMFT expertise
- **Computational cost**: DMFT is expensive
- **Documentation**: Limited
- **Installation**: Complex (DFT+DMFT stack)

## Comparison with Other Codes
- **vs exchanges**: Jx_DMFT includes DMFT, exchanges is DFT-only
- **vs TB2J**: Jx_DMFT has DMFT, TB2J uses torque method
- **vs SPR-KKR**: Jx_DMFT is DFT+DMFT, SPR-KKR is KKR
- **Unique strength**: Exchange parameters with DMFT for correlated magnets, beyond-DFT magnetism

## Application Areas

### Strongly Correlated Magnets:
- Transition metal oxides
- Rare-earth compounds
- Heavy fermion systems
- Mott insulators

### DMFT Magnetism:
- Correlation-enhanced exchange
- Temperature-dependent magnetism
- Orbital-selective magnetism
- Spin-orbit coupling effects

### Magnon Spectroscopy:
- Correlated magnon dispersion
- Damping from DMFT
- Comparison with INS
- Temperature-dependent spectra

## Best Practices

### DMFT Setup:
- Use converged DMFT self-energy
- Validate against DFT results
- Test impurity solver convergence
- Compare with experimental spectra

### Exchange Calculation:
- Include sufficient neighbor shells
- Test k-point convergence
- Validate against known systems
- Compare DFT vs DMFT exchange

## Community and Support
- Open source on GitHub
- Developed at KAIST
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/KAIST-ELST/Jx_DMFT
2. Related publications from KAIST ELST group

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Exchange parameters with DMFT for correlated magnets, beyond-DFT magnetism
