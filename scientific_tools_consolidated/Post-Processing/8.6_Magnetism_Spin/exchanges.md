# exchanges

## Official Resources
- Source Repository: https://github.com/dkorotin/exchanges
- Documentation: Included in repository
- License: Open source

## Overview
**exchanges** is a Fortran code for calculating Heisenberg exchange parameters for magnetic compounds using the Green's function formalism within Density Functional Theory. It interfaces with Quantum ESPRESSO to compute exchange coupling constants (Jij) from first principles.

**Scientific domain**: Magnetic exchange parameters, Heisenberg model from DFT  
**Target user community**: Researchers computing magnetic exchange coupling constants from DFT for use in spin models

## Theoretical Methods
- Green's function formalism for exchange parameters
- Magnetic force theorem (Lichtenstein formula)
- Heisenberg exchange coupling (Jij)
- Dzyaloshinskii-Moriya interaction (Dij)
- Self-consistent spin-polarized DFT
- Multiple-scattering theory
- Korringa-Kohn-Rostoker (KKR) approach

## Capabilities (CRITICAL)
- Heisenberg exchange parameter (Jij) calculation
- Dzyaloshinskii-Moriya interaction (Dij) calculation
- Green's function approach
- QE interface for DFT
- Model Hamiltonian generation
- Magnon dispersion from Jij
- Critical temperature estimation
- Multiple neighbor shell calculation

**Sources**: GitHub repository, J. Phys.: Condens. Matter

## Key Strengths

### Green's Function Method:
- Accurate Jij from DFT
- Magnetic force theorem
- Systematic convergence
- Well-established formalism

### QE Integration:
- Direct interface with Quantum ESPRESSO
- Same pseudopotentials
- Consistent calculation flow
- Automated workflow

### Comprehensive Exchange:
- Isotropic exchange (Jij)
- Anisotropic exchange
- DM interaction (Dij)
- Multiple neighbor shells

## Inputs & Outputs
- **Input formats**:
  - QE output files
  - Exchange calculation parameters
  - k-point specifications
  
- **Output data types**:
  - Exchange parameters (Jij)
  - DM vectors (Dij)
  - Magnon dispersion
  - Critical temperature estimates

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: DFT backend
- **Spin dynamics codes**: Parameter input
- **Fortran**: Core computation

## Performance Characteristics
- **Speed**: Fast (post-processing of QE data)
- **Accuracy**: Good (force theorem)
- **System size**: Limited by QE
- **Memory**: Moderate

## Computational Cost
- **Jij calculation**: Minutes to hours
- **QE pre-requisite**: Hours (separate)
- **Typical**: Moderate

## Limitations & Known Constraints
- **QE only**: No VASP or other code support
- **Force theorem**: Approximate (perturbative)
- **Collinear reference**: Non-collinear may be less accurate
- **Documentation**: Limited
- **Fortran**: Less accessible than Python

## Comparison with Other Codes
- **vs TB2J**: exchanges uses Green's function, TB2J uses torque method
- **vs Jx_DMFT**: exchanges is DFT-only, Jx_DMFT includes DMFT
- **vs SPR-KKR**: exchanges uses QE, SPR-KKR is KKR code
- **Unique strength**: Green's function exchange parameters from QE, Lichtenstein formula

## Application Areas

### Magnetic Materials:
- Ferromagnets (Fe, Co, Ni)
- Antiferromagnets
- Ferrimagnets
- Heusler alloys

### Critical Temperature:
- Tc estimation from Jij
- Mean-field and RPA
- Comparison with experiment
- Composition dependence

### Spin Model Parameters:
- Input for Monte Carlo
- Input for spin dynamics
- Input for SpinW
- Input for UppASD

## Best Practices

### QE Calculation:
- Use well-converged SCF
- Adequate k-point density
- Include spin-orbit for DM
- Use consistent pseudopotentials

### Jij Convergence:
- Test convergence with neighbor shells
- Check k-point convergence in Green's function
- Validate against known systems
- Compare with experimental Tc

## Community and Support
- Open source on GitHub
- Developed by D. Korotin
- Research code
- Limited documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/dkorotin/exchanges
2. M. I. Katsnelson et al., Phys. Rev. B 61, 15522 (2000) (Lichtenstein formula)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: Included in repository
- Active development: Research code
- Specialized strength: Green's function Heisenberg exchange parameters from QE, Lichtenstein formula
