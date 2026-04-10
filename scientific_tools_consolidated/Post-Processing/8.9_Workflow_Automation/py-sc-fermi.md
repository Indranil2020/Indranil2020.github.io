# py-sc-fermi

## Official Resources
- Source Repository: https://github.com/bjmorgan/py-sc-fermi
- Documentation: https://py-sc-fermi.readthedocs.io/
- PyPI: https://pypi.org/project/py-sc-fermi/
- License: MIT License

## Overview
**py-sc-fermi** is a materials modelling code for calculating self-consistent Fermi energies and defect concentrations under thermodynamic equilibrium (or quasi-equilibrium) given defect formation energies. It determines the equilibrium Fermi level and defect/carrier concentrations from DFT defect calculations.

**Scientific domain**: Self-consistent defect thermodynamics, Fermi level determination  
**Target user community**: Researchers studying defect concentrations and Fermi levels in semiconductors and insulators

## Theoretical Methods
- Self-consistent Fermi level calculation
- Defect concentration thermodynamics
- Charge neutrality condition
- Mass-action law for defects
- Temperature-dependent concentrations
- Quasi-equilibrium defect chemistry
- Formation energy from DFT

## Capabilities (CRITICAL)
- Self-consistent Fermi energy calculation
- Defect concentration at given temperature
- Carrier concentration (n, p)
- Temperature-dependent defect concentrations
- Annealing temperature effects
- Multiple defect species
- Formation energy input from DFT
- Chemical potential constraints

**Sources**: GitHub repository

## Key Strengths

### Self-Consistent Solution:
- Solves charge neutrality self-consistently
- Fermi level determined by defect concentrations
- No a priori Fermi level assumption
- Physically consistent treatment

### Temperature Dependence:
- Defect concentrations vs temperature
- Carrier concentrations vs temperature
- Fermi level vs temperature
- Annealing simulation

### Comprehensive Defect Chemistry:
- Multiple defect species simultaneously
- Competing defects
- Charge state transitions
- Doping effects

## Inputs & Outputs
- **Input formats**:
  - Defect formation energies
  - Entropy data (optional)
  - Temperature range
  
- **Output data types**:
  - Self-consistent Fermi level
  - Defect concentrations
  - Carrier concentrations
  - Temperature-dependent plots

## Interfaces & Ecosystem
- **doped**: Compatible input format
- **pymatgen-analysis-defects**: Compatible
- **VASP**: DFT backend (via other tools)
- **Python**: Scripting

## Performance Characteristics
- **Speed**: Instant (analytical solution)
- **Accuracy**: Depends on input formation energies
- **System size**: Any number of defects
- **Memory**: Low

## Computational Cost
- **Calculation**: Seconds
- **DFT pre-requisite**: Hours (separate)
- **Typical**: Very efficient

## Limitations & Known Constraints
- **Equilibrium only**: No kinetic effects
- **Requires formation energies**: From DFT (separate)
- **No finite-size corrections**: Uses pre-corrected energies
- **Point defects only**: No extended defects

## Comparison with Other Codes
- **vs doped**: py-sc-fermi does Fermi-level SCF, doped does full workflow
- **vs pymatgen-analysis-defects**: py-sc-fermi focuses on concentrations
- **vs PyCDT**: py-sc-fermi is concentration-focused, PyCDT is broader (unmaintained)
- **Unique strength**: Self-consistent Fermi level and defect concentration calculation, temperature-dependent defect chemistry

## Application Areas

### Semiconductor Defects:
- Native defect concentrations
- Dopant activation
- Compensation mechanisms
- Carrier concentration prediction

### Thermoelectric Materials:
- Defect engineering
- Carrier optimization
- Doping strategies
- Concentration vs temperature

### Photovoltaics:
- Defect tolerance
- Fermi level engineering
- Carrier lifetime implications
- Stability assessment

### Battery Materials:
- Cation disorder
- Oxygen vacancy concentrations
- Electrolyte decomposition products
- Degradation mechanisms

## Best Practices

### Input Data:
- Use well-corrected formation energies
- Include all relevant defects
- Consider all charge states
- Use consistent DFT settings

### Temperature Range:
- Cover relevant temperature range
- Check convergence of SCF
- Compare with experimental data
- Consider annealing effects

## Community and Support
- Open source (MIT License)
- PyPI installation available
- ReadTheDocs documentation
- Developed by B. J. Morgan
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/bjmorgan/py-sc-fermi
2. Documentation: https://py-sc-fermi.readthedocs.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE (ReadTheDocs)
- PyPI: AVAILABLE
- Active development: Ongoing
- Specialized strength: Self-consistent Fermi level and defect concentration calculation, temperature-dependent defect chemistry
