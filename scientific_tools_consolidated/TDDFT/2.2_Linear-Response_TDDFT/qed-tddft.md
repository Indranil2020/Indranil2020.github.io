# qed-tddft

## Official Resources
- Homepage: https://github.com/cc-ats/qed-tddft
- Source Repository: https://github.com/cc-ats/qed-tddft
- Documentation: README with examples
- License: Open Source

## Overview
qed-tddft is a specialized Python package for Quantum-Electrodynamical Time-Dependent Density Functional Theory (QED-TDDFT) within Gaussian atomic basis sets. It enables the simulation of molecules strongly coupled to quantized electromagnetic field modes in optical cavities, capturing light-matter interactions at the quantum level. Built on top of PySCF, it provides a framework for cavity QED calculations in molecular systems.

**Scientific domain**: Cavity QED, polaritonic chemistry, strong light-matter coupling  
**Target user community**: Researchers in quantum optics, polaritonic chemistry, and cavity-modified molecular properties

## Theoretical Methods
- Time-Dependent Density Functional Theory (TDDFT)
- Quantum Electrodynamics (QED) coupling
- Pauli-Fierz Hamiltonian (TDDFT-PF)
- Gaussian atomic basis sets
- Cavity photon modes
- Light-matter coupling tensors
- Analytic energy gradients

## Capabilities
- Ground-state DFT with cavity coupling
- QED-TDDFT excited states
- Polaritonic state calculations
- Cavity frequency specification
- Cavity mode direction control
- Multiple excited state roots
- Analytic gradients for geometry optimization
- Integration with PySCF workflows

## Key Strengths

### Cavity QED Implementation:
- Full Pauli-Fierz Hamiltonian
- Multiple cavity modes
- Tunable coupling strength
- Arbitrary photon frequencies

### PySCF Integration:
- Leverages PySCF infrastructure
- Compatible with all PySCF basis sets
- Uses PySCF SCF methods
- Standard Python workflow

### Published Methodology:
- J. Chem. Phys. 155, 064107 (2021)
- J. Chem. Phys. 156, 124104 (2022)
- Peer-reviewed implementation

### Gradient Capability:
- Analytic energy gradients
- Geometry optimization in cavities
- Polaritonic potential energy surfaces

## Inputs & Outputs
- **Input formats**:
  - PySCF molecule objects
  - Standard basis set specifications
  - Cavity frequency arrays (NumPy)
  - Cavity mode vectors (NumPy)
  - XC functional specification
  
- **Output data types**:
  - Polaritonic excitation energies
  - Oscillator strengths
  - Cavity-matter coupling analysis
  - Gradients for optimization

## Interfaces & Ecosystem
- **Core dependency**:
  - PySCF (required)
  - NumPy for array operations

- **Workflow integration**:
  - Standard PySCF RKS/UKS objects
  - Any PySCF-supported XC functional
  - Any PySCF-supported basis set

## Example Usage
```python
from pyscf import gto, scf
import qed

mol = gto.Mole()
mol.atom = '''H 0 0 0; H 0 0 0.74'''
mol.basis = 'cc-pVDZ'
mol.build()

mf = scf.RKS(mol)
mf.xc = "b3lyp"
mf.kernel()

cavity_freq = numpy.asarray([0.200])
cavity_mode = numpy.asarray([[0.001, 0.0, 0.0]])

cav_model = qed.PF(mf, cavity_mode=cavity_mode, cavity_freq=cavity_freq)
td = qed.TDDFT(mf, cav_obj=cav_model)
td.nroots = 5
td.kernel()
```

## Performance Characteristics
- **Speed**: Depends on PySCF performance
- **System size**: Medium molecules (standard TDDFT limits)
- **Memory**: Standard TDDFT memory requirements
- **Scalability**: Single-node calculations

## Limitations & Known Constraints
- **Cavity modes**: Single or few modes typical
- **Coupling regime**: Strong coupling focus
- **Relativistic**: Non-relativistic only
- **Periodic**: Molecular systems only
- **Platform**: Requires PySCF installation

## Comparison with Other Codes
- **vs standard TDDFT**: Adds cavity QED coupling
- **vs OpenMolcas QED-CASSCF**: Different theory level (DFT vs multi-reference)
- **vs Molpro cavity**: Open-source alternative
- **Unique strength**: Gaussian basis QED-TDDFT with gradients

## Application Areas

### Polaritonic Chemistry:
- Cavity-modified reaction rates
- Polaritonic potential energy surfaces
- Ground state modification under strong coupling

### Optical Cavities:
- Molecule-cavity interactions
- Purcell effect simulations
- Cavity-induced energy splittings

### Spectroscopy:
- Modified absorption spectra
- Cavity-dressed molecular states
- Light-matter hybridization

## Best Practices
- Start with known cavity frequencies
- Test coupling strength convergence
- Compare with cavity-free TDDFT
- Use appropriate XC functionals

## Community and Support
- Open-source on GitHub (cc-ats organization)
- Published methodology with references
- Python 100%
- Academic development

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/cc-ats/qed-tddft
2. J. Yang et al., J. Chem. Phys. 155, 064107 (2021)
3. J. Yang et al., J. Chem. Phys. 156, 124104 (2022)

**Confidence**: VERIFIED - Published methodology with active GitHub

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Documentation: README with examples
- Academic citations: 2 J. Chem. Phys. papers
- Purpose: Research (cavity QED)
- Language: Python 100%
