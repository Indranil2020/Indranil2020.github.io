# kdotp-symmetry

## Official Resources
- Homepage: https://kdotp-symmetry.greschd.ch/
- GitHub: https://github.com/greschd/kdotp-symmetry
- Documentation: https://kdotp-symmetry.greschd.ch/
- PyPI: https://pypi.org/project/kdotp-symmetry/
- Publication: D. Gresch et al., Phys. Rev. Materials 2, 103805 (2018)
- License: Apache License 2.0

## Overview
kdotp-symmetry is a Python tool for calculating the general form of a k·p Hamiltonian under given symmetry constraints. It automatically derives the allowed terms in an effective Hamiltonian near a high-symmetry k-point by analyzing the little group symmetry operations and basis function transformations.

**Scientific domain**: k·p perturbation theory, effective mass theory, band structure modeling
**Target user community**: Researchers developing effective models for semiconductors and topological materials

## Theoretical Methods
- k·p perturbation theory
- Little group symmetry analysis
- Invariant theory for Hamiltonians
- Representation theory
- Löwdin partitioning basis

## Capabilities (CRITICAL)
- **Symmetry-constrained k·p**: Derive allowed Hamiltonian terms
- **Arbitrary Order**: k-expansion to any order
- **General Symmetries**: Works with any little group
- **Basis Functions**: Custom orbital/spin basis support
- **Spin-Orbit Coupling**: SOC-compatible formulation
- **Symbolic Output**: Mathematica-compatible expressions

**Sources**: kdotp-symmetry documentation, Phys. Rev. Materials publication

## Key Strengths

### Automatic Derivation:
- No manual symmetry analysis needed
- Systematic term generation
- Complete basis for allowed terms
- Verified by symmetry

### Flexible Framework:
- Any little group
- Arbitrary k-expansion order
- Custom basis functions
- SOC inclusion

### Integration:
- Works with Z2Pack
- TBmodels compatible
- Symbolic output
- Python ecosystem

## Inputs & Outputs
- **Input formats**:
  - Symmetry operations (rotation matrices)
  - Basis function representations
  - Expansion order specification
  
- **Output data types**:
  - Symbolic Hamiltonian terms
  - Coefficient matrices
  - Basis function list
  - Mathematica expressions

## Installation
```bash
pip install kdotp-symmetry
```

## Usage Examples
```python
import kdotp_symmetry as kp
import numpy as np

# Define symmetry operations (e.g., C3 rotation)
symmetries = [
    kp.SymmetryOperation(
        rotation_matrix=np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]]),
        repr_matrix=np.array([[np.exp(-1j*np.pi/3), 0], [0, np.exp(1j*np.pi/3)]])
    )
]

# Generate k.p Hamiltonian to second order
basis = kp.monomial_basis(dim=2)  # 2D k-space
hamiltonian = kp.symmetric_hamiltonian(
    symmetries=symmetries,
    kp_variable='k',
    order=2
)
print(hamiltonian)
```

## Performance Characteristics
- **Speed**: Fast symbolic computation
- **Accuracy**: Exact symmetry constraints
- **Scalability**: Higher orders computationally heavier

## Limitations & Known Constraints
- **Symbolic focus**: Not for numerical band structure
- **Input format**: Requires explicit symmetry matrices
- **Learning curve**: Representation theory knowledge helpful

## Comparison with Other Tools
- **vs kdotp-generator**: kdotp-symmetry original, kdotp-generator extends
- **vs MagneticKP**: kdotp-symmetry non-magnetic focus
- **vs manual derivation**: Automated and systematic
- **Unique strength**: Rigorous symmetry-based k·p derivation

## Application Areas
- Semiconductor effective mass models
- Topological semimetal models
- Weyl/Dirac point Hamiltonians
- Quantum well structures
- Band edge modeling

## Best Practices
- Verify symmetry operation matrices
- Check basis function representations
- Validate against known models
- Use appropriate expansion order

## Community and Support
- GitHub issue tracker
- Published methodology
- Z2Pack ecosystem

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/greschd/kdotp-symmetry
2. D. Gresch et al., Phys. Rev. Materials 2, 103805 (2018)

**Confidence**: VERIFIED - Published in Phys. Rev. Materials

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, Apache-2.0)
- Academic citations: Phys. Rev. Materials
- Active development: Maintained
