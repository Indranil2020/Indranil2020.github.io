# kdotp-generator

## Official Resources
- GitHub: https://github.com/yjiang-iop/kdotp-generator
- Publication: Y. Jiang et al., related publications
- License: Check repository

## Overview
kdotp-generator is a tool for automatically generating k·p effective Hamiltonians, extending the kdotp-symmetry package. It includes support for magnetic space groups and provides pre-computed results for many magnetic and non-magnetic space groups, enabling rapid model development.

**Scientific domain**: k·p perturbation theory, magnetic materials, effective Hamiltonians
**Target user community**: Researchers studying magnetic topological materials and developing effective models

## Theoretical Methods
- k·p perturbation theory
- Magnetic space group symmetry
- Time-reversal symmetry handling
- Magnetic little group analysis
- Invariant subspace method

## Capabilities (CRITICAL)
- **Magnetic Groups**: Full magnetic space group support
- **Pre-computed Database**: Results for many groups available
- **kdotp-symmetry Extension**: Built on proven framework
- **Time-Reversal**: Proper TR symmetry handling
- **Automated Generation**: Minimal user input needed
- **Multiple Representations**: Various basis choices

**Sources**: GitHub repository, related publications

## Key Strengths

### Magnetic Support:
- All 1651 magnetic space groups
- Magnetic little groups
- Time-reversal combinations
- Antiferromagnetic systems

### Pre-computed Results:
- Database of common cases
- Quick lookup
- Validated results
- Reference implementations

### User-Friendly:
- Simple input format
- Automated workflow
- Example library
- Documentation

## Inputs & Outputs
- **Input formats**:
  - Space group number
  - Magnetic space group type
  - k-point specification
  - Basis information
  
- **Output data types**:
  - k·p Hamiltonian terms
  - Symmetry-allowed coefficients
  - Model parameters
  - Symbolic expressions

## Installation
```bash
git clone https://github.com/yjiang-iop/kdotp-generator.git
cd kdotp-generator
pip install -e .
```

## Usage Examples
```python
from kdotp_generator import generate_kp_hamiltonian

# Generate k.p Hamiltonian for magnetic space group
hamiltonian = generate_kp_hamiltonian(
    magnetic_space_group=62.441,  # Pnma with AFM ordering
    kpoint="Gamma",
    order=2
)

# Access pre-computed results
from kdotp_generator.database import get_precomputed
result = get_precomputed(msg_number=62, kpoint="X")
```

## Performance Characteristics
- **Speed**: Fast with pre-computed database
- **Coverage**: Extensive magnetic group support
- **Accuracy**: Symmetry-exact results

## Limitations & Known Constraints
- **kdotp-symmetry dependency**: Requires base package
- **Documentation**: May be limited
- **Magnetic focus**: Specialized for magnetic systems

## Comparison with Other Tools
- **vs kdotp-symmetry**: kdotp-generator adds magnetic groups
- **vs MagneticKP**: Different implementations, similar goals
- **Unique strength**: Magnetic space group database

## Application Areas
- Antiferromagnetic topological materials
- Magnetic Weyl semimetals
- Altermagnets
- Magnetic semiconductors
- Spin-orbit coupled magnets

## Best Practices
- Verify magnetic space group assignment
- Check basis conventions
- Validate against DFT bands
- Use appropriate order expansion

## Community and Support
- GitHub repository
- Academic development
- Related publications

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/yjiang-iop/kdotp-generator

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Source code: OPEN
- Method: Magnetic k·p generation
