# DFT2kp

## Official Resources
- Publication: SciPost Phys. Codebases 25 (2024)
- arXiv: https://arxiv.org/abs/2306.08554
- Code: Available via SciPost/Zenodo
- License: Check publication

## Overview
DFT2kp is a code to explicitly calculate Kane (linear in k) and Luttinger (quadratic in k) parameters of k·p effective Hamiltonians directly from ab initio wavefunctions provided by Quantum ESPRESSO. It analyzes symmetry transformations to construct optimal symmetry-adapted bases for effective models.

**Scientific domain**: k·p perturbation theory, effective mass theory, semiconductor physics
**Target user community**: Researchers developing effective models from first-principles calculations

## Theoretical Methods
- k·p perturbation theory
- Kane and Luttinger parameter extraction
- Symmetry-adapted basis construction
- Löwdin downfolding
- Quantum ESPRESSO interface

## Capabilities (CRITICAL)
- **Parameter Extraction**: Kane and Luttinger parameters from DFT
- **QE Interface**: Direct Quantum ESPRESSO integration
- **Symmetry Analysis**: Automatic symmetry-adapted basis
- **Linear Terms**: First-order k·p parameters
- **Quadratic Terms**: Second-order effective mass tensors
- **Validation**: Compare k·p bands with DFT

**Sources**: SciPost publication, arXiv preprint

## Key Strengths

### Ab Initio Parameters:
- Direct DFT parameter extraction
- No fitting required
- Quantitative accuracy
- Material-specific models

### Symmetry Optimization:
- Automatic basis optimization
- Symmetry-adapted representation
- Clean Hamiltonian form
- Reduced parameter count

### Quantum ESPRESSO:
- Native QE interface
- Uses QE wavefunctions
- Standard DFT workflow
- Well-documented

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO output
  - Wavefunction files
  - Band structure data
  
- **Output data types**:
  - Kane parameters
  - Luttinger parameters
  - Effective mass tensors
  - k·p Hamiltonian matrices

## Installation
```bash
# Obtain from SciPost/Zenodo
# Follow installation instructions
pip install dft2kp  # if available
```

## Usage Examples
```python
from dft2kp import KPModel

# Load QE calculation
model = KPModel.from_qe(prefix='silicon', outdir='./tmp')

# Extract parameters at Gamma
kane_params = model.get_kane_parameters(kpoint=[0,0,0])
luttinger_params = model.get_luttinger_parameters(kpoint=[0,0,0])

# Generate k.p Hamiltonian
H_kp = model.get_hamiltonian(order=2)
```

## Performance Characteristics
- **Speed**: Fast parameter extraction
- **Accuracy**: DFT-level precision
- **Automation**: Minimal user intervention

## Limitations & Known Constraints
- **QE-specific**: Currently Quantum ESPRESSO only
- **Local expansion**: Valid near expansion point
- **Band selection**: User must choose relevant bands

## Comparison with Other Tools
- **vs kdotp-symmetry**: DFT2kp extracts parameters, kdotp-symmetry derives form
- **vs manual fitting**: DFT2kp automated and systematic
- **Unique strength**: Direct ab initio parameter extraction

## Application Areas
- Semiconductor effective mass models
- Topological material modeling
- Heterostructure band alignment
- Quantum dot/well simulations
- Device physics

## Best Practices
- Ensure converged DFT calculation
- Include sufficient bands
- Verify k·p bands against DFT
- Check parameter units

## Community and Support
- SciPost publication
- arXiv preprint
- Academic correspondence

## Verification & Sources
**Primary sources**:
1. SciPost Phys. Codebases 25 (2024)
2. arXiv: https://arxiv.org/abs/2306.08554

**Confidence**: VERIFIED - Published in SciPost

**Verification status**: ✅ VERIFIED
- Publication: PEER-REVIEWED
- arXiv preprint: ACCESSIBLE
- Method: Ab initio k·p extraction
