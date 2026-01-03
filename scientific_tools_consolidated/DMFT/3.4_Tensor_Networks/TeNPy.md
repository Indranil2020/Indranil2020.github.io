# TeNPy (Tensor Network Python)

## Official Resources
- Homepage: https://tenpy.readthedocs.io/
- Documentation: https://tenpy.readthedocs.io/en/latest/
- Source Repository: https://github.com/tenpy/tenpy
- License: GNU General Public License v3.0

## Overview
TeNPy is a comprehensive Python library for tensor network algorithms, particularly focused on one-dimensional quantum systems and DMRG (Density Matrix Renormalization Group). Developed with emphasis on accessibility and ease of use, TeNPy provides a pure Python implementation of state-of-the-art tensor network methods with excellent documentation and a clean, modular design. The library is ideal for researchers and students learning tensor network methods while maintaining production-quality capabilities.

**Scientific domain**: Tensor networks, DMRG, 1D quantum systems  
**Target user community**: Python users, researchers, educators, tensor network learners

## Theoretical Methods
- Density Matrix Renormalization Group (DMRG)
- Matrix Product States (MPS)
- Matrix Product Operators (MPO)
- Time evolution (TEBD, TDVP)
- Infinite-size systems (iDMRG)
- Finite-temperature MPS
- Purification methods
- Variational algorithms
- Entanglement measures

## Capabilities (CRITICAL)
**Category**: Open-source tensor network library (Python)
- DMRG (finite and infinite)
- Ground and excited states
- Time evolution (TEBD, TDVP)
- Finite-temperature methods
- Built-in models
- Custom Hamiltonians
- Correlation functions
- Entanglement analysis
- Conserved quantum numbers
- Pure Python implementation
- NumPy/SciPy based
- Jupyter-friendly
- Production quality
- Educational value

**Sources**: Official documentation, GitHub, publications

## Key Strengths

### Python-Native:
- Pure Python
- NumPy/SciPy ecosystem
- Easy to learn
- Jupyter integration
- Readable code

### Excellent Documentation:
- Comprehensive tutorials
- User guide
- Algorithm explanations
- Example notebooks
- API reference

### Modular Design:
- Clean architecture
- Extensible
- Object-oriented
- Well-tested
- Maintainable

### Educational Value:
- Learning tensor networks
- Teaching tool
- Clear implementations
- Algorithm transparency
- Student-friendly

## Inputs & Outputs
- **Input formats**:
  - Python code
  - Model definitions
  - Parameters via dictionaries
  - HDF5 loading
  
- **Output data types**:
  - MPS objects
  - NumPy arrays
  - Observables
  - Correlation functions
  - HDF5 exports

## Interfaces & Ecosystem

### Python Scientific Stack:
- NumPy arrays
- SciPy sparse matrices
- Matplotlib visualization
- Jupyter notebooks
- Pandas (analysis)

### Data Management:
- HDF5 support
- Pickle serialization
- Parameter management
- Result organization

## Workflow and Usage

### Installation:
```bash
# Via pip
pip install physics-tenpy

# From source
git clone https://github.com/tenpy/tenpy.git
cd tenpy
pip install -e .
```

### Basic DMRG:
```python
import tenpy
from tenpy.networks.mps import MPS
from tenpy.models.tf_ising import TFIChain
from tenpy.algorithms import dmrg

# Define model
model_params = {
    'L': 20,
    'J': 1.0,
    'g': 1.5,
    'bc_MPS': 'finite'
}
M = TFIChain(model_params)

# Initial state
psi = MPS.from_product_state(M.lat.mps_sites(), [0]*20, bc='finite')

# DMRG
dmrg_params = {
    'trunc_params': {'chi_max': 100, 'svd_min': 1.e-10},
    'max_sweeps': 10,
    'verbose': 1
}
info = dmrg.run(psi, M, dmrg_params)

print("Energy:", info['E'])
```

### Custom Hamiltonian:
```python
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import SpinHalfSite

class MyModel(CouplingMPOModel):
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.5)
        
        # Nearest-neighbor
        for u in range(len(self.lat.unit_cell)):
            self.add_coupling(J, u, 'Sz', u, 'Sz', 1)
            self.add_onsite(g, u, 'Sx')
```

### Time Evolution:
```python
from tenpy.algorithms import tebd

# TEBD engine
tebd_params = {
    'N_steps': 10,
    'dt': 0.1,
    'order': 2,
    'trunc_params': {'chi_max': 50}
}

engine = tebd.TEBDEngine(psi, M, tebd_params)
engine.run()

# Measure observables over time
Sz = psi.expectation_value('Sz')
```

### Infinite System:
```python
from tenpy.algorithms import dmrg

# Infinite DMRG
model_params = {'L': 2, 'bc_MPS': 'infinite'}
M = TFIChain(model_params)

psi = MPS.from_product_state(M.lat.mps_sites(), [0,1], bc='infinite')
info = dmrg.run(psi, M, dmrg_params)
```

## Advanced Features

### Conserved Quantum Numbers:
- U(1) charge
- Z2 symmetries
- Custom conserved quantities
- Efficiency improvements
- Automatic handling

### Correlation Functions:
- Two-point correlations
- Multi-point functions
- Connected correlations
- Structure factors
- Analysis tools

### Entanglement:
- Entanglement entropy
- Entanglement spectrum
- Mutual information
- Topological invariants
- Quantum information

### Purification:
- Finite-temperature MPS
- Thermal states
- Mixed states
- Density matrices
- Temperature evolution

## Performance Characteristics
- **Speed**: Python (moderate)
- **Accuracy**: DMRG precision
- **System size**: 100s sites (1D)
- **Purpose**: Research and education
- **Typical**: Desktop to small HPC

## Computational Cost
- Python overhead
- NumPy efficiency
- Moderate performance
- Good for most applications
- Not fastest but accessible

## Limitations & Known Constraints
- **Performance**: Python slower than C++
- **1D focus**: Primarily 1D systems
- **2D**: Limited 2D capability
- **HPC**: Not optimized for supercomputers
- **System size**: Moderate vs C++ codes

## Comparison with Other Tensor Network Codes
- **vs ITensor**: TeNPy pure Python, ITensor C++/Julia
- **vs Block**: TeNPy accessible, Block specialized
- **Unique strength**: Pure Python, excellent docs, educational value, Jupyter-friendly, ease of learning, clean code

## Application Areas

### 1D Quantum Systems:
- Spin chains
- 1D fermions
- Bosonic systems
- Quantum wires
- Topological phases

### Research:
- Tensor network methods
- Algorithm development
- Many-body physics
- Quantum simulation
- Method testing

### Education:
- Learning DMRG
- Teaching tensor networks
- Student projects
- Algorithm understanding
- Interactive demonstrations

### Production:
- 1D calculations
- Moderate systems
- Rapid prototyping
- Python workflows

## Best Practices

### Python Usage:
- Jupyter notebooks
- Interactive development
- NumPy efficiency
- Documentation reading
- Example following

### DMRG Parameters:
- Appropriate chi_max
- Convergence criteria
- Truncation errors
- Sweeping strategy
- Validation

### Model Building:
- Use built-in models when possible
- Extend CouplingMPOModel
- Test on small systems
- Clear parameter definitions

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Active development
- Issue tracking
- User forum
- Responsive developers
- Growing community

## Educational Resources
- Excellent documentation
- Tutorial notebooks
- User guide
- Algorithm descriptions
- Example gallery
- API reference
- Physics background

## Development
- Academic project (Germany)
- Active development
- Community contributions
- Regular releases
- Feature additions
- Bug fixes
- Modern Python practices

## Research Impact
TeNPy has become a popular choice for tensor network calculations in Python, particularly valued for teaching and learning tensor network methods while maintaining research-quality capabilities.

## Verification & Sources
**Primary sources**:
1. Homepage: https://tenpy.readthedocs.io/
2. GitHub: https://github.com/tenpy/tenpy
3. Publications: SciPost Phys. Lect. Notes 5 (2018)

**Secondary sources**:
1. Tensor network literature
2. User publications
3. DMRG tutorials

**Confidence**: VERIFIED - Python tensor network library

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Open-source tensor network library
- Status: Actively developed
- Community: Growing, international
- Specialized strength: Pure Python tensor networks, DMRG implementation, excellent documentation, educational value, Jupyter-friendly, clean architecture, 1D quantum systems, accessible learning curve, production-capable, comprehensive tutorials, NumPy-based, modular design
