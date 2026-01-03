# PyQMC (Python Quantum Monte Carlo)

## Official Resources
- Homepage: https://github.com/WagnerGroup/pyqmc
- Documentation: https://pyqmc.readthedocs.io/
- Source Repository: https://github.com/WagnerGroup/pyqmc
- License: MIT License

## Overview
PyQMC is a modern Python implementation of real-space quantum Monte Carlo methods for electronic structure calculations. Developed by the Wagner group, PyQMC emphasizes ease of use, rapid prototyping, and integration with the Python scientific ecosystem. The code implements VMC and DMC with a clean, object-oriented design that makes it accessible for both research and educational purposes while maintaining production capability.

**Scientific domain**: Quantum Monte Carlo, electronic structure, Python scientific computing  
**Target user community**: Python users, researchers, educators, method developers

## Theoretical Methods
- Variational Monte Carlo (VMC)
- Diffusion Monte Carlo (DMC)
- Slater-Jastrow wavefunctions
- Wavefunction optimization
- Real-space methods
- Pseudopotentials
- Periodic boundary conditions

## Capabilities (CRITICAL)
**Category**: Open-source QMC library (Python)
- VMC and DMC methods
- Pure Python implementation
- NumPy/SciPy based
- PySCF integration
- Slater-Jastrow wavefunctions
- Wavefunction optimization
- Molecules and solids
- Periodic systems
- Energy and forces
- Jupyter notebook friendly
- Educational applications
- Research and production

**Sources**: GitHub repository, documentation

## Key Strengths

### Python-Native:
- Pure Python
- NumPy/SciPy ecosystem
- Jupyter integration
- Readable code
- Easy modification

### Ease of Use:
- Clean API
- User-friendly
- Rapid prototyping
- Interactive workflows
- Educational value

### PySCF Integration:
- Trial wavefunctions from PySCF
- Python quantum chemistry
- Seamless workflow
- Modern tools

### Modern Design:
- Object-oriented
- Well-documented
- Unit tested
- Active development
- Community-driven

## Inputs & Outputs
- **Input formats**:
  - Python objects
  - PySCF wavefunctions
  - NumPy arrays
  - Structured data
  
- **Output data types**:
  - Python objects
  - NumPy arrays
  - Pandas DataFrames
  - HDF5 files
  - Visualization-ready

## Interfaces & Ecosystem

### PySCF:
- Native integration
- Trial wavefunctions
- Molecular orbitals
- Python workflow

### Python Scientific Stack:
- NumPy arrays
- Pandas DataFrames
- Matplotlib visualization
- Jupyter notebooks
- Modern Python tools

## Workflow and Usage

### Installation:
```bash
# Via pip
pip install pyqmc

# From source
git clone https://github.com/WagnerGroup/pyqmc.git
cd pyqmc
pip install -e .
```

### Basic VMC:
```python
import pyqmc
from pyscf import gto, scf

# Define molecule in PySCF
mol = gto.M(atom='H 0 0 0; H 0 0 1.4', basis='ccpvdz')
mf = scf.RHF(mol).run()

# Create QMC wavefunction from PySCF
wf = pyqmc.slater_jastrow(mol, mf)

# Run VMC
configs = pyqmc.initial_guess(mol, 1000)
df, configs = pyqmc.vmc(wf, configs, nsteps=1000)

print("VMC energy:", df['total'].mean())
```

### Wavefunction Optimization:
```python
# Optimize Jastrow parameters
wf, to_opt = pyqmc.default_sj(mol, mf)

# Variance minimization
pyqmc.line_minimization(
    wf, configs, to_opt,
    nsteps=1000
)
```

### DMC Calculation:
```python
# Run DMC
dmc_data, configs = pyqmc.dmc(
    wf, configs,
    nsteps=5000,
    tstep=0.01,
    accumulators={'energy': pyqmc.EnergyAccumulator(mol)}
)

print("DMC energy:", dmc_data['energytotal'].mean())
```

### Jupyter Workflow:
```python
# Interactive analysis
import pandas as pd
import matplotlib.pyplot as plt

# Visualize energy convergence
df.plot(y='total')
plt.ylabel('Energy (Ha)')
plt.xlabel('Step')
plt.show()

# Statistical analysis
print(df.describe())
```

## Advanced Features

### Custom Wavefunctions:
- Extensible design
- Custom trial functions
- Python implementation
- Research development

### Analysis Tools:
- Built-in accumulators
- Statistical analysis
- Visualization helpers
- Pandas integration

### Periodic Systems:
- PBC support
- Crystal calculations
- k-points
- Solid-state physics

## Performance Characteristics
- **Speed**: Python (slower than compiled)
- **Usability**: Excellent
- **Purpose**: Research, education, prototyping
- **Typical**: Desktop calculations, small to medium systems

## Computational Cost
- Python overhead
- Moderate performance
- Not HPC-optimized
- Best for development/education
- Small to medium production

## Limitations & Known Constraints
- **Performance**: Python overhead
- **System size**: Limited vs compiled codes
- **HPC**: Not optimized for supercomputers
- **Production**: Best for moderate systems
- **Speed**: Slower than QMCPACK/CASINO

## Comparison with Other QMC Codes
- **vs QMCPACK**: PyQMC Python/easy, QMCPACK HPC/fast
- **vs CASINO**: PyQMC accessible, CASINO feature-rich
- **Unique strength**: Python-native, ease of use, PySCF integration, Jupyter workflows, rapid prototyping, educational value

## Application Areas

### Rapid Prototyping:
- Method development
- Algorithm testing
- Research ideas
- Custom features
- Python workflows

### Education:
- Teaching QMC
- Learning tool
- Interactive demonstrations
- Student projects
- Code understanding

### Small Production:
- Molecules
- Small systems
- Validation studies
- Benchmarks
- Python-integrated workflows

## Best Practices

### Python Usage:
- Jupyter notebooks
- Interactive development
- NumPy efficiency
- Vectorization
- Profiling

### PySCF Integration:
- Quality trial wavefunctions
- Basis set choices
- SCF convergence
- Molecular orbitals

### Research:
- Rapid iteration
- Custom modifications
- Testing ideas
- Validation
- Prototyping

## Community and Support
- Open-source (MIT)
- Wagner Group (UIUC)
- GitHub repository
- Issue tracking
- Active development
- Growing community
- Educational focus

## Educational Resources
- Comprehensive docs
- Tutorials
- Jupyter notebooks
- Example scripts
- API reference
- Python QMC teaching

## Development
- University of Illinois
- Wagner research group
- Active development
- Community contributions
- Modern Python practices
- Regular updates

## Research Impact
PyQMC enables accessible QMC method development and education, bringing modern Python workflows to quantum Monte Carlo and lowering barriers to QMC research and learning.

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/WagnerGroup/pyqmc
2. Documentation: https://pyqmc.readthedocs.io/
3. Publications: J. Chem. Phys. 153, 174111 (2020)

**Secondary sources**:
1. QMC literature
2. PySCF ecosystem
3. User publications

**Confidence**: VERIFIED - Python QMC library

**Verification status**: ✅ VERIFIED
- GitHub: ACCESSIBLE
- License: MIT (open-source)
- **Category**: Open-source QMC library
- Status: Actively developed
- Institution: UIUC (Wagner Group)
- Specialized strength: Python-native quantum Monte Carlo, PySCF integration, user-friendly API, Jupyter workflows, educational value, rapid prototyping, modern design, NumPy-based, interactive development, accessible QMC for Python ecosystem
