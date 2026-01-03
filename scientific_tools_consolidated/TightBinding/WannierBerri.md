# WannierBerri

## Official Resources
- Homepage: https://wannier-berri.org/
- Documentation: https://wannier-berri.org/
- Source Repository: https://github.com/stepan-tsirkin/wannier-berri
- License: GNU General Public License v2.0

## Overview
WannierBerri is a Python code for calculating Berry curvature-related properties of materials using Wannier tight-binding models. Developed by Stepan Tsirkin, WannierBerri provides efficient and accurate calculations of anomalous Hall conductivity, orbital magnetization, shift currents, and other Berry phase phenomena from first principles. The code uses adaptive mesh refinement and modern algorithms to achieve high precision in Berry curvature integration.

**Scientific domain**: Berry phase physics, transport properties, topological materials  
**Target user community**: Wannier users, Berry physics researchers, transport calculations

## Theoretical Methods
- Berry curvature calculations
- Anomalous Hall effect
- Orbital magnetization
- Berry curvature dipole
- Shift currents
- Nonlinear Hall effect
- Gyrotropic effects
- Modern theory of polarization
- Adaptive refinement

## Capabilities (CRITICAL)
**Category**: Open-source Berry phase calculation tool
- Berry curvature integration
- Anomalous Hall conductivity (AHC)
- Orbital magnetization
- Berry curvature dipole
- Shift current
- Injection current
- Nonlinear Hall conductivity
- Gyrotropic magnetic effect (GME)
- Second harmonic generation
- Fermi surface properties
- Spin Hall conductivity
- Wannier interpolation
- Adaptive mesh refinement
- Python implementation
- Production quality

**Sources**: Official website, GitHub, publications

## Key Strengths

### Berry Phase Specialist:
- Comprehensive Berry properties
- High-precision integration
- Adaptive algorithms
- Modern theory implementation
- Production quality

### Python-Based:
- User-friendly
- NumPy/SciPy ecosystem
- Easy integration
- Extensible
- Rapid development

### Adaptive Refinement:
- Intelligent k-mesh
- Accuracy control
- Efficient computation
- Peak detection
- Convergence guaranteed

### Wannier90 Integration:
- Direct tb.dat input
- Standard workflow
- ab-initio to Berry physics
- Seamless pipeline

## Inputs & Outputs
- **Input formats**:
  - Wannier90 tb.dat files
  - Python configuration
  - System parameters
  
- **Output data types**:
  - Transport tensors
  - Berry curvature maps
  - Fermi surface data
  - Tabulated results
  - NumPy arrays
  - Publication plots

## Interfaces & Ecosystem

### Wannier90:
- Native tb.dat reading
- Tight-binding input
- Standard integration
- DFT workflow

### Python Scientific:
- NumPy arrays
- Matplotlib visualization
- Jupyter notebooks
- Data analysis tools

## Workflow and Usage

### Installation:
```bash
# Via pip
pip install wannierberri

# From source
git clone https://github.com/stepan-tsirkin/wannier-berri.git
cd wannier-berri
pip install -e .
```

### Basic Usage:
```python
import wannierberri as wberri
import numpy as np

# Load tight-binding model
system = wberri.System_tb(
    tb_file='wannier90_tb.dat',
    getAA=True  # Get Berry connection
)

# Setup k-mesh
grid = wberri.Grid(system, NK=[10,10,10])

# Calculate anomalous Hall conductivity
ahc = wberri.tabulate(
    system,
    grid=grid,
    quantities=['ahc'],
    fout_name='ahc'
)

print("AHC (S/cm):", ahc.data['ahc'])
```

### Adaptive Refinement:
```python
# Adaptive grid for high precision
ahc_adaptive = wberri.integrate(
    system,
    grid=wberri.Grid(system, NK=[10,10,10]),
    Efermi=np.linspace(-0.5, 0.5, 101),
    quantities=['ahc'],
    adpt_num_iter=10,  # Adaptive iterations
    fout_name='ahc_adaptive'
)
```

### Fermi Surface:
```python
# Calculate Fermi surface properties
fs = wberri.tabulate(
    system,
    grid=grid,
    quantities=['dos', 'cumdos'],
    Efermi=np.linspace(-2, 2, 201)
)

# Plot
import matplotlib.pyplot as plt
plt.plot(fs.Efermi, fs.data['dos'])
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
plt.show()
```

### Multiple Properties:
```python
# Calculate multiple Berry properties
results = wberri.integrate(
    system,
    grid=grid,
    Efermi=0.0,
    quantities=[
        'ahc',           # Anomalous Hall conductivity
        'morb',          # Orbital magnetization  
        'berry_dipole',  # Berry curvature dipole
        'spin_current'   # Spin Hall conductivity
    ]
)
```

## Advanced Features

### Adaptive K-Mesh:
- Automatic refinement
- Peak detection
- Error control
- Convergence criteria
- Efficiency optimization

### Symmetry:
- Symmetry exploitation
- Reduced k-mesh
- Faster calculations
- Automatic detection

### Parallel Execution:
- Multiprocessing support
- Efficient k-point distribution
- Scalable calculations
- Production HPC

### Fermi Surface Analysis:
- DOS calculations
- Fermi surface visualization
- k-resolved properties
- Integration over FS

## Performance Characteristics
- **Speed**: Efficient Python/NumPy
- **Accuracy**: Adaptive high precision
- **System size**: Any (uses TB model)
- **Purpose**: Berry phase calculations
- **Typical**: Minutes to hours

## Computational Cost
- Post-Wannier90 processing
- k-mesh dependent
- Adaptive refinement efficient
- Parallel capable
- Production suitable

## Limitations & Known Constraints
- **Requires Wannier90**: Tight-binding input
- **Python speed**: Slower than Fortran
- **TB quality**: Depends on MLWFs
- **Adaptive**: May need tuning
- **Learning curve**: Berry phase physics

## Comparison with Other Tools
- **vs BoltzWann**: WannierBerri Berry focus, BoltzWann Boltzmann
- **vs WannierTools**: WannierBerri transport, WannierTools topology
- **Unique strength**: Adaptive Berry curvature integration, comprehensive Berry properties, Python ease-of-use

## Application Areas

### Anomalous Hall Effect:
- Intrinsic AHC
- Ferromagnets
- Magnetic topological materials
- Berry curvature mechanism
- Experimental comparison

### Orbital Magnetization:
- Magnetic moments
- Topological materials
- Berry phase contribution
- Modern theory

### Nonlinear Transport:
- Berry curvature dipole
- Nonlinear Hall effect
- Shift currents
- Second-order response
- Photocurrents

### Materials Screening:
- High-throughput AHC
- Property prediction
- Materials discovery
- Database generation

## Best Practices

### Wannier Input:
- Quality tight-binding
- Validated band structure
- Appropriate energy windows
- Converged MLWFs

### k-Mesh Convergence:
- Start coarse, refine
- Use adaptive refinement
- Check convergence
- Balance accuracy/cost

### Berry Properties:
- Physical interpretation
- Symmetry checks
- Energy window selection
- Comparison with experiment

## Community and Support
- Open-source (GPL v2)
- GitHub repository
- Active development
- Issue tracking
- Documentation
- User community
- Publications

## Educational Resources
- Comprehensive documentation
- Tutorials
- Example notebooks
- API reference
- Berry phase theory
- Publication list

## Development
- Stepan Tsirkin (lead, Zurich)
- Active development
- Regular updates
- Community contributions
- Feature additions
- Python ecosystem

## Research Impact
WannierBerri enables high-precision calculations of Berry curvature-related properties, advancing understanding of topological transport, anomalous Hall effects, and nonlinear optical phenomena in materials.

## Verification & Sources
**Primary sources**:
1. Homepage: https://wannier-berri.org/
2. GitHub: https://github.com/stepan-tsirkin/wannier-berri
3. Publications: npj Comput. Mater. 7, 33 (2021)

**Secondary sources**:
1. Berry phase literature
2. Anomalous Hall effect papers
3. User publications

**Confidence**: VERIFIED - Berry phase calculation tool

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v2 (open-source)
- **Category**: Open-source Berry phase tool
- Status: Actively developed
- Specialized strength: Adaptive Berry curvature integration, anomalous Hall conductivity, orbital magnetization, Berry curvature dipole, nonlinear transport properties, Python-based, Wannier90 integration, high-precision algorithms, comprehensive Berry phase calculations, production quality
