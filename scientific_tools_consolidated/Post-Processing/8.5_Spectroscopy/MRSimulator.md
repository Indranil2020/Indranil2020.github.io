# MRSimulator

## Official Resources
- Homepage: https://mrsimulator.readthedocs.io/
- GitHub: https://github.com/deepanshs/mrsimulator
- Documentation: https://mrsimulator.readthedocs.io/
- Publication: D. Srivastava et al., J. Chem. Phys. 161, 212501 (2024)
- License: BSD 3-Clause License

## Overview
MRSimulator is a fast, versatile, and open-source Python package for simulating one- and higher-dimensional solid-state NMR spectra. It supports static, magic-angle spinning (MAS), and variable-angle spinning (VAS) conditions for nuclei experiencing chemical shift and quadrupolar coupling interactions.

**Scientific domain**: Solid-state NMR spectroscopy simulation
**Target user community**: NMR spectroscopists and materials scientists

## Theoretical Methods
- Chemical shift (nuclear shielding) tensors
- Quadrupolar coupling interactions
- Magic-angle spinning (MAS)
- Variable-angle spinning (VAS)
- Static powder patterns
- Multi-dimensional spectra

## Capabilities (CRITICAL)
- **1D/2D NMR**: One and two-dimensional spectra
- **MAS/VAS/Static**: Multiple spinning conditions
- **Quadrupolar**: First and second-order effects
- **Chemical Shift**: Full tensor treatment
- **Efficient**: OpenCL GPU acceleration
- **Python API**: Scriptable interface
- **Visualization**: Built-in plotting

**Sources**: MRSimulator documentation, JCP publication

## Key Strengths

### Performance:
- GPU acceleration via OpenCL
- Fast powder averaging
- Efficient algorithms
- Large-scale simulations

### Versatility:
- Multiple spinning modes
- Various interactions
- Multi-dimensional
- Extensible framework

### Python Ecosystem:
- NumPy/SciPy based
- Jupyter compatible
- Modern API
- Active development

## Inputs & Outputs
- **Input formats**:
  - Python objects
  - JSON/YAML parameters
  - CIF structures (via integration)
  
- **Output data types**:
  - Simulated spectra
  - NumPy arrays
  - Matplotlib figures
  - CSDM data format

## Installation
```bash
pip install mrsimulator
# With GPU support
pip install mrsimulator[gpu]
```

## Usage Examples
```python
from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.method.lib import BlochDecaySpectrum

# Create spin system
site = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.0,
    shielding_symmetric={"zeta": 59.8, "eta": 0.62}
)
spin_system = SpinSystem(sites=[site])

# Set up MAS experiment
method = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=9.4,
    rotor_frequency=5000
)

# Run simulation
sim = Simulator(spin_systems=[spin_system], methods=[method])
sim.run()
```

## Performance Characteristics
- **Speed**: GPU-accelerated
- **Memory**: Efficient algorithms
- **Scalability**: Large systems supported

## Limitations & Known Constraints
- **Solid-state focus**: Not for solution NMR
- **Spin interactions**: Limited to specific types
- **Learning curve**: NMR concepts required
- **GPU setup**: OpenCL configuration needed

## Comparison with Other Tools
- **vs SIMPSON**: MRSimulator Python-native, SIMPSON older
- **vs SpinEvolution**: Different implementation
- **vs DMFit**: MRSimulator open-source, scriptable
- **Unique strength**: Python ecosystem, GPU acceleration

## Application Areas
- Solid-state NMR analysis
- Materials characterization
- Glasses and disordered systems
- Pharmaceuticals
- Battery materials

## Best Practices
- Validate parameters against experiments
- Use appropriate powder averaging
- Check convergence of integration
- Compare with experimental spectra

## Community and Support
- GitHub repository
- ReadTheDocs documentation
- BSD licensed
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/deepanshs/mrsimulator
2. D. Srivastava et al., J. Chem. Phys. 161, 212501 (2024)
3. Documentation: https://mrsimulator.readthedocs.io/

**Confidence**: VERIFIED - Published in J. Chem. Phys.

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (BSD-3)
- Academic citations: Growing
- Active development: Maintained
