# PAOFLOW (Projections of Atomic Orbitals FLOW)

## Official Resources
- Homepage: http://www.paoflow.org/
- Documentation: http://www.paoflow.org/
- Source Repository: https://github.com/marcobn/PAOFLOW
- License: GNU General Public License v3.0

## Overview
PAOFLOW is a Python-based post-processing tool for electronic structure calculations, designed to compute transport, topological, and optical properties from DFT calculations using tight-binding models constructed with atomic orbital projections. Developed primarily at the University of North Texas, PAOFLOW provides automated workflows for property calculations from first principles with emphasis on ease of use and comprehensive analysis capabilities.

**Scientific domain**: Electronic transport, topological properties, DFT post-processing  
**Target user community**: DFT users, transport calculations, materials properties

## Theoretical Methods
- Projections of atomic orbitals
- Tight-binding from DFT
- Boltzmann transport
- Berry phase properties
- Topological invariants
- Optical conductivity
- Spin textures

## Capabilities (CRITICAL)
**Category**: Open-source DFT post-processing tool
- DFT to TB conversion (atomic projections)
- Band structure interpolation
- Boltzmann transport (σ, S, κₑ)
- Berry curvature and AHC
- Topological invariants (Z2, Chern)
- Optical properties
- Spin Hall conductivity
- Spin textures
- DFT interface (Quantum ESPRESSO, others)
- Python implementation
- Automated workflows
- Production quality

**Sources**: Official website, GitHub, publications

## Key Strengths

### Comprehensive Properties:
- Transport coefficients
- Topological invariants
- Optical properties
- Berry phase phenomena
- Spin properties
- All-in-one tool

### Automated Workflows:
- DFT to properties
- Minimal user intervention
- Standard pipelines
- Production ready
- User-friendly

### Python-Based:
- Accessible interface
- NumPy/SciPy
- Visualization tools
- Extensible
- Modern design

### Atomic Projections:
- PAO-based TB construction
- DFT integration
- No Wannier90 required
- Alternative approach
- Flexible

## Inputs & Outputs
- **Input formats**:
  - DFT outputs (QE primarily)
  - Atomic orbital projections
  - Configuration files
  
- **Output data types**:
  - Transport coefficients
  - Topological invariants
  - Band structures
  - Optical conductivity
  - Spin textures
  - Berry curvature maps

## Interfaces & Ecosystem

### DFT Codes:
- Quantum ESPRESSO (primary)
- Extensible to others
- Atomic projections
- Standard workflow

### Visualization:
- Built-in plotting
- Matplotlib integration
- Publication-ready
- Property maps

## Workflow and Usage

### Installation:
```bash
# Via pip
pip install paoflow

# From source
git clone https://github.com/marcobn/PAOFLOW.git
cd PAOFLOW
python setup.py install
```

### Basic Workflow:
```python
from PAOFLOW import PAOFLOW

# Initialize
paoflow = PAOFLOW(
    savedir='./paoflow_output',
    verbose=True
)

# Read DFT data (Quantum ESPRESSO)
paoflow.read_atomic_proj_QE()

# Build Hamiltonian
paoflow.projectability()
paoflow.pao_hamiltonian()

# Band structure
paoflow.bands(
    ibrav=2,
    nk=1000,
    band_path=['L', 'G', 'X']
)

# Transport properties
paoflow.transport(
    tmin=100,   # K
    tmax=800,
    tstep=50
)

# Topological properties
paoflow.z2_pack()
paoflow.Berry_curvature()
paoflow.anomalous_Hall()

# Optical properties
paoflow.optical_conductivity()

# Generate output
paoflow.finish_execution()
```

### Quantum ESPRESSO Interface:
```bash
# 1. SCF calculation
pw.x < scf.in > scf.out

# 2. NSCF with projections
pw.x < nscf.in > nscf.out

# 3. Projections
projwfc.x < proj.in > proj.out

# 4. PAOFLOW processing
python paoflow_script.py
```

## Advanced Features

### Transport:
- Electrical conductivity
- Seebeck coefficient
- Electronic thermal conductivity
- Power factor
- Temperature dependence
- Boltzmann equation

### Topological:
- Z2 invariants
- Chern numbers
- Berry curvature
- Anomalous Hall conductivity
- Topological characterization

### Optical:
- Optical conductivity
- Dielectric function
- Absorption
- Interband transitions
- Frequency-dependent

### Spin:
- Spin Hall conductivity
- Spin textures
- Spin-orbit effects
- Rashba/Dresselhaus

## Performance Characteristics
- **Speed**: Fast post-processing
- **Accuracy**: DFT quality
- **Purpose**: Comprehensive properties
- **Typical**: Minutes to hours

## Computational Cost
- Post-DFT processing
- DFT calculation dominant
- Property calculations fast
- Production capable
- Automated

## Limitations & Known Constraints
- **Requires DFT**: Post-processing tool
- **PAO projections**: Alternative to Wannier
- **QE focus**: Primarily Quantum ESPRESSO
- **Python speed**: Moderate for large systems
- **Documentation**: Growing

## Comparison with Other Tools
- **vs Wannier90 ecosystem**: PAOFLOW PAO-based, W90 MLWF-based
- **vs BoltzTraP**: PAOFLOW comprehensive, BoltzTraP transport specialist
- **vs WannierBerri**: Different projection approach
- **Unique strength**: All-in-one DFT post-processing, automated workflows, PAO projections

## Application Areas

### Thermoelectrics:
- Transport coefficients
- Material screening
- Doping optimization
- ZT estimation

### Topological Materials:
- Topological classification
- Berry phase properties
- Anomalous Hall
- Material discovery

### Optical Properties:
- Absorption spectra
- Optical conductivity
- Dielectric response
- Spectroscopy theory

### Spintronics:
- Spin Hall effect
- Spin textures
- SOC materials
- Device applications

## Best Practices

### DFT Input:
- Quality SCF/NSCF
- Appropriate projections
- Sufficient k-points
- Converged calculations

### Workflows:
- Follow examples
- Standard pipelines
- Property selection
- Validation

### Analysis:
- Visualize results
- Physical interpretation
- Compare with experiment
- Systematic studies

## Community and Support
- Open-source (GPL v3)
- University of North Texas
- GitHub repository
- Documentation website
- Active development
- User community

## Educational Resources
- Official documentation
- Tutorial examples
- Example gallery
- Publication list
- Workflow guides

## Development
- Marco Buongiorno Nardelli (lead)
- University of North Texas
- Active development
- Regular updates
- Community contributions

## Research Impact
PAOFLOW provides comprehensive electronic structure property calculations from DFT, enabling systematic materials screening for transport, topological, and optical applications.

## Verification & Sources
**Primary sources**:
1. Homepage: http://www.paoflow.org/
2. GitHub: https://github.com/marcobn/PAOFLOW
3. Publications: Comp. Phys. Comm. 235, 415 (2019)

**Secondary sources**:
1. User publications
2. Materials property papers

**Confidence**: VERIFIED - DFT post-processing tool

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: DFT post-processing tool
- Status: Actively developed
- Institution: University of North Texas
- Specialized strength: Comprehensive property calculations from DFT using atomic orbital projections, transport (Boltzmann), topological invariants, optical properties, spin textures, automated workflows, Python-based, all-in-one post-processing, production quality, alternative to Wannier-based approaches
