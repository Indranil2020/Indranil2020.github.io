# JARVIS-Tools

## Official Resources
- **Homepage**: https://pages.nist.gov/jarvis/
- **GitHub**: https://github.com/usnistgov/jarvis
- **Documentation**: https://jarvis-tools.readthedocs.io/
- **PyPI**: https://pypi.org/project/jarvis-tools/
- **Web Portal**: https://jarvis.nist.gov/
- **License**: NIST Software License

## Overview
JARVIS-Tools (Joint Automated Repository for Various Integrated Simulations) is a comprehensive open-source Python package for atomistic data-driven materials design developed at NIST. It provides tools for setting up calculations, post-processing, analysis, visualization, and machine learning across multiple simulation codes. The package is tightly integrated with the JARVIS databases containing >100,000 materials.

**Scientific domain**: High-throughput materials science, DFT post-processing, machine learning, database-driven research
**Target user community**: Materials scientists, computational researchers, ML practitioners in materials science

## Theoretical Background
JARVIS-Tools interfaces with multiple levels of theory:
- DFT calculations (VASP, QE, Wien2k)
- Tight-binding (Wannier90, WTBH)
- Classical MD (LAMMPS)
- Machine learning potentials
- Property predictions from descriptors

## Capabilities (CRITICAL)
- **Workflow Automation**: VASP, QE, Wien2k, WTBH, Wannier90, LAMMPS
- **Band Structure**: Electronic band analysis and plotting
- **DOS/PDOS**: Density of states processing
- **Database Access**: JARVIS-DFT (>75,000 materials), JARVIS-FF
- **Machine Learning**: ALIGNN, CGCNN, descriptors
- **Structure Analysis**: Symmetry, defects, surfaces
- **Property Prediction**: Formation energy, band gap, elastic constants
- **Visualization**: Structure and property plotting

## Key Strengths

### Multi-Code Support:
- VASP (comprehensive)
- Quantum ESPRESSO
- Wien2k
- BoltzTraP
- Wannier90
- LAMMPS
- GPAW
- WTBH (Wannier tight-binding)

### Database Integration:
- JARVIS-DFT: 75,000+ 3D materials
- JARVIS-2D: 2D materials database
- JARVIS-FF: Force field database
- JARVIS-ML: ML model repository
- REST API access

### Machine Learning:
- ALIGNN (graph neural network)
- CGCNN integration
- Descriptor generation
- Pre-trained models
- Property prediction

### High-Throughput:
- Automated workflow generation
- Error handling
- Job management
- Result parsing

## Inputs & Outputs
- **Input formats**:
  - POSCAR, CIF, XYZ structures
  - VASP input/output files
  - QE input/output files
  - Database queries (materials ID)
  
- **Output data types**:
  - Band structures, DOS
  - Elastic tensors
  - Optical properties
  - ML predictions
  - JSON/CSV exports

## Interfaces & Ecosystem
- **Python integration**:
  - NumPy, Pandas for data
  - Matplotlib for visualization
  - PyTorch for ML models
  - ASE compatibility
  
- **Framework compatibility**:
  - Jupyter notebooks
  - REST API
  - Command-line tools
  - Web interface

## Installation
```bash
pip install jarvis-tools
```

With ML dependencies:
```bash
pip install jarvis-tools[ai]
```

## Usage Examples
```python
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import data
from jarvis.tasks.vasp.vasp import VaspJob

# Load structure from database
dft_3d = data("dft_3d")
atoms = Atoms.from_dict(dft_3d[0]["atoms"])

# Get band structure
from jarvis.io.vasp.outputs import Vasprun
vrun = Vasprun("vasprun.xml")
bands = vrun.get_bandstructure()

# ML prediction with ALIGNN
from jarvis.core.atoms import Atoms
from alignn.pretrained import get_figshare_model
model = get_figshare_model("jv_formation_energy_peratom_alignn")
prediction = model.predict(atoms)
```

## Performance Characteristics
- **Speed**: Efficient file parsing
- **Scalability**: High-throughput ready
- **ML**: GPU-accelerated models
- **Database**: Fast API queries

## Limitations & Known Constraints
- **Learning curve**: Many features require exploration
- **Documentation**: Extensive but distributed
- **VASP-centric**: Best support for VASP
- **Dependencies**: ML features need PyTorch

## Comparison with Other Tools
- **vs pymatgen**: JARVIS has ML focus, database integration
- **vs ASE**: JARVIS specialized for high-throughput, ML
- **vs atomate**: Different workflow philosophy
- **Unique strength**: NIST databases, ALIGNN ML, comprehensive ecosystem

## Application Areas
- High-throughput materials screening
- Machine learning for materials
- Database-driven discovery
- Property prediction
- Workflow automation
- 2D materials research
- Defect calculations

## Best Practices
- Use database for initial screening
- Leverage pre-trained ML models
- Validate ML predictions with DFT
- Use workflows for systematic studies
- Cite JARVIS papers appropriately

## Community and Support
- GitHub issue tracker
- NIST support
- YouTube tutorials
- Active development
- Regular database updates

## Verification & Sources
**Primary sources**:
1. Official website: https://pages.nist.gov/jarvis/
2. GitHub: https://github.com/usnistgov/jarvis
3. K. Choudhary et al., npj Comput. Mater. 6, 173 (2020)
4. JARVIS-DFT database papers

**Confidence**: VERIFIED - NIST-developed, peer-reviewed publications

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub)
- Developer: NIST (K. Choudhary et al.)
- Academic citations: >500 citations
- Active development: Regular releases, database updates
