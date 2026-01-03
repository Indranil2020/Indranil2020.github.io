# Z2Pack

## Official Resources
- Homepage: https://z2pack.ethz.ch/
- Documentation: https://z2pack.greschd.ch/
- Source Repository: https://github.com/Z2PackDev/Z2Pack
- License: GNU General Public License v3.0

## Overview
Z2Pack is a Python tool for calculating topological invariants and Berry phases using ab-initio or tight-binding methods. Developed by Dominik Gresch, Alexey Soluyanov, and collaborators at ETH Zurich, Z2Pack automates the calculation of Z2 invariants, Chern numbers, and Weyl point chiralities through Wilson loop and Wannier charge center methods. The code emphasizes automation, convergence control, and integration with first-principles calculations.

**Scientific domain**: Topological invariants, Berry phase, Z2 classification  
**Target user community**: Topological materials researchers, ab-initio users, Python workflows

## Theoretical Methods
- Wilson loop calculations
- Wannier charge centers (WCC)
- Z2 topological invariants
- Chern number calculations
- Weyl point chirality
- Berry phase
- Automated convergence
- Surface calculations (BZ)

## Capabilities (CRITICAL)
**Category**: Open-source topological invariant calculator
- Z2 invariant calculation
- Chern number determination
- Weyl chirality
- Wilson loops
- Wannier charge centers
- Automated convergence
- ab-initio integration (VASP, Quantum ESPRESSO, ABINIT, Wannier90)
- Tight-binding models
- Parallel execution
- HDF5 output
- Visualization tools
- Production quality

**Sources**: Official website, documentation, publications

## Key Strengths

### Automation:
- Automatic convergence
- Adaptive refinement
- Error control
- Minimal user intervention
- Production workflows

### Z2 Specialist:
- Rigorous Z2 calculation
- Wilson loop method
- WCC tracking
- Gap detection
- Robust algorithms

### Universal Interface:
- ab-initio codes (VASP, QE, ABINIT)
- Wannier90 models
- Tight-binding (TBmodels)
- Extensible systems
- Python integration

### Python Ecosystem:
- Clean design
- NumPy/SciPy
- Matplotlib plots
- Jupyter workflows
- Modern tools

## Inputs & Outputs
- **Input formats**:
  - Python configuration
  - DFT interface scripts
  - Wannier90 models
  - TBmodels objects
  
- **Output data types**:
  - Z2 invariants
  - Chern numbers
  - WCC evolution
  - Convergence data
  - HDF5 archives
  - Visualization plots

## Interfaces & Ecosystem

### ab-initio Codes:
- VASP (vasp.py interface)
- Quantum ESPRESSO
- ABINIT
- Wannier90
- Extensible framework

### TBmodels:
- Native integration
- Tight-binding calculations
- Fast prototyping
- Model testing

### Visualization:
- Matplotlib output
- WCC plots
- Convergence graphs
- Publication-ready

## Workflow and Usage

### Installation:
```bash
pip install z2pack
```

### Tight-Binding Example:
```python
import z2pack
import tbmodels

# Load TB model
model = tbmodels.Model.from_hr_file('wannier90_hr.dat')

# Create Z2Pack system
system = z2pack.tb.System(model)

# Calculate Z2 for k_z=0 surface
result = z2pack.surface.run(
    system=system,
    surface=lambda s, t: [s/2, t, 0],
    save_file='z2_result.json'
)

# Get Z2 invariant
print("Z2 invariant:", result.z2)
```

### VASP Integration:
```python
import z2pack

# Define VASP system
system = z2pack.fp.System(
    input_files=['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS'],
    kpt_fct=z2pack.fp.kpoint.vasp,
    kpt_path='KPOINTS',
    command='mpirun -n 16 vasp_std'
)

# Run Z2 calculation
result = z2pack.surface.run(
    system=system,
    surface=lambda s, t: [s/2, t, 0.0]
)
```

### Chern Number:
```python
# 2D system - Chern number
result = z2pack.surface.run(
    system=system,
    surface=lambda s, t: [s, t, 0]
)

print("Chern number:", result.chern)
```

### Weyl Points:
```python
# Weyl point chirality
result = z2pack.surface.run(
    system=system,
    surface=lambda s, t: [0.1, s, t],
    pos_tol=1e-3
)

# Check for Weyl point
if result.is_gapped:
    print("Gapped")
else:
    print("Weyl point detected")
    print("Chirality:", result.chern)
```

### Custom Convergence:
```python
# Fine control over convergence
result = z2pack.surface.run(
    system=system,
    surface=lambda s, t: [s/2, t, 0],
    num_lines=11,
    min_neighbour_dist=0.01,
    iterator=range(8, 27, 2),
    pos_tol=1e-3
)
```

## Advanced Features

### Adaptive Refinement:
- Automatic line placement
- Gap-based refinement
- Convergence criteria
- Minimal k-points
- Efficiency optimization

### Parallel Execution:
- Multiple DFT calls
- Parallel k-points
- HPC integration
- Production scaling

### Data Management:
- HDF5 storage
- Result archiving
- Restart capability
- Analysis tools

### Visualization:
- WCC plots
- Gap evolution
- Convergence monitoring
- Publication plots

## Performance Characteristics
- **Speed**: DFT-limited for ab-initio
- **Accuracy**: Controlled convergence
- **Automation**: High
- **Purpose**: Topological invariants
- **Typical**: Hours to days (ab-initio)

## Computational Cost
- DFT calculations dominant
- Adaptive minimizes k-points
- Tight-binding fast
- Convergence-dependent
- Production capable

## Limitations & Known Constraints
- **ab-initio**: Requires DFT code
- **Learning curve**: Moderate
- **System setup**: Code-specific scripts
- **2D focus**: Surface calculations
- **Documentation**: Good but technical

## Comparison with Other Topological Tools
- **vs WannierTools**: Z2Pack automated Z2, WannierTools comprehensive
- **vs Wannier90**: Z2Pack topology calculator, Wannier90 MLWF builder
- **Unique strength**: Automated Z2 invariants, adaptive convergence, Wilson loops, ab-initio integration

## Application Areas

### Topological Insulators:
- 3D TI characterization
- Z2 classification
- (ν0;ν1ν2ν3) calculation
- Material screening
- High-throughput

### Weyl Semimetals:
- Weyl point detection
- Chirality determination
- Fermi arc connection
- Node characterization

### Chern Insulators:
- 2D topological systems
- Chern number calculation
- Quantum Hall physics
- QAHE materials

### Materials Discovery:
- Automated screening
- Database generation
- Topological classification
- Phase identification

## Best Practices

### DFT Setup:
- Sufficient k-points
- Converged wavefunctions
- Appropriate functionals
- SOC when needed

### Convergence:
- Start with coarse
- Monitor WCC evolution
- Check gap closure
- Validate results

### ab-initio:
- Test with tight-binding first
- Validate interface
- Parallel execution
- Error handling

## Community and Support
- Open-source (GPL v3)
- ETH Zurich development
- GitHub repository
- Documentation
- Active maintenance
- User community

## Educational Resources
- Comprehensive documentation
- Tutorials
- Example gallery
- API reference
- Scientific background
- Workflow guides

## Development
- Dominik Gresch (lead)
- Alexey Soluyanov group (ETH Zurich)
- Active development
- Regular updates
- Community contributions
- TBmodels integration

## Research Impact
Z2Pack has enabled systematic topological characterization of materials from first principles, cited in numerous publications on topological insulators and Weyl semimetals.

## Verification & Sources
**Primary sources**:
1. Homepage: https://z2pack.ethz.ch/
2. Documentation: https://z2pack.greschd.ch/
3. GitHub: https://github.com/Z2PackDev/Z2Pack
4. Publications: Phys. Rev. B 95, 075146 (2017)

**Secondary sources**:
1. Topological materials papers
2. User publications
3. Z2 invariant literature

**Confidence**: VERIFIED - Topological invariant tool

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- GitHub: ACCESSIBLE
- License: GPL v3 (open-source)
- **Category**: Topological invariant calculator
- Status: Actively developed
- Institution: ETH Zurich
- Specialized strength: Automated Z2 topological invariant calculation, Wilson loops, Wannier charge centers, adaptive convergence, ab-initio integration, Weyl chirality, Chern numbers, Python-based, TBmodels integration, production quality, minimal user intervention
