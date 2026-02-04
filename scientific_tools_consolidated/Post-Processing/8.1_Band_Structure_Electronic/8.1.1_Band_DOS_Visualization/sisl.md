# sisl

## Official Resources
- **Homepage**: https://zerothi.github.io/sisl/
- **GitHub**: https://github.com/zerothi/sisl
- **Documentation**: https://zerothi.github.io/sisl/
- **PyPI**: https://pypi.org/project/sisl/
- **Zenodo**: https://doi.org/10.5281/zenodo.597181
- **License**: Mozilla Public License 2.0

## Overview
sisl (Scientific Python toolbox for large-scale tight-binding and DFT/NEGF calculations) is a comprehensive Python library for electronic structure analysis. Originally developed to handle SIESTA/TranSIESTA output, it has grown into a versatile tool supporting multiple DFT codes. sisl excels at creating, manipulating, and analyzing tight-binding Hamiltonians and performing large-scale transport calculations.

**Scientific domain**: Electronic structure, tight-binding, quantum transport, NEGF
**Target user community**: Researchers working with SIESTA, transport calculations, tight-binding models, and large-scale electronic structure

## Theoretical Background
sisl implements:
- Tight-binding Hamiltonians: H = Σ t_ij |i⟩⟨j|
- Overlap matrices for non-orthogonal bases
- Green's function methods: G(E) = (E·S - H)⁻¹
- NEGF for quantum transport
- Density matrix and charge analysis
- Band structure from Hamiltonian diagonalization

## Capabilities (CRITICAL)
- **Tight-Binding**: Create orthogonal/non-orthogonal TB Hamiltonians
- **Electronic Structure**: Band structure, DOS, PDOS analysis
- **NEGF Transport**: Non-equilibrium Green's function calculations
- **File I/O**: Read/write multiple DFT code formats
- **Geometry**: Structure manipulation, supercells, defects
- **Sparse Matrices**: Efficient handling of large systems
- **Wannier90**: Interface for Wannier function analysis
- **Visualization**: Geometry and orbital plotting

## Key Strengths

### Multi-Code Support:
- SIESTA / TranSIESTA (native, comprehensive)
- VASP (WAVECAR, vasprun.xml)
- GULP
- BigDFT
- Wannier90
- OpenMX
- ScaleUp

### Tight-Binding Framework:
- Build custom TB models
- Non-orthogonal basis support
- Sparse matrix efficiency
- Large-scale systems (millions of orbitals)

### Transport Calculations:
- NEGF implementation
- Transmission functions
- Current calculations
- Multi-terminal devices

### Geometry Tools:
- Supercell generation
- Defect creation
- Structure manipulation
- Coordinate transformations

## Inputs & Outputs
- **Input formats**:
  - SIESTA: .fdf, .XV, .DM, .TSHS, .TBT.nc
  - VASP: POSCAR, WAVECAR, vasprun.xml
  - Wannier90: .win, _hr.dat
  - Generic: XYZ, CIF, PDB
  
- **Output data types**:
  - Band structures
  - DOS/PDOS
  - Transmission spectra
  - Charge densities
  - Hamiltonian matrices

## Interfaces & Ecosystem
- **Python integration**:
  - NumPy/SciPy for numerical operations
  - Matplotlib for visualization
  - xarray for multi-dimensional data
  - Sparse matrices (scipy.sparse)
  
- **Framework compatibility**:
  - ASE (Atomic Simulation Environment)
  - pymatgen structures
  - Jupyter notebooks

## Installation
```bash
pip install sisl
conda install -c conda-forge sisl
```

With optional dependencies:
```bash
pip install sisl[analysis,viz]
```

## Usage Examples
```python
import sisl

# Read SIESTA output
geometry = sisl.get_sile("RUN.fdf").read_geometry()
H = sisl.get_sile("RUN.TSHS").read_hamiltonian()

# Calculate band structure
band = sisl.BandStructure(H, [[0,0,0], [0.5,0,0], [0.5,0.5,0]], 100)
eigs = band.eigh()

# Calculate DOS
E = np.linspace(-5, 5, 1000)
dos = H.DOS(E)

# Build custom tight-binding model
H_tb = sisl.Hamiltonian(geometry)
H_tb[0, 1] = -2.7  # Hopping parameter
```

## Performance Characteristics
- **Speed**: Optimized sparse matrix operations
- **Memory**: Efficient for large systems
- **Scalability**: Handles millions of orbitals
- **Parallelization**: NumPy/SciPy parallelism

## Limitations & Known Constraints
- **Learning curve**: Rich API requires time to master
- **Documentation**: Extensive but can be complex
- **SIESTA-centric**: Best support for SIESTA family
- **Visualization**: Basic plotting, external tools recommended

## Comparison with Other Tools
- **vs pymatgen**: sisl specialized for TB/transport, pymatgen broader
- **vs ASE**: sisl has TB/NEGF focus, ASE is general calculator
- **vs PythTB**: sisl more comprehensive, PythTB simpler TB
- **Unique strength**: SIESTA/TranSIESTA integration, NEGF, large-scale TB

## Application Areas
- Quantum transport in nanostructures
- Molecular electronics
- 2D material devices
- Defect calculations
- Tight-binding model development
- Large-scale electronic structure
- Wannier function analysis

## Best Practices
- Use sparse matrices for large systems
- Leverage sisl's geometry tools for structure manipulation
- Use xarray for multi-dimensional data analysis
- Check convergence with k-point sampling
- Validate TB models against DFT

## Community and Support
- GitHub issue tracker
- Active development
- Comprehensive tutorials
- Zenodo DOI for citation

## Verification & Sources
**Primary sources**:
1. Official documentation: https://zerothi.github.io/sisl/
2. GitHub repository: https://github.com/zerothi/sisl
3. Zenodo DOI: 10.5281/zenodo.597181
4. N. Papior et al., Comput. Phys. Commun. 212, 8 (2017)

**Confidence**: VERIFIED - Well-established tool with academic publications

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, MPL 2.0)
- Developer: Nick Papior (zerothi)
- Academic citations: Published in CPC
- Active development: Regular releases, 1000+ commits
