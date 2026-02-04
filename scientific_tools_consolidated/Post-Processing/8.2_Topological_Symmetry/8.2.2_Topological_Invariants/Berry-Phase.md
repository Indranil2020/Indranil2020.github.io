# berry

## Official Resources
- Homepage: https://ricardoribeiro-2020.github.io/berry/
- GitHub: https://github.com/ricardoribeiro-2020/berry
- Documentation: https://ricardoribeiro-2020.github.io/berry/
- Publication: R. Ribeiro-Palau et al., Comput. Phys. Commun. 267, 108064 (2021)
- License: MIT License

## Overview
berry is a Python code for the differentiation of Bloch wavefunctions from DFT calculations. It calculates Berry connections, Berry curvature, and topological invariants (Chern numbers, Z2) by unwinding the phase of the Bloch states. It also computes nonlinear optical properties like second harmonic generation (SHG) and anomalous Hall conductivity.

**Scientific domain**: Berry phase physics, topological invariants, nonlinear optics
**Target user community**: Researchers studying topological properties and nonlinear optical responses

## Theoretical Methods
- Berry connection calculation
- Berry curvature from Bloch states
- Phase unwinding algorithms
- Modern theory of polarization
- Wannier charge center evolution
- Nonlinear optical response theory

## Capabilities (CRITICAL)
- **Berry Connection**: Calculation from DFT wavefunctions
- **Berry Curvature**: Full Brillouin zone mapping
- **Chern Numbers**: Topological invariant computation
- **Z2 Invariant**: Time-reversal protected topology
- **SHG Conductivity**: Second harmonic generation
- **Anomalous Hall**: Intrinsic Hall conductivity
- **QE Interface**: Quantum ESPRESSO support
- **Wannier90 Interface**: Wannier interpolation support

**Sources**: berry documentation, CPC publication

## Key Strengths

### Comprehensive Berry Physics:
- Full Berry connection/curvature
- Multiple invariant types
- Optical response properties
- Unified framework

### DFT Integration:
- Quantum ESPRESSO interface
- Wannier90 compatibility
- Standard file formats
- Automated workflows

### Nonlinear Optics:
- SHG conductivity
- Shift current
- Injection current
- Optical applications

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO output
  - Wannier90 files
  - Band structure data
  
- **Output data types**:
  - Berry curvature maps
  - Topological invariants
  - Optical conductivities
  - Hall conductivity

## Installation
```bash
pip install berry
# Or from source
git clone https://github.com/ricardoribeiro-2020/berry.git
cd berry
pip install -e .
```

## Usage Examples
```python
import berry

# Load DFT calculation
calc = berry.Calculation.from_qe(prefix='material')

# Calculate Berry curvature
bc = calc.berry_curvature(kpoints, bands)

# Calculate Chern number
chern = calc.chern_number(occupied_bands)

# Calculate SHG conductivity
shg = calc.shg_conductivity(omega)
```

## Performance Characteristics
- **Speed**: Efficient wavefunction processing
- **Accuracy**: Phase unwinding precision
- **Scalability**: Handles dense k-meshes

## Limitations & Known Constraints
- **DFT codes**: Currently QE and Wannier90
- **Phase unwinding**: Requires careful handling
- **Memory**: Large for fine k-meshes
- **Documentation**: Evolving

## Comparison with Other Tools
- **vs BerryPI**: berry more comprehensive, BerryPI VASP-focused
- **vs WannierBerri**: Different approaches to Berry physics
- **vs Z2Pack**: berry includes optical properties
- **Unique strength**: Combined topological and optical properties

## Application Areas
- Topological material characterization
- Nonlinear optical materials
- Anomalous Hall effect studies
- Ferroelectric polarization
- Photovoltaic materials (shift current)

## Best Practices
- Ensure converged DFT calculation
- Use dense k-mesh for curvature
- Check phase continuity
- Validate against known materials

## Community and Support
- GitHub repository
- CPC publication
- Active development
- Academic support

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ricardoribeiro-2020/berry
2. R. Ribeiro-Palau et al., Comput. Phys. Commun. 267, 108064 (2021)
3. Documentation: https://ricardoribeiro-2020.github.io/berry/

**Note**: Berry phase calculations are also implemented natively in major DFT codes:
- VASP: `LBERRY = .TRUE.`
- Quantum ESPRESSO: `bp` calculation
- ABINIT: Berry phase polarization
- Wannier90: Berry curvature via interpolation

**Confidence**: VERIFIED - Published in CPC

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: AVAILABLE
- Source code: OPEN (GitHub, MIT)
- Academic citations: CPC publication
- Active development: Maintained
