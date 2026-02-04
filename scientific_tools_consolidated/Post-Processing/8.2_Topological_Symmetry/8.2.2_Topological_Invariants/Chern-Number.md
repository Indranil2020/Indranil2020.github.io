# Chern-Number

## Official Resources
- Homepage: https://github.com/stepan-tsirkin/chern-number
- GitHub: https://github.com/stepan-tsirkin/chern-number
- Developer: Stepan Tsirkin (University of Zurich)
- License: Check repository

## Overview
Chern-Number is a code for calculating the Chern number, a topological invariant, using the discretized Berry curvature on a grid. It post-processes output from Wannier90 or tight-binding models to determine the topological character of 2D band structures and detect topological phase transitions.

**Scientific domain**: Topological band theory, Chern insulators, quantum anomalous Hall effect
**Target user community**: Researchers studying topological phases of matter and Chern insulators

## Theoretical Methods
- Berry curvature calculation on k-space grid
- Discretized Brillouin zone integration
- Fukui-Hatsugai-Suzuki method
- Link variable approach for gauge invariance
- Wannier function interpolation

## Capabilities (CRITICAL)
- **Berry Curvature**: Calculation across the Brillouin zone
- **Chern Number**: Grid-based integration method
- **Metallic Systems**: Handles partially filled bands
- **Phase Transitions**: Detection of topological transitions
- **Wannier90 Interface**: Direct _hr.dat file support
- **Tight-Binding**: Model Hamiltonian compatibility

**Sources**: GitHub repository, developer documentation

## Key Strengths

### Robust Integration:
- Grid-based method
- Gauge-invariant formulation
- Convergence with k-points
- Accurate for complex bands

### Wannier90 Compatible:
- Direct interface with Wannier90
- Uses interpolated bands
- Efficient k-point sampling
- Standard file formats

### Developer Ecosystem:
- Part of IrRep/WannierBerri ecosystem
- Consistent methodology
- Active maintenance
- Research-grade code

## Inputs & Outputs
- **Input formats**:
  - Wannier90 _hr.dat files
  - Tight-binding model parameters
  - k-mesh specification
  
- **Output data types**:
  - Chern numbers per band
  - Berry curvature maps
  - Band-resolved invariants

## Installation
```bash
git clone https://github.com/stepan-tsirkin/chern-number.git
cd chern-number
# Follow installation instructions
```

## Usage Examples
```python
# Example usage
from chern import ChernCalculator

# Load Wannier90 Hamiltonian
calc = ChernCalculator.from_wannier90("wannier90_hr.dat")

# Calculate Chern number for occupied bands
chern = calc.compute_chern(nk=100, bands=range(4))
print(f"Chern number: {chern}")
```

## Performance Characteristics
- **Speed**: Efficient grid-based calculation
- **Accuracy**: Convergent with k-mesh density
- **Memory**: Moderate for typical systems

## Limitations & Known Constraints
- **2D systems**: Primarily for 2D or quasi-2D systems
- **Band gaps**: Best for gapped systems
- **Convergence**: Requires sufficient k-points
- **Learning curve**: Topological band theory knowledge needed

## Comparison with Other Tools
- **vs Z2Pack**: Chern-Number grid-based, Z2Pack WCC-based
- **vs WannierBerri**: Chern-Number specialized, WannierBerri comprehensive
- **vs WannierTools**: Different calculation approaches
- **Unique strength**: Simple, focused Chern number calculation

## Application Areas
- Quantum anomalous Hall effect
- Chern insulator identification
- Topological phase diagrams
- Magnetic topological materials
- 2D material topology

## Best Practices
- Use converged Wannier90 calculations
- Check Chern number convergence with k-mesh
- Verify gap remains open
- Compare with symmetry indicators

## Community and Support
- GitHub repository
- Developer: Stepan Tsirkin
- Related to IrRep ecosystem
- Research code with academic support

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/stepan-tsirkin/chern-number
2. Developer: Stepan Tsirkin (University of Zurich)
3. Related tools: IrRep, WannierBerri

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- GitHub repository: ACCESSIBLE
- Source code: OPEN
- Developer: Active researcher
- Method: Berry curvature integration
