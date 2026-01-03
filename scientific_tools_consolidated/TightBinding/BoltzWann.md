# BoltzWann

## Official Resources
- Homepage: Part of Wannier90 ecosystem
- Documentation: Wannier90 documentation
- Source Repository: https://github.com/wannier-developers/wannier90 (integrated)
- License: GNU General Public License v2.0

## Overview
BoltzWann is a module within the Wannier90 package for calculating Boltzmann transport properties using Wannier interpolation. The tool computes electrical conductivity, Seebeck coefficient, electronic thermal conductivity, and other transport coefficients from first principles via Wannier tight-binding models. BoltzWann implements the Boltzmann transport equation in the relaxation time approximation.

**Scientific domain**: Boltzmann transport, thermoelectrics, electronic transport  
**Target user community**: Wannier90 users, transport properties, thermoelectrics

## Theoretical Methods
- Boltzmann transport equation
- Relaxation time approximation (RTA)
- Wannier interpolation
- Band structure integration
- Transport distribution function
- Temperature-dependent transport

## Capabilities (CRITICAL)
**Category**: Wannier90 transport module
- Electrical conductivity
- Seebeck coefficient
- Electronic thermal conductivity
- Power factor
- Figure of merit (ZT estimation)
- Temperature-dependent properties
- Wannier interpolation
- Fine k-mesh integration
- Chemical potential variation
- Integrated with Wannier90

**Sources**: Wannier90 documentation

## Key Strengths

### Wannier90 Integration:
- Native module
- Direct hr.dat usage
- Standard workflow
- Seamless integration
- DFT to transport

### Transport Properties:
- Comprehensive coefficients
- Thermoelectric properties
- Temperature dependence
- Chemical potential scans
- Production quality

### Interpolation:
- Fine k-mesh from coarse DFT
- Accurate integration
- Efficient computation
- Wannier quality

## Inputs & Outputs
- **Input formats**:
  - Wannier90 hr.dat
  - BoltzWann input section in .win
  - Temperature range
  - Chemical potential range
  
- **Output data types**:
  - Conductivity tensors
  - Seebeck coefficient
  - Thermal conductivity
  - Power factor
  - Transport distribution
  - Tabulated results

## Workflow and Usage

### Wannier90 Input (.win):
```
# After standard Wannier90 parameters
boltzwann = true

# BoltzWann specific
boltz_kmesh = 50 50 50
boltz_relax_time = 1e-14
boltz_2d_dir = 3

# Temperature range (K)
boltz_temp = 100 200 300 400 500

# Chemical potential range (eV)
boltz_mu_min = -0.5
boltz_mu_max = 0.5
boltz_mu_step = 0.01
```

### Run Wannier90 with BoltzWann:
```bash
# After standard Wannier90 workflow
wannier90.x silicon

# BoltzWann output in silicon_boltzdist.dat
```

### Analyze Results:
```python
import numpy as np
import matplotlib.pyplot as plt

# Load BoltzWann output
data = np.loadtxt('silicon_boltzdist.dat')
mu = data[:, 0]  # Chemical potential
sigma = data[:, 1]  # Conductivity
seebeck = data[:, 2]  # Seebeck coefficient

# Plot Seebeck coefficient
plt.plot(mu, seebeck)
plt.xlabel('Chemical potential (eV)')
plt.ylabel('Seebeck coefficient (μV/K)')
plt.show()
```

## Advanced Features

### Transport Coefficients:
- σ (conductivity tensor)
- S (Seebeck coefficient)
- κₑ (electronic thermal conductivity)
- Power factor (S²σ)
- Lorenz number

### Integration:
- Fine k-mesh via interpolation
- Fermi surface integration
- Temperature dependence
- Chemical potential variation
- Convergence control

### 2D Systems:
- Layered materials
- Reduced dimensions
- Direction specification
- Quasi-2D transport

## Performance Characteristics
- **Speed**: Post-Wannier90, fast
- **Accuracy**: Wannier quality
- **Purpose**: Transport properties
- **Typical**: Minutes

## Computational Cost
- Post-processing (after Wannier90)
- k-mesh dependent
- Temperature/μ scan cost
- Efficient interpolation
- Production capable

## Limitations & Known Constraints
- **RTA**: Relaxation time approximation
- **τ parameter**: Requires relaxation time input
- **No phonons**: Electronic only
- **Requires Wannier90**: Part of ecosystem
- **Constant τ**: Simplified scattering

## Comparison with Other Transport Tools
- **vs WannierBerri**: BoltzWann Boltzmann, WannierBerri Berry
- **vs BoltzTraP**: BoltzWann Wannier-based, BoltzTraP DFT-direct
- **Unique strength**: Wannier interpolation, Wannier90 integration, thermoelectrics

## Application Areas

### Thermoelectrics:
- Seebeck coefficient
- Power factor
- ZT estimation
- Material screening
- Doping optimization

### Electronic Transport:
- Conductivity tensors
- Carrier mobility
- Temperature dependence
- Chemical potential effects

### Materials Screening:
- High-throughput
- Transport databases
- Property prediction
- Optimization studies

## Best Practices

### Wannier Functions:
- Quality MLWFs from Wannier90
- Converged tight-binding
- Appropriate bands
- Validated model

### k-Mesh:
- Dense BoltzWann mesh
- Convergence testing
- Fermi surface resolution
- Integration accuracy

### Relaxation Time:
- Appropriate τ values
- Literature/experiment
- Material-dependent
- Constant τ approximation

## Community and Support
- Part of Wannier90
- Wannier90 community
- Documentation
- Mailing list
- Established tool

## Educational Resources
- Wannier90 user guide
- BoltzWann tutorial
- Example calculations
- Transport theory background

## Development
- Wannier90 developers
- International collaboration
- Integrated module
- Regular updates with Wannier90

## Research Impact
BoltzWann enables first-principles thermoelectric property calculations, widely used for screening thermoelectric materials and understanding electronic transport from ab-initio.

## Verification & Sources
**Primary sources**:
1. Wannier90: https://wannier.org/
2. Wannier90 documentation
3. GitHub: https://github.com/wannier-developers/wannier90

**Secondary sources**:
1. Boltzmann transport literature
2. Thermoelectric papers

**Confidence**: VERIFIED - Wannier90 module

**Verification status**: ✅ VERIFIED
- **Category**: Wannier90 transport module
- Status: Part of Wannier90
- **Note**: BoltzWann is an integrated module within Wannier90 for Boltzmann transport calculations. Computes electrical conductivity, Seebeck coefficient, and electronic thermal conductivity using Wannier interpolation. Relaxation time approximation. Essential for first-principles thermoelectric property prediction. Part of standard Wannier90 distribution.
