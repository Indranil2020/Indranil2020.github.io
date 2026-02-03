# ShengBTE

## Official Resources
- Homepage: http://www.shengbte.org/
- Documentation: http://www.shengbte.org/index.php/usage
- Source Repository: http://www.shengbte.org/downloads
- License: GNU General Public License v3.0

## Overview
ShengBTE is a software package for computing lattice thermal conductivity of crystalline materials from first principles by solving the Boltzmann transport equation (BTE) for phonons. It uses harmonic and anharmonic interatomic force constants from DFT calculations.

**Scientific domain**: Thermal transport, thermoelectrics, phonon physics  
**Target user community**: Researchers studying thermal conductivity and phonon-mediated properties

## Theoretical Methods
- Boltzmann transport equation (BTE) for phonons
- Iterative and relaxation time approximation solutions
- Third-order anharmonic force constants
- Three-phonon scattering processes
- Normal and Umklapp processes
- Isotope scattering
- Boundary scattering
- Finite-size effects

## Capabilities (CRITICAL)
- Lattice thermal conductivity from first principles
- Temperature-dependent thermal conductivity
- Phonon relaxation times
- Phonon mean free paths
- Phonon group velocities
- Mode-resolved thermal conductivity contributions
- Spectral thermal conductivity
- Cumulative thermal conductivity
- Grain boundary scattering effects
- Nanostructure size effects
- Isotope disorder scattering
- Anisotropic thermal conductivity tensors

**Sources**: Official ShengBTE documentation, cited in 6/7 source lists

## Key Strengths
- **Production quality**: Widely used and validated
- **Full BTE solution**: Both RTA and iterative methods
- **Comprehensive**: Mode-resolved and spectral analysis
- **Well-documented**: Established methodology and publications

## Inputs & Outputs
- **Input formats**:
  - CONTROL file (calculation parameters)
  - FORCE_CONSTANTS_2ND (harmonic)
  - FORCE_CONSTANTS_3RD (anharmonic)
  - From thirdorder.py or other phonon codes
  
- **Output data types**:
  - BTE.kappa (thermal conductivity tensor)
  - BTE.KappaTensor (full output)
  - BTE.ReciprocalLattice
  - Mode-resolved properties
  - Cumulative thermal conductivity data

## Interfaces & Ecosystem
- **Force constant generation**:
  - thirdorder.py - automated 3rd order FC calculation
  - Phonopy for 2nd order
  - Interfaces with VASP, Quantum ESPRESSO, etc.
  
- **Post-processing**:
  - Python scripts for analysis
  - Plotting utilities
  - cumulative_kappa.py for analysis

## Workflow and Usage

### Typical ShengBTE Workflow:
```bash
# 1. Generate 2nd order force constants (phonopy)
phonopy -d --dim="2 2 2" -c POSCAR
# Run DFT on displaced structures
phonopy --fc vasprun.xml

# 2. Generate 3rd order force constants (thirdorder.py)
thirdorder.py sow POSCAR
# Run DFT on displaced structures
thirdorder.py reap POSCAR

# 3. Run ShengBTE
ShengBTE
```

## Advanced Features
- **Iterative BTE**: Full solution beyond RTA for accurate transport
- **Spectral analysis**: Frequency-resolved thermal conductivity
- **Cumulative functions**: Mean free path accumulation
- **Size effects**: Grain boundary and nanostructure scattering
- **Isotope scattering**: Natural isotope disorder effects

## Performance Characteristics
- **Computational cost**: Moderate for BTE solution
- **Scalability**: Handles standard q-point grids efficiently
- **Parallelization**: MPI support for large calculations
- **Typical runtime**: Hours to days depending on system and grid

## Computational Cost
- Force constant calculations (DFT): Dominant cost
- ShengBTE BTE solution: Hours to days
- Iterative BTE more expensive than RTA
- Dense q-grids increase cost significantly

## Limitations & Known Constraints
- **Requires force constants**: Needs 2nd and 3rd order from DFT
- **Computational cost**: 3rd order FC calculations expensive
- **Supercell size**: Large supercells needed for accurate 3rd order FC
- **Memory**: Can be intensive for large systems
- **Classical approximation**: Uses classical phonon occupations (fixable)
- **Three-phonon only**: Higher-order processes not included
- **Learning curve**: Requires understanding of phonon BTE
- **Platform**: Linux/Unix

## Verification & Sources
**Primary sources**:
1. Official website: http://www.shengbte.org/
2. Documentation: http://www.shengbte.org/index.php/usage
3. W. Li et al., Comput. Phys. Commun. 185, 1747 (2014) - ShengBTE paper
4. thirdorder.py companion code

**Secondary sources**:
1. ShengBTE tutorials and examples
2. Published thermal conductivity calculations
3. Thermoelectric material studies
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: ACCESSIBLE
- Source code: Available (download from website)
- Community support: Active (email, publications)
- Academic citations: >500
- Standard tool: Widely used for thermal conductivity

## Best Practices
- Converge supercell size for force constants
- Systematic q-point grid convergence
- Test RTA vs iterative BTE
- Validate against experimental data
- Appropriate cutoff distances for 3rd order IFCs

## Community and Support
- Open-source (GPL v3)
- Active user community
- Email support and publications
- Extensive literature using ShengBTE
- Workshop materials available

## Development
- Wu Li and colleagues
- Established and well-maintained
- Regular updates
- Standard tool in thermal transport community

## Research Impact
ShengBTE is a standard tool for first-principles lattice thermal conductivity calculations, widely cited in thermoelectric and thermal transport literature with >500 citations.
