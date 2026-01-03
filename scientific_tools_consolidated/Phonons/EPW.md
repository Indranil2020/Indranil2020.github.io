# EPW (Electron-Phonon Wannier)

## Official Resources
- Homepage: https://epw-code.org/
- Documentation: https://docs.epw-code.org/
- Source Repository: https://github.com/EPW-code/EPW (also part of QE distribution)
- License: GNU General Public License v2.0

## Overview
EPW is a highly specialized code for calculating electron-phonon coupling and related properties from first principles using maximally localized Wannier functions. Part of the Quantum ESPRESSO distribution, EPW uses Wannier interpolation to achieve extremely efficient calculations of electron-phonon matrix elements on ultra-dense k and q-point grids, enabling accurate predictions of superconductivity, electrical transport, optical properties, and carrier mobilities.

**Scientific domain**: Electron-phonon coupling, superconductivity, transport properties  
**Target user community**: Superconductivity researchers, transport properties, materials science

## Theoretical Methods
- Electron-phonon coupling via Wannier interpolation
- Wannier function representation (from Wannier90)
- Eliashberg theory of superconductivity
- Migdal-Eliashberg equations
- Boltzmann transport for electrons and phonons
- Polar correction (Fröhlich interaction)
- Phonon-assisted optical absorption
- Temperature-dependent band structures
- Carrier mobility calculations
- Ultra-dense k/q-mesh interpolation

## Capabilities (CRITICAL)
**Category**: Open-source electron-phonon code
- Electron-phonon matrix elements calculation
- Superconducting critical temperature (Tc)
- Eliashberg spectral function α²F(ω)
- Electron-phonon coupling constant λ
- Phonon linewidths from electron-phonon coupling
- Temperature-dependent electron self-energies
- Carrier mobility (electrons and holes)
- Electrical conductivity with electron-phonon scattering
- Seebeck coefficient
- Phonon-assisted optical absorption
- Polar materials (Fröhlich interaction)
- Anisotropic superconducting gaps
- Real-axis Eliashberg equations
- Ultra-fast Wannier interpolation
- Integration with Quantum ESPRESSO and Wannier90
- Production quality

**Sources**: Official EPW documentation, publications

## Key Strengths

### Wannier Interpolation:
- Ultra-dense k/q grids (millions of points)
- Orders of magnitude faster than direct DFT
- Smooth interpolation
- Production efficiency
- Accurate convergence

### Superconductivity:
- Full Eliashberg theory
- Anisotropic gaps
- Critical temperature prediction
- Spectral functions
- Research and prediction

### Transport Properties:
- Carrier mobility from first principles
- Electron-phonon scattering
- Temperature dependence
- Anisotropic tensors
- Thermoelectric properties

### QE Integration:
- Seamless workflow
- Standard DFT input
- Wannier90 connection
- Production pipeline

## Inputs & Outputs
- **Input formats**:
  - Quantum ESPRESSO DFT output
  - Wannier90 Wannier functions
  - EPW input file (epw.in)
  - Phonon calculations from QE
  
- **Output data types**:
  - Electron-phonon matrix elements
  - Eliashberg spectral function
  - Superconducting Tc
  - Carrier mobility
  - Transport coefficients
  - Self-energies
  - Spectral functions

## Interfaces & Ecosystem

### Quantum ESPRESSO:
- Native integration
- Part of QE distribution
- Standard workflow
- pw.x and ph.x input

### Wannier90:
- Essential for Wannier functions
- Wannier interpolation
- MLWFs for electrons
- Standard pipeline

### Workflow:
- QE SCF/NSCF → QE phonons → Wannier90 → EPW

## Workflow and Usage

### Typical EPW Workflow:

```bash
# 1. DFT calculation (Quantum ESPRESSO)
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out

# 2. Phonon calculation on coarse grid
ph.x < ph.in > ph.out

# 3. Wannierization (Wannier90)
wannier90.x -pp silicon
pw2wannier90.x < pw2wan.in
wannier90.x silicon

# 4. EPW calculation
epw.x < epw.in > epw.out
```

### EPW Input Example (epw.in):
```fortran
&inputepw
  prefix = 'silicon'
  amass(1) = 28.0855
  outdir = './'
  
  elph = .true.
  kmaps = .false.
  epbwrite = .true.
  epbread = .false.
  
  etf_mem = 1
  
  nbndsub = 8
  
  wannierize = .false.
  num_iter = 500
  iprint = 2
  
  ephwrite = .true.
  
  fsthick = 10.0
  eptemp = 300
  degaussw = 0.05
  
  dvscf_dir = './save'
  
  nkf1 = 40
  nkf2 = 40
  nkf3 = 40
  
  nqf1 = 40
  nqf2 = 40
  nqf3 = 40
  
  nk1 = 8
  nk2 = 8
  nk3 = 8
  
  nq1 = 4
  nq2 = 4
  nq3 = 4
/
```

## Advanced Features

### Superconductivity:
- Migdal-Eliashberg equations
- Anisotropic Eliashberg (on Fermi surface)
- Real-axis solutions
- Critical temperature
- Gap functions

### Transport:
- Carrier mobility (iterative BTE)
- Relaxation time approximation
- Electron-phonon scattering rates
- Temperature-dependent
- Anisotropic tensors

### Polar Materials:
- Long-range Fröhlich interaction
- Polar optical phonons
- Quadrupole corrections
- LO-TO splitting

### Optical Properties:
- Phonon-assisted absorption
- Indirect transitions
- Temperature-dependent spectra

## Performance Characteristics
- **Speed**: Very fast via Wannier interpolation
- **Accuracy**: High-quality convergence
- **k/q-mesh**: Millions of points feasible
- **Purpose**: Electron-phonon specialist
- **Typical**: Hours to days (depending on convergence)

## Computational Cost
- DFT/phonon calculations most expensive
- EPW interpolation very efficient
- Dense k/q-grid feasible
- Convergence testing important
- Production capable

## Limitations & Known Constraints
- **Requires QE+Wannier90**: Part of workflow
- **Wannier functions**: Quality critical
- **Phonon calculations**: Expensive for large systems
- **Learning curve**: Steep, complex workflow
- **Convergence**: Many parameters to test
- **Polar materials**: Requires special treatment
- **Documentation**: Comprehensive but technical
- **Imaginary modes**: Can cause issues

## Comparison with Other Codes
- **Unique capability**: Wannier interpolation for electron-phonon
- **vs direct DFT**: EPW orders of magnitude faster
- **QE integration**: Standard workflow
- **Gold standard**: For superconductivity/transport from first principles

## Application Areas

### Superconductivity:
- Tc prediction
- Conventional superconductors
- Anisotropic gaps
- Eliashberg theory
- Material discovery

### Transport:
- Carrier mobility
- Semiconductors
- Thermoelectrics
- Electron-phonon scattering
- Temperature dependence

### Optical Properties:
- Indirect absorption
- Phonon-assisted processes
- Temperature effects
- Spectroscopy theory

## Best Practices

### Workflow:
- Quality DFT convergence
- Good Wannier functions
- Phonon convergence
- Dense k/q-grid testing
- Systematic convergence

### Wannier Functions:
- Appropriate projections
- Check spreads
- Band structure validation
- Energy window selection

### Convergence:
- k-point and q-point grids
- Broadening parameters
- Temperature convergence
- Fermi surface sampling

## Community and Support
- Open-source (GPL v2)
- Part of Quantum ESPRESSO
- Active development
- Mailing list
- Workshops and schools
- Comprehensive documentation
- Tutorial materials

## Educational Resources
- EPW documentation
- QE schools tutorials
- Hands-on workshops
- Example calculations
- Publication list
- Theory background

## Development
- Samuel Poncé (lead, Oxford)
- Quantum ESPRESSO developers
- International collaboration
- Active development
- Regular updates
- Feature additions

## Research Impact
EPW is the standard tool for first-principles electron-phonon calculations, enabling accurate predictions of superconducting properties and carrier mobilities, cited in hundreds of publications.

## Verification & Sources
**Primary sources**:
1. Homepage: https://epw-code.org/
2. Documentation: https://docs.epw-code.org/
3. GitHub: https://github.com/EPW-code/EPW
4. Publications: Comp. Phys. Comm. 209, 116 (2016); Rev. Mod. Phys. 89, 015003 (2017)

**Secondary sources**:
1. Quantum ESPRESSO distribution
2. Superconductivity literature
3. Transport property papers
4. User publications (hundreds)

**Confidence**: VERIFIED - Standard electron-phonon code

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub, part of QE)
- License: GPL v2 (open-source)
- **Category**: Open-source electron-phonon code
- Status: Actively developed
- Community: Large, international
- Specialized strength: Wannier interpolation for electron-phonon coupling, superconducting Tc prediction, carrier mobility calculations, Eliashberg theory, ultra-dense k/q-mesh capability, Quantum ESPRESSO integration, production quality, gold standard for superconductivity and transport from first principles
