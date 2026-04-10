# ppafm

## Official Resources
- Source Repository: https://github.com/Probe-Particle/ppafm
- Documentation: https://probe-particle.github.io/ppafm/
- PyPI: https://pypi.org/project/ppafm/
- License: Open source (Apache-2.0)

## Overview
**ppafm** (Probe-Particle AFM) is a simple and efficient simulation software for high-resolution atomic force microscopy (HR-AFM) and other scanning probe microscopy (SPM) techniques with sub-molecular resolution. It simulates the deflection of a probe particle (typically CO or Xe) attached to the tip, enabling realistic AFM, STM, IETS, and TERS simulations.

**Scientific domain**: Atomic force microscopy, scanning probe microscopy, surface science  
**Target user community**: Researchers simulating and interpreting high-resolution AFM and SPM experiments

## Theoretical Methods
- Probe-particle model (CO/Xe tip functionalization)
- Classical force field for tip-sample interaction
- Lennard-Jones potentials
- Point-charge electrostatics
- Hartree potential from DFT
- Tersoff-Hamann STM approximation
- Inelastic tunneling (IETS)
- Tip-enhanced Raman spectroscopy (TERS)

## Capabilities (CRITICAL)
- High-resolution AFM image simulation
- STM image simulation
- IETS (inelastic electron tunneling spectroscopy)
- TERS (tip-enhanced Raman spectroscopy)
- Kelvin probe force microscopy (KPFM)
- Sub-molecular resolution imaging
- CO/Xe tip functionalization
- Multiple DFT code interfaces
- 3D force map calculation
- Frequency shift calculation

**Sources**: GitHub repository, Comput. Phys. Commun. 305, 109341 (2024)

## Key Strengths

### Probe-Particle Model:
- Realistic tip functionalization (CO, Xe, etc.)
- Sub-molecular resolution
- Efficient classical simulation
- Captures key experimental features
- Well-validated against experiment

### Multi-Mode SPM:
- AFM, STM, IETS, TERS, KPFM
- Comprehensive SPM simulation
- Consistent model across modes
- Direct comparison with experiment

### DFT Integration:
- VASP, QE, FHI-aims, CP2K, GPAW
- Reads DFT-calculated densities
- Hartree potential for electrostatics
- Orbital data for STM

### Efficient:
- Fast classical simulation
- GPU acceleration available
- Python API
- PyPI installation

## Inputs & Outputs
- **Input formats**:
  - DFT charge density (VASP LOCPOT, QE rho)
  - DFT Hartree potential
  - DFT orbitals for STM
  - Force field parameters
  - Probe-particle configuration
  
- **Output data types**:
  - AFM images (frequency shift maps)
  - STM images
  - IETS spectra and maps
  - TERS spectra
  - KPFM images
  - 3D force maps

## Interfaces & Ecosystem
- **PPSTM**: STM/STS companion code
- **VASP**: DFT input
- **Quantum ESPRESSO**: DFT input
- **FHI-aims**: DFT input
- **CP2K**: DFT input
- **GPAW**: DFT input
- **Python**: Scripting and visualization

## Performance Characteristics
- **Speed**: Very fast (seconds per image)
- **Accuracy**: Good qualitative agreement with experiment
- **System size**: Hundreds of atoms
- **Memory**: Low
- **GPU**: Available for acceleration

## Computational Cost
- **AFM image**: Seconds to minutes
- **3D force map**: Minutes
- **STM image**: Seconds
- **Typical**: Very efficient

## Limitations & Known Constraints
- **Classical model**: Not fully quantum mechanical
- **Force field**: Simplified tip-sample interaction
- **No full NEGF**: Approximate tunneling
- **Parameter dependence**: Results depend on probe-particle stiffness
- **No chemical bond formation**: Cannot simulate bond-breaking

## Comparison with Other Codes
- **vs PPSTM**: ppafm focuses on AFM, PPSTM on STM/STS
- **vs cp2k-spm-tools**: ppafm is classical model, cp2k-spm is DFT-based
- **vs full DFT AFM**: ppafm is much faster, less accurate
- **Unique strength**: Efficient HR-AFM simulation with probe-particle model, multi-mode SPM, multi-DFT-code interface

## Application Areas

### On-Surface Molecular Imaging:
- Organic molecules on surfaces
- Bond-resolved AFM
- Molecular structure determination
- Intermolecular bonds

### 2D Materials:
- Graphene, hBN, TMDs
- Moiré patterns
- Defect characterization
- Edge structure

### Tip Functionalization:
- CO-tip AFM
- Xe-tip AFM
- Cl-tip imaging
- Tip-induced contrast

### Surface Science:
- Adsorption geometry
- Surface reconstruction
- Charge distribution (KPFM)
- Vibrational mapping (IETS/TERS)

## Best Practices

### Probe-Particle Parameters:
- Calibrate stiffness (k ~ 0.5 N/m for CO)
- Test tip-sample distance
- Compare with experimental contrast
- Consider lateral force effects

### DFT Input:
- Use well-converged Hartree potential
- Adequate vacuum for surface
- Include enough atoms for force field
- Check electrostatic accuracy

### Image Interpretation:
- Consider both AFM and STM
- Compare with experimental resolution
- Account for thermal drift
- Validate with known systems

## Community and Support
- Open source (Apache-2.0)
- PyPI installation available
- Active development (Probe-Particle team)
- Published in Comput. Phys. Commun. (2024)
- Used by major SPM groups worldwide
- Tutorial examples provided

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/Probe-Particle/ppafm
2. N. Oinonen et al., Comput. Phys. Commun. 305, 109341 (2024)
3. P. Hapala et al., Phys. Rev. B 90, 085421 (2014)
4. PyPI: https://pypi.org/project/ppafm/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Documentation: ACCESSIBLE
- PyPI: AVAILABLE
- Community support: Active (SPM community)
- Academic citations: >500 (method papers)
- Active development: Ongoing
- Specialized strength: Efficient HR-AFM/SPM simulation with probe-particle model, multi-mode SPM, multi-DFT-code interface
