# easyunfold

## Official Resources
- **Homepage**: https://smtg-bham.github.io/easyunfold/
- **GitHub**: https://github.com/SMTG-Bham/easyunfold
- **Documentation**: https://smtg-bham.github.io/easyunfold/
- **PyPI**: https://pypi.org/project/easyunfold/
- **JOSS Paper**: https://joss.theoj.org/papers/10.21105/joss.05974
- **License**: MIT License

## Overview
easyunfold is a Python package for band structure unfolding, making it easy to obtain effective band structures (EBS) from supercell calculations. It properly accounts for symmetry-breaking by sampling necessary additional k-points and averaging spectral weights appropriately. The package supports multiple DFT codes and provides publication-ready visualizations.

**Scientific domain**: Band structure unfolding, supercell calculations, defect physics, alloy band structures
**Target user community**: Researchers studying defects, alloys, interfaces, and disordered systems using supercell methods

## Theoretical Background
Band unfolding recovers the primitive cell band structure from supercell calculations:
- Spectral weight: W(k,E) = |⟨ψ_K|ψ_k⟩|²
- Maps supercell K-points back to primitive cell k-points
- Accounts for zone folding in supercells
- Handles symmetry-breaking from defects/disorder
- Based on Popescu & Zunger methodology (Phys. Rev. B 85, 085201)

## Capabilities (CRITICAL)
- **K-point Generation**: Generate k-points for supercell BZ sampling
- **Band Unfolding**: Extract effective band structure from supercells
- **Symmetry Handling**: Proper treatment of symmetry-breaking
- **Spectral Weights**: Calculate and visualize unfolding weights
- **Multi-code Support**: VASP, CASTEP, Quantum ESPRESSO
- **Visualization**: Publication-ready unfolded band plots
- **Automation**: Streamlined workflow from generation to plotting

## Key Strengths

### Symmetry-Aware Unfolding:
- Properly samples additional k-points for broken symmetry
- Averages spectral weights correctly
- Handles partial symmetry breaking
- Robust for defect calculations

### Multi-Code Support:
- VASP (primary support)
- CASTEP
- Quantum ESPRESSO
- Extensible architecture

### User-Friendly Workflow:
- Simple two-step process (generate → unfold)
- Command-line interface
- Python API for scripting
- Minimal user input required

### Publication Quality:
- Matplotlib-based visualization
- Customizable color maps
- Spectral function plotting
- Export to various formats

## Inputs & Outputs
- **Input formats**:
  - Primitive cell structure (POSCAR, CIF)
  - Supercell structure
  - DFT wavefunction files (WAVECAR for VASP)
  - K-point paths
  
- **Output data types**:
  - Unfolded band structure data
  - Spectral weights
  - Matplotlib figures
  - JSON data files

## Interfaces & Ecosystem
- **Python integration**:
  - NumPy for numerical operations
  - Matplotlib for visualization
  - pymatgen for structure handling
  - spglib for symmetry analysis
  
- **DFT code interfaces**:
  - VASP (via PyVaspwfc)
  - CASTEP
  - Quantum ESPRESSO

## Installation
```bash
pip install easyunfold
```

With visualization extras:
```bash
pip install easyunfold[plotting]
```

## Usage Examples
```bash
# Step 1: Generate k-points for supercell
easyunfold generate primitive_POSCAR supercell_POSCAR --kpoints KPOINTS_band

# Step 2: Run DFT calculation with generated k-points

# Step 3: Unfold and plot
easyunfold unfold --data-file easyunfold.json WAVECAR
```

Python API:
```python
from easyunfold.unfold import UnfoldKSet

# Load and unfold
unfoldset = UnfoldKSet.from_file("easyunfold.json")
unfoldset.get_spectral_weights(wavecar="WAVECAR")
unfoldset.plot_spectral_function()
```

## Performance Characteristics
- **Speed**: Efficient wavefunction reading
- **Memory**: Handles large WAVECAR files
- **Scalability**: Suitable for large supercells
- **Accuracy**: Proper symmetry averaging

## Limitations & Known Constraints
- **Wavefunction required**: Needs WAVECAR or equivalent
- **Memory**: Large supercells need significant RAM
- **k-point density**: Dense sampling increases cost
- **Spin-orbit**: Limited SOC support in some codes

## Comparison with Other Tools
- **vs BandUP**: easyunfold is Python-native, easier installation
- **vs fold2Bloch**: easyunfold has better symmetry handling
- **vs VaspBandUnfolding**: easyunfold is more automated
- **Unique strength**: Symmetry-aware k-point sampling, multi-code support

## Application Areas
- Defect band structures
- Alloy electronic structure
- Interface/heterostructure bands
- Disordered system analysis
- Dopant level identification
- Band gap engineering studies

## Best Practices
- Use dense k-point sampling for smooth spectra
- Verify primitive-supercell relationship
- Check spectral weight convergence
- Use appropriate energy broadening
- Compare with pristine band structure

## Community and Support
- GitHub issue tracker
- JOSS publication for citation
- Active development by SMTG group
- Documentation with tutorials

## Verification & Sources
**Primary sources**:
1. Official documentation: https://smtg-bham.github.io/easyunfold/
2. GitHub repository: https://github.com/SMTG-Bham/easyunfold
3. JOSS paper: https://joss.theoj.org/papers/10.21105/joss.05974
4. Popescu & Zunger, Phys. Rev. B 85, 085201 (2012) - Theory

**Confidence**: VERIFIED - Published in JOSS, active development

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE
- Source code: OPEN (GitHub, MIT)
- Developer: SMTG Birmingham
- Academic citations: JOSS publication
- Active development: Regular releases
