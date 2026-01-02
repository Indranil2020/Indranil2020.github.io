# xTB

## Official Resources
- Homepage: https://github.com/grimme-lab/xtb
- Documentation: https://xtb-docs.readthedocs.io/
- Source Repository: https://github.com/grimme-lab/xtb
- License: GNU Lesser General Public License v3.0

## Overview
xTB (extended tight-binding) is a semiempirical quantum chemistry package implementing various tight-binding methods with parametrizations ranging from GFN0-xTB to GFN2-xTB. It is extremely fast, robust, and covers the entire periodic table, making it ideal for large-scale screening, conformer searches, and preliminary geometry optimizations.

**Scientific domain**: Computational chemistry, conformational analysis, screening, high-throughput calculations  
**Target user community**: Chemists needing fast quantum calculations for large molecules or extensive screening

## Theoretical Methods
- GFN0-xTB: Minimal basis tight-binding
- GFN1-xTB: Generalized Force Field for xTB version 1
- GFN2-xTB: Improved parametrization with H-bonding
- GFN-FF: Force field approximation (ultra-fast)
- Implicit solvation models (GBSA, ALPB)
- D3 and D4 dispersion corrections
- Geometry, Frequency, Noncovalent interactions (GFN methods)
- Periodic boundary conditions
- External electric fields

## Capabilities (CRITICAL)
- Geometry optimization (molecules and periodic systems)
- Single-point energy calculations
- Molecular dynamics (NVE, NVT, NVT ensembles)
- Metadynamics for conformational sampling
- Vibrational frequencies and IR/Raman spectra
- Thermochemistry (enthalpies, entropies, Gibbs energies)
- Transition state searches
- Reaction path calculations
- Conformational searches and ensemble generation
- Solvation free energies
- pKa value estimation
- Protonation site prediction
- Noncovalent interaction analysis
- Molecular properties (dipole, polarizability)
- Band structure and DOS for periodic systems
- Systems up to 10,000+ atoms routinely
- Extremely fast (seconds to minutes for typical molecules)

**Sources**: Official xTB documentation (https://xtb-docs.readthedocs.io/), cited in 6/7 source lists

## Key Advantages

### Speed:
- Orders of magnitude faster than DFT
- Covers entire periodic table (H-Rn)
- Typical small molecule: <1 second
- Protein (1000 atoms): minutes
- Enables extensive conformer generation

### Robustness:
- Stable SCF convergence
- Minimal failures compared to other semi-empirical methods
- Handles challenging systems (metal complexes, radicals)
- Reasonable geometries across chemical space

### Accuracy:
- Good geometries (RMSD ~0.1 Å for organic molecules)
- Reasonable energetics for conformers
- Reliable thermochemistry
- Better than most force fields, faster than DFT

### Coverage:
- All elements H through Rn
- Organic, inorganic, organometallic
- Main group and transition metals
- No parameterization gaps

## GFN Method Comparison

### GFN0-xTB:
- Fastest variant
- Minimal basis
- Good for initial screening
- Less accurate than GFN1/GFN2

### GFN1-xTB:
- Balanced speed and accuracy
- Good for general applications
- Reliable energies and geometries

### GFN2-xTB:
- Most accurate variant
- Improved H-bonding description
- Better thermochemistry
- Recommended for most applications

### GFN-FF:
- Ultra-fast force field mode
- For very large systems or long MD
- Quantum-informed force field
- Smooth potential energy surfaces

## Inputs & Outputs
- **Input formats**:
  - XYZ coordinates (primary)
  - SDF molecular files
  - Turbomole coord format
  - POSCAR for periodic systems
  - Command-line driven
  
- **Output data types**:
  - Energies and gradients
  - Optimized geometries (XYZ, Turbomole)
  - Vibrational frequencies
  - Thermochemical data
  - Charges (Mulliken, CM5, EEQ)
  - Molecular orbitals
  - Property files

## Interfaces & Ecosystem
- **CREST integration**:
  - Conformer-Rotamer Ensemble Sampling Tool
  - Automated conformer searches
  - Extensive sampling capabilities
  
- **Python interfaces**:
  - ASE calculator for xTB
  - Direct Python bindings (xtb-python)
  - Workflow automation
  
- **Standalone utilities**:
  - crest for conformer generation
  - xtb4stda for excited states preparation
  
- **Workflow tools**:
  - Compatible with standard quantum chemistry workflows
  - Pre-optimizer for DFT calculations

## Workflow and Usage

### Typical Workflows:

#### 1. Conformer Search:
```bash
# Generate conformers with CREST
crest molecule.xyz --gfn2 --alpb water

# xTB optimization
xtb molecule.xyz --opt --gfn2 --alpb water
```

#### 2. Geometry Optimization:
```bash
# Single optimization
xtb molecule.xyz --opt tight --gfn2

# With solvation
xtb molecule.xyz --opt --gbsa water
```

#### 3. Thermochemistry:
```bash
# Frequency calculation
xtb molecule.xyz --hess --gfn2

# Extract thermochemical data from output
```

#### 4. High-Throughput Screening:
```bash
# Loop over many structures
for file in *.xyz; do
  xtb $file --gfn2 --chrg 0 >> results.txt
done
```

## Advanced Features

### Conformational Sampling (CREST):
- Metadynamics-based sampling
- iMTD-GC algorithm
- Automated conformer generation
- Degenerate conformer identification
- Boltzmann weighting

### Solvation Models:
- GBSA: Generalized Born with surface area
- ALPB: Analytical Linearized Poisson-Boltzmann
- Wide range of solvents parameterized
- Solvation free energies

### Property Calculations:
- Fukui functions
- Electron localization function (ELF)
- Partial charges (multiple schemes)
- Bond orders
- Noncovalent interaction index

### Periodic Systems:
- Crystal structure optimization
- Band structures
- Density of states
- Bulk modulus calculations

## Performance Characteristics
- **Speed**: 10-1000x faster than DFT
- **Scaling**: Near-linear for many operations
- **Memory**: Very low; minimal RAM requirements
- **Parallelization**: OpenMP support
- **Typical times**:
  - Small molecule (20 atoms): <1 second
  - Medium molecule (100 atoms): 1-10 seconds
  - Protein fragment (500 atoms): 1-5 minutes

## Limitations & Known Constraints
- **Semiempirical accuracy**: Not as accurate as high-level ab initio
- **Energetics**: Barriers and reaction energies less reliable
- **Excited states**: Not directly calculated (use sTDA separately)
- **Strong correlation**: Not designed for multireference systems
- **Absolute energies**: Relative energies more reliable
- **Metal complexes**: Good but not always quantitative
- **Learning curve**: Low; command-line interface straightforward
- **Documentation**: Good and accessible
- **Platform**: Linux, macOS, Windows (via WSL)

## Comparison with Other Methods
- **vs DFT**: 100-1000x faster, lower accuracy
- **vs DFTB+**: xTB no parameter files needed, better coverage
- **vs PM6/PM7**: xTB more robust, broader coverage
- **vs Force Fields**: More accurate, still very fast
- **Sweet spot**: Preliminary optimizations, conformer searches, screening

## Application Areas

### Drug Discovery:
- Conformer generation for docking
- Ligand preparation
- pKa prediction
- Tautomer enumeration

### Chemical Reaction Screening:
- Reaction pathway exploration
- Barrier estimation
- Product prediction

### Materials Science:
- Crystal structure prediction
- MOF screening
- Supramolecular assemblies

### Methodology:
- Pre-optimization for DFT
- Initial guess generation
- Filtering for expensive calculations

## Integration in Computational Workflows

### As Pre-optimizer:
1. xTB geometry optimization (fast)
2. DFT single-point or optimization (accurate)
3. Post-processing

### For Ensemble Generation:
1. CREST conformer search (extensive)
2. xTB re-optimization and ranking
3. DFT refinement of top conformers

### High-Throughput Screening:
1. Generate large library
2. xTB property calculation (fast)
3. ML model training or direct selection
4. DFT validation

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/grimme-lab/xtb
2. Documentation: https://xtb-docs.readthedocs.io/
3. C. Bannwarth et al., J. Chem. Theory Comput. 15, 1652 (2019) - GFN2-xTB
4. S. Grimme et al., J. Chem. Theory Comput. 13, 1989 (2017) - GFN1-xTB
5. P. Pracht et al., Phys. Chem. Chem. Phys. 22, 7169 (2020) - CREST

**Secondary sources**:
1. xTB tutorials and examples
2. CREST documentation
3. Published applications in conformer generation
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE (GitHub)
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, LGPL v3)
- Community support: Very active (GitHub issues, discussions)
- Academic citations: >500 (GFN papers)
- Active development: Regular releases, continuous improvements
- Benchmark validation: Extensive benchmarks published
- Wide adoption: Standard tool for conformer generation and pre-optimization
