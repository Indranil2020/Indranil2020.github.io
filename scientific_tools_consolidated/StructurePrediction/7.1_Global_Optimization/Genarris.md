# Genarris

## Overview
Genarris is an open-source Python package for generating random molecular crystal structures with physical constraints. It serves as a structure generator for seeding crystal structure prediction algorithms and training machine learning models.

## Theoretical Basis
- Random structure generation with symmetry constraints
- Wyckoff position occupation (general and special)
- Physical constraints (van der Waals radii, density)
- Space group compatibility with molecular symmetry
- Volume estimation from molecular properties

## Key Capabilities
- Random molecular crystal structure generation
- All compatible space groups supported
- General and special Wyckoff positions
- Physical constraint enforcement
- Parallel structure generation

**Sources**: Comp. Phys. Comm. 250, 107170 (2020)

## Key Strengths

### Structure Generation:
- All 230 space groups
- Special Wyckoff positions
- Flexible molecules supported

### Physical Constraints:
- Van der Waals overlap prevention
- Density constraints
- Volume estimation

### Flexibility:
- Python-based
- Easy integration
- Customizable constraints

## Inputs & Outputs
- **Input formats**: Molecular geometry (xyz, mol2), configuration file
- **Output data types**: Crystal structures (cif, json), generation logs

## Interfaces & Ecosystem
- **Structure formats**: CIF, JSON, Pymatgen
- **CSP codes**: Can seed USPEX, AIRSS, GAtor
- **ML**: Training data generation

## Workflow and Usage
1. Prepare molecular geometry
2. Configure generation parameters
3. Run: `python -m Genarris`
4. Filter generated structures
5. Use for CSP or ML training

## Performance Characteristics
- Fast structure generation
- Parallelizable
- Efficient constraint checking

## Computational Cost
- Minimal (no energy calculations)
- Scales with number of structures
- Parallel generation supported

## Best Practices
- Use appropriate van der Waals radii
- Generate diverse space groups
- Filter by density constraints
- Validate molecular geometry

## Limitations & Known Constraints
- Generation only (no optimization)
- Requires external CSP for ranking
- Molecular crystals focus

## Application Areas
- Molecular crystal CSP seeding
- ML training data generation
- Polymorph screening
- Crystal engineering

## Comparison with Other Codes
- **vs PyXtal**: Genarris molecular-focused, PyXtal more general
- **vs GAtor**: Genarris generation only, GAtor includes optimization
- **vs RandSpg**: Similar purpose, different implementation
- **Unique strength**: Molecular crystals, special Wyckoff positions, CSP seeding

## Community and Support
- Open-source (MIT License)
- GitHub repository
- Active development (Marom group)
- Published methodology

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/timcrose/Genarris
2. Publication: Comp. Phys. Comm. 250, 107170 (2020)
3. arXiv: https://arxiv.org/abs/1909.10629

**Secondary sources**:
1. Genarris documentation
2. Molecular CSP literature

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Documentation: AVAILABLE
- Source: OPEN (GitHub, MIT)
- Development: ACTIVE
- Applications: Molecular crystal generation, CSP seeding
