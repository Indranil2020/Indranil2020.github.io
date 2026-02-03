# RandSpg

## Overview
RandSpg is a random symmetric crystal structure generator that creates structures using Wyckoff positions. It is part of the XtalOpt ecosystem and generates structures with specified space group symmetry.

## Theoretical Basis
- Space group symmetry
- Wyckoff position occupation
- Random coordinate generation
- Symmetry-constrained structures
- Volume estimation

## Key Capabilities
- Random structure generation with symmetry
- All 230 space groups supported
- Wyckoff position handling
- Integration with XtalOpt
- Fast structure generation

**Sources**: XtalOpt documentation, GitHub repository

## Key Strengths

### Symmetry:
- All space groups
- Proper Wyckoff handling
- Symmetry preservation

### Speed:
- Fast generation
- No energy calculations
- Efficient sampling

### Integration:
- XtalOpt ecosystem
- Standard output formats

## Inputs & Outputs
- **Input formats**: Composition, space group, constraints
- **Output data types**: Crystal structures (POSCAR, CIF)

## Interfaces & Ecosystem
- **XtalOpt**: Primary integration
- **Structure formats**: VASP, CIF
- **CSP codes**: Seeding for optimization

## Workflow and Usage
1. Specify composition and space group
2. Set volume/distance constraints
3. Generate random structures
4. Use for CSP seeding or ML training

## Performance Characteristics
- Very fast generation
- No computational bottleneck
- Scales with structure complexity

## Computational Cost
- Minimal (generation only)
- No energy calculations
- Fast constraint checking

## Best Practices
- Use appropriate volume estimates
- Generate diverse space groups
- Apply distance constraints
- Filter unphysical structures

## Limitations & Known Constraints
- Generation only (no optimization)
- Requires external ranking
- May generate unstable structures

## Application Areas
- CSP algorithm seeding
- ML training data
- Structure space exploration
- Symmetry-constrained generation

## Comparison with Other Codes
- **vs PyXtal**: Similar purpose, different implementation
- **vs Genarris**: RandSpg inorganic-focused
- **Unique strength**: XtalOpt integration, Wyckoff handling

## Community and Support
- Open-source (GitHub)
- XtalOpt community
- Documentation available

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/xtalopt/randSpg
2. XtalOpt: https://xtalopt.github.io/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: MAINTAINED
- Applications: Random structure generation
