# CRYSOL

## Official Resources
- Homepage: https://www.embl-hamburg.de/biosaxs/software.html
- Documentation: https://www.embl-hamburg.de/biosaxs/manuals/crysol.html
- Source Repository: Part of ATSAS package (Academic license)
- License: Proprietary (Free for academic use)

## Overview
CRYSOL is a program for evaluating the solution scattering (SAXS/SANS) from macromolecules with known atomic structure and fitting it to experimental data. It uses multipole expansion of the scattering amplitudes to calculate the spherically averaged scattering pattern, taking into account the hydration shell and solvent density. It is part of the ATSAS software suite.

**Scientific domain**: Small-angle scattering (SAXS/SANS), structural biology, solution scattering  
**Target user community**: Structural biologists, biophysicists

## Theoretical Methods
- Multipole expansion of scattering intensity
- Spherical harmonics
- Hydration shell modeling (boundary layer)
- Excluded volume calculation
- Fitting to experimental data (chi-square minimization)

## Capabilities (CRITICAL)
- Calculation of theoretical scattering curves from PDB structures
- Fitting of theoretical curves to experimental SAXS/SANS data
- Estimation of hydration shell contrast and volume
- Handling of mixtures
- Prediction of Radius of Gyration (Rg) and Forward Scattering (I(0))

**Sources**: CRYSOL documentation, J. Appl. Cryst. 28, 768 (1995)

## Key Strengths

### Fast Algorithm:
- Multipole expansion
- Efficient calculation
- Quick fitting
- Large proteins

### Hydration Modeling:
- Boundary layer
- Contrast variation
- Solvent effects
- Realistic model

### ATSAS Integration:
- Complete suite
- GUI available
- Standard tool
- Active support

## Inputs & Outputs
- **Input formats**: PDB file (atomic coordinates), experimental data file (.dat)
- **Output data types**: .fit (fitted curve), .log (parameters), .alm (multipole coefficients)

## Interfaces & Ecosystem
- **ATSAS**: Core component of the ATSAS suite
- **Primus**: GUI for analysis
- **Python**: ATSAS provides Python bindings (sas7bdat)

## Workflow and Usage
1. Prepare PDB structure.
2. Run CRYSOL: `crysol structure.pdb experimental_data.dat`
3. Check fit quality (Chi-square).
4. Adjust parameters (contrast, hydration) if necessary.

## Performance Characteristics
- Fast multipole algorithm
- Efficient for standard proteins
- Scales with number of atoms and multipole order

## Limitations & Known Constraints
- **Atomic model needed**: Requires PDB structure
- **Spherical averaging**: Loses orientation info
- **Hydration parameters**: May need adjustment
- **Proprietary**: Academic license required

## Comparison with Other Tools
- **vs FoXS**: Similar capabilities, different algorithms
- **vs SASSIE**: CRYSOL faster, SASSIE more flexible
- **vs experimental**: Validation tool
- **Unique strength**: Standard SAXS fitting tool

## Application Areas
- Protein structure validation
- Conformational analysis in solution
- Complex formation studies
- Unfolding/folding monitoring

## Best Practices
- Use high-quality PDB structures
- Check hydration shell parameters
- Validate chi-square values
- Compare multiple conformations

## Community and Support
- Developed by EMBL Hamburg (Svergun Group)
- Standard tool in SAXS community
- Annual courses and active support

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.embl-hamburg.de/biosaxs/software.html
2. Publication: D. I. Svergun et al., J. Appl. Cryst. 28, 768 (1995)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: ACADEMIC (ATSAS)
- Development: ACTIVE (EMBL)
- Applications: SAXS simulation, fitting, protein structure
