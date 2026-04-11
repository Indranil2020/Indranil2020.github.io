# OpenKIM

## Official Resources
- Website: https://openkim.org/
- Source Repository: https://github.com/openkim
- License: Open source (various per model)

## Overview
**OpenKIM** (Knowledgebase of Interatomic Models) is a standardized framework for testing, archiving, and distributing interatomic potentials. It provides the KIM API for seamless integration with LAMMPS, ASE, and other simulation codes, with automatic verification tests.

**Scientific domain**: Standardized interatomic model repository and API  
**Target user community**: Researchers finding, testing, and using interatomic potentials

## Theoretical Methods
- KIM API standard
- Model verification tests
- Predictive capability assessment
- Standardized model format
- Cross-code compatibility

## Capabilities (CRITICAL)
- 1000+ interatomic models
- KIM API (LAMMPS, ASE, etc.)
- Automatic verification tests
- Model comparison tools
- KLIFF fitting framework
- KUSP deployment utility

**Sources**: https://openkim.org/

## Key Strengths

### Repository:
- 1000+ models
- EAM, MEAM, Tersoff, ReaxFF, MLIPs
- Verified models
- Metadata and citations

### Standards:
- KIM API
- Verification tests
- Reproducibility
- Cross-code compatibility

### Tools:
- KLIFF (fitting)
- KUSP (deployment)
- kim-query (data access)
- Model comparison

## Inputs & Outputs
- **Input formats**: KIM model specifications
- **Output data types**: Verified potentials, test results

## Interfaces & Ecosystem
- **LAMMPS**: KIM pair_style
- **ASE**: KIM calculator
- **GULP**: Integration
- **Python/C**: API

## Performance Characteristics
- **Speed**: Model-dependent
- **Accuracy**: Model-dependent
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Model access**: Free
- **Verification**: Minutes
- **MD**: Model-dependent

## Limitations & Known Constraints
- **Model quality varies**: Not all equally accurate
- **KIM API required**: For integration
- **Limited MLIPs**: Mostly classical FFs
- **Registration**: Required for publishing

## Comparison with Other Codes
- **vs NIST potentials**: OpenKIM is more comprehensive
- **vs LAMMPS potential list**: OpenKIM has verification
- **vs KLIFF**: OpenKIM is repository, KLIFF is fitting
- **Unique strength**: Standardized repository with 1000+ verified interatomic models and KIM API

## Application Areas

### Potential Discovery:
- Find potentials for elements
- Compare model accuracy
- Access verified models
- Citation tracking

### Model Development:
- Publish to OpenKIM
- Run verification tests
- KIM API integration
- Community contribution

## Best Practices
- Check verification tests before use
- Compare multiple models
- Use KIM API for consistency
- Cite model and OpenKIM

## Community and Support
- Open source
- NSF-funded
- University of Minnesota
- Active community

## Verification & Sources
**Primary sources**:
1. Website: https://openkim.org/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACCESSIBLE
- Specialized strength: Standardized repository with 1000+ verified interatomic models and KIM API
