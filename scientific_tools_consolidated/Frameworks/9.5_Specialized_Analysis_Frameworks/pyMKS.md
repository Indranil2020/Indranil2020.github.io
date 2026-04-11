# pyMKS

## Official Resources
- Source Repository: https://github.com/materialsinnovation/pymks
- Documentation: https://pymks.readthedocs.io/
- PyPI: https://pypi.org/project/pymks/
- License: Open source (MIT)

## Overview
**pyMKS** (Materials Knowledge System) is a Python framework for materials data analytics using the Materials Knowledge System paradigm. It provides tools for microstructure quantification, localization relationships, and homogenization/linkage using 2-point statistics and MKS regression.

**Scientific domain**: Materials Knowledge System, microstructure quantification, 2-point statistics  
**Target user community**: Researchers analyzing microstructure-property relationships in materials science

## Theoretical Methods
- 2-point spatial correlations (statistics)
- Materials Knowledge System (MKS) regression
- Microstructure quantification
- Localization relationships
- Homogenization linkages
- Discrete Fourier transform for correlations

## Capabilities (CRITICAL)
- 2-point spatial correlation calculation
- MKS regression for localization
- Microstructure generation and sampling
- PCA of microstructure statistics
- Homogenization linkages
- Discrete Fourier transform acceleration

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### Microstructure Analysis:
- 2-point statistics
- Microstructure quantification
- PCA dimensionality reduction
- Statistical representation

### MKS Framework:
- Localization (micro→local property)
- Homogenization (micro→effective property)
- Regression-based linkages
- Efficient computation via FFT

### Data-Driven:
- No physics simulation needed
- Purely data-driven
- Statistical learning
- Scalable

## Inputs & Outputs
- **Input formats**: Microstructure images, 2D/3D grids
- **Output data types**: Spatial correlations, property predictions, MKS coefficients

## Interfaces & Ecosystem
- **NumPy**: Computation
- **scikit-learn**: ML models
- **matplotlib**: Visualization
- **Python**: Core language

## Performance Characteristics
- **Speed**: Fast (FFT-based)
- **Accuracy**: Data-dependent
- **System size**: Microstructure grids
- **Memory**: Moderate

## Computational Cost
- **2-point statistics**: Seconds to minutes
- **MKS regression**: Minutes
- **No DFT needed**: Image/data-based

## Limitations & Known Constraints
- **Microstructure focus**: Not for electronic structure
- **Grid-based**: Requires discretized microstructure
- **Linear assumption**: Basic MKS is linear
- **Data quantity**: Needs sufficient samples

## Comparison with Other Codes
- **vs pymatgen**: pyMKS is microstructure, pymatgen is atomic
- **vs matminer**: pyMKS is spatial statistics, matminer is descriptors
- **vs DREAM3D**: pyMKS is Python, DREAM3D is C++ with GUI
- **Unique strength**: Materials Knowledge System with 2-point spatial correlations and MKS regression for microstructure-property linkages

## Application Areas

### Microstructure-Property:
- Structure-property linkages
- Homogenization models
- Localization predictions
- Process-structure-property

### Materials Design:
- Microstructure optimization
- Inverse design
- Statistical learning
- High-throughput screening

## Best Practices

### Data:
- Use sufficient microstructure samples
- Check spatial correlation convergence
- Validate with known systems
- Use PCA for dimensionality reduction

## Community and Support
- Open source (MIT)
- PyPI installable
- ReadTheDocs documentation
- Materials Innovation team
- Published in Computational Materials Science

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/materialsinnovation/pymks

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Materials Knowledge System with 2-point spatial correlations and MKS regression
