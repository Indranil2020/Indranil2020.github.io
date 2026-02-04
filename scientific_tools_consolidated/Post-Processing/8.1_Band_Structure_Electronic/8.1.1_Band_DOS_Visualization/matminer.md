# matminer

## Official Resources
- **Homepage**: https://hackingmaterials.lbl.gov/matminer/
- **GitHub**: https://github.com/hackingmaterials/matminer
- **Documentation**: https://hackingmaterials.lbl.gov/matminer/
- **PyPI**: https://pypi.org/project/matminer/
- **Publication**: L. Ward et al., Comput. Mater. Sci. 152, 60 (2018)
- **License**: BSD 3-Clause

## Overview
matminer is a comprehensive Python library for data mining in materials science, developed at Lawrence Berkeley National Laboratory. It provides tools for extracting features (descriptors) from materials data including compositions, structures, band structures, and density of states for machine learning applications.

**Scientific domain**: Materials informatics, machine learning, feature engineering
**Target user community**: Materials scientists using ML, data-driven materials discovery researchers

## Theoretical Background
matminer implements feature extraction based on:
- Compositional descriptors (Magpie, element properties)
- Structural descriptors (Voronoi, coordination)
- Electronic structure features (band gap, DOS moments)
- Site-based features (local environment)

## Capabilities (CRITICAL)
- **Band Structure Features**: Extract ML features from band structures
- **DOS Features**: Density of states featurization (moments, fingerprints)
- **Composition Features**: Elemental property descriptors (Magpie, etc.)
- **Structure Features**: Crystal structure descriptors (Voronoi, RDF)
- **Data Retrieval**: Access Materials Project, AFLOW, Citrine databases
- **ML Pipeline**: Scikit-learn compatible featurizers

## Key Strengths

### Comprehensive Featurizers:
- 50+ built-in featurizers
- Composition-based (Magpie, element stats)
- Structure-based (Voronoi, coordination)
- DOS-based (moments, fingerprints)
- Band structure-based (gap, effective mass)

### Database Integration:
- Materials Project API
- AFLOW database
- Citrine Informatics
- MPDS (Materials Platform)

### ML Integration:
- Scikit-learn compatible
- Pandas DataFrame workflow
- Feature selection tools
- Cross-validation support

## Inputs & Outputs
- **Input formats**:
  - pymatgen Composition objects
  - pymatgen Structure objects
  - pymatgen BandStructure objects
  - pymatgen DOS objects
  
- **Output data types**:
  - Pandas DataFrames
  - NumPy arrays
  - Feature vectors for ML

## Installation
```bash
pip install matminer
```

## Usage Examples
```python
from matminer.featurizers.composition import ElementProperty
from matminer.featurizers.dos import DOSFeaturizer
from matminer.featurizers.bandstructure import BandFeaturizer
import pandas as pd

# Composition features
ep = ElementProperty.from_preset("magpie")
df = pd.DataFrame({"composition": ["Fe2O3", "NaCl"]})
df = ep.featurize_dataframe(df, "composition")

# DOS features
dos_feat = DOSFeaturizer()
dos_features = dos_feat.featurize(dos_object)

# Band structure features
band_feat = BandFeaturizer()
band_features = band_feat.featurize(band_structure)
```

## Performance Characteristics
- **Speed**: Efficient vectorized operations
- **Scalability**: Handles large datasets
- **Memory**: Optimized for batch processing

## Limitations & Known Constraints
- **pymatgen dependency**: Requires pymatgen objects
- **API limits**: Database queries may have rate limits
- **Feature selection**: Many features may need pruning

## Comparison with Other Tools
- **vs JARVIS-Tools**: matminer focused on featurization
- **vs dscribe**: Different descriptor sets
- **Unique strength**: Comprehensive featurizer library, database integration

## Application Areas
- Materials property prediction
- High-throughput screening
- Structure-property relationships
- Band gap prediction
- Stability prediction

## Best Practices
- Use preset featurizers for standard applications
- Perform feature selection for large feature sets
- Validate features on held-out data
- Cite appropriate featurizer papers

## Verification & Sources
**Primary sources**:
1. Documentation: https://hackingmaterials.lbl.gov/matminer/
2. GitHub: https://github.com/hackingmaterials/matminer
3. L. Ward et al., Comput. Mater. Sci. 152, 60 (2018)

**Confidence**: VERIFIED - Published, widely used

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, BSD 3-Clause)
- Developer: Hacking Materials (LBNL)
- Academic citations: >500 citations
