# dpdata

## Official Resources
- Source Repository: https://github.com/deepmodeling/dpdata
- Documentation: https://docs.deepmodeling.com/projects/dpdata/
- PyPI: https://pypi.org/project/dpdata/
- License: Open source (LGPL-3.0)

## Overview
**dpdata** is a Python package for manipulating atomistic data from various computational chemistry software. It provides unified data conversion between VASP, QE, Gaussian, LAMMPS, ABINIT, CP2K, and many other codes, with DeepMD-kit integration for ML potential training data.

**Scientific domain**: Atomistic data format conversion, ML training data preparation  
**Target user community**: Researchers converting between DFT/MD code formats and preparing ML training data

## Theoretical Methods
- Multi-format data conversion
- Structure and trajectory handling
- ML training data formatting
- DeepMD-kit integration
- Force/energy data management

## Capabilities (CRITICAL)
- Read/write 20+ code formats (VASP, QE, Gaussian, LAMMPS, etc.)
- Structure conversion between codes
- Trajectory handling
- DeepMD-kit training data format
- Batch data processing
- Type map management

**Sources**: GitHub repository, documentation

## Key Strengths

### Multi-Code:
- VASP, QE, Gaussian, LAMMPS, ABINIT, CP2K, SIESTA, etc.
- Unified API across formats
- Structure and trajectory
- Force and energy data

### ML Integration:
- DeepMD-kit training format
- Training/test data splitting
- Type map handling
- Batch conversion

### Efficient:
- Batch processing
- Minimal memory
- Fast I/O
- Command-line interface

## Inputs & Outputs
- **Input formats**: 20+ DFT/MD code formats
- **Output data types**: Converted structures, DeepMD data

## Interfaces & Ecosystem
- **DeepMD-kit**: ML potential training
- **ASE**: Structure handling
- **NumPy**: Computation
- **Python**: Core language

## Performance Characteristics
- **Speed**: Fast (I/O)
- **System size**: Any
- **Memory**: Low to moderate

## Computational Cost
- **Conversion**: Seconds
- **No DFT needed**: Data manipulation only

## Limitations & Known Constraints
- **I/O only**: No analysis features
- **DeepMD focus**: Optimized for DeepMD-kit
- **Format limitations**: Some formats partially supported
- **Documentation**: Could be more extensive

## Comparison with Other Codes
- **vs ASE I/O**: dpdata supports more formats, DeepMD integration
- **vs pymatgen I/O**: dpdata is format conversion, pymatgen is analysis
- **vs cif2cell**: dpdata is multi-format, cif2cell is CIF to DFT input
- **Unique strength**: Unified multi-format atomistic data conversion with DeepMD-kit integration for ML training

## Application Areas

### Data Conversion:
- VASP to QE and vice versa
- DFT output to ML training data
- Trajectory format conversion
- Batch data processing

### ML Potential Training:
- DeepMD-kit data preparation
- Training/test splitting
- Multi-code data aggregation
- Type map management

## Best Practices

### Conversion:
- Check format compatibility
- Validate converted structures
- Use batch mode for large datasets
- Preserve type maps

## Community and Support
- Open source (LGPL-3.0)
- PyPI installable
- DeepModeling community
- Active development

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/deepmodeling/dpdata

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Unified multi-format atomistic data conversion with DeepMD-kit integration
