# sisl

## Overview
**sisl** is a Python toolbox for electronic structure calculations, providing capabilities for post-analysis and large-scale tight-binding DFT/NEGF calculations. It was developed to handle, manipulate, and analyze output from DFT programs.

## Official Resources
- **GitHub**: https://github.com/zerothi/sisl
- **Documentation**: https://zerothi.github.io/sisl/
- **PyPI**: https://pypi.org/project/sisl/

## Capabilities
- **Tight-Binding**: Create orthogonal/non-orthogonal TB matrices
- **Electronic Structure**: Band structure and DOS analysis
- **NEGF**: Non-equilibrium Green's function support
- **File I/O**: Read/write multiple DFT code formats
- **Geometry**: Structure manipulation and analysis

## Key Features
- SIESTA/TranSIESTA native support
- VASP, GULP, BigDFT, ScaleUp compatibility
- Wannier90 interface
- Large-scale calculations
- Sparse matrix handling

## Supported Codes
- SIESTA / TranSIESTA
- VASP
- GULP
- BigDFT
- Wannier90
- OpenMX

## Inputs & Outputs
- **Inputs**: Various DFT output formats
- **Outputs**: Processed data, TB models, plots

## Installation
```bash
pip install sisl
conda install -c conda-forge sisl
```

## Verification & Sources
- **Status**: ✅ VERIFIED
- **Confidence**: VERIFIED
- **Developer**: Nick Papior (zerothi)
