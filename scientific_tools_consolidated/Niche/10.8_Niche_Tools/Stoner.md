# Stoner (Stoner Python Code)

## Official Resources
- Homepage: https://stoner.readthedocs.io/
- Documentation: https://stoner.readthedocs.io/
- Source Repository: https://github.com/PhysicsStoner/Stoner
- License: GPL v3

## Overview
The Stoner package is a Python library containing a collection of utilities for managing and analyzing experimental data, particularly in physics. It provides a `Data` class (subclassing Pandas DataFrame) with metadata handling and specialized plotting/analysis methods useful for experimentalists (e.g., hysteresis loops, IV curves).

**Scientific domain**: Experimental physics data analysis  
**Target user community**: Experimental physicists (Leeds University origin)

## Capabilities (CRITICAL)
- **Data Loading**: Readers for various instrument formats (VSM, resistivity, etc.).
- **Metadata**: Dictionary-like metadata associated with data.
- **Analysis**: Curve fitting, smoothing, differentiation.
- **Plotting**: Publication-quality plotting wrappers around Matplotlib.

**Sources**: Stoner documentation

## Inputs & Outputs
- **Input formats**: CSV, instrument text files
- **Output data types**: Plots, processed files

## Interfaces & Ecosystem
- **Pandas/Matplotlib/Scipy**: Built on top of the standard stack.

## Workflow and Usage
1. `d = Stoner.Data.load('file.dat')`
2. `d.plot(x='Field', y='Moment')`
3. `d.analyze(...)`

## Performance Characteristics
- Python convenience wrapper.

## Application Areas
- Magnetism (VSM data)
- Transport measurements
- General lab data analysis

## Community and Support
- Developed at University of Leeds (Cond. Matt. group).

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/PhysicsStoner/Stoner

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Experimental data analysis
