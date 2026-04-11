# Matbench

## Official Resources
- Homepage: https://matbench.materialsproject.org/
- Documentation: https://hackingmaterials.lbl.gov/matbench/
- Source Repository: https://github.com/materialsproject/matbench
- License: MIT License

## Overview
Matbench is an automated leaderboard and benchmark suite for materials science machine learning. It consists of a curated set of 13 materials datasets (covering properties like band gap, formation energy, elastic moduli) and a python package to simplify testing and submission. It aims to standardize the comparison of ML algorithms in materials science.

**Scientific domain**: Machine learning benchmarking  
**Target user community**: ML researchers, materials informaticians

## Capabilities (CRITICAL)
- **Datasets**: 13 diverse datasets (experimental and computational).
- **Leaderboard**: Online ranking of algorithms (RF, Graph Networks, etc.).
- **Python API**: `MatbenchBenchmark` class to automate cross-validation and scoring.
- **Metrics**: MAE, RMSE, etc.

**Sources**: Matbench website, npj Comput. Mater. 6, 181 (2020)

## Inputs & Outputs
- **Input formats**: Algorithms (Python functions)
- **Output data types**: Benchmark scores

## Interfaces & Ecosystem
- **matminer**: Often used for baselines.
- **scikit-learn**: Compatible.

## Workflow and Usage
1. `mb = MatbenchBenchmark()`
2. `for task in mb.tasks: task.record(predictions)`
3. `mb.to_file("results.json.gz")`
4. Submit to leaderboard.

## Performance Characteristics
- Standardized nested cross-validation ensures fair comparison.

## Application Areas
- Validating new ML models
- Comparing descriptors

## Community and Support
- Materials Project (Dunn, Jain, et al.)
- Active community challenges

## Verification & Sources
**Primary sources**:
1. Homepage: https://matbench.materialsproject.org/
2. Publication: A. Dunn et al., npj Comput. Mater. 6, 181 (2020)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN
- Development: ACTIVE
- Applications: ML benchmarking
