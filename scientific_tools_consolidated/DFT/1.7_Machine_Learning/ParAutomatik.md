# ParAutomatik

## Official Resources
- Homepage: https://github.com/Teoroo-CMC/ParAutomatik
- Source Repository: https://github.com/Teoroo-CMC/ParAutomatik
- License: MIT License

## Overview
ParAutomatik is a cutting-edge workflow automation tool designed for the parameterization of Density Functional Tight Binding (DFTB) models using Machine Learning. It addresses one of the biggest bottlenecks in semi-empirical methods: the difficulty of creating accurate parameters. By combining high-throughput DFT calculations with Neural Network training, ParAutomatik automates the fitting of repulsive potentials and electronic parameters, significantly accelerating the development of transferable DFTB models.

**Scientific domain**: Machine Learning, Parameterization, Tight-Binding
**Target user community**: Method developers, Researchers needing custom potentials for specific chemistries

## Theoretical Methods
- **SCC-DFTB**: Self-Consistent Charge Density Functional Tight Binding.
- **Machine Learning**: Neural Networks (PyTorch) for regression.
- **Active Learning**: Iterative improvement of training sets.
- **Repulsive Potential Fitting**: Fitting the difference between DFT and electronic DFTB energy ($E_{rep} = E_{DFT} - E_{elec}$).
- **Spline Interpolation**: Generation of traditional spline-based repulsive potentials.

## Capabilities (CRITICAL)
- **Dataset Generation**: Automated sampling of geometries (dimers, trimers, clusters).
- **Reference Calculation**: Automatic execution of DFT reference runs (via ASE calculators).
- **ML Fitting**: Training of NNs to predict repulsive energies.
- **Validation**: Automatic benchmarking against DFT forces and energies.
- **Output Generation**: Production of `.skf` files compatible with DFTB+.

## Key Strengths

### Automation:
- Replaces manual "by-eye" fitting or simple polynomial fits with robust ML workflows.
- Handles the complexity of multi-element parameterization consistency.

### Physics-Informed ML:
- Uses Neural Networks to learn the complex environment dependence of repulsion, or to generate improved 2-body splines.

### Reproducibility:
- Defines the parameterization protocol as code, making the provenance of potentials clear.

## Inputs & Outputs
- **Inputs**:
  - List of elements (e.g., `["C", "H", "O"]`).
  - Reference method definition (e.g., PBE/def2-TZVP).
  - Configuration settings (cutoff distances, NN architecture).
- **Outputs**:
  - `dataset.db`: ASE database of training structures.
  - `model.pt`: Trained PyTorch model.
  - `*-*.skf`: Final Slater-Koster files ready for DFTB+.
  - `report.pdf`: Validation report showing RMS errors.

## Interfaces & Ecosystem
- **ASE**: Built entirely around the Atomic Simulation Environment.
- **DFTB+**: The primary target engine for the resulting parameters.
- **PyTorch**: The ML backend.
- **dftbpara**: Compatible with/Alternative to other fitting tools.

## Advanced Features
- **Delta Learning**: Predicting the correction to a baseline model.
- **Force Matching**: Training on forces (gradients) as well as energies for better dynamics stability.

## Performance Characteristics
- **Speed**: Training takes minutes to hours (GPU supported); Reference data generation is the bottleneck.
- **Scalability**: Can parameterize complex multi-element sets if training data is sufficient.

## Computational Cost
- **Moderate**: dominated by the cost of the *ab initio* reference calculations.

## Limitations & Known Constraints
- **Two-Body Limit**: Standard `.skf` format is strictly two-body; if ParAutomatik is used to generate these, it compresses many-body physics into 2-body terms (loss of accuracy). (Note: Newer versions may support Many-Body Repulsion).
- **Data Hungry**: Quality depends entirely on the coverage of the training set.

## Comparison with Other Codes
- **vs TBFIT**: TBFIT fits electronic band structures (SK integrals); ParAutomatik typically fits the *repulsive* part (total energy/forces). They are complementary.
- **vs Hotbit**: Hotbit is a code that *can* fit parameters manually/scripted; ParAutomatik is a dedicated ML workflow.
- **Unique strength**: Bringing modern ML workflows to the tedious task of repulsive potential fitting.

## Application Areas
- **MOFs/COFs**: Creating parameters for novel porous materials.
- **Catalysis**: Tuning parameters for specific metal-organic interfaces.
- **High-Pressure**: Fitting potentials for extreme conditions where standard sets fail.

## Best Practices
- **Span Config Space**: Ensure training data covers the distances/angles seen in production runs.
- **Check Limits**: Verify behavior at short range (avoid hole collapse) and long range (smooth cutoff).
- **Validate**: Always run a test calculation on a system *not* in the training set (e.g., a crystal bulk modulus).

## Community and Support
- **GitHub**: Developed by the Teoroo-CMC group (Calcula).
- **Documentation**: Jupyter notebook tutorials available.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/Teoroo-CMC/ParAutomatik
2. Publications by the Teoroo group (e.g., on DFTB parameterization).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (MIT)
- Functionality: Functional ML workflow.
