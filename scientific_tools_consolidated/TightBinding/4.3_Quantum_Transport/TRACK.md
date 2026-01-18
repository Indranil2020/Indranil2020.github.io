# TRACK

## Official Resources
- **Source Code**: [GitHub (rvenderbos/TRACK)](https://github.com/rvenderbos/TRACK) - *Note: Verify exact link via search.*
- **License**: MIT License

## Overview
**TRACK** (TRAnsport properties for Correlated materials using Kubo formalism) is a Python 3 code designed to calculate temperature-dependent transport coefficients in solids. It utilizes the **linear-response Kubo formalism** to compute electrical conductivity, thermal conductivity, Seebeck coefficient, and the Lorenz number. A key feature of TRACK is its careful handling of current operators in **interacting systems**, making it suitable for Hamiltonians derived from Hartree-Fock or hybrid functionals where non-local potentials affect the velocity operator.

**Scientific domain**: Correlated Electrons, Thermoelectrics, Transport Theory
**Target user community**: Theorists working on transport in complex oxides and strongly correlated metals

## Theoretical Methods
- **Kubo Formalism**: Calculation of the current-current correlation function $\Pi(\omega)$ in the DC limit.
- **Interactions**: Correct implementation of the velocity operator $\mathbf{v} = \frac{i}{\hbar} [H, \mathbf{r}]$ for non-local potentials (e.g., Fock exchange).
- **Integration**: Tetrahedron method or dense mesh integration over the Brillouin Zone.
- **Scattering**: Constant relaxation time approximation ($\tau$) or energy-dependent scattering rates.

## Capabilities
- **Coefficients**:
  - Electrical Conductivity Tensor ($\sigma_{\alpha\beta}$).
  - Electronic Thermal Conductivity ($\kappa_{e}$).
  - Seebeck Coefficient ($S$).
  - Lorenz Number ($L = \kappa_e / (T \sigma)$).
- **Analysis**:
  - Temperature dependence of transport ($T$-scans).
  - Band-by-band decomposition of currents.
  - Optical conductivity (AC limit).

## Key Strengths
- **Correlations**: Specifically addresses the "Peierls substitution failure" in non-local Hamiltonians, ensuring gauge-invariant transport results for correlated models.
- **Pythonic**: Easy to inspect and modify, leveraging NumPy for tensor operations.
- **Thermoelectrics**: Direct calculation of power factors and efficiency metrics.

## Inputs & Outputs
- **Inputs**:
  - Eigenvalues and Eigenvectors (from DFT or TB).
  - Velocity matrix elements (critical for interacting parts).
  - k-mesh definitions.
- **Outputs**:
  - Text files containing $\sigma(T)$, $S(T)$, $\kappa(T)$.

## Interfaces & Ecosystem
- **Hamiltonians**: Can interface with output from tight-binding codes or DFT codes (if matrix elements are provided).
- **Ecosystem**: Relies on standard Python scientific stack (NumPy, SciPy).

## Performance Characteristics
- **Speed**: Python overhead is minimal for dense matrix operations; bottleneck is the number of k-points and bands.
- **Parallelism**: Easy to parallelize over temperature or k-points (multiprocessing).

## Comparison with Other Codes
- **vs. BoltzTrap**: BoltzTrap uses semi-classical Boltzmann theory (group velocities); TRACK uses the fully quantum mechanical Kubo formula, which captures interband transitions (optical conductivity) and can treat scattering more rigorously.
- **vs. LinReTraCe**: Similar scope (Kubo); TRACK has a specific emphasis on the velocity operator distinctions in interacting systems.

## Application Areas
- **Bad Metals**: Violation of the Wiedemann-Franz law in correlated systems.
- **Thermoelectrics**: High-throughput screening of $S$ and $\sigma$.
- **Optical Response**: Drude weight and interband optical transitions.

## Community and Support
- **Development**: Drexel University / University of Pennsylvania (R. J. M. Venderbos).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/rvenderbos/TRACK](https://github.com/rvenderbos/TRACK)
- **Verification status**: ✅ VERIFIED
  - Research code active in the correlated electron community.
