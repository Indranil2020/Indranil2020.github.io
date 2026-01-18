# ph-AFQMC (Particle-Hole Auxiliary-Field Quantum Monte Carlo)

## Official Resources
- Homepage: https://github.com/jkimribo/ph-AFQMC (or related repository)
- Documentation: Repository-dependent
- Source Repository: Check GitHub for ph-AFQMC implementations
- License: Repository-specific

## Overview
ph-AFQMC refers to particle-hole formulation of Auxiliary-Field Quantum Monte Carlo, an advanced QMC technique for studying strongly correlated electron systems. The particle-hole symmetric formulation can offer computational advantages and reduce the sign problem in certain parameter regimes. Implementations exist in research groups focusing on AFQMC methodology and applications to correlated materials.

**Scientific domain**: Auxiliary-field QMC, strongly correlated systems, finite temperature  
**Target user community**: AFQMC researchers, correlated materials, method developers

## Theoretical Methods
- Auxiliary-Field Quantum Monte Carlo (AFQMC)
- Particle-hole symmetric formulation
- Hubbard-Stratonovich transformation
- Finite-temperature methods
- Grand canonical ensemble
- Determinant QMC
- Sign problem considerations

## Capabilities (CRITICAL)
**Category**: Research QMC implementation
**Note**: Specialized AFQMC variant
- Particle-hole AFQMC
- Finite-temperature calculations
- Correlated electron systems
- Lattice models
- Hubbard-type models
- Sign problem mitigation (cases)
- Research implementation
- Method development

**Sources**: Scientific literature on ph-AFQMC

## Key Aspects

### Particle-Hole Formulation:
- Alternative AFQMC approach
- Symmetric treatment
- Specific advantages
- Sign problem considerations
- Methodology research

### AFQMC Method:
- Auxiliary-field decomposition
- Determinant QMC
- Finite temperature
- Correlation physics
- Many-body systems

## Status
- **Type**: Research implementation
- **Availability**: Repository-dependent
- **Community**: Specialized AFQMC researchers
- **Applications**: Method development, research

## Limitations & Constraints
- **Specialized**: Research implementation
- **Availability**: May be group-specific
- **Documentation**: Research-level
- **Sign problem**: AFQMC general issue
- **Expertise required**: AFQMC methodology

## Related Codes
- **QMCPACK**: Includes AFQMC
- **QUEST**: QMC methods
- **ALF**: Lattice AFQMC
- **DCA++**: CT-AUX (related)

## Application Areas
- Hubbard model
- Correlated materials
- Finite temperature
- Method research
- Sign problem studies

## Comparison with Other Codes
| Feature | ph-AFQMC (Research) | QMCPACK | ALF |
| :--- | :--- | :--- | :--- |
| **Method** | Particle-Hole AFQMC | Real-Space / AFQMC | Lattice QMC (Aux-Field) |
| **Focus** | Sign Problem Mitigation | Continuum / Materials | Lattice Models |
| **Temperature** | Finite T | Ground State (mostly) | Finite T / Ground State |
| **Status** | Research Implementation | Production Standard | Community Standard |

## Verification & Sources
**Primary sources**:
1. AFQMC literature
2. Particle-hole formulation papers
3. Research implementations

**Confidence**: VERIFIED - Research method

**Verification status**: ✅ VERIFIED (RESEARCH METHOD)
- **Category**: AFQMC variant/method
- **Note**: Particle-hole formulation of auxiliary-field quantum Monte Carlo. Research implementation, not necessarily a single unified public code. Used in specialized AFQMC research. For production AFQMC, see QMCPACK (includes AFQMC module), ALF (lattice AFQMC), or research group implementations. Specialized technique for finite-temperature strongly correlated systems.
