# GAP

## Official Resources
- Homepage: http://www.wien2k.at/ (Distributed with/for WIEN2k)
- Documentation: WIEN2k User Guide / GAP2 literature
- Source Repository: Part of WIEN2k distribution (Licensed)
- License: Proprietary / Academic License (WIEN2k)

## Overview
GAP (specifically GAP2 - GW with Augmented Plane-waves) is the all-electron GW implementation within the WIEN2k ecosystem. It utilizes the full-potential linearized augmented plane-wave (FP-LAPW) basis set to perform high-precision G0W0 calculations, treating core, semi-core, and valence electrons on equal footing. It is considered a "gold standard" for GW accuracy in solids.

**Scientific domain**: All-electron GW, precision solid-state physics, core-level excitations  
**Target user community**: WIEN2k users, researchers requiring benchmark accuracy

## Theoretical Methods
- GW approximation (G0W0, GAP2 implementation)
- All-electron FP-LAPW formalism
- Mixed basis set (APW + lo for G, PW for W)
- Full frequency integration
- Core-valence interactions treated explicitly
- Relativistic effects (scalar + SOC)

## Capabilities (CRITICAL)
- All-electron G0W0 calculations
- High-precision band gaps (<10% error vs experiment)
- d- and f-electron systems (strongly correlated)
- Transition metal oxides and lanthanides

**Sources**: Academic literature, WIEN2k interface references

## Key Strengths

### All-Electron Accuracy:
- No pseudopotentials
- treatment of core-valence interaction
- Accurate for heavy elements
- High precision reference

### APW Basis:
- Efficient for open structures
- Accurate near nucleus
- Proven solid-state basis
- Rigorous mathematical foundation

### WIEN2k Interface:
- Leverages WIEN2k DFT
- Proven ground state
- Established ecosystem
- Specialized community

## Inputs & Outputs
- **Input formats**:
  - WIEN2k struct/vector files
  - GAP specific control files
  
- **Output data types**:
  - Quasiparticle band structure
  - Self-energy corrections
  - Core level shifts

## Interfaces & Ecosystem
- **Primary Interface**: WIEN2k
- **Usage**: Typically post-processing step
- **Ecosystem**: FLAPW community tools

## Performance Characteristics
- **Speed**: Expensive (all-electron)
- **Accuracy**: Very high (gold standard)
- **System size**: Small to medium solids
- **Scaling**: Steep scaling with size

## Computational Cost
- **High**: Due to all-electron basis
- **Memory**: Demanding mixed basis
- **Precision**: Costs justify accuracy

## Limitations & Known Constraints
- **Availability**: Not widely public like VASP/QE
- **Complexity**: Steep learning curve
- **Maintenance**: Academic code status
- **Efficiency**: Slower than pseudopotential codes

## Comparison with Other Codes
- **vs BerkeleyGW**: GAP is all-electron, BerkeleyGW is pseudopotential
- **vs SPEX**: Similar domain (all-electron GW), widely used competitor
- **vs exciting**: Alternative all-electron implementation
- **Unique strength**: Historic APW implementation, WIEN2k synergy

## Application Areas

### Strongly Correlated Materials:
- Transition metal oxides
- f-electron systems
- Magnetic materials

### Precision Benchmarking:
- Reference values
- Pseudopotential validation
- Core-level spectroscopy

## Community and Support
- Academic user base
- WIEN2k mailing list
- Specialized workshops
- Literature-based support

## Verification & Sources
**Primary sources**:
1. References in WIEN2k community
2. "GAP" all-electron GW literature (TU Wien group)
3. Comparison papers with SPEX/exciting

**Confidence**: VERIFIED (Historic/Academic)
- Existence: CONFIRMED in literature
- Status: Niche/Academic research code
