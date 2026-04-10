# ESI-3D

## Official Resources
- Resource Listing: https://quantchemdev.github.io/resources.html
- Historical Announcement and Manual Links: https://server.ccl.net//chemistry/resources/messages/2006/06/06.007-dir/
- Maintainer Contact: ematito@dipc.org

## Overview
ESI-3D is a program for calculating electron sharing and delocalization descriptors used in chemical bonding and aromaticity analysis. It evaluates two-center and multicenter electron sharing indices and a range of aromaticity indicators from atomic overlap matrices derived from Hilbert-space or real-space partitions such as QTAIM.

**Scientific domain**: Electron sharing indices, delocalization analysis, aromaticity  
**Target user community**: Theoretical chemists studying bonding delocalization and aromaticity using QTAIM or Hilbert-space partitions

## Theoretical Methods
- Electron sharing indices (ESI)
- Delocalization indices
- Multicenter electron sharing analysis
- Aromaticity indicators including FLU, PDI, HOMA, MCI, Iring, AV1245, and AVmin
- Use of atomic overlap matrices from Hilbert-space or real-space partitions

## Capabilities (CRITICAL)
- Calculates two-center and multicenter electron sharing indices
- Computes multiple aromaticity indicators from a single workflow
- Accepts atomic overlap matrices derived from QTAIM or Hilbert-space partitions
- Relevant to bonding interpretation in conjugated, aromatic, and multicenter systems
- Publicly described in group resources, with program access available via contact with the maintainer

**Sources**: Quantum Chemistry Development Group resources page and archived CCL program announcement

## Key Strengths

### Delocalization-Focused Analysis:
- Directly targets electron sharing and delocalization
- Useful complement to QTAIM critical-point analysis
- Covers both bonding and aromaticity descriptors

### Rich Aromaticity Toolkit:
- FLU and PDI
- HOMA and MCI
- Iring and related multicenter measures
- Suitable for comparative aromaticity studies

### Methodological Flexibility:
- Compatible with Hilbert-space and real-space partitioning frameworks
- Relevant to users working with QTAIM-derived quantities
- Long-standing specialized code in the electron-delocalization community

## Inputs & Outputs
- **Input formats**:
  - Atomic overlap matrices from Hilbert-space partitions
  - Atomic overlap matrices derived from real-space partitions such as QTAIM

- **Output data types**:
  - Electron sharing indices
  - Delocalization measures
  - Aromaticity indicators

## Workflow and Usage
1. Generate the required atomic overlap matrices from the chosen partition scheme.
2. Run ESI-3D on the prepared data.
3. Extract two-center, multicenter, and aromaticity descriptors.
4. Compare the resulting delocalization metrics across molecules or fragments.

## Performance Characteristics
- Specialized descriptor workflow for delocalization and aromaticity analysis
- Best suited to targeted bonding studies rather than general-purpose electronic-structure post-processing
- Longstanding niche tool with focused methodology

## Limitations & Known Constraints
- **Distribution model**: Current version is obtained through maintainer contact rather than a mainstream package manager
- **Specialized inputs**: Requires precomputed atomic overlap matrices from compatible partitioning schemes
- **Niche scope**: Focused on delocalization and aromaticity rather than broad wavefunction analysis

## Comparison with Other Tools
- **vs EDDB**: ESI-3D focuses on electron sharing and aromaticity indices derived from overlap information, while EDDB emphasizes delocalized electron density descriptors
- **vs QTAIM suites**: ESI-3D complements topological analysis with quantitative delocalization metrics
- **Unique strength**: Broad set of electron sharing and aromaticity indices within a specialized framework

## Application Areas
- Aromaticity quantification
- Delocalization analysis
- Multicenter bonding studies
- Comparative bonding descriptors for conjugated systems

## Community and Support
- Listed in the Quantum Chemistry Development Group resources
- Maintainer contact provided for current versions
- Historical public announcement and manual links preserved in CCL archives

## Verification & Sources
**Primary sources**:
1. Resources page: https://quantchemdev.github.io/resources.html
2. Archived program announcement: https://server.ccl.net//chemistry/resources/messages/2006/06/06.007-dir/
3. Maintainer contact listed in the resources page

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Public description: ACCESSIBLE
- Historical program links: AVAILABLE via archive
- Current access path: CONTACT WITH MAINTAINER
- Primary use case: Electron sharing and aromaticity analysis
