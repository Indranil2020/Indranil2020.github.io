# TCSP

## Overview
TCSP (Template-based Crystal Structure Prediction) is a fast and accurate algorithm for crystal structure prediction using template structures from known materials databases.

## Theoretical Basis
- Template-based prediction
- Element substitution
- Structure similarity matching
- Database-driven approach
- Fast screening without DFT

## Key Capabilities
- Template-based CSP
- Web server available
- Fast prediction
- Wide chemical space
- No DFT required

**Sources**: Inorg. Chem. 61, 8431 (2022), arXiv:2111.14049

## Key Strengths

### Speed:
- No DFT calculations
- Fast template matching
- Large-scale screening

### Accessibility:
- Web server available
- Easy to use
- No installation needed

### Coverage:
- Large template database
- Multiple crystal systems
- Wide applicability

## Inputs & Outputs
- **Input formats**: Chemical composition
- **Output data types**: Predicted structures, template matches

## Interfaces & Ecosystem
- **Web server**: Online interface
- **Databases**: Template database
- **Validation**: DFT codes

## Workflow and Usage
1. Input chemical composition
2. Search template database
3. Match and substitute elements
4. Generate candidate structures
5. Validate with DFT (optional)

## Performance Characteristics
- Very fast prediction
- Database-limited
- Good for known structure types

## Computational Cost
- Minimal (database search)
- No DFT required
- Fast screening

## Best Practices
- Use for initial screening
- Validate top predictions
- Consider multiple templates
- Check structural validity

## Limitations & Known Constraints
- Template-dependent
- May miss novel structures
- Limited to known structure types

## Application Areas
- High-throughput screening
- Initial structure guessing
- Materials discovery
- Database expansion

## Comparison with Other Codes
- **vs CSPML**: Similar approach, different implementation
- **vs USPEX**: TCSP faster, less general
- **Unique strength**: Web server, template database, fast

## Community and Support
- Web server available
- Academic development
- Published methodology

## Verification & Sources
**Primary sources**:
1. Publication: Inorg. Chem. 61, 8431 (2022)
2. arXiv: https://arxiv.org/abs/2111.14049
3. TCSP 2.0: arXiv:2503.23183

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (web server)
- Documentation: AVAILABLE
- Development: ACTIVE
- Applications: Template-based CSP
