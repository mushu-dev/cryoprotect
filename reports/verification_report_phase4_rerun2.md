# Database Population Verification Report

Generated: 2025-05-01T17:39:38.800030

## Overall Status: FAIL

## Summary

| Metric | Value | Requirement | Status |
|--------|-------|-------------|--------|
| Total molecules | 786 | ≥ 5,000 | FAIL |
| Reference compounds | 0/9 | All complete | FAIL |
| Property completeness | 0.0% | ≥ 90% | FAIL |
| Average query time | 414.37 ms | < 50 ms | FAIL |

## Molecule Counts

- **Total molecules:** 786
- **With PubChem CID:** 693
- **With ChEMBL ID:** 73
- **With cross-references:** 0
- **With properties:** 519

## Reference Compounds

- **Total reference compounds:** 9
- **Found reference compounds:** 9
- **Complete reference compounds:** 0
- **Incomplete reference compounds:** 9
- **Missing reference compounds:** 0

### Incomplete Reference Compounds

| ChEMBL ID | Name | Missing Properties |
|-----------|------|-------------------|
| CHEMBL1098659 | CHEMBL1098659 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL1487 | CHEMBL1487 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL262548 | CHEMBL262548 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL388978 | CHEMBL388978 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL500033 | CHEMBL500033 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL6196 | CHEMBL6196 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL66195 | CHEMBL66195 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL6752 | CHEMBL6752 | logP, h_bond_donors, h_bond_acceptors |
| CHEMBL967 | CHEMBL967 | logP, h_bond_donors, h_bond_acceptors |

## Property Completeness

- **Total molecules with properties:** 0
- **Molecules with complete properties:** 0
- **Molecules with incomplete properties:** 0
- **Property completeness percentage:** 0.0%

## Query Performance

- **Overall average query time:** 414.37 ms
- **Performance acceptable:** No

### Individual Query Performance

| Query | Average (ms) | Min (ms) | Max (ms) |
|-------|-------------|----------|----------|
| Fetch single molecule by ID | 404.16 | 381.14 | 432.95 |
| Fetch molecule with properties | 445.26 | 397.84 | 544.16 |
| Search molecules by name | 399.73 | 365.18 | 422.38 |
| Count molecules by property value range | 408.31 | 373.54 | 449.03 |

