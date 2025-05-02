# Database Verification Report Template

## Overview

This report provides a comprehensive verification of the database population status, including molecule counts, reference compounds, property completeness, cross-references, and query performance.

## Database Information

- **Adapter Type**: [Adapter Type]
- **Connection Status**: [Connection Status]

## Test Results

### Molecule Counts

Status: **[Status]**

- Total molecules: [Count]
- PubChem molecules: [Count]
- ChEMBL molecules: [Count]
- Cross-referenced molecules: [Count]

### Reference Compounds

Status: **[Status]**

- Found: [Count] / [Expected Count]
- Missing: [List of missing compounds if any]

### Property Completeness

Status: **[Status]**

- Total molecules checked: [Count]
- Complete properties: [Count] ([Percentage]%)
- With logP: [Count] ([Percentage]%)
- With H-bond donors: [Count] ([Percentage]%)
- With H-bond acceptors: [Count] ([Percentage]%)

#### Molecules with Missing Properties (Sample)

| ID | ChEMBL ID | Missing Properties |
|----|-----------|-------------------|
| [ID] | [ChEMBL ID] | [Properties] |
| [ID] | [ChEMBL ID] | [Properties] |
| [ID] | [ChEMBL ID] | [Properties] |

### Cross References

Status: **[Status]**

- Cross-referenced molecules: [Count]

### Query Performance

Status: **[Status]**

- Query time: [Time]ms
- Result count: [Count]

## Detailed Diagnostics

### Molecule Sources

| Source | Count |
|--------|-------|
| [Source] | [Count] |
| [Source] | [Count] |

### Recently Added Molecules

| ID | Name | ChEMBL ID | PubChem CID | Created At |
|----|------|-----------|-------------|------------|
| [ID] | [Name] | [ChEMBL ID] | [PubChem CID] | [Timestamp] |
| [ID] | [Name] | [ChEMBL ID] | [PubChem CID] | [Timestamp] |

### Cross-Referenced Molecule Samples

| ID | Name | ChEMBL ID | PubChem CID |
|----|------|-----------|-------------|
| [ID] | [Name] | [ChEMBL ID] | [PubChem CID] |
| [ID] | [Name] | [ChEMBL ID] | [PubChem CID] |

### Database Indexes

| Table | Index | Column |
|-------|-------|--------|
| [Table] | [Index] | [Column] |
| [Table] | [Index] | [Column] |

## Execution Information

- Total execution time: [Time] seconds
- Timestamp: [Timestamp]