# CryoProtect v2 Data Integrity Framework

## Overview

This document establishes a comprehensive data integrity framework for CryoProtect v2, with specific focus on the ChEMBL integration but extending to all data sources. The framework addresses validation, verification, and ongoing monitoring of data quality throughout the system.

## Core Principles

1. **Source Truth**: Every data point must be traceable to an authoritative source
2. **Multi-Source Validation**: Critical data should be validated against multiple sources when possible
3. **Computational Verification**: Property values should be verified by computational methods where applicable
4. **Continuous Monitoring**: Data quality must be continuously monitored, not just during import
5. **Transparent Reporting**: Data quality metrics must be visible to users and administrators
6. **Automated Remediation**: Common data issues should have automated remediation paths
7. **Manual Review Process**: Complex data issues must have a defined manual review workflow

## Data Quality Dimensions

### 1. Accuracy

**Definition**: The degree to which data correctly reflects the real-world entity or event it represents.

**Measurement Approaches**:
- Compare property values against authoritative sources (ChEMBL, PubChem)
- Cross-reference with computational methods (RDKit)
- Validate against experimental results when available

**Key Metrics**:
- Percentage of properties within acceptable tolerance of authoritative values
- Number of conflicts between different authoritative sources
- Root mean square error between predicted and experimental values

### 2. Completeness

**Definition**: The degree to which all required data is present and in usable form.

**Measurement Approaches**:
- Track presence of required properties for each molecule
- Monitor coverage of optional properties
- Assess completeness of related data (mixtures, experiments)

**Key Metrics**:
- Percentage of molecules with complete required properties
- Coverage percentage for optional properties
- Completeness trend over time

### 3. Consistency

**Definition**: The degree to which data is consistent within the database and across the system.

**Measurement Approaches**:
- Validate referential integrity across related tables
- Check for logical consistency in related properties
- Ensure property values conform to chemical principles

**Key Metrics**:
- Number of referential integrity violations
- Percentage of properties with logical inconsistencies
- Violation of chemical principles (e.g., impossible LogP values)

### 4. Timeliness

**Definition**: The degree to which data is up-to-date and reflects current knowledge.

**Measurement Approaches**:
- Track data source version and update dates
- Monitor refresh frequency
- Identify outdated property values

**Key Metrics**:
- Average age of property values
- Percentage of data updated within target timeframe
- Number of properties superseded by newer information

### 5. Uniqueness

**Definition**: The degree to which entities are represented once in the database without duplication.

**Measurement Approaches**:
- Check for duplicate molecule records
- Identify redundant property values
- Validate uniqueness constraints

**Key Metrics**:
- Number of duplicate molecule records
- Percentage of redundant property values
- Uniqueness constraint violations

## ChEMBL-Specific Data Integrity Measures

### Property Validation Strategy

#### LogP Validation
LogP values require special attention due to previous discrepancies:

1. **Multi-source Validation**:
   - Compare ChEMBL LogP values with PubChem values
   - Calculate LogP using RDKit computational methods
   - Compare with experimental LogP values where available

2. **Tolerance Definitions**:
   - Primary tolerance: ±0.1 between authoritative sources
   - Secondary tolerance: ±0.3 between calculated and authoritative values
   - Experimental tolerance: ±0.5 between experimental and authoritative values

3. **Resolution Approach**:
   - Within primary tolerance: Accept ChEMBL value
   - Outside primary tolerance: Flag for review and use weighted average
   - Outside secondary tolerance: Manual review required

#### Molecular Structure Validation

1. **InChI/InChIKey Consistency**:
   - Validate InChIKey derivation from InChI
   - Cross-reference InChIKey between ChEMBL and PubChem
   - Verify structure consistency using RDKit

2. **SMILES Canonicalization**:
   - Convert all SMILES to canonical form
   - Verify SMILES against InChI representation
   - Ensure SMILES generates consistent molecular graphs

3. **Structure-Property Consistency**:
   - Validate molecular weight against structure
   - Check for consistency between structure and reported properties
   - Verify atom counts and compositions

### Reference Compound Validation

For critical reference compounds (e.g., CHEMBL25, CHEMBL1, etc.):

1. **Enhanced Validation Process**:
   - Cross-reference against multiple authoritative sources
   - Manual verification of all properties
   - Computational verification of all calculable properties

2. **Property Completeness Requirements**:
   - 100% completion of all standard properties
   - Experimental values for critical properties
   - Complete citation information

3. **Versioning and Audit Trail**:
   - Track all changes to reference compounds
   - Document validation history
   - Maintain persistent identifiers

## Implementation: Multi-tiered Validation System

### Tier 1: Import-time Validation

```python
def validate_molecule_import(molecule_data, authoritative_sources=None):
    """
    Validate molecule data during import process.
    
    Args:
        molecule_data: Dictionary of molecule properties
        authoritative_sources: Optional dictionary of values from other sources
        
    Returns:
        Dictionary with validation results and flags
    """
    validation = {
        "passed": True,
        "warnings": [],
        "errors": [],
        "property_validations": {}
    }
    
    # Validate structural identifiers
    if not validate_structure_consistency(molecule_data):
        validation["passed"] = False
        validation["errors"].append("Structural identifier inconsistency")
    
    # Validate individual properties
    for prop, value in molecule_data["properties"].items():
        prop_validation = validate_property(
            prop, value, molecule_data, authoritative_sources
        )
        validation["property_validations"][prop] = prop_validation
        
        if prop_validation["severity"] == "error":
            validation["passed"] = False
            validation["errors"].append(f"{prop}: {prop_validation['message']}")
        elif prop_validation["severity"] == "warning":
            validation["warnings"].append(f"{prop}: {prop_validation['message']}")
    
    # Check completeness
    completeness = check_completeness(molecule_data)
    if completeness < 0.8:  # 80% threshold
        validation["warnings"].append(f"Low completeness score: {completeness}")
    
    return validation
```

### Tier 2: Post-import Verification

```sql
-- Example SQL verification query for LogP consistency
WITH property_sources AS (
    SELECT 
        m.id AS molecule_id,
        m.name,
        m.chembl_id,
        mp_chembl.numeric_value AS chembl_logp,
        mp_pubchem.numeric_value AS pubchem_logp,
        mp_rdkit.numeric_value AS calculated_logp
    FROM 
        molecules m
    LEFT JOIN 
        molecular_properties mp_chembl ON m.id = mp_chembl.molecule_id 
        AND mp_chembl.property_type_id = (SELECT id FROM property_types WHERE name = 'LogP' AND data_source LIKE '%ChEMBL%')
    LEFT JOIN 
        molecular_properties mp_pubchem ON m.id = mp_pubchem.molecule_id 
        AND mp_pubchem.property_type_id = (SELECT id FROM property_types WHERE name = 'LogP' AND data_source LIKE '%PubChem%')
    LEFT JOIN 
        molecular_properties mp_rdkit ON m.id = mp_rdkit.molecule_id 
        AND mp_rdkit.property_type_id = (SELECT id FROM property_types WHERE name = 'LogP' AND data_source LIKE '%RDKit%')
    WHERE 
        m.chembl_id IS NOT NULL
)
SELECT 
    molecule_id,
    name,
    chembl_id,
    chembl_logp,
    pubchem_logp,
    calculated_logp,
    ABS(COALESCE(chembl_logp, 0) - COALESCE(pubchem_logp, 0)) AS chembl_pubchem_diff,
    ABS(COALESCE(chembl_logp, 0) - COALESCE(calculated_logp, 0)) AS chembl_calculated_diff,
    CASE 
        WHEN ABS(COALESCE(chembl_logp, 0) - COALESCE(pubchem_logp, 0)) > 0.1 THEN 'Primary Tolerance Exceeded'
        WHEN ABS(COALESCE(chembl_logp, 0) - COALESCE(calculated_logp, 0)) > 0.3 THEN 'Secondary Tolerance Exceeded'
        ELSE 'Within Tolerance'
    END AS validation_status
FROM 
    property_sources
WHERE 
    ABS(COALESCE(chembl_logp, 0) - COALESCE(pubchem_logp, 0)) > 0.1
    OR ABS(COALESCE(chembl_logp, 0) - COALESCE(calculated_logp, 0)) > 0.3
ORDER BY 
    chembl_pubchem_diff DESC;
```

### Tier 3: Continuous Monitoring

```python
def run_data_integrity_monitoring():
    """
    Scheduled job to continuously monitor data integrity.
    """
    # Check for new issues
    new_issues = identify_new_data_issues()
    
    # Check previously identified issues for changes
    resolved_issues = check_issue_resolutions()
    
    # Update data quality metrics
    update_data_quality_dashboard()
    
    # Send notifications if necessary
    if critical_issues_present(new_issues):
        send_alert("Critical data integrity issues detected")
    
    # Generate monitoring report
    generate_monitoring_report(new_issues, resolved_issues)
    
    # Schedule remediation for automated fixes
    schedule_automated_remediation()
    
    # Update issue tracking system
    update_issue_tracker(new_issues, resolved_issues)
    
    return {
        "run_timestamp": datetime.now().isoformat(),
        "new_issues": len(new_issues),
        "resolved_issues": len(resolved_issues),
        "data_quality_score": calculate_overall_quality_score()
    }
```

## Data Remediation Workflows

### Automated Remediation

For issues that can be automatically resolved:

1. **Property Normalization**:
   - Automatically convert units to standard form
   - Apply standard precision rounding
   - Format string values consistently

2. **Missing Value Imputation**:
   - Calculate missing values where possible
   - Use predictive models for estimation
   - Apply defaults for non-critical properties

3. **Duplicate Resolution**:
   - Apply merge rules for duplicate molecules
   - Consolidate properties from multiple sources
   - Maintain provenance information

### Manual Review Process

For issues requiring human judgment:

1. **Issue Classification**:
   - Categorize by severity and impact
   - Assign priority based on usage patterns
   - Group related issues

2. **Review Assignment**:
   - Route to appropriate subject matter expert
   - Provide context and recommendations
   - Set resolution timeframes

3. **Resolution Workflow**:
   - Document decision rationale
   - Apply changes with audit trail
   - Verify fix effectiveness
   - Update knowledge base

## Data Quality Reporting

### Dashboard Metrics

1. **Overall Quality Score**: Composite metric of all quality dimensions
2. **Dimension Scores**: Individual scores for accuracy, completeness, etc.
3. **Issue Trends**: Number of issues over time by category
4. **Resolution Rate**: Percentage of issues resolved within target timeframe
5. **Source Reliability**: Quality score by data source

### User-Facing Quality Indicators

1. **Property Confidence Indicators**:
   - Visual indicators of property value confidence
   - Tooltips explaining validation sources
   - Links to authoritative sources

2. **Data Source Attribution**:
   - Clear indication of data provenance
   - Multiple source acknowledgment
   - Last verification timestamp

3. **Completeness Indicators**:
   - Visual representation of profile completeness
   - Missing property indicators
   - Suggestions for additional data sources

## Integration with Development Workflow

### Development Guidelines

1. **New Property Addition**:
   - Validation requirements definition
   - Authoritative source identification
   - Tolerance specification
   - Verification query implementation

2. **Data Source Integration**:
   - Source reliability assessment
   - Conflict resolution strategy
   - Import validation implementation
   - Cross-reference mapping

3. **Schema Evolution**:
   - Impact assessment on validation processes
   - Migration validation procedures
   - Backward compatibility verification
   - Data quality regression testing

### Testing Requirements

1. **Validation Logic Tests**:
   - Unit tests for all validation functions
   - Test cases covering edge conditions
   - Performance testing with large datasets

2. **Verification Query Tests**:
   - Correctness verification
   - Performance optimization
   - Result reproducibility

3. **End-to-End Data Quality Tests**:
   - Import-verification integration
   - Full pipeline validation
   - Cross-source consistency checks

## Summary: Implementation Roadmap

### Phase 1: Foundation (Week 1-2)
- Implement basic validation for critical properties (LogP, MW)
- Create initial verification queries
- Set up issue tracking system

### Phase 2: Comprehensive Validation (Week 3-4)
- Extend validation to all property types
- Implement multi-source validation
- Create automated test suite

### Phase 3: Monitoring System (Week 5-6)
- Implement continuous monitoring
- Create data quality dashboard
- Set up alerting system

### Phase 4: Remediation Workflows (Week 7-8)
- Implement automated remediation for common issues
- Create manual review workflow
- Develop knowledge base for recurring issues

### Phase 5: User Integration (Week 9-10)
- Implement user-facing quality indicators
- Create data quality API
- Develop user documentation

## Conclusion

This data integrity framework provides a comprehensive approach to ensuring the highest quality data in the CryoProtect v2 system. By implementing multi-tiered validation, continuous monitoring, and effective remediation workflows, we can maintain data integrity throughout the lifecycle of the system, even as the database scales to include thousands of compounds and properties.

The framework pays special attention to the issues encountered with LogP values and other critical properties while establishing a generalizable approach that can be applied to all current and future data in the system.