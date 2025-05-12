#!/bin/bash
# Script to create GitHub issues for ChEMBL import improvements

# Ensure the GitHub CLI is installed and authenticated
if ! command -v gh &> /dev/null; then
    echo "GitHub CLI not found. Please install it first."
    exit 1
fi

# Check if we're in a GitHub repository
if ! gh repo view &> /dev/null; then
    echo "Not in a GitHub repository or gh not authenticated."
    exit 1
fi

# Function to create an issue
create_issue() {
    local title="$1"
    local body="$2"
    local labels="$3"
    
    echo "Creating issue: $title"
    gh issue create --title "$title" --body "$body" --label "$labels"
}

# Issue 1: Fix molecular property calculation
create_issue \
    "Fix molecular property calculation for ChEMBL imported molecules" \
    "## Background
Based on our recent verification of the ChEMBL import data, we discovered that while 95.2% of ChEMBL molecules have some properties calculated, all of them are missing one or more critical properties.

## Requirements
- Fix the property calculation component to ensure all standard properties are computed:
  - LogP
  - TPSA
  - Molecular Weight
  - Heavy Atom Count
  - H-Bond Donor Count
  - H-Bond Acceptor Count
  - Rotatable Bond Count
  - Ring Count
- Implement a reprocessing script to calculate missing properties for existing molecules
- Add logging to track property calculation success/failure

## Technical Details
- The verification report is available at \`reports/chembl_data_verification_*.json\`
- The molecular_properties table includes partial data for most molecules
- Property values may be null in some cases even when records exist

## Acceptance Criteria
- All ChEMBL molecules should have complete property data
- New imports should automatically calculate all properties
- A repair script should fix existing data
- Verification script should confirm 100% property coverage" \
    "bug,enhancement,chembl,molecular-properties"

# Issue 2: Implement ChEMBL-to-PubChem cross-reference resolution
create_issue \
    "Implement ChEMBL-to-PubChem cross-reference resolution" \
    "## Background
Our verification found that only 4.8% of ChEMBL molecules have corresponding PubChem CIDs. This missing cross-reference significantly limits our ability to combine data from both sources.

## Requirements
- Implement a robust ChEMBL-to-PubChem identifier mapping mechanism
- Populate missing PubChem CIDs for existing ChEMBL molecules
- Ensure new ChEMBL imports automatically resolve PubChem identifiers

## Technical Approaches
1. Use the PubChem PUG-REST API to lookup ChEMBL IDs
2. Create a local mapping table from publicly available cross-reference data
3. Use InChIKey as a common identifier for matching

## Acceptance Criteria
- At least 90% of ChEMBL molecules should have corresponding PubChem CIDs
- The cross-reference process should be automatically applied during import
- Implement fallback mechanisms for when direct mapping fails
- Include detailed logging of the resolution process" \
    "enhancement,chembl,pubchem,data-integration"

# Issue 3: Create data quality validation checks for molecule imports
create_issue \
    "Create data quality validation checks for molecule imports" \
    "## Background
Our verification of imported chemical data revealed quality issues including missing properties and cross-references. We need systematic validation during the import process.

## Requirements
- Implement a comprehensive validation framework for molecule imports that checks:
  - Structural integrity (valid SMILES, InChI)
  - Property calculation completeness
  - Cross-reference resolution success
  - Duplicate detection
- Create informative error messages for failed validations
- Generate data quality reports after batch imports

## Technical Specifications
- Create a validation module that can be integrated into all import pipelines
- Define validation levels: critical (blocking) and warning (non-blocking)
- Implement retry mechanisms for transient failures
- Store validation results for auditing

## Acceptance Criteria
- All new molecule imports pass critical validation checks
- Failed imports are properly logged with detailed error messages
- Batch imports include a data quality summary
- Validation framework can be extended for new data sources" \
    "enhancement,data-quality,chembl,pubchem"

# Issue 4: Optimize molecular property storage using JSONB
create_issue \
    "Optimize molecular property storage using JSONB fields" \
    "## Background
The current database schema uses separate records in the molecular_properties table for each property, which increases query complexity and reduces performance. The molecules table already has a JSONB 'properties' field that's currently unused.

## Requirements
- Migrate and normalize property data to use the JSONB 'properties' field
- Create database functions to efficiently query and filter by properties in JSONB
- Implement indexing strategies for JSONB property queries
- Update import pipelines to store properties in the JSONB format

## Technical Specifications
- Define a standard schema for the JSONB properties structure
- Create migration scripts to convert existing property data
- Add GIN or BTREE indexes for common property queries
- Update application code to use the new property access patterns

## Acceptance Criteria
- Property queries show measurable performance improvement
- All existing property data is preserved during migration
- New imports automatically use the JSONB storage format
- Documentation is updated to reflect the new approach" \
    "enhancement,performance,database,molecular-properties"

echo "Issues created successfully!"