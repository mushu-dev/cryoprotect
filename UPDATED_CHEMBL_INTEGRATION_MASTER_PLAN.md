# ChEMBL Integration Master Plan for CryoProtect v2

## Executive Summary

This plan outlines our approach to integrating ChEMBL data into the CryoProtect v2 database. By leveraging our successful PubChem integration architecture and adapting it for ChEMBL's specific requirements, we'll populate our Supabase database with high-quality scientific data essential for cryoprotectant research.

## Key Components

1. **Purpose**: Populate the CryoProtect database with scientifically accurate cryoprotectant data from ChEMBL
2. **Implementation Approach**: MCP-based direct SQL execution for database operations
3. **Execution Method**: Delegated implementation through ROO code agents with clear task boundaries
4. **Success Criteria**: 1,000+ molecules, 10,000+ properties, verified scientific accuracy

## Implementation Strategy

### Phase 1: Foundation (ROO Agent Task)

**Task 1.1: Implement Robust ChEMBL Client Wrapper**
- File: `/chembl/client.py`
- Purpose: Create a resilient wrapper around the official ChEMBL client
- Implementation Details:
  - Sunday/Monday rate limiting adjustments
  - Multi-level caching (memory and disk-based)
  - Circuit breaker pattern for failure protection
  - Transparent fallback to cached data

**Task 1.2: Set Up Cache Management System**
- File: `/chembl/cache.py`
- Purpose: Optimize API usage with proper caching
- Implementation Details:
  - TTL-based cache invalidation
  - JSON serialization for property storage
  - Memory cache with LRU eviction policy

**Task 1.3: Configure Logging System**
- File: `/chembl/logging.py`
- Purpose: Standardize logging across ChEMBL operations
- Implementation Details:
  - Structured JSON logging
  - Progress tracking
  - Error categorization
  - Summary report generation

### Phase 2: Data Pipeline (ROO Agent Task)

**Task 2.1: Implement Data Acquisition**
- File: `/ChEMBL_Integrated_Import.py`
- Purpose: Fetch relevant cryoprotectant data from ChEMBL
- Implementation Details:
  - Retrieve standard reference compounds (glycerol, DMSO, etc.)
  - Execute search for potential cryoprotectants
  - Implement checkpointing for resumable operations
  - Batch processing with progress tracking

**Task 2.2: Implement Data Transformation**
- File: `/ChEMBL_Integrated_Import.py`
- Purpose: Convert ChEMBL data to CryoProtect schema
- Implementation Details:
  - Transform molecule core data (smiles, inchi, formula, etc.)
  - Extract and normalize properties
  - Map property types to CryoProtect schema
  - Fallback handling for missing property types

**Task 2.3: Implement Database Operations**
- File: `/ChEMBL_Integrated_Import.py`
- Purpose: Insert data into Supabase using MCP
- Implementation Details:
  - Batch insertion with transaction support
  - SQL injection protection
  - Conflict resolution for duplicate molecules
  - Proper error handling and reporting

### Phase 3: Verification & Optimization (ROO Agent Task)

**Task 3.1: Implement Verification System**
- File: `/verify_chembl_data.py`
- Purpose: Verify imported data for accuracy and completeness
- Implementation Details:
  - Database count verification
  - Property distribution analysis
  - Sample molecule validation
  - Data quality reports

**Task 3.2: Implement Property Reconciliation**
- File: `/reconcile_chembl_properties.py`
- Purpose: Ensure consistency between ChEMBL and existing properties
- Implementation Details:
  - Identify equivalent properties
  - Standardize property names and units
  - Merge duplicate property types
  - Update molecule scores based on new properties

**Task 3.3: Implement Performance Optimizations**
- Purpose: Optimize database performance for the increased data volume
- Implementation Details:
  - Add appropriate indexes
  - Optimize query plans
  - Configure connection pooling
  - Update RLS policies for performance

## Integration with PubChem Data

1. **Shared Code Reuse**
   - Reuse progress tracking dashboard from PubChem implementation
   - Adapt batch processing framework from PubChem
   - Leverage existing SQL insertion templates

2. **Data Consolidation**
   - Ensure no duplicates between PubChem and ChEMBL data
   - Standardize property naming and units across both sources
   - Add source tracking for all molecules (PubChem vs ChEMBL)

3. **Combined Scoring System**
   - Update scoring algorithm to work consistently across both data sources
   - Reconcile property weights for uniform evaluation
   - Preserve source information for traceability

## Execution Plan

1. **Preparation (Day 1)**
   - Task: Review existing ChEMBL code and identify reusable components
   - Agent Type: PM/Design 
   - Output: Technical specification with file/line references

2. **Foundation Implementation (Day 1-2)**
   - Task: Implement core ChEMBL client wrapper
   - Agent Type: RooCode
   - Output: Functional client with caching and rate limiting

3. **Data Pipeline Implementation (Day 2-3)**
   - Task: Implement data acquisition, transformation, and storage
   - Agent Type: RooCode
   - Output: Working pipeline with checkpoint system

4. **Verification System (Day 3-4)**
   - Task: Implement verification and reconciliation
   - Agent Type: RooCode
   - Output: Validation scripts and repair utilities

5. **Testing & Optimization (Day 4-5)**
   - Task: Test and optimize the complete system
   - Agent Type: RooCode
   - Output: Performance report and optimization patches

## ROO Agent Guidelines

### Implementation Principles

1. **Focus on Implementation Only**
   - Implement functionality exactly as specified
   - Do not attempt to redesign architecture
   - Defer complex decisions back to PM

2. **Strict File Boundaries**
   - Modify only assigned files
   - Do not create new files without approval
   - Follow proper import structure

3. **Follow Established Patterns**
   - Use consistent error handling
   - Follow existing code style
   - Maintain coherent logging approach

### ROO Agent Task Structure

For each ROO agent task:

1. **Input**
   - Specific file to modify
   - Line number ranges for changes
   - References to example code

2. **Expected Output**
   - Modified file(s) implementing the specified functionality
   - Test verification results
   - No unnecessary explanations or designs

3. **Communication Protocol**
   - Report only implementation status
   - Flag specific technical issues
   - Request clarification on specific implementation details only

## Technical Resources

- **ChEMBL API Documentation**: https://chembl.gitbook.io/chembl-interface-documentation/web-services
- **Official Python Client**: https://github.com/chembl/chembl_webresource_client
- **Reference Implementation**: `/PubChem_CryoProtectants_Supabase_Enhanced_MCP.py`

## Success Metrics

1. **Data Volume**
   - At least 1,000 molecules loaded from ChEMBL
   - At least 10,000 molecular properties
   - Complete property coverage for reference compounds

2. **Performance**
   - < 500ms query time for molecule retrieval
   - < 1s for complex property queries
   - < 30 minutes for full import process

3. **Reliability**
   - 99% success rate for ChEMBL API requests
   - Zero data corruption events
   - Resumable operation after interruptions

## Appendix: Implementation Example

Example of delegating implementation to ROO agent:

```
Task: Implement transform_chembl_to_molecule function
File: /ChEMBL_Integrated_Import.py
Lines: 570-615
Reference: See PubChem_CryoProtectants_Supabase_Enhanced_MCP.py lines 766-834

The function should transform ChEMBL molecule data to match our database schema:
- Extract molecule name, formulas, identifiers
- Handle missing values gracefully
- Include data source attribution
- Support user tracking

Implementation must match the function signature and behavior in the reference file.
```

This structured approach ensures ROO agents can focus entirely on implementation without needing to spend tokens on design decisions or explanations.