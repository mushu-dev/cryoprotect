# ROO ChEMBL Remediation Task Breakdown

## Task Summary
The ChEMBL remediation process is currently blocked by three technical issues related to MCP tool limitations. This task requires implementing fixes to the orchestration script to address these issues and successfully complete the ChEMBL data remediation process.

## Technical Blockers Identified

| Blocker ID | Issue Description | Error | Status |
|------------|-------------------|-------|--------|
| #1 | MCP DDL Execution | `ALTER TABLE query does not return tuples` | Blocked-Debug |
| #2 | PostgreSQL Type Compatibility | JSON incompatible with `sql_identifier` type | Blocked-Debug |
| #3 | Script Argument Format | Incorrect argument names with underscores | Blocked-Debug |

## Task Breakdown

### Task 1: Fix MCP DDL Execution Issues

**Objective**: Modify the orchestration script to handle DDL statements through MCP

**Approach**:
1. Create the SQL function wrapper pattern for DDL operations
2. Split migration into individual statements
3. Handle each statement appropriately based on type
4. Implement proper error handling and verification

**Expected Output**:
- Updated `apply_migration` function that can successfully execute DDL operations via MCP
- Migration SQL with compatibility modifications

**Time Estimate**: 2-3 hours

### Task 2: Fix PostgreSQL Type Compatibility Issues

**Objective**: Ensure all database query results are JSON-compatible for MCP

**Approach**:
1. Identify all queries that return PostgreSQL-specific types
2. Update these queries to explicitly cast columns to text
3. Add appropriate column aliases to maintain field names
4. Verify query results work with MCP's JSON serialization

**Expected Output**:
- Updated verification queries that explicitly cast PostgreSQL types to text
- Clean JSON-compatible results for all MCP tool operations

**Time Estimate**: 1-2 hours

### Task 3: Fix Script Argument Format Issues

**Objective**: Update the ChEMBL import function to use correct argument naming convention

**Approach**:
1. Review current argument formatting in orchestration script
2. Update all command arguments to use hyphens instead of underscores
3. Verify argument compatibility with ChEMBL_Integrated_Import.py
4. Implement robust error handling for subprocess execution

**Expected Output**:
- Updated `run_import` function with correct argument format
- Clean execution of ChEMBL_Integrated_Import.py through subprocess

**Time Estimate**: 1 hour

### Task 4: Implement Consolidated Fix and Testing

**Objective**: Create a single fixed orchestration script and test the end-to-end process

**Approach**:
1. Consolidate all fixes into chembl_remediation_main_fixed.py
2. Test each component individually
3. Perform a full end-to-end test of the remediation process
4. Add proper logging and error reporting

**Expected Output**:
- Complete working implementation of chembl_remediation_main_fixed.py
- Verification that each step in the remediation process executes successfully
- Detailed log output for debugging and verification

**Time Estimate**: 2-3 hours

## Implementation Plan

1. **Code Review & Setup** (30 mins)
   - Review existing chembl_remediation_main.py
   - Create chembl_remediation_main_fixed.py
   - Set up test environment and logging

2. **MCP DDL Execution Fix** (2 hours)
   - Implement SQL wrapper functions for DDL
   - Modify migration SQL for compatibility
   - Test SQL execution with MCP

3. **Type Compatibility Fix** (1 hour)
   - Identify all type issues
   - Add explicit casts to all queries
   - Test verification queries

4. **Argument Format Fix** (1 hour)
   - Update subprocess command arguments
   - Test ChEMBL import execution

5. **Integration Testing** (2 hours)
   - Test full remediation process
   - Verify all steps complete successfully
   - Update project_state.json status

6. **Documentation & Handoff** (1 hour)
   - Update documentation
   - Create troubleshooting guide
   - Prepare handoff notes

## Success Criteria

1. All SQL operations execute successfully through MCP
2. ChEMBL ID column is added to molecules table
3. ChEMBL data is imported with correct arguments
4. Verification step confirms successful schema changes
5. project_state.json status is updated from "Blocked-Debug" to "Completed"

## Additional Context

The DDL execution issue is a common limitation with ORM and database abstraction layers. Our solution involves wrapping DDL statements in functions that return tuples, making them compatible with MCP. The PostgreSQL type compatibility issue is addressed by explicitly casting database-specific types to text to ensure JSON compatibility.

The team should:

1. Review the CHEMBL_REMEDIATION_ACTION_PLAN.md for detailed technical information
2. Use the MCP_DDL_EXECUTION_GUIDE.md for best practices on executing DDL via MCP
3. Test each fix individually before performing full end-to-end testing

## Expected Timeline

Total estimated time: 7-8 hours

| Task | Estimated Time | Dependencies |
|------|----------------|--------------|
| Code Review & Setup | 30 mins | None |
| MCP DDL Execution Fix | 2 hours | Code Review |
| Type Compatibility Fix | 1 hour | Code Review |
| Argument Format Fix | 1 hour | Code Review |
| Integration Testing | 2 hours | All fixes |
| Documentation & Handoff | 1 hour | Integration Testing |

## Assigned To

- Primary: ROO Agent #1 (SQL & DDL Expert)
- Support: ROO Agent #2 (Python & Subprocess Expert)
- Reviewer: ROO Master Orchestrator