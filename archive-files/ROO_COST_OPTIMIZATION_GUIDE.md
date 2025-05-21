# ROO Agent Cost Optimization Guide

This guide provides strategies to minimize API calls and reduce token usage when working with ROO agents on the ChEMBL integration project.

## Core Optimization Principles

1. **Implementation-Only Focus**: Agents should only implement, never explain or think aloud
2. **File Boundaries**: Strict scope control to prevent unnecessary exploration
3. **Reference Reuse**: Maximize code reuse from existing implementations
4. **Resource Bundling**: Combine related resources in single files
5. **Incremental Testing**: Test small units before proceeding to next task

## Prompt Optimization

### Token-Efficient Prompts

**Before (Wasteful)**:
```
Please implement a direct Supabase connection with proper connection pooling. 
When designing the connection pool, consider threading issues, connection
timeout handling, and proper resource management. Be sure to explain your
implementation approach and why you chose specific patterns.
```

**After (Efficient)**:
```
Implement supabase_direct.py:50-150
- Add connection pooling class (see connection_pool_wrapper.py:10-85)
- Validate credentials from environment variables
- Include execute_query and execute_batch methods
- Use singleton pattern from reference implementation
```

### Include Only Essential Context

**DO**:
- Provide exact file:line references
- Include minimal code snippets for reference
- Link to existing implementations
- Specify exact interface requirements

**DON'T**:
- Include full file contents
- Ask for explanations or justifications
- Request multiple implementation options
- Include theoretical background

## Implementation Instructions

### 1. Resource Consolidation

1. **Combine Reference Material**
   - Create condensed reference files that agents can access
   - Extract only relevant patterns from existing code
   - Place reference implementation in `/reference/chembl_patterns.py`

2. **Create Template Files**
   - Provide skeleton implementations to reduce exploration
   - Include import statements and function signatures
   - Define core interfaces for implementation

### 2. Execution Optimization

1. **Implementation Batching**
   - Group related functions in single implementation tasks
   - Focus on file-level completeness to reduce context switching
   - Prioritize building foundational components first

2. **Sequential Dependencies**
   - Define clear implementation order to prevent rework
   - Ensure dependent components are implemented before consumers
   - Establish clear interfaces between components

### 3. Resource-Efficient Testing

1. **Minimal Verification**
   - Include specific test cases in the prompt
   - Focus on key functionality only
   - Use assertions instead of print statements

2. **Staged Testing**
   - Test core functions before integration
   - Verify integration points separately
   - Only run end-to-end tests after all components work

## ROO Agent Task Structure

### Task Format

```
TASK: [Specific implementation task]
FILE: [Absolute file path]
LINES: [Line range]
REFERENCE: [Reference file:line range]

REQUIREMENTS:
1. [Specific requirement 1]
2. [Specific requirement 2]
3. [Specific requirement 3]

INTERFACE:
[Function/class signatures to implement]

TEST:
[Specific test case to verify implementation]
```

### Example Task

```
TASK: Implement direct Supabase connection pool
FILE: /mnt/c/Users/1edwa/Documents/CryoProtect v2/supabase_direct.py
LINES: 1-100
REFERENCE: connection_pool_wrapper.py:10-85

REQUIREMENTS:
1. Use psycopg2.pool.ThreadedConnectionPool
2. Read credentials from environment variables
3. Implement singleton pattern
4. Add get_connection and release_connection methods

INTERFACE:
class SupabaseDirectConnection:
    @classmethod
    def get_instance(cls)
    def execute_query(self, query, params=None)
    def close_all(self)

TEST:
db = SupabaseDirectConnection.get_instance()
result = db.execute_query("SELECT 1 as test")
assert result[0]['test'] == 1
```

## Toolkit File Optimization

Let's optimize our three existing guidance files to minimize token usage:

### 1. Integration Guide

- Reduce to core implementation patterns
- Remove redundant explanations
- Focus on code samples and interfaces
- Consolidate related sections

### 2. Task Breakdown

- Use strict tabular format
- Remove narrative elements
- Include only implementation details
- Link directly to reference code

### 3. Template Structure

- Minimize boilerplate
- Focus on interface definitions
- Reduce example complexity
- Include only essential parameters

## Communication Protocol

### Agent Output Requirements

1. **Implementation Only**:
   ```
   IMPLEMENTED: [task description]
   FILE: [file path]
   LINES: [lines implemented]
   TESTS: PASS
   ```

2. **Issue Reporting**:
   ```
   ISSUE: [concise issue description]
   FILE: [file path]
   LINE: [problematic line]
   CONTEXT: [minimal context]
   ```

3. **Completion Confirmation**:
   ```
   TASK COMPLETE
   REQUIREMENTS MET: [Yes/No]
   VERIFICATION: [verification results]
   ```

### PM Input Format

1. **Task Assignment**:
   ```
   ASSIGN: [task description]
   AGENT: [agent name]
   FILE: [file path]
   REFERENCE: [reference file:lines]
   ```

2. **Requirement Clarification**:
   ```
   CLARIFY: [specific question]
   CONTEXT: [minimal necessary context]
   ```

3. **Feedback**:
   ```
   FEEDBACK: [pass/fail]
   ISSUE: [specific issue]
   FIX: [specific correction needed]
   ```

## Implementation Reference Bundle

For the ChEMBL integration, create a single reference file with these components:

```python
# chembl_reference_patterns.py
"""
Reference implementation patterns for ChEMBL integration.
Each section contains minimal, essential code patterns.
"""

# 1. Direct Connection Pattern
def get_connection():
    """Connection acquisition pattern."""
    conn = connection_pool.getconn()
    return conn

# 2. Query Execution Pattern
def execute_query(query, params=None):
    """Query execution pattern with proper resource management."""
    conn = get_connection()
    try:
        with conn.cursor() as cursor:
            cursor.execute(query, params)
            result = cursor.fetchall() if cursor.description else None
            conn.commit()
            return result
    except Exception as e:
        conn.rollback()
        raise
    finally:
        release_connection(conn)

# 3. Batch Processing Pattern
def process_batch(batch_items, batch_size=100):
    """Batch processing pattern with chunking."""
    results = []
    for i in range(0, len(batch_items), batch_size):
        chunk = batch_items[i:i+batch_size]
        # Process chunk
        results.extend(process_chunk(chunk))
    return results

# 4. Data Transformation Pattern
def transform_item(source_item):
    """Data transformation pattern."""
    return {
        "target_field1": source_item.get("source_field1"),
        "target_field2": source_item.get("source_field2"),
        # Add remaining field mappings
    }

# 5. Progress Tracking Pattern
def update_progress(current, total):
    """Simple progress tracking pattern."""
    percent = (current / total) * 100
    return {
        "current": current,
        "total": total,
        "percent": percent,
        "remaining": total - current
    }
```

## Final ROO Optimization Checklist

Before sending tasks to ROO agents, verify:

1. ✅ Task includes only implementation instructions
2. ✅ Reference code is specific and minimal
3. ✅ File paths and line numbers are exact
4. ✅ No requests for explanations or reasoning
5. ✅ Interface requirements are complete
6. ✅ Success criteria are objective and measurable
7. ✅ Dependencies are clearly identified
8. ✅ Test cases validate only required functionality
9. ✅ Communication protocol is token-efficient
10. ✅ Implementation sequence minimizes rework