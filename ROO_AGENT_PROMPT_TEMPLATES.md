# ROO Agent Prompt Templates

This guide provides standardized templates for prompting ROO agents to implement specific tasks. These templates focus on clarity, specificity, and implementation guidance without excessive explanation or thinking.

## Design Philosophy

1. **Implementation Focus**: Prompts should focus on "what" and "how" to implement, not "why"
2. **File Specificity**: Always include precise file paths and line numbers
3. **Reference Examples**: Provide reference implementations from existing code
4. **Clear Boundaries**: Define exactly what the agent should modify
5. **Success Criteria**: Specify concrete validation metrics

## Template 1: File Creation Task

```
# Task: [Brief Task Description]

## Implementation Details
- Create file: [Absolute file path]
- Purpose: [1-2 sentence functional description]
- Dependencies: [List required imports]

## Core Functionality
Implement these key components:
1. [Component 1 name]
   - Lines: [Line range, e.g., 1-50]
   - Function: [Specific function(s) to implement]
   - Reference: [Reference file:line range]

2. [Component 2 name]
   - Lines: [Line range]
   - Function: [Specific function(s) to implement]
   - Reference: [Reference file:line range]

## Interface Requirements
The file must expose these functions/classes:
- `[function_name]([parameters])`: [Brief description of expected behavior]
- `class [ClassName]`: [Brief description of class purpose]

## Example Usage
```python
# Example code showing how the implementation will be used
from [module] import [component]
result = [component].[method]([sample_input])
assert result == [expected_output]
```

## Success Criteria
- [Specific measurable outcome 1]
- [Specific measurable outcome 2]
- [Specific measurable outcome 3]
```

## Template 2: File Modification Task

```
# Task: [Brief Task Description]

## Implementation Details
- Modify file: [Absolute file path]
- Lines to update: [Line range, e.g., 120-158]
- Purpose: [1-2 sentence functional description]

## Required Changes

1. [Change 1 description]
   - From:
   ```python
   # Current code snippet to be replaced
   ```
   - To:
   ```python
   # New implementation code
   ```
   - Reference: [Reference file:line range]

2. [Change 2 description]
   - From:
   ```python
   # Current code snippet to be replaced
   ```
   - To:
   ```python
   # New implementation code
   ```
   - Reference: [Reference file:line range]

## Key Requirements
- [Specific requirement 1]
- [Specific requirement 2]
- [Specific requirement 3]

## Testing Instructions
```python
# Test code to verify changes work correctly
```

## Success Criteria
- [Specific measurable outcome 1]
- [Specific measurable outcome 2]
- [Specific measurable outcome 3]
```

## Template 3: Function Implementation Task

```
# Task: Implement [function_name] in [file_path]

## Implementation Details
- Function: `[function_name]([parameters])`
- Location: [file_path:line_range]
- Purpose: [1-2 sentence functional description]

## Function Signature
```python
def [function_name]([parameter1]: [type], [parameter2]: [type]) -> [return_type]:
    """
    [Function documentation]
    
    Args:
        [parameter1]: [Parameter description]
        [parameter2]: [Parameter description]
        
    Returns:
        [Return value description]
        
    Raises:
        [Exception1]: [When this exception is raised]
        [Exception2]: [When this exception is raised]
    """
    # Implementation goes here
```

## Implementation Requirements
1. [Specific requirement 1]
2. [Specific requirement 2]
3. [Specific requirement 3]

## Reference Implementation
See [reference_file:line_range] for similar functionality:
```python
# Relevant code snippet from reference file
```

## Test Cases
The function should handle these cases:
1. [Test case 1 description]
   - Input: [Input values]
   - Expected output: [Expected output]
2. [Test case 2 description]
   - Input: [Input values]
   - Expected output: [Expected output]

## Success Criteria
- [Specific measurable outcome 1]
- [Specific measurable outcome 2]
- [Specific measurable outcome 3]
```

## Template 4: Bug Fix Task

```
# Task: Fix [brief bug description] in [file_path]

## Bug Details
- File: [file_path]
- Lines: [line_range]
- Issue: [Concise description of the bug]
- Impact: [What functionality is affected]

## Current Behavior
```python
# Problematic code snippet
```

## Expected Behavior
The code should [description of correct behavior]

## Fix Requirements
1. [Specific requirement 1]
2. [Specific requirement 2]
3. [Specific requirement 3]

## Testing Instructions
```python
# Test code to verify the fix works correctly
```

## Success Criteria
- [Specific measurable outcome 1]
- [Specific measurable outcome 2]
- [Specific measurable outcome 3]
```

## Template 5: Integration Task

```
# Task: Integrate [component1] with [component2]

## Integration Details
- Primary File: [file_path1]
- Secondary File: [file_path2]
- Integration Points: [line_ranges in both files]

## Required Changes

1. [Change 1 description]
   - File: [file_path1]
   - Lines: [line_range]
   - Implementation:
   ```python
   # Implementation code
   ```

2. [Change 2 description]
   - File: [file_path2]
   - Lines: [line_range]
   - Implementation:
   ```python
   # Implementation code
   ```

## Interface Requirements
- [component1] must expose [specific interface]
- [component2] must consume [specific interface]
- Error handling must include [specific scenarios]

## Testing Instructions
```python
# Test code to verify integration works correctly
```

## Success Criteria
- [Specific measurable outcome 1]
- [Specific measurable outcome 2]
- [Specific measurable outcome 3]
```

## Example: Real-World Prompt

Here's an example of a well-structured prompt for Task 1 in our ChEMBL integration:

```
# Task: Implement Direct Supabase Connection

## Implementation Details
- Create file: /mnt/c/Users/1edwa/Documents/CryoProtect v2/supabase_direct.py
- Purpose: Provide a direct PostgreSQL connection to Supabase for improved performance
- Dependencies: psycopg2, threading, logging

## Core Functionality
Implement these key components:

1. Connection Pool Management
   - Lines: 1-100
   - Class: SupabaseDirectConnection with singleton pattern
   - Reference: connection_pool_wrapper.py:10-85
   - Features:
     - Thread-safe connection pooling
     - Environment variable configuration
     - Connection health monitoring

2. SQL Execution Functions
   - Lines: 100-200
   - Functions: execute_query, execute_batch
   - Reference: PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834
   - Features:
     - Parameterized queries
     - Result dictionary conversion
     - Error handling with rollback

3. Transaction Management
   - Lines: 200-250
   - Methods for explicit transaction control
   - Reference: ChEMBL_Integrated_Import.py:257-339

## Interface Requirements
The file must expose these functions/classes:
- `SupabaseDirectConnection.get_instance()`: Returns singleton instance
- `execute_query(query, params=None)`: Executes single query
- `execute_batch(queries, transaction=True)`: Executes multiple queries
- `close_all()`: Properly closes all connections

## Example Usage
```python
# Example usage
from supabase_direct import SupabaseDirectConnection

# Get connection instance
db = SupabaseDirectConnection.get_instance()

# Execute query
result = db.execute_query("SELECT current_setting('role'), current_user;")
print(result)

# Execute batch in transaction
results = db.execute_batch([
    "BEGIN;",
    "INSERT INTO test_table (name) VALUES ('test');",
    "SELECT * FROM test_table;",
    "COMMIT;"
])
```

## Success Criteria
- Connection pool establishes and maintains connections to Supabase
- Queries execute successfully with proper parameter handling
- Transactions commit or rollback appropriately
- Connection resources are properly cleaned up
- Environment variables are used for credentials (never hardcoded)
```

## Anti-Patterns to Avoid

1. **Vague Boundaries**:
   - ❌ "Update the ChEMBL client to be more efficient"
   - ✅ "Update cache expiration logic in chembl/client.py:162-173"

2. **Missing References**:
   - ❌ "Implement SQL execution like we did before"
   - ✅ "Implement SQL execution based on PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834"

3. **Unclear Success Metrics**:
   - ❌ "Make sure it works well"
   - ✅ "Verify connections complete in <100ms and throughput exceeds 50 queries/second"

4. **Design Discussions**:
   - ❌ "Consider different approaches to connection pooling..."
   - ✅ "Implement connection pooling exactly as specified in reference implementation"

5. **Excessive Explanation Requests**:
   - ❌ "Explain how your implementation works and why you chose it"
   - ✅ "Implement according to specifications and verify with provided test cases"