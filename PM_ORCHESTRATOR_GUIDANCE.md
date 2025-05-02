# Project Manager Guidance: ChEMBL Integration

## Orchestration Strategy

As Project Manager, I'm providing high-level guidance for your ChEMBL integration tasks. This document focuses on **implementation efficiency** to reduce API token costs while maintaining quality.

### Core Principles

1. **Implementation-only focus**: Have ROO agents implement code directly without explanations
2. **Reference existing code**: Leverage our successful PubChem integration patterns
3. **Direct SQL execution**: Implement direct Supabase connection instead of MCP for better performance
4. **Clear boundaries**: Assign specific file sections with precise line references
5. **Sequential dependencies**: Build foundation first, then dependent components

## Key Components and Implementation Order

1. **Direct Supabase Connection (First)**
   - Create `supabase_direct.py` based on `connection_pool_wrapper.py`
   - Implement connection pooling, transaction handling, and query execution
   - Reference: Look at our PubChem implementation with MCP for SQL patterns

2. **Enhanced ChEMBL Client**
   - Modify `chembl/client.py` to improve rate limiting for Monday restrictions
   - Enhance caching and error handling for robustness
   - Reference: Our PubChem client implementation

3. **Data Transformation**
   - Update `ChEMBL_Integrated_Import.py` transformation functions
   - Map ChEMBL properties to our database schema
   - Reference: Our successful PubChem data transformation approach

4. **Database Operations Migration**
   - Replace MCP calls with direct connection in `ChEMBL_Integrated_Import.py`
   - Implement batch processing for efficiency
   - Reference: Our PubChem implementation's batch insertion pattern

5. **Progress Tracking**
   - Reuse existing progress tracking dashboard with ChEMBL adaptations
   - Ensure checkpoint-based resumable operations
   - Reference: PubChem implementation's progress tracking

## Agent Instruction Optimization

When instructing agents, follow these patterns to reduce token usage:

1. **Direct Pointers**: Use specific file:line references
   ```
   Implement in supabase_direct.py:50-80
   ```

2. **Minimal Context**: Only include essential reference information
   ```
   Reference connection_pool_wrapper.py:10-85 for pattern
   ```

3. **Clear Interfaces**: Define expected function signatures
   ```
   def execute_query(query, params=None) -> List[Dict]:
   ```

4. **Objective Criteria**: Provide specific success metrics
   ```
   Verify with: assert result[0]['test'] == 1
   ```

5. **Minimal Reporting**: Request only essential status information
   ```
   Report: Implemented connection pool ✓ | test passed ✓
   ```

## Integration Steps

Guide the orchestrator to execute these steps in order:

1. **Set up direct Supabase connection** (Foundation)
2. **Improve ChEMBL client robustness** (Data Source)
3. **Implement data transformation logic** (Processing)
4. **Update database operations** (Storage)
5. **Add progress tracking** (Monitoring)
6. **Verify implementation** (Quality)

For each step, provide:
- Specific file paths and line numbers to modify
- Reference implementation examples from existing code
- Clear success criteria for verification

## Communication Protocol

When instructing the orchestrator, use this concise format:

```
IMPLEMENT: [specific functionality]
FILE: [file path]
LINES: [line range]
REFERENCE: [existing implementation to follow]
VERIFY: [concrete success criteria]
```

Emphasize that agents should focus purely on implementation, not explanation or design discussions.

## Cost-Saving Strategies

- Have agents read reference code ONCE and implement without re-reading
- Focus on implementation efficiency rather than explanation
- Batch-assign related functions within the same file
- Process tasks sequentially with clear dependency chain
- Provide interface definitions upfront to reduce back-and-forth

## Success Criteria

The complete integration should deliver:
- Direct Supabase connection with improved performance over MCP
- ChEMBL data integration with 1,000+ compounds
- Property mapping and normalization
- Implementation following established patterns
- Robustness with error handling and recovery

Guide the orchestrator to deliver these outcomes through efficient, focused implementation tasks.