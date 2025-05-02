# Optimized Agent Task Structure

This document outlines a refined approach to task delegation for the Roo agent team, focusing on smaller, more specialized task chunks that optimize context usage and improve implementation efficiency.

## Core Principles

1. **Specialized Agent Roles**: Each agent focuses on a specific expertise area
2. **Minimalist Context**: Tasks are scoped to minimize required context
3. **Clear Boundaries**: Explicit file and line references for all tasks
4. **Independent Execution**: Tasks designed for parallel implementation when possible
5. **Verification First**: Built-in validation for each task

## Agent Role Framework

### 1. Database Specialists

- **Database Schema Expert**: Focuses on table structure and relationships
- **Query Optimization Expert**: Specializes in performance tuning
- **Data Integrity Expert**: Ensures data consistency and validation
- **Migration Specialist**: Handles schema evolution

### 2. API Specialists

- **Client Implementation Expert**: Creates resilient API clients
- **Rate Limiting Specialist**: Optimizes API usage patterns
- **Caching Strategist**: Implements multi-level caching
- **Response Parser**: Transforms API responses to internal models

### 3. Data Processing Specialists

- **Transformation Expert**: Converts between data formats
- **Validation Specialist**: Ensures data quality
- **Pipeline Architect**: Creates efficient data processing flows
- **Batch Processing Expert**: Optimizes bulk operations

### 4. Security Specialists

- **Authentication Expert**: Implements secure auth flows
- **RLS Policy Specialist**: Creates row-level security policies
- **API Security Expert**: Ensures endpoint security
- **Token Management Expert**: Handles secure token operations

### 5. Testing Specialists

- **Unit Test Expert**: Creates focused component tests
- **Integration Test Specialist**: Tests component interactions
- **Performance Test Engineer**: Validates system performance
- **Mocking Specialist**: Creates test fixtures and mocks

### 6. Orchestration Specialists

- **Process Coordinator**: Manages execution flow
- **Error Recovery Expert**: Handles failure scenarios
- **Monitoring Specialist**: Tracks and reports progress
- **Resource Manager**: Optimizes resource utilization

## Task Structure Template

Each task should follow this template format for clarity and consistent implementation:

```
Task: [Task ID]: [Brief descriptive name]
Specialist: [Specific agent role required]
Files: 
  - [file_path:line_range]
  - [additional_files as needed]
  
Reference Implementation:
  - [reference_file:line_range]
  
Dependencies:
  - [dependent_task_ids]
  
Implementation Steps:
  1. [Step 1 with specific details]
  2. [Step 2 with specific details]
  3. [...]
  
Interface:
```python
# Expected function/class signatures
def example_function(param1, param2):
    """Function docstring with clear purpose"""
    pass
```

Verification:
```python
# Test code to verify implementation
result = example_function("test", 123)
assert result == expected_value
```

Success Criteria:
  - [Measurable outcome 1]
  - [Measurable outcome 2]
  - [...]
```

## Example: Database Population Task Breakdown

### Original Large Task

```
TASK: Execute PubChem data import leveraging Sunday's higher rate limits
FILES:
- PubChem_CryoProtectants_Supabase_Enhanced.py:1-489
- CID-Synonym-curated
- logs/

IMPLEMENTATION:
1. Run PubChem_CryoProtectants_Supabase_Enhanced.py with:
   - batch_size=100 (increased from default 50)
   - PUBCHEM_API_DELAY=0.15 (optimal for Sunday limits)
   - checkpoint enabled for resumability
2. Monitor progress and API response times
3. Adjust API delay parameter as needed based on rate limit responses
4. Target 5,000+ compounds with property profiles

SUCCESS CRITERIA:
- At least 5,000 molecules imported from PubChem
- Complete property profiles for >90% of molecules
- Checkpoint files regularly updated
- Final import report shows <2% rate limit errors
```

### Optimized Task Breakdown

```
Task: DB-PUB-1: Implement Adaptive Rate Limiter
Specialist: Rate Limiting Specialist
Files:
  - pubchem/rate_limiter.py:1-80
  
Reference Implementation:
  - PubChem_CryoProtectants_Supabase_Enhanced.py:55-75
  
Dependencies:
  - None
  
Implementation Steps:
  1. Create AdaptiveRateLimiter class with day-of-week awareness
  2. Implement exponential backoff for 429 responses
  3. Add dynamic delay adjustment based on response headers
  4. Include logging of rate limit encounters
  
Interface:
```python
class AdaptiveRateLimiter:
    def __init__(self, initial_delay=0.2, max_retries=5):
        """Initialize rate limiter with configurable parameters"""
        pass
        
    def wait(self):
        """Wait appropriate time before next request"""
        pass
        
    def adjust_for_response(self, status_code, headers):
        """Adjust delay based on response"""
        pass
```

Verification:
```python
limiter = AdaptiveRateLimiter(initial_delay=0.1)
limiter.wait()  # Should pause briefly
limiter.adjust_for_response(429, {"Retry-After": "2"})
# Should increase delay and implement backoff
```

Success Criteria:
  - Correctly implements day-of-week based delays
  - Adapts to 429 responses with proper backoff
  - Logs rate limit encounters with timestamps
  - Successfully recovers from rate limiting
```

```
Task: DB-PUB-2: Implement Checkpoint Manager
Specialist: Data Pipeline Specialist
Files:
  - pubchem/checkpoint.py:1-75
  
Reference Implementation:
  - PubChem_CryoProtectants_Supabase_Enhanced.py:120-160
  
Dependencies:
  - None
  
Implementation Steps:
  1. Create CheckpointManager class for saving/loading state
  2. Implement JSON serialization with metadata
  3. Add rotation of checkpoint files (keep last 3)
  4. Include progress statistics in checkpoint data
  
Interface:
```python
class CheckpointManager:
    def __init__(self, checkpoint_dir, prefix="import"):
        """Initialize checkpoint manager"""
        pass
        
    def save_checkpoint(self, state_data):
        """Save current state to checkpoint file"""
        pass
        
    def load_latest_checkpoint(self):
        """Load the most recent checkpoint"""
        pass
        
    def get_progress_stats(self):
        """Get statistics about progress"""
        pass
```

Verification:
```python
manager = CheckpointManager("./checkpoints")
test_data = {"processed": 100, "current_batch": 2}
manager.save_checkpoint(test_data)
loaded = manager.load_latest_checkpoint()
assert loaded["processed"] == 100
```

Success Criteria:
  - Successfully saves and loads checkpoint data
  - Maintains checkpoint history (last 3 files)
  - Includes timestamp and metadata
  - Recovers correctly after process restart
```

```
Task: DB-PUB-3: Implement PubChem API Client
Specialist: API Client Specialist
Files:
  - pubchem/client.py:1-120
  
Reference Implementation:
  - PubChem_CryoProtectants_Supabase_Enhanced.py:200-280
  
Dependencies:
  - DB-PUB-1 (Rate Limiter)
  
Implementation Steps:
  1. Create PubChemClient class with configurable endpoints
  2. Implement methods for compound retrieval by CID
  3. Add methods for property retrieval
  4. Integrate with rate limiter
  
Interface:
```python
class PubChemClient:
    def __init__(self, rate_limiter=None):
        """Initialize PubChem client with optional rate limiter"""
        pass
        
    def get_compound(self, cid):
        """Get compound data by CID"""
        pass
        
    def get_properties(self, cid, property_list):
        """Get specific properties for compound"""
        pass
```

Verification:
```python
client = PubChemClient()
compound = client.get_compound("2244")
assert compound["CID"] == "2244"
assert "IUPACName" in compound
```

Success Criteria:
  - Successfully retrieves compound data
  - Handles API errors gracefully
  - Respects rate limits
  - Returns data in consistent format
```

## Benefits of This Approach

1. **Context Efficiency**: Each specialist only needs to understand their specific domain, reducing context consumption
2. **Parallel Execution**: Independent tasks can be executed simultaneously by different agents
3. **Fault Isolation**: Issues in one task don't impact others
4. **Clearer Accountability**: Each task has a single responsible specialist
5. **Improved Quality**: Specialists develop deeper expertise in their focus areas
6. **Easier Testing**: Smaller tasks are more thoroughly testable
7. **Flexible Scheduling**: Tasks can be prioritized and scheduled independently

## Implementation Guidelines

1. **Start with Infrastructure**: Begin with foundation tasks that others depend on
2. **Prioritize API Tasks**: Focus on Sunday optimization for PubChem tasks
3. **Implement Verification First**: Create validation before implementation
4. **Use Reference Code**: Leverage existing implementations as templates
5. **Monitor Progress**: Track completion status of each specialized task
6. **Optimize Critical Path**: Focus resources on tasks blocking others

## Orchestration Process

The Project Manager should:

1. **Identify Specialists**: Assign specific roles to available agents
2. **Schedule Critical Path**: Prioritize tasks on the critical path
3. **Monitor Dependencies**: Ensure dependent tasks are completed in order
4. **Review Interfaces**: Verify that component interfaces are compatible
5. **Track Progress**: Maintain a dashboard of task completion status
6. **Adjust Dynamically**: Reallocate resources based on progress

## Conclusion

By breaking down tasks into smaller, specialist-focused chunks, we can dramatically improve the efficiency of our agent team. This approach optimizes context usage, allows for parallel execution, and leverages specialized expertise for higher-quality implementations.