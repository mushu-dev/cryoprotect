# Agent Task Guide for CryoProtect v2

This guide provides essential information for Roo agents working on the CryoProtect v2 project. It explains how to approach tasks efficiently to minimize token usage while maintaining high-quality implementation.

## Task Structure Principles

1. **Micro-Task Architecture**: Each task should modify <100 lines of code
2. **Single Responsibility**: Focus on one clear purpose per task
3. **Minimal Context**: Require minimal file context to understand the task
4. **Reference Implementation**: Include explicit examples to follow
5. **Clear Interfaces**: Define precise input/output contracts

## Agent Roles and Specializations

### Database Specialists
- **Error Handler**: Implements resilient error handling systems
- **State Manager**: Creates checkpoint and resumption systems
- **Data Transformer**: Converts between data formats

### Parallel Processing Specialists
- **Worker Designer**: Creates worker thread implementations
- **Workload Manager**: Designs task distribution systems
- **Result Processor**: Builds result collection and aggregation

### Chemical Data Specialists
- **Property Mapper**: Maps chemical properties between systems
- **Structure Normalizer**: Standardizes chemical representations
- **Metadata Handler**: Adds consistent metadata to entities

### Execution Specialists
- **Runner Builder**: Creates main execution scripts
- **Progress Tracker**: Implements real-time monitoring
- **Verification Expert**: Builds data quality verification

## Task Context Optimization

To minimize token usage, follow these best practices:

1. **Focus on Changed Lines**: Only consider the specific lines being modified
2. **Reference Existing Patterns**: Study reference implementations first
3. **Interface-First Design**: Focus on interfaces before implementation details
4. **Targeted File Reading**: Only load the specific files/lines needed
5. **Skip Boilerplate**: Don't waste tokens regenerating standard boilerplate

## Task Implementation Template

When implementing a task, follow this template structure:

```
TASK: [TASK-ID]: [Brief Name]

CONTEXT:
This task implements [specific functionality] for the [component].
It addresses the need for [problem statement].

FILES:
- Primary: [file_path:line_range]
- Reference: [reference_file:line_range]

IMPLEMENTATION STEPS:
1. [Step 1]
2. [Step 2]
3. [Step 3]

INTERFACE:
```python
# Key functions/classes with clear signature and docstrings
def example_function(param1: str, param2: int) -> bool:
    """
    Purpose of this function with param descriptions
    """
    pass
```

PATTERNS TO FOLLOW:
- Pattern 1: [description] - see [reference_file:line_range]
- Pattern 2: [description] - see [reference_file:line_range]

SUCCESS CRITERIA:
- [Verifiable outcome 1]
- [Verifiable outcome 2]
- [Verifiable outcome 3]
```

## Example: Optimized Task Implementation

Here's an example of a well-structured task:

```
TASK: DB-ERR-1: Error Classification System

CONTEXT:
This task implements the error classification system for ChEMBL API integration.
It addresses the need for categorizing errors and determining recovery strategies.

FILES:
- Primary: chembl/error_handler.py:1-50
- Reference: pubchem/rate_limiter.py:55-75 (for error patterns)

IMPLEMENTATION STEPS:
1. Create ErrorCategory enum for error types
2. Implement classify_error() function
3. Implement recovery_strategy() function

INTERFACE:
```python
class ErrorCategory(Enum):
    API_RATE_LIMIT = "rate_limit"
    CONNECTION_ERROR = "connection"
    DATA_VALIDATION = "validation"
    TRANSFORMATION = "transformation"
    UNKNOWN = "unknown"

def classify_error(exception: Exception) -> ErrorCategory:
    """Classify an exception into an error category"""
    pass

def recovery_strategy(category: ErrorCategory) -> Dict[str, Any]:
    """Determine recovery strategy for error category"""
    pass
```

PATTERNS TO FOLLOW:
- Error classification: see pubchem/rate_limiter.py:55-60
- Recovery strategies: see PubChem_CryoProtectants_Supabase_Enhanced.py:180-200

SUCCESS CRITERIA:
- Correctly classifies all common ChEMBL API errors
- Provides appropriate recovery strategies
- Handles unknown errors gracefully
```

## File Structure Reference

Key project files for quick reference:

```
chembl/
  |-- __init__.py
  |-- client.py         # ChEMBL client wrapper
  |-- cache.py          # Caching system
  |-- error_handler.py  # Error handling (to be implemented)
  |-- checkpoint.py     # Checkpoint system (to be implemented)
  |-- logging.py        # Logging configuration

pubchem/
  |-- __init__.py
  |-- client.py         # Reference: PubChem client
  |-- cache.py          # Reference: Caching implementation
  |-- rate_limiter.py   # Reference: Rate limiting logic

ChEMBL_Integrated_Import.py  # Main ChEMBL import script
PubChem_CryoProtectants_Supabase_Enhanced.py  # Reference implementation
```

## Database Schema Quick Reference

**Molecules Table**
- `id`: UUID PRIMARY KEY
- `name`: VARCHAR NOT NULL
- `formula`: VARCHAR
- `molecular_weight`: DOUBLE PRECISION
- `smiles`: VARCHAR
- `inchi`: VARCHAR
- `inchi_key`: VARCHAR
- `chembl_id`: VARCHAR
- `pubchem_cid`: VARCHAR
- `data_source`: VARCHAR
- `created_at`: TIMESTAMPTZ DEFAULT now()
- `updated_at`: TIMESTAMPTZ DEFAULT now()

**Molecular Properties Table**
- `id`: UUID PRIMARY KEY
- `molecule_id`: UUID REFERENCES molecules(id)
- `property_name`: VARCHAR NOT NULL
- `property_type`: VARCHAR NOT NULL
- `value`: DOUBLE PRECISION
- `unit`: VARCHAR
- `source`: VARCHAR
- `confidence`: DOUBLE PRECISION
- `created_at`: TIMESTAMPTZ DEFAULT now()
- `updated_at`: TIMESTAMPTZ DEFAULT now()

## Common Implementation Patterns

### Error Handling Pattern
```python
try:
    result = do_something_risky()
    return result
except RateLimitError as e:
    # Implement backoff and retry
    logger.warning(f"Rate limit hit: {e}, backing off for {backoff_time}s")
    time.sleep(backoff_time)
    return handle_rate_limit(e)
except ConnectionError as e:
    # Network error recovery
    logger.error(f"Connection error: {e}")
    if retries < max_retries:
        return retry_with_backoff(retries + 1)
    else:
        raise MaxRetriesExceeded(f"Failed after {max_retries} attempts") from e
except Exception as e:
    # Catch-all error handling
    logger.exception(f"Unexpected error: {e}")
    raise
```

### Checkpoint Pattern
```python
def save_checkpoint(state):
    """Save checkpoint atomically to prevent corruption"""
    temp_file = f"{checkpoint_path}.tmp"
    with open(temp_file, 'w') as f:
        json.dump(state, f, indent=2)
    os.replace(temp_file, checkpoint_path)  # Atomic operation

def load_checkpoint():
    """Load most recent checkpoint file"""
    if not os.path.exists(checkpoint_path):
        return create_initial_state()
    
    try:
        with open(checkpoint_path, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Error loading checkpoint: {e}")
        # Try backup file if main is corrupted
        if os.path.exists(f"{checkpoint_path}.bak"):
            with open(f"{checkpoint_path}.bak", 'r') as f:
                return json.load(f)
        return create_initial_state()
```

### Worker Pattern
```python
class Worker:
    def __init__(self, worker_id, task_queue, result_queue):
        self.worker_id = worker_id
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.running = False
        
    def run(self):
        self.running = True
        while self.running:
            try:
                task = self.task_queue.get(timeout=1)
                if task is None:  # Poison pill
                    break
                result = self.process_task(task)
                self.result_queue.put(result)
            except Empty:
                continue  # Queue is empty, try again
            except Exception as e:
                # Handle error but don't crash worker
                self.result_queue.put({
                    "status": "error",
                    "task_id": task.get("id"),
                    "error": str(e)
                })
                
    def process_task(self, task):
        # Implement specific task processing logic
        pass
        
    def stop(self):
        self.running = False
```

## Verification Best Practices

1. **Component Testing**: Test each small component in isolation
2. **Parameter Validation**: Add defensive parameter checking
3. **Error Scenario Testing**: Test failure modes explicitly
4. **Reference Validation**: Compare results with expected outcomes
5. **State Verification**: Confirm system state is as expected after operations

## Token Optimization Tips

1. **Minimize File Reading**: Only read what's necessary
2. **Leverage References**: Study reference implementations first
3. **Focus on Interfaces**: Understand interfaces before implementation details
4. **Template Reuse**: Use standardized patterns when possible
5. **Implementation-First**: Focus on implementation over explanation

By following these guidelines, you'll be able to implement high-quality code while minimizing token usage and maximizing efficiency.