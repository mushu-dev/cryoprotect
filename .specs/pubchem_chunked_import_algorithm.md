# PubChem Chunked Processing Algorithm Specification

**Spec for:** ROO_PUBCHEM_IMPORT_OPTIMIZATION.md Task 2  
**Status:** Draft  
**Author:** Solution Architect  
**Date:** 2025-04-27

---

## 1. Overview

This document specifies the design for a robust, adaptive chunked processing algorithm for PubChem data import, integrating smart chunking, adaptive scheduling, checkpointing, and a circuit breaker pattern. The goal is to maximize throughput and resilience under severe PubChem API rate limiting, as outlined in ROO_PUBCHEM_IMPORT_OPTIMIZATION.md Task 2.

---

## 2. Objectives

- **Process large CID lists in optimal-sized, dynamically adjusted chunks**
- **Adapt chunk size based on recent API response times and error rates**
- **Checkpoint progress for resumability and fault tolerance**
- **Integrate with a circuit breaker to avoid repeated API failures**
- **Support smart backoff and retry mechanisms**
- **Enable detailed logging and statistics for monitoring**

---

## 3. Architecture & Component Interactions

### 3.1. Main Components

- **Chunk Generator:** Splits CID list into chunks, adapts chunk size based on feedback.
- **Chunk Processor:** Processes each chunk, handles API calls, errors, and retries.
- **Checkpoint Manager:** Saves and loads progress to allow resuming after interruption.
- **Circuit Breaker:** Monitors API failures, blocks requests when threshold exceeded, and recovers after cooldown.
- **Scheduler:** Manages timing and sequencing of chunk processing, including delays and backoff.

### 3.2. Data Flow Diagram (Mermaid)

```mermaid
flowchart TD
    A[CID List] --> B[Chunk Generator]
    B --> C[Chunk Processor]
    C --> D[PubChem Client (with Circuit Breaker)]
    C --> E[Checkpoint Manager]
    D -->|API Response| C
    C -->|Success/Failure Stats| B
    C -->|Progress| E
    C -->|Error| F[Circuit Breaker]
    F -->|Open/Close| D
```

---

## 4. Algorithm Details

### 4.1. Smart Chunking & Adaptive Sizing

- **Initial chunk size**: Configurable (default: 100)
- **Min/max chunk size**: Configurable (e.g., min=10, max=200)
- **Adaptation logic**:
  - After each chunk, record response times and error rates.
  - If recent average response time is low and error rate is low, increase chunk size (up to max).
  - If response time increases or error rate rises, decrease chunk size (down to min).
  - Use a sliding window of recent N chunks for adaptation.

#### Pseudocode

```python
def get_optimal_chunk_size(recent_response_times, recent_error_rates, current_chunk_size):
    if mean(recent_error_rates) > error_threshold:
        return max(current_chunk_size // 2, min_chunk_size)
    elif mean(recent_response_times) < fast_response_threshold:
        return min(current_chunk_size + chunk_increment, max_chunk_size)
    else:
        return current_chunk_size
```

### 4.2. Chunk Processing & Scheduling

- **For each chunk**:
  - Wait for any required delay (rate limiting, backoff).
  - Process all CIDs in the chunk using the PubChem client.
  - On API failure, use retry_with_backoff and circuit breaker.
  - On repeated failures, mark chunk as failed and checkpoint progress.
  - Log detailed stats for each chunk (success, skip, error, timing).

### 4.3. Checkpointing

- **Checkpoint file**: JSON or similar, stores processed CIDs and current state.
- **On startup**: Load checkpoint if present, resume from last successful chunk.
- **On each chunk completion**: Save checkpoint with processed CIDs and stats.
- **On interruption or crash**: Resume from last checkpoint.

### 4.4. Circuit Breaker Integration

- **CircuitBreaker class**: Used as a decorator for API calls.
- **Failure threshold**: Configurable (e.g., 5 consecutive failures).
- **Recovery timeout**: Configurable (e.g., 30 seconds).
- **When open**: All API calls fail fast; chunk processor waits and retries after cooldown.
- **Integration**: Chunk processor should check circuit state before processing; if open, wait and retry.

### 4.5. Backoff and Retry

- **retry_with_backoff**: Used for transient API errors (e.g., 503).
- **Exponential backoff**: Delay increases with each retry, up to a max.
- **Reset on success**.

---

## 5. API Contracts

### 5.1. Chunked Processing Functions

```python
def generate_chunks(cids, initial_chunk_size=100):
    # Yields chunks of CIDs, adapts size based on feedback

def process_chunk(chunk, session, delay):
    # Processes a single chunk, returns stats

def get_optimal_chunk_size(recent_response_times):
    # Returns next chunk size based on recent performance

def save_checkpoint(processed_cids):
    # Saves progress to checkpoint file

def load_checkpoint():
    # Loads progress from checkpoint file
```

### 5.2. Circuit Breaker Components

```python
class CircuitBreaker:
    def __init__(self, failure_threshold=5, reset_timeout=30)
    def __call__(self, func)
    def record_success(self)
    def record_failure(self)
    def reset(self)

def exponential_backoff(retry_count, base_delay=0.5, max_delay=60)
```

---

## 6. Configuration

- All parameters (chunk size, thresholds, delays, circuit breaker settings) should be configurable via config.py (see ROO_PUBCHEM_IMPORT_OPTIMIZATION.md Integration I2).
- Logging should use the enhanced logging system (see Integration I4).

---

## 7. Acceptance Criteria

- Successfully processes 5,000+ compounds in manageable, adaptive chunks.
- Dynamically adjusts chunk size based on API response times and error rates.
- Resumes processing from last successful point after interruptions.
- Maintains optimal throughput without triggering API rate limits or repeated failures.
- Circuit breaker prevents repeated API failures and recovers gracefully.
- Detailed logs and statistics are generated for each chunk and overall run.

---

## 8. References

- [ROO_PUBCHEM_IMPORT_OPTIMIZATION.md](../ROO_PUBCHEM_IMPORT_OPTIMIZATION.md)
- [pubchem/client.py](../pubchem/client.py)
- [pubchem/cache.py](../pubchem/cache.py)
- [config.py](../config.py)
- [logging_enhanced.py](../logging_enhanced.py)

---

## 9. Implementation Notes

- If import_pubchem_data_chunked.py does not exist, create it as a new module.
- Ensure all new functions/classes are unit tested.
- Use dependency injection for configuration and logging to facilitate testing and future enhancements.