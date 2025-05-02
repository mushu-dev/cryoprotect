# ChEMBL-Supabase Preparation & Verification Specification

## Objective
Establish a robust, auditable, and recoverable foundation for large-scale ChEMBL-driven population of the CryoProtect v2 Supabase database. This includes configuration verification, database resource checks, checkpointing, error logging, and operational safeguards for high-volume data ingestion.

---

## 1. Configuration Verification

### 1.1. Environment Variables & Config
- **Required:** `SUPABASE_URL`, `SUPABASE_KEY`, `CHEMBL_API_DELAY`, `CHEMBL_ID_FILE` (optional), and any batch/checkpoint parameters.
- **Validation:**
  - All required variables must be present at runtime; fail-fast if missing.
  - Use the hierarchical config system (`config.py`, `.env`) for all config access.
  - Log all config values (except secrets) at startup for traceability.
  - Validate types and formats (URLs, numbers, file paths).

### 1.2. API Connectivity
- **Supabase:** Attempt a test query (e.g., list tables) at startup; log and abort on failure.
- **ChEMBL:** Attempt a test search (e.g., for "glycerol"); log and abort on failure.

---

## 2. Database Resource Verification

### 2.1. Schema & Table Checks
- **Pre-population audit:** Run `supabase_database_audit.py` or equivalent before any data operation.
- **Checks:**
  - All required tables exist: `molecule`, `molecular_property`, `mixture`, `mixture_component`, `experiment`, `experiment_property`, `prediction`, `calculation_method`, etc.
  - Table relationships (foreign keys) match expected schema.
  - Log a summary report; abort if any critical resource is missing.

### 2.2. Quota & Access
- **Check:** Database quotas (row, storage) and user permissions.
- **Log:** Any warnings about approaching limits.

---

## 3. Checkpointing & Error Logging

### 3.1. Checkpointing
- **Mechanism:** 
  - Use a JSON or CSV checkpoint file to record progress (e.g., last processed ChEMBL ID, batch number, timestamp).
  - Update checkpoint after each batch or N molecules (parameterizable).
  - On resume, load checkpoint and skip completed work.
- **Location:** Store in a dedicated `checkpoints/` directory, with timestamped filenames for traceability.
- **Format Example:**
  ```json
  {
    "last_chembl_id": "CHEMBL12345",
    "batch": 7,
    "timestamp": "2025-04-26T00:00:00Z",
    "total_processed": 3500
  }
  ```

### 3.2. Error Logging
- **Mechanism:** 
  - Log all errors (API, DB, validation) to a dedicated error log file (e.g., `skipped_chembl_molecules.log`).
  - Include context: ChEMBL ID, error type, message, stack trace if available.
  - For skipped molecules, log reason and relevant data.
- **Structure:** Use structured logging (JSON lines or CSV) for easy post-hoc analysis.

---

## 4. Operational Safeguards

### 4.1. Rate Limiting & Batch Processing
- **ChEMBL API:** Enforce rate limiting (parameterizable, default ≤5 requests/sec).
- **Batch Size:** Parameterizable; default to 100–500 molecules per batch.
- **Retries:** Implement exponential backoff for transient errors (API, DB).
- **Memory Management:** Monitor memory usage; avoid loading all data at once.

### 4.2. Progress Reporting & Output Standardization
- **Progress:** Log progress and ETA after each batch (to file and console).
- **Output:** Standardize all logs and outputs for downstream analysis and reproducibility.

---

## 5. Acceptance Criteria

- All required config and credentials are validated at startup; failures are logged and abort execution.
- Database schema and table relationships are verified before data population; failures abort execution.
- Checkpointing is robust, supports resumability, and is updated after each batch.
- All errors and skipped records are logged with full context.
- Rate limiting, batch processing, and memory safeguards are enforced.
- Progress and summary reports are logged in a standardized format.

---

## 6. References

- `ChEMBL_CryoProtectants_Supabase.py` (API config, checkpointing, batch logic)
- `supabase_adapter.py` (response normalization, error handling)
- `supabase_database_audit.py` (schema verification)
- `.env.template`, `config.py` (config system)
- `checkpoints/`, `logs/` directories

---

## 7. Mermaid Diagram: Preparation Workflow

```mermaid
flowchart TD
    A[Start: Load Config] --> B{Validate Config}
    B -- OK --> C[Test Supabase Connection]
    B -- Fail --> Z[Abort & Log Error]
    C -- OK --> D[Test ChEMBL API]
    C -- Fail --> Z
    D -- OK --> E[Run Database Audit]
    D -- Fail --> Z
    E -- OK --> F{Checkpoint Exists?}
    E -- Fail --> Z
    F -- Yes --> G[Resume from Checkpoint]
    F -- No --> H[Start New Population]
    G --> I[Batch Processing Loop]
    H --> I
    I --> J[Update Checkpoint]
    I --> K[Log Errors/Skipped]
    I --> L[Log Progress]
    I --> M{More Data?}
    M -- Yes --> I
    M -- No --> N[Summary Report]
    N --> O[End]