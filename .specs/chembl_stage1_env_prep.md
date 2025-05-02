# ChEMBL Integration Stage 1: Environment Preparation Specification

## Objective
Establish a robust, secure, and resilient environment for ChEMBL data import into CryoProtect v2, covering:
- Secure ChEMBL API credential management
- Centralized, structured logging
- Checkpointing for resumable, fault-tolerant processing

---

## 1. Credential Management

- **Configuration Source:** All credentials and config values are loaded from environment variables via the config class system (`BaseConfig` and its subclasses in `config.py`).
- **ChEMBL API Credentials:**
  - If the ChEMBL API requires authentication, add a `CHEMBL_API_KEY` field to `BaseConfig` and ensure it is loaded from the environment (e.g., `os.getenv('CHEMBL_API_KEY')`).
  - Document required environment variables in `.env.template` and README.
  - Validate presence of required credentials in `validate_config()`.
- **Supabase/MCP Credentials:** Already managed via existing config fields.

---

## 2. Logging Architecture

- **Progress Logging:**
  - Use a dedicated logger (`progress_logger`) for structured progress reporting.
  - Log to both file (e.g., `logs/chembl_progress.jsonl`) and console.
  - Use JSON format for log entries to facilitate parsing and monitoring.
- **Error Logging:**
  - Capture and log all errors with structured context (error type, message, stack trace, source).
  - Log both to file (e.g., `logs/chembl_errors.jsonl`) and console.
- **Logger Initialization:**
  - Ensure log directory exists at startup.
  - Avoid duplicate log entries by managing handlers carefully.
- **Log Rotation:** (Optional, for large imports) Consider rotating logs or archiving after import.

---

## 3. Checkpointing System

- **Purpose:** Enable robust, resumable ChEMBL data import by recording progress at regular intervals.
- **Checkpoint Data:**
  - Store the last successfully processed ChEMBL ID or batch index.
  - Store in a dedicated checkpoint file (e.g., `cache/chembl/checkpoint.json`).
  - Include timestamp, batch size, and any relevant metadata.
- **Checkpointing Logic:**
  - After each successful batch, update the checkpoint file.
  - On startup, check for an existing checkpoint and resume from the last recorded position.
  - Integrate with the caching system (`ChEMBLCache`) for efficient lookups and to avoid redundant API calls.
- **Failure Recovery:**
  - On error, log the failure and ensure checkpoint reflects the last successful state.
  - Support manual reset/override of checkpoint for reprocessing if needed.

---

## 4. Implementation Guidance

- **Reference Files:**
  - `config.py` for configuration and credential management.
  - `ChEMBL_CryoProtectants_Supabase.py` for logging and checkpoint integration.
  - `chembl/cache.py` for cache and checkpoint storage.
- **Testing:**
  - Validate that missing/invalid credentials cause a clear error and abort.
  - Simulate interruption and verify that import resumes from the last checkpoint.
  - Confirm logs are written in the correct format and location.

---

## Acceptance Criteria

- All required credentials are loaded securely from environment variables and validated at startup.
- Progress and error logs are written in structured (JSON) format to both file and console.
- Checkpointing system records and restores progress, enabling resumable import.
- All design elements are documented and referenced in implementation tasks.