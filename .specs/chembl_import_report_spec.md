# ChEMBL Import Report Specification

## Purpose
Define the structure, required analyses, and acceptance criteria for the comprehensive ChEMBL import report (`reports/chembl_import_report.md`). The report provides a data quality assessment, security verification, and actionable recommendations for each ChEMBL import run.

## Report Structure

### 1. Executive Summary
- Brief overview of the import run (date, duration, batch size, dry run/production).
- High-level results: number of compounds processed, inserted, skipped, errors.

### 2. Data Quality Assessment
- **Completeness:** Compare total_fetched vs. total_processed vs. inserted vs. skipped.
- **Property Coverage:** Number of properties inserted; highlight any missing or anomalous property types.
- **Error Analysis:** Summarize error count and types (if any). If error log is empty, state "No errors detected."
- **Skipped Records:** If any, provide reasons and representative examples.
- **Import Arguments:** Document key parameters (limit, batch_size, dry_run, etc.).

### 3. Security Verification
- **RLS Protocol:** Confirm that Row Level Security (RLS) was verified and remediated before/after import, referencing the decorator and audit log.
- **Policy Coverage:** State that all required RLS policies for `molecule` and `molecular_property` tables were present and correct, per `.specs/rls_restoration_protocol.md`.
- **Audit Logging:** Reference the audit log location and summarize any security-related events (if present).
- **References:** Link to `docs/rls_verification_guide.md` and `.specs/rls_restoration_protocol.md`.

### 4. Actionable Recommendations
- **Data Quality:** Recommend actions if any errors, skipped records, or property anomalies are detected.
- **Security:** Recommend any further RLS or policy improvements if issues were found.
- **Process Improvements:** Suggest improvements to import parameters, logging, or validation as needed.

### 5. Supporting Evidence & References
- List and link to all relevant logs, specs, and code (e.g., `logs/chembl_summary.json`, `ChEMBL_Integrated_Import.py`, RLS protocol docs).

## Inputs
- `logs/chembl_summary.json`
- `logs/chembl_integrated_import_progress.jsonl` (if present)
- `logs/chembl_integrated_import_errors.jsonl` (if present)
- `docs/rls_verification_guide.md`
- `.specs/rls_restoration_protocol.md`
- `.specs/chembl_property_types_remediation.md` (for property type context)
- `ChEMBL_Integrated_Import.py`

## Outputs
- `reports/chembl_import_report.md` (comprehensive, human-readable report)

## Acceptance Criteria
- Report includes all sections above, even if some (e.g., errors) state "None detected."
- Data quality and security verification are clearly documented and referenced.
- Actionable recommendations are specific and justified by findings.
- All referenced files and logs are linked or cited.
- Report is suitable for both technical and non-technical stakeholders.

## References
- `.specs/rls_restoration_protocol.md`
- `docs/rls_verification_guide.md`
- `.specs/chembl_property_types_remediation.md`
- `reports/chembl_property_types_verification.md`
- `reports/chembl_molecules_schema_inspection.md`
- `ChEMBL_Integrated_Import.py`