# RLS Restoration & Verification Protocol for `molecule` and `molecular_property` Tables

## Objective

Guarantee that Row Level Security (RLS) and all required policies for the `molecule` and `molecular_property` tables are always enabled and correct after any import or DDL operation, minimizing risk of accidental exposure or privilege escalation.

---

## 1. RLS State & Policy Verification

**Tables Covered:**
- `molecule`
- `molecular_property`

**Verification Steps:**
1. **RLS Enabled:** Confirm `ALTER TABLE ... ENABLE ROW LEVEL SECURITY` is active for both tables.
2. **Required Policies Present:** Confirm all policies from `migrations/006_rls_policies.sql` and `migrations/007_service_role_rls.sql` are present and enabled:
   - Project member SELECT/INSERT/UPDATE/DELETE policies.
   - Service role INSERT policy.
3. **No Unauthorized Policies:** Ensure no policies exist that would allow broader access than intended.
4. **Policy Logic Integrity:** Optionally, verify policy definitions match the canonical SQL (hash or text compare).

---

## 2. Python Decorator/Wrapper: `@ensure_rls_restored`

**Purpose:**  
Guarantee that, before and after any import or DDL operation, RLS is enabled and all required policies are present for the two tables.

**Behavior:**
- On function entry:
  - Verify RLS and policies as above.
  - If any are missing/incorrect, auto-remediate (re-apply canonical SQL).
  - Log/audit any changes or discrepancies.
- On function exit (including exceptions):
  - Re-verify and auto-remediate as above.
  - Log/audit as above.

**Usage Example:**
```python
from rls_utils import ensure_rls_restored

@ensure_rls_restored
def import_compounds_to_database(...):
    ...
```

---

## 3. Logging & Auditing

- All RLS state changes, discrepancies, and remediations must be logged to a dedicated audit log (e.g., `logs/rls_audit.log`).
- Log entries should include timestamp, table, action (verify/remediate), details of the change, and user/script context.

---

## 4. Acceptance Criteria

- RLS and all required policies are always enabled and correct for `molecule` and `molecular_property` after any import or DDL operation.
- Decorator/wrapper is applied to all relevant import/DDL entry points.
- All RLS state changes and remediations are auditable.
- Unit and integration tests verify:
  - RLS is not accidentally disabled.
  - Policies are not missing or altered.
  - Remediation occurs if drift is detected.

---

## 5. References

- `migrations/006_rls_policies.sql`
- `migrations/007_service_role_rls.sql`
- `ChEMBL_Integrated_Import.py`
- `use_mcp_tool.py`