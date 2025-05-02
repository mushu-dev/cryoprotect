# CryoProtect v2 Database Remediation Verification Report

Generated: 2025-04-18 12:54:00

## Overall Status: ⚠️ WARNING

## Summary

| Verification | Status |
|-------------|--------|
| Schema Standardization | ✅ SUCCESS |
| Foreign Key Constraints | ⚠️ WARNING |
| Rls Policies | ⚠️ WARNING |
| Api Endpoints | ✅ SUCCESS |
| Performance Benchmarks | ❌ FAILED |
| Data Integrity | ✅ SUCCESS |

## Detailed Results

### Schema Standardization

Status: SUCCESS

Message: All table names are properly standardized to plural form

### Foreign Key Constraints

Status: WARNING

Message: Found 2 unexpected foreign key relationships

Missing relationships:

- projects.team_id -> teams.id

Unexpected relationships:

- mixture_components.experiment_id -> experiments.id
- experiments.project_id -> projects.id

### RLS Policies

Status: WARNING

Message: Found 2 tables without owner access policy

Tables without owner policy:

- calculation_methods
- property_types

### API Endpoints

Status: SUCCESS

Message: All API endpoints are working correctly

Success rate: 100.0%

### Performance Benchmarks

Status: FAILED

Message: Performance benchmark success rate is only 60.0%

Success rate: 60.0%

Benchmark results:

| Benchmark | Avg Time (s) | Expected Max Time (s) | Status |
|-----------|--------------|----------------------|--------|
| List all molecules | 0.045 | 0.500 | ✅ |
| Get molecule with properties | 0.087 | 0.500 | ✅ |
| Get mixture with components | 0.092 | 0.500 | ✅ |
| Search molecules by name | 0.638 | 0.500 | ❌ |
| Complex join with RLS | 1.256 | 1.000 | ❌ |

### Data Integrity

Status: SUCCESS

Message: All data integrity checks passed

Orphaned record checks:

| Check | Count | Status |
|-------|-------|--------|
| Mixture components without valid mixture | 0 | ✅ |
| Mixture components without valid molecule | 0 | ✅ |
| Molecular properties without valid molecule | 0 | ✅ |
| Predictions without valid mixture | 0 | ✅ |
| Experiments without valid mixture | 0 | ✅ |

Empty table checks:

| Check | Count | Status |
|-------|-------|--------|
| Empty molecules table | 1245 | ✅ |
| Empty property_types table | 32 | ✅ |
| Empty calculation_methods table | 8 | ✅ |