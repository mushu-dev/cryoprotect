# Phase 2 Implementation Plan: Minimal Viable Deployment

This document provides a detailed implementation plan for Phase 2 of the CryoProtect v2 project, focusing on achieving a minimal viable deployment with essential maintenance utilities, predictive models, CI/CD pipeline, and monitoring.

## Current Status

Phase 1 (Essential Core Functionality) has been successfully completed:
- Lab verification workflow is fully implemented
- Authentication requirements have been removed from core pages
- API implementation has been standardized and finalized

## Phase 2 Objectives (3 Weeks)

### Task 1: Essential Maintenance Utilities (Week 1)

#### 1.1: Complete Foreign Key Relationship Fixes

**Issue:** Some database relationships have inconsistencies that could affect data integrity.

**Implementation:**
1. Analyze remaining foreign key issues by inspecting models and database schema
2. Create new migration script to fix the remaining foreign key relationships
3. Focus on fixing relationships between:
   - Molecules and Mixtures
   - Experiments and Verifications
   - Users and shared content

**Files to modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/017_fix_foreign_keys.sql` (new file)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/fix_foreign_key_relationships.py` (update)

#### 1.2: Complete RLS Implementation Tools

**Issue:** Certain tables are missing Row Level Security (RLS) policies, which affects data access control.

**Implementation:**
1. Audit all database tables for missing RLS policies
2. Implement proper RLS policies for all tables
3. Create service role bypass policies for admin functions
4. Add verification tools to check RLS effectiveness

**Files to modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/018_complete_rls_policies.sql` (new file)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/apply_missing_rls_policies.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/verify_rls_effectiveness.py` (update)

#### 1.3: Database Health Check Utilities

**Issue:** The system lacks comprehensive database health monitoring.

**Implementation:**
1. Create a database health check module with:
   - Schema validation (table structure, column types)
   - Data integrity checks (orphaned records, invalid data)
   - Performance diagnostics (slow queries, missing indexes)
   - Regular health report generation

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/database/utils/health_check.py` (new file)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/system_resources.py` (update for API endpoint)

#### 1.4: Basic Backup/Restore Tools

**Issue:** The backup system needs completion for operational readiness.

**Implementation:**
1. Complete the scheduled backup system
   - Implement rotation policies for backups
   - Add compression for storage efficiency
   - Create metadata for backup tracking
2. Create easy-to-use restore scripts
   - Point-in-time recovery options
   - Partial restore capabilities
   - Verification steps after restore

**Files to modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/backup/backup_manager.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/backup/restore_backup.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/storage_manager.py` (update)

### Task 2: Predictive Models Completion (Week 2)

#### 2.1: Core Prediction Algorithms

**Issue:** Some core prediction algorithms for key properties need completion and optimization.

**Implementation:**
1. Complete implementation of remaining prediction algorithms:
   - Glass transition temperature (Tg) model
   - Ice recrystallization inhibition (IRI) model
   - Cell viability prediction model
2. Optimize model parameter selection for all models
3. Implement model versioning system

**Files to modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models.py` (update)

#### 2.2: Validation Framework

**Issue:** Need a robust way to validate predictive model results against experimental data.

**Implementation:**
1. Create a cross-validation framework for predictive models
2. Implement metrics for model evaluation:
   - Mean Absolute Error (MAE)
   - Root Mean Square Error (RMSE)
   - R² coefficient of determination
3. Create visualization tools for model validation
4. Add automated validation against test datasets

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/validation_framework.py` (new file)

#### 2.3: Visualization Components

**Issue:** Need better visualization of model results for scientific interpretation.

**Implementation:**
1. Enhance the visualization of model results:
   - Predicted vs. actual value scatter plots
   - Residual analysis charts
   - Feature importance visualizations
2. Add comparison charts between predicted and measured values
3. Implement uncertainty visualization (error bars, confidence bands)

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/static/js/predictive-models.js` (create)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/templates/predictive_models.html` (create)

#### 2.4: Model Training/Evaluation

**Issue:** Need a complete pipeline for model training and evaluation.

**Implementation:**
1. Complete model training pipeline:
   - Feature extraction from molecule properties
   - Hyperparameter optimization
   - Model persistence and loading
2. Add model evaluation metrics dashboard
3. Implement model versioning for reproducibility
4. Create automated retraining schedule

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models_resources.py` (create)

### Task 3: Basic CI/CD Pipeline (Week 3, First Half)

#### 3.1: GitHub Actions Workflow

**Issue:** Need automated CI/CD pipeline for reliable deployments.

**Implementation:**
1. Complete `.github/workflows/ci-cd.yml`:
   - Test jobs for Python, JavaScript, and integration tests
   - Build job for Docker image creation and testing
   - Deployment jobs for staging and production
   - Security scanning for vulnerabilities
2. Set up proper secrets management for CI/CD

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/.github/workflows/ci-cd.yml` (update)

#### 3.2: Automated Testing

**Issue:** Need comprehensive automated testing as part of CI/CD.

**Implementation:**
1. Configure test automation with pytest:
   - Unit tests for all components
   - Integration tests for end-to-end functionality
   - API tests for endpoint verification
2. Set up code coverage reporting
3. Create test fixtures for reproducible testing
4. Implement end-to-end testing with Selenium

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/run_tests.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/conftest.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/test_predictive_models.py` (create)

#### 3.3: Deployment Scripts

**Issue:** Need reliable deployment scripts for different environments.

**Implementation:**
1. Create scripts for deployment to different environments:
   - Staging environment deployment script
   - Production environment deployment script
   - Blue/Green deployment implementation
2. Implement database migration as part of deployment
3. Add rollback capability in case of deployment failure

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/deploy-blue.sh` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/deploy-green.sh` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/rollback.sh` (update)

#### 3.4: Environment Configuration

**Issue:** Need standardized environment configuration for different deployment targets.

**Implementation:**
1. Set up proper environment configuration:
   - Development environment settings
   - Staging environment settings
   - Production environment settings
2. Implement secure handling of environment variables
3. Create a template for environment configuration
4. Implement Docker secrets for sensitive configuration

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/.env.template` (create)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/config_production.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/config_staging.py` (update)

### Task 4: Minimal Monitoring (Week 3, Second Half)

#### 4.1: Error Logging

**Issue:** Need enhanced error logging for production deployments.

**Implementation:**
1. Enhance the existing logging system:
   - Structured logging with JSON format
   - Context-rich error messages
   - Proper error categorization
2. Implement log aggregation
3. Create log retention policies
4. Set up log rotation

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/logging_enhanced.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/observability.py` (update)

#### 4.2: Performance Metrics

**Issue:** Need performance monitoring for system optimization.

**Implementation:**
1. Implement performance monitoring for:
   - API response times by endpoint
   - Database query performance
   - System resource usage (CPU, memory, disk)
2. Create metrics storage system
3. Implement metrics visualization dashboard
4. Set up trending analysis for performance degradation detection

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/prometheus_metrics.py` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/grafana/dashboards/api_performance.json` (create)

#### 4.3: Health Checks

**Issue:** Need comprehensive health check system for uptime monitoring.

**Implementation:**
1. Create a comprehensive health check system:
   - System component health endpoints
   - Database connection checking
   - External service dependency verification
2. Implement monitoring for all critical components
3. Add a health check dashboard
4. Create health check API endpoint

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/check-health.sh` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/system_resources.py` (update)

#### 4.4: Basic Alerts

**Issue:** Need alerting system for critical errors and performance issues.

**Implementation:**
1. Set up alerting for critical errors:
   - Configure email alerts
   - Set up Slack notifications
   - Create on-call rotation
2. Configure performance threshold alerts
3. Implement notification channels
4. Create alert escalation policies

**Files to create/modify:**
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/alertmanager/alertmanager.yml` (update)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/prometheus/rules/api_alerts.yml` (create)
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/prometheus/rules/system_alerts.yml` (create)

## Success Criteria

Phase 2 will be considered complete when:

### Essential Maintenance Utilities
- All database foreign key relationships are properly implemented
- RLS policies are in place for all tables
- Database health check utilities are operational
- Backup/restore tools are functioning reliably

### Predictive Models
- Core prediction algorithms are fully implemented
- Validation framework is operational
- Visualization components are working
- Model training/evaluation is functional

### Basic CI/CD Pipeline
- GitHub Actions workflow is operational
- Automated testing is configured and running
- Deployment scripts are working for all environments
- Environment configuration is standardized and secure

### Minimal Monitoring
- Error logging is enhanced with structured outputs
- Performance metrics are being collected
- Health checks are operational for all components
- Basic alerting is configured for critical issues

## Testing Strategy

1. **Unit Testing**: Test individual components in isolation
2. **Integration Testing**: Test interactions between components
3. **End-to-End Testing**: Test complete workflows
4. **Performance Testing**: Test system under load
5. **Validation Testing**: Validate predictive model outputs

## Timeline

| Week | Tasks |
|------|-------|
| 1 | Essential Maintenance Utilities |
| 2 | Predictive Models Completion |
| 3 | Basic CI/CD Pipeline + Minimal Monitoring |