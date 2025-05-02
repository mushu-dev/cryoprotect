# Roo Code PM Implementation Directive: Phase 2 - Minimal Viable Deployment

## Overview

This directive provides detailed implementation instructions for Phase 2 of the CryoProtect v2 project. Phase 2 focuses on achieving a Minimal Viable Deployment with four key components:

1. **Essential Maintenance Utilities**
2. **Predictive Models Completion**
3. **Basic CI/CD Pipeline**
4. **Minimal Monitoring**

The Phase 2 implementation builds on our successful completion of Phase 1, which delivered the core Lab Verification workflow.

## Implementation Strategy

For optimal efficiency, use the following implementation strategy:

1. **Line References**: All implementation tasks should include specific file paths and line numbers where changes are needed, following the pattern established in Phase 1
2. **Component-Based Approach**: Implement each major component in sequence to ensure a cohesive solution
3. **Leverage Existing Patterns**: Follow established code patterns in the codebase for consistency
4. **Test-Driven Development**: Create or update tests for each component before implementation

## Task 1: Essential Maintenance Utilities (1 Week)

### 1.1: Foreign Key Relationship Fixes

**Objective**: Complete all database relationship fixes to ensure data integrity.

**Implementation**:
1. Implement `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/017_fix_foreign_keys.sql` with:
   - Identification and cleanup of orphaned records
   - Addition of missing foreign key constraints
   - Support for cascade deletion where appropriate

2. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/fix_foreign_key_relationships.py` to properly apply the migration

3. Focus on relationships between:
   - Molecules and Mixtures
   - Experiments and Verifications
   - Users and shared content

**Success Criteria**:
- Run `python verify_database_integrity.py` with no errors
- Verify through SQL queries that no orphaned records exist

### 1.2: RLS Implementation Tools

**Objective**: Ensure all tables have proper Row Level Security policies.

**Implementation**:
1. Create `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/018_complete_rls_policies.sql` that:
   - Defines standard policies for SELECT, INSERT, UPDATE, DELETE operations
   - Implements service role bypass for admin functions
   - Applies policies to all tables consistently

2. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/apply_missing_rls_policies.py` to apply these policies

3. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/verify_rls_effectiveness.py` to test policy effectiveness

**Success Criteria**:
- All tables have RLS enabled
- All tables have appropriate policies for each operation
- Service role can bypass RLS for admin functions
- Regular users can only access their own data

### 1.3: Database Health Check Utilities

**Objective**: Create comprehensive database health monitoring.

**Implementation**:
1. Create `/mnt/c/Users/1edwa/Documents/CryoProtect v2/database/utils/health_check.py` with:
   - Schema validation (table structure, column types)
   - Data integrity checks (orphaned records, invalid data)
   - Performance diagnostics (slow queries, missing indexes)
   - Regular health report generation

2. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/system_resources.py` to add API endpoints for:
   - `/api/v1/system/database/health` - Get health check report
   - `/api/v1/system/database/fix` - Apply automatic fixes

**Success Criteria**:
- Health check runs without errors
- Health check identifies common database issues
- API endpoints return proper responses

### 1.4: Basic Backup/Restore Tools

**Objective**: Complete the backup system for operational readiness.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/backup/backup_manager.py` with:
   - Rotation policies for backups
   - Compression for storage efficiency
   - Metadata for backup tracking

2. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/backup/restore_backup.py` with:
   - Point-in-time recovery options
   - Partial restore capabilities
   - Verification steps after restore

3. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/storage_manager.py` with:
   - Remote storage options
   - Secure transfer mechanisms
   - Backup integrity verification

**Success Criteria**:
- Successful backup creation with compression
- Successful restoration from backup
- Backup rotation policy works as expected
- Backup metadata is correctly stored and retrieved

## Task 2: Predictive Models Completion (1 Week)

### 2.1: Core Prediction Algorithms

**Objective**: Complete all core prediction algorithms for cryoprotectant properties.

**Implementation**:
1. Enhance `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models.py` with:
   - Glass transition temperature (Tg) model (lines 1190-1240)
   - Ice recrystallization inhibition (IRI) model (new method)
   - Cell viability prediction model (new method)

2. Implement model versioning system in the same file

**Success Criteria**:
- All prediction algorithms work correctly
- Models can be versioned and retrieved by version
- Predictions have appropriate confidence intervals

### 2.2: Validation Framework

**Objective**: Create a robust validation framework for predictive models.

**Implementation**:
1. Create `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/validation_framework.py` with:
   - Cross-validation functionality
   - Model evaluation metrics
   - Comparison between predictions and experimental data

2. Update the `evaluate()` method in `PredictiveModel` class to use this framework

**Success Criteria**:
- Cross-validation returns appropriate metrics
- Validation results can be visualized
- Framework handles different model types consistently

### 2.3: Visualization Components

**Objective**: Enhance visualization of model results for better interpretation.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/static/js/predictive-models.js` to add:
   - Predicted vs. actual scatter plots
   - Residual analysis charts
   - Feature importance visualizations
   - Uncertainty visualization (error bars, confidence bands)

2. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/templates/predictive_models.html` to include:
   - New visualization containers
   - Tabs for different visualization types
   - Interactive controls for visualization options

**Success Criteria**:
- All visualizations render correctly
- Visualizations are interactive and responsive
- Visualizations accurately represent model data

### 2.4: Model Training/Evaluation

**Objective**: Complete the end-to-end model training and evaluation pipeline.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models.py` to enhance:
   - Feature extraction (lines 130-190)
   - Hyperparameter optimization (lines 538-570)
   - Model persistence and loading (lines 572-640)

2. Create `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/predictive_models_resources.py` with:
   - API endpoints for model training
   - API endpoints for model evaluation
   - API endpoints for model comparison

**Success Criteria**:
- Models can be trained from the UI
- Models can be evaluated with test data
- Models persist correctly between application restarts

## Task 3: Basic CI/CD Pipeline (0.5 Week)

### 3.1: GitHub Actions Workflow

**Objective**: Set up automated CI/CD pipeline for reliable deployments.

**Implementation**:
1. Complete `.github/workflows/ci-cd.yml` with:
   - Test jobs for Python, JavaScript, and integration tests
   - Build job for Docker image creation and testing
   - Deployment jobs for staging and production
   - Security scanning for vulnerabilities

**Success Criteria**:
- Workflow runs successfully on push
- Tests run and failures are reported
- Docker image is built and tested
- Security scanning identifies vulnerabilities

### 3.2: Automated Testing

**Objective**: Enhance automated testing for reliable CI/CD.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/run_tests.py` to support:
   - Parallel test execution
   - Test categorization (unit, integration, etc.)
   - Coverage reporting

2. Create new test files:
   - `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/test_predictive_models.py`
   - `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/test_database_health.py`
   - `/mnt/c/Users/1edwa/Documents/CryoProtect v2/tests/test_backup_restore.py`

**Success Criteria**:
- All tests run successfully
- Code coverage meets minimum threshold (70%)
- Test results are reported in CI pipeline

### 3.3: Deployment Scripts

**Objective**: Create reliable deployment scripts for different environments.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/deploy-blue.sh` and `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/deploy-green.sh` with:
   - Environment-specific deployment steps
   - Database migration execution
   - Service restart procedures
   - Health check verification

2. Create `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/rollback.sh` for deployment failure recovery

**Success Criteria**:
- Blue/green deployment works correctly
- Database migrations are applied during deployment
- Rollback restores previous working state
- Health checks verify successful deployment

### 3.4: Environment Configuration

**Objective**: Standardize environment configuration for different deployment targets.

**Implementation**:
1. Create `.env.template` with all required environment variables
2. Update `config_production.py` and `config_staging.py` to load from environment
3. Implement Docker secrets handling in the application

**Success Criteria**:
- Application loads configuration from environment
- Sensitive values are stored securely
- Different environments use appropriate configurations

## Task 4: Minimal Monitoring (0.5 Week)

### 4.1: Error Logging

**Objective**: Enhance error logging for production deployments.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/logging_enhanced.py` with:
   - Structured logging with JSON format
   - Context-rich error messages
   - Log rotation and retention

2. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/observability.py` to:
   - Add request tracing with correlation IDs
   - Capture request context in logs
   - Track request durations

**Success Criteria**:
- Logs contain structured data in JSON format
- Error logs include full context (user, request, etc.)
- Logs are properly rotated and archived

### 4.2: Performance Metrics

**Objective**: Implement performance monitoring for system optimization.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/prometheus_metrics.py` to track:
   - API response times by endpoint
   - Database query performance
   - System resource usage (CPU, memory, disk)

2. Create Grafana dashboards in `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/grafana/dashboards/`

**Success Criteria**:
- Metrics are collected correctly
- Dashboards display real-time metrics
- Performance bottlenecks can be identified

### 4.3: Health Checks

**Objective**: Create comprehensive health check system for uptime monitoring.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/check-health.sh` to check:
   - System components health
   - Database connection
   - External service dependencies

2. Add API endpoint for health status in `/mnt/c/Users/1edwa/Documents/CryoProtect v2/api/system_resources.py`

**Success Criteria**:
- Health check script returns accurate status
- API endpoint returns health status in standard format
- Health checks can identify common failure modes

### 4.4: Basic Alerts

**Objective**: Set up alerting for critical errors and performance issues.

**Implementation**:
1. Update `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/alertmanager/alertmanager.yml` with:
   - Email alert configuration
   - Slack notification configuration
   - Alert routing rules

2. Create Prometheus alert rules in:
   - `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/prometheus/rules/api_alerts.yml`
   - `/mnt/c/Users/1edwa/Documents/CryoProtect v2/monitoring/prometheus/rules/system_alerts.yml`

**Success Criteria**:
- Alerts fire when conditions are met
- Notifications are delivered to configured channels
- Alert thresholds are appropriate (not too noisy)

## Advanced Techniques for Implementation

### Optimize API Calls

When implementing Phase 2, leverage advanced techniques to make your API calls more efficient:

1. **Batched Operations**: Use BatchTool wherever possible to perform multiple database operations in a single request
2. **Connection Pooling**: Leverage connection pooling for database operations to reduce overhead
3. **Caching**: Implement appropriate caching for frequently accessed data
4. **Optimized Queries**: Use query optimization techniques to minimize database load

### Leverage Chart.js for Visualizations

When implementing visualizations, use Chart.js effectively:

1. **Responsive Design**: Ensure all charts are responsive and mobile-friendly
2. **Interactive Elements**: Add tooltips, zooming, and panning for better user experience
3. **Consistent Styling**: Use the established color palette for visual consistency
4. **Accessibility**: Ensure all visualizations are accessible with proper ARIA attributes

### Implement Effective Workflow Patterns

For CI/CD implementation:

1. **Trunk-Based Development**: Configure workflow for trunk-based development
2. **Feature Flags**: Implement feature flags for controlled rollout
3. **Atomic Deployments**: Ensure deployments are atomic (all-or-nothing)
4. **Canary Deployments**: Add support for canary releases for risk mitigation

## Testing Strategy

For Phase 2, follow this testing strategy:

1. **Unit Tests**: Test individual components in isolation
2. **Integration Tests**: Test interactions between components
3. **End-to-End Tests**: Test complete workflows
4. **Performance Tests**: Test system under load
5. **Validation Tests**: Validate predictive model outputs

Ensure all tests are automated and included in the CI pipeline.

## Final Deliverables

At the completion of Phase 2, the following must be delivered:

1. All code implemented according to specifications
2. Comprehensive test suite with at least 70% code coverage
3. Documentation for all new functionality
4. Migration scripts for database changes
5. Configuration templates for all environments
6. Detailed report on Phase 2 completion

## Security Considerations

Throughout Phase 2 implementation, ensure:

1. No sensitive information is hardcoded or committed to the repository
2. All user inputs are properly validated and sanitized
3. Database queries are protected against SQL injection
4. API endpoints have appropriate authentication and authorization

## Collaboration Strategy

For optimal collaboration:

1. **Daily Updates**: Provide daily updates on task progress
2. **Clear Documentation**: Document all implementation decisions
3. **Code Reviews**: Request code reviews for critical components
4. **Issue Tracking**: Track all issues and blockers in the issue tracker

## Next Steps After Phase 2

After completing Phase 2, we will proceed to Phase 3: Essential Documentation, which will focus on:

1. README Standardization
2. API Documentation
3. User Documentation
4. Deployment Guides

## Timeline

| Week | Tasks |
|------|-------|
| 1 | Essential Maintenance Utilities |
| 2 | Predictive Models Completion |
| 3 | Basic CI/CD Pipeline + Minimal Monitoring |