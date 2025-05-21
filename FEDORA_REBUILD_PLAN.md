# CryoProtect Fedora Rebuild Plan

This document outlines the comprehensive plan for rebuilding and optimizing CryoProtect on Fedora Linux with incremental testing.

> **IMPORTANT UPDATE**: We are now taking a more practical approach to database operations. See [PRACTICAL_DATABASE_APPROACH.md](PRACTICAL_DATABASE_APPROACH.md) for details on our streamlined database strategy that avoids unnecessary complexity.

## Completed Tasks

- ✅ **Verify Fedora environment readiness and dependencies**
  - Created verification scripts (`verify_fedora_setup.sh`, `check_rdkit.py`)
  - Confirmed Python 3.x, SELinux, and core dependencies are working
  - Verified RDKit installation and proper detection

- ✅ **Configure SELinux contexts for application directories**
  - Applied proper contexts for application files and directories
  - Verified file access permissions for the Python environment
  - Documented SELinux configuration approach

- ✅ **Test Podman installation and compatibility**
  - Verified Podman is installed and functioning
  - Created minimal container configuration
  - Checked volume mounting with SELinux contexts

- ✅ **Configure network settings for proper connectivity**
  - Set up proper network configuration for outbound connections
  - Configured firewalld for required ports
  - Verified DNS resolution and connectivity

- ✅ **Test database connection with Supabase in Fedora**
  - Created simplified test application to verify connectivity
  - Successfully connected to Supabase REST API
  - Documented database schema in DATABASE_SCHEMA_SUMMARY.md

## Phase 1: Database and Infrastructure (High Priority)

### 1.1 Implement and test database migrations and schema updates

- **Analysis**:
  - Review existing migrations in `/migrations` directory
  - Examine current database schema from DATABASE_SCHEMA_SUMMARY.md
  - Identify any missing tables or schema updates

- **Implementation**:
  - Update migration scripts for Fedora compatibility
  - Run migration scripts incrementally with validation between steps
  - Verify schema matches expected structure

- **Testing**:
  - Create verification script to validate database schema
  - Test data insertion and retrieval for each table
  - Document any schema differences from original implementation

### 1.2 Set up connection pooling for improved performance

- **Analysis**:
  - Review CONNECTION_POOL_OPTIMIZATION_GUIDE.md
  - Identify current connection management in `database/connection_manager.py`
  - Determine optimal pool size and timeout settings for Fedora

- **Implementation**:
  - Implement optimized connection pool configuration
  - Set up retry logic and connection validation
  - Add monitoring for connection pool health

- **Testing**:
  - Create stress test script to verify connection pool behavior
  - Measure performance with different pool configurations
  - Document optimal settings for Fedora environment

### 1.3 Implement and test RLS policies for enhanced security

- **Analysis**:
  - Review RLS_OPTIMIZATION_GUIDE.md
  - Examine existing RLS policies from migrations
  - Identify any missing or suboptimal policies

- **Implementation**:
  - Apply RLS policies from migrations
  - Optimize RLS policies using security definer functions
  - Implement row-level security for all tables

- **Testing**:
  - Create verification script to test RLS policy effectiveness
  - Test access control for different user roles
  - Document RLS implementation and any Fedora-specific considerations

### 1.4 Create comprehensive test suite for Fedora environment

- **Design**:
  - Define test categories (unit, integration, system)
  - Identify critical paths and functionality to test
  - Design Fedora-specific test cases

- **Implementation**:
  - Create test script framework with setup and teardown
  - Implement tests for core functionality
  - Add Fedora-specific tests for SELinux, Podman, and systemd integration

- **Automation**:
  - Set up test automation scripts
  - Create reporting and result aggregation
  - Document test procedures and expected results

## Phase 2: Core Application (High Priority)

### 2.1 Build and test core API endpoints for molecular property calculation

- **Analysis**:
  - Review existing API endpoints in `/api` directory
  - Identify core molecular property calculation endpoints
  - Check RDKit integration points

- **Implementation**:
  - Adapt API endpoint code for Fedora compatibility
  - Implement core API endpoints incrementally
  - Ensure RDKit integration works properly

- **Testing**:
  - Create test suite for API endpoints
  - Verify property calculations match expected results
  - Document any differences from original implementation

## Phase 3: Performance Optimization (Medium Priority)

### 3.1 Optimize database queries and implement performance indexes

- **Analysis**:
  - Review database query patterns
  - Identify slow or inefficient queries
  - Analyze index requirements

- **Implementation**:
  - Apply performance indexes from migrations
  - Optimize query patterns for PostgreSQL
  - Implement query caching where appropriate

- **Testing**:
  - Measure query performance before and after optimizations
  - Verify index effectiveness
  - Document performance improvements

### 3.2 Set up service role authentication for production environment

- **Analysis**:
  - Review SERVICE_ROLE_OPTIMIZATION_GUIDE.md
  - Identify authentication flow and requirements
  - Determine proper JWT implementation

- **Implementation**:
  - Set up service role authentication
  - Implement JWT validation
  - Configure secure token management

- **Testing**:
  - Test authentication flow with different roles
  - Verify token security and validation
  - Document authentication setup and configuration

### 3.3 Test and optimize RDKit visualization on Fedora

- **Analysis**:
  - Identify visualization requirements
  - Check RDKit visualization dependencies
  - Review X11/Wayland considerations

- **Implementation**:
  - Install required visualization libraries
  - Optimize RDKit visualization for Fedora
  - Implement proper error handling for missing dependencies

- **Testing**:
  - Test visualization with different molecule types
  - Verify image generation and display
  - Document any Fedora-specific considerations

## Phase 4: Containerization and Deployment (Medium Priority)

### 4.1 Create Podman containers for development and testing

- **Analysis**:
  - Review Dockerfile and docker-compose configuration
  - Identify required modifications for Podman
  - Plan container organization

- **Implementation**:
  - Create Podman-compatible container definitions
  - Set up volume mounts with proper SELinux contexts
  - Configure networking for containers

- **Testing**:
  - Test container startup and operation
  - Verify application functionality in containers
  - Document Podman container setup

### 4.2 Set up monitoring and logging infrastructure for Fedora environment

- **Analysis**:
  - Review existing monitoring configuration
  - Identify Fedora-specific logging considerations
  - Determine monitoring requirements

- **Implementation**:
  - Configure centralized logging
  - Set up Prometheus metrics collection
  - Implement health check endpoints

- **Testing**:
  - Verify log collection and rotation
  - Test monitoring alerts and notifications
  - Document monitoring setup and configuration

### 4.3 Set up and test backup and restore functionality

- **Analysis**:
  - Review backup requirements
  - Identify critical data for backup
  - Determine backup strategy

- **Implementation**:
  - Configure database backup procedures
  - Set up file backup for critical configurations
  - Implement restoration testing

- **Testing**:
  - Test backup creation and verification
  - Verify restoration process
  - Document backup and restore procedures

## Phase 5: Additional Features (Medium Priority)

### 5.1 Implement toxicity data schema and optimizations

- **Analysis**:
  - Review TOXICITY_OPTIMIZATION_GUIDE.md
  - Identify toxicity data requirements
  - Determine schema optimization opportunities

- **Implementation**:
  - Apply toxicity schema migrations
  - Implement optimized queries for toxicity data
  - Set up data import for toxicity information

- **Testing**:
  - Test toxicity data storage and retrieval
  - Verify query performance
  - Document toxicity data implementation

## Phase 6: System Integration (Low Priority)

### 6.1 Create systemd service files for automatic startup

- **Analysis**:
  - Determine service requirements
  - Identify dependencies and ordering
  - Plan systemd integration

- **Implementation**:
  - Create systemd service definitions
  - Configure proper user and permissions
  - Set up automatic startup and restart

- **Testing**:
  - Test service startup and shutdown
  - Verify automatic restart on failure
  - Document systemd configuration

### 6.2 Implement firewalld rules for enhanced security

- **Analysis**:
  - Identify required network access
  - Determine port requirements
  - Plan firewalld configuration

- **Implementation**:
  - Configure firewalld zones
  - Set up port openings for required services
  - Implement access control rules

- **Testing**:
  - Test connectivity with firewall enabled
  - Verify proper access control
  - Document firewall configuration

## Phase 7: Documentation and Knowledge Transfer

### 7.1 Document Fedora-specific installation and configuration steps

- **Analysis**:
  - Identify all Fedora-specific configurations
  - Compile installation procedures
  - Determine troubleshooting guidance

- **Implementation**:
  - Create comprehensive installation guide
  - Document configuration options
  - Develop troubleshooting guide

- **Validation**:
  - Verify documentation accuracy
  - Test installation following documentation
  - Update based on feedback

## Incremental Testing Approach

Each phase will follow this testing approach:

1. **Unit Testing**: Test individual components in isolation
2. **Integration Testing**: Test component interactions
3. **System Testing**: Test end-to-end functionality
4. **Performance Testing**: Measure performance metrics
5. **Security Testing**: Verify security controls

After each phase, a comprehensive report will be generated documenting:
- Completed tasks
- Test results
- Issues encountered and resolutions
- Performance metrics
- Next steps

## Working Environment

Development will proceed in the following environment:

- **OS**: Fedora Linux 42
- **Python**: 3.13.3
- **Database**: PostgreSQL via Supabase
- **SELinux**: Enforcing mode
- **Container Engine**: Podman
- **IDE**: Cursor

## Timeline and Priorities

Tasks are prioritized as follows:
- **High Priority**: Critical for basic functionality
- **Medium Priority**: Important for production readiness
- **Low Priority**: Enhances functionality and maintainability

Focus on completing Phase 1 and Phase 2 before moving to other phases, but maintain flexibility to address blocking issues in any phase.

## Verification

The rebuild progress will be verified using:
- Automated test suite execution
- Manual verification of key functionality
- Performance benchmarking
- Security assessment

Updates to the verification tools will be made as needed to ensure comprehensive coverage.

## Acceptance Criteria

The Fedora rebuild will be considered successful when:

1. All high-priority tasks are completed
2. The application passes all tests
3. Performance meets or exceeds the original implementation
4. All security controls are properly implemented
5. Documentation is complete and verified

## Logging Progress

Progress will be tracked in:
- The GitHub issue system
- Task-specific documentation
- Test reports and metrics
- The project todo list