# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository. As the Project Manager (PM), you coordinate the overall project and delegate tasks to the Roo agent team, which uses a minimalist orchestration approach based on the rooroo framework.

## Project Management & Delegation Structure

### Hierarchy and Responsibilities
- **Claude Code CLI (You)**: Top-level Project Manager that oversees the entire project
  - Responsible for high-level planning, assessment, and coordination
  - Delegates complex tasks to the Roo agent team
  - Handles simple tasks directly when appropriate
  - Tracks overall project progress and priorities
  - Makes strategic decisions about resource allocation

- **Roo Agent Team**: Implementation team with specialized roles
  - Coordinated by the Master Orchestrator agent in Roo
  - Consists of specialized "Swiss Army Knife" agents
  - Works on specific implementation tasks based on your delegation
  - Reports progress back to you for integration into overall project plan

### Task Delegation Protocol
- Assess task complexity before delegation
- For simple tasks: Handle directly with Claude Code CLI
- For complex tasks: Create delegation plan for Roo agent team
- Include specific file:line references in all delegations (example.py:45-72)
- Provide clear success criteria for each delegated task
- Request structured reports for completed Roo agent work

### Interaction Pattern
1. PM (you) assesses project needs and creates plans
2. PM delegates appropriate tasks to Roo agent team
3. Roo Master Orchestrator distributes work to specialized agents
4. Specialized agents complete work and report to Master Orchestrator
5. Master Orchestrator synthesizes results and reports back to PM
6. PM integrates results into overall project

## MCP Server Integration

### Context7 MCP Server
- **Purpose**: Provides up-to-date documentation access and context-aware assistance
- **Usage**: Use for library version issues, API documentation, and syntax updates
- **Commands**: 
  - `generate_documentation`: Create documentation for specific features
  - `extract_concepts`: Extract key concepts from documentation
  - `create_diagram`: Generate visual representations
  - `summarize_document`: Create summaries for large documentation files
- **Best Practices**: 
  - Check Context7 when encountering versioning issues
  - Use for real-time documentation validation
  - Leverage for API specification updates

### Supabase MCP Server
- **Purpose**: Direct integration with Supabase database operations
- **Usage**: Use for database queries, schema management, and RLS policy testing
- **Commands**:
  - `execute_sql`: Run SQL queries directly
  - `list_tables`: View database schema
  - `get_project`: Access project configuration
  - `apply_migration`: Apply database migrations
- **Best Practices**:
  - Test RLS policies using execute_sql
  - Verify database changes before implementing in code
  - Use for real-time database state inspection

## Workflow Optimization Guidelines

### API Cost Management
- Prefer delta context mode: Only analyze changed files when possible
- Use chunking: Focus on one directory/component at a time
- Minimize document generation: Don't create unnecessary documentation files
- Use precise file:line references for targeted changes
- Optimize token usage by delegating appropriate tasks to Roo
- Leverage MCP servers for documentation and database queries to reduce token usage

### Task Management Approach
- Handle simple tasks directly through Claude Code CLI
- Delegate complex tasks to Roo agent team through Master Orchestrator
- Provide detailed instructions with file paths and line numbers
- Store task plans in existing issue tracking system, not as separate files
- Track completion status of each delegated task
- Use Context7 for documentation-related tasks instead of regenerating documentation

### File Management
- Prune temporary files after task completion
- Avoid creating new versions of implementation plans; update existing ones
- Use explicit file paths when referring to code that needs modification
- Reuse existing documentation rather than regenerating
- Maintain clean directory structure according to project standards

### Thinking Modes
- Use "think" for simple tasks
- Use "think hard" for complex architecture decisions
- Use "think harder" for challenging debugging
- Use "ultrathink" for system-wide refactoring

## Build/Run/Test Commands

- **Run application**: `python app.py` (Windows: `run_app.bat`)
- **Run all tests**: `python tests/run_tests.py` or `run_tests.bat`/`run_tests.sh`
- **Run single test**: `python tests/run_tests.py -t test_file_name.py` 
- **Specific test pattern**: `python -m unittest tests/test_specific_file.py::TestClass::test_method`
- **API tests only**: `python tests/run_supabase_api_tests.py`
- **Set up environment**: `conda env create -f environment.yml`

## Code Style Guidelines

- **Docstrings**: Triple-quoted docstrings for all modules, classes, and functions
- **Imports**: Group imports (stdlib, third-party, local) with empty line between groups
- **Formatting**: 4-space indentation, 79-character line length
- **Typing**: Use type hints for function arguments and return values
- **Naming**: snake_case for functions/variables, PascalCase for classes, UPPER_CASE for constants
- **Error handling**: Use try/except with specific exception types and proper logging
- **Logging**: Use the logging module instead of print statements
- **Documentation**: Include detailed docstrings with parameters and return values
- **Testing**: Write unit tests with clear test cases and assertions

## Project Completion Plan

### Current Phase: Technical Foundation
- ✅ Code Cleanup (Phase 1.1)
- ✅ Testing Infrastructure (Phase 1.2)
- 🔄 Database Architecture (Phase 1.3)
- 🔄 Authentication System (Phase 1.4)

### Priority Tasks
1. Optimize RLS implementation for complex queries
2. Stress test and fine-tune connection pooling
3. Implement robust migration framework
4. Verify data integrity across all tables

### Phase 1: Technical Foundation

#### Phase 1.3: Database Architecture
- Optimize RLS implementation for complex queries
- Stress test and fine-tune connection pooling
- Implement robust migration framework
- Verify data integrity across all tables

#### Phase 1.4: Authentication System
- Replace service role workaround with proper implementation
- Enhance user session handling
- Implement secure token management
- Add proper role-based access controls

### Phase 2: Feature Completion

#### Phase 2.1: API Layer Completion
- Finish implementing remaining API endpoints
- Standardize error handling across endpoints
- Implement rate limiting for production readiness
- Add comprehensive API documentation

#### Phase 2.2: Core Functionality
- Complete predictive models implementation
- Finalize protocol designer functionality
- Enhance export/sharing capabilities
- Implement integration with external systems

#### Phase 2.3: User Interface
- Improve UI responsiveness for all screen sizes
- Complete molecular visualization features
- Implement accessibility standards
- Enhance user experience workflows

### Phase 3: Production Readiness

#### Phase 3.1: Deployment Infrastructure
- Complete CI/CD pipeline with test automation
- Optimize Docker configuration for production
- Standardize environment configuration
- Implement blue/green deployment strategy

#### Phase 3.2: Monitoring and Maintenance
- Implement centralized logging system
- Set up performance monitoring and alerting
- Configure scheduled backups
- Create maintenance runbooks

#### Phase 3.3: Security
- Conduct comprehensive security audit
- Add scanning for vulnerable dependencies
- Enhance data encryption for sensitive information
- Implement security best practices

### Phase 4: Documentation and Knowledge Transfer

#### Phase 4.1: Documentation
- Complete user documentation
- Create operations guide for administrators
- Finalize developer documentation
- Update API documentation

#### Phase 4.2: Knowledge Transfer
- Create onboarding materials
- Conduct handover sessions
- Document known issues and workarounds
- Create video tutorials for key workflows