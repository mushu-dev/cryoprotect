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

### Context7 MCP Server ✅
- **Purpose**: Provides up-to-date documentation access and context-aware assistance
- **Status**: Functioning correctly
- **Usage**: Use for library version issues, API documentation, and syntax updates
- **Commands**: 
  - `resolve-library-id`: Find the exact Context7-compatible library ID
  - `get-library-docs`: Retrieve comprehensive documentation including code examples
- **Test Results**: Successfully resolved and retrieved documentation for React, including:
  - Modern React patterns with hooks
  - Functional components
  - State management techniques
  - Context API usage
  - Client component patterns
- **Best Practices**: 
  - Check Context7 when encountering versioning issues
  - Use for real-time documentation validation
  - Leverage for API specification updates

### Supabase MCP Server ✅
- **Purpose**: Direct integration with Supabase database operations
- **Status**: Functioning correctly
- **Usage**: Use for database queries, schema management, and RLS policy testing
- **Commands**:
  - `execute_sql`: Run SQL queries directly
  - `list_tables`: View database schema
  - `get_project`: Access project configuration
  - `apply_migration`: Apply database migrations
  - `list_organizations`: List available Supabase organizations
- **Test Results**: Successfully listed available organizations
- **Best Practices**:
  - Test RLS policies using execute_sql
  - Verify database changes before implementing in code
  - Use for real-time database state inspection

### Sequential-thinking MCP Server ✅
- **Purpose**: Provides structured, step-by-step reasoning for complex problems
- **Status**: Functioning correctly
- **Usage**: Use for breaking down complex problems, planning, and tool recommendation
- **Commands**:
  - This is an integrated thinking process rather than individual commands
- **Test Results**: Successfully demonstrated:
  - Progressive, step-by-step reasoning
  - Tool recommendation with confidence scores and rationale
  - Tracking of executed steps and results verification
  - Structured problem-solving approach
- **Best Practices**:
  - Use for complex algorithmic problems
  - Leverage for planning multi-step implementations
  - Apply when debugging requires systematic reasoning

### Playwright MCP Server ✅
- **Purpose**: Browser automation for testing and web scraping
- **Status**: Fully functional with containerized solution
- **Usage**: Use for UI testing, web scraping, and generating screenshots
- **Linux Compatibility**: On Fedora/newer Linux distributions, our containerized solution ensures full compatibility:
  ```bash
  # The MCP command is only accessible through Claude Code CLI
  # But you can use our container solution for testing:
  ./mcp-playwright-final.sh browser_navigate https://example.com
  ./mcp-playwright-final.sh browser_take_screenshot https://example.com screenshot.png
  ./mcp-playwright-final.sh browser_snapshot https://example.com
  ```
- **Commands**:
  - `browser_navigate`: Navigate to a URL
  - `browser_take_screenshot`: Capture screenshots of pages
  - `browser_click`: Interact with page elements
  - `browser_type`: Input text into forms
  - `browser_snapshot`: Get accessibility snapshot of a page
- **Implementation**:
  We've created a robust containerized solution that works on any Linux system:
  - Uses official Microsoft Playwright Docker image
  - No system-wide changes required
  - Handles all library compatibility issues automatically
  - Provides JSON-formatted output for easy parsing
  - Supports all MCP Playwright commands
- **Multiple Solutions Available**:
  - **mcp-playwright-final.sh**: Production-ready solution using Docker/Podman container
  - **mcp-playwright-direct.sh**: Alternative implementation with direct container access
  - **mcp-cli.sh**: Original symbolic link approach for systems without containers
- **Best Practices**:
  - Use the containerized solution for maximum compatibility
  - Automated tests ensure reliability across environments
  - The integration is transparent to Claude Code MCP
  - Container starts on-demand and has minimal resource usage
  - Use browser snapshot for interactive analysis rather than screenshots

## MCP Testing and Verification

### Testing Tools
- **Automated test script**: `test_mcp_services.js` - tests Playwright compatibility and library setup
- **Manual test script**: `test_playwright.js` - simple Playwright browser test
- **MCP integration test**: `test_mcp_playwright.js` - tests the compatibility layer with MCP-style operations 
- **Compatibility wrapper**: `mcp-cli.sh` - wrapper script for running commands with proper library paths
- **Library setup script**: `setup_mcp_playwright.sh` - sets up the symbolic links needed for compatibility

### How to Verify MCP Services
1. Run the automated test scripts:
   ```bash
   # Test the containerized solution (recommended)
   node test-final-playwright.js
   
   # Test direct Playwright solution 
   node test-direct-playwright.js
   
   # Test original symbolic link solution
   ./mcp-cli.sh node test_mcp_services.js
   ```
2. Check test results in the test output
3. For manual verification of specific services within Claude Code, use the Task tool to test:
   - Context7: Resolve documentation for a common library
   - Supabase: List organizations or tables
   - Sequential-thinking: Run a simple reasoning task
   - Playwright: Try browser navigation and screenshots

### Testing Each MCP Playwright Command
You can test each MCP Playwright command individually:

```bash
# Check status of Playwright container
./mcp-playwright-final.sh status

# Navigate to a URL and print information
./mcp-playwright-final.sh browser_navigate https://example.com

# Take a screenshot of a website
./mcp-playwright-final.sh browser_take_screenshot https://example.com test.png

# Get accessibility snapshot of a website
./mcp-playwright-final.sh browser_snapshot https://example.com
```

### Handling Issues
- **Container issues**: 
  - The solution automatically pulls the container image if needed
  - Make sure podman/docker is installed on your system
  - Verify you have permissions to create containers

- **Compatibility issues**:
  - The container approach eliminates most compatibility issues
  - If using symbolic links, run `./setup_mcp_playwright.sh` to recreate links
  - For browser install errors, manually run `npx playwright install` in the container

- **MCP service unavailability**:
  - Verify status with `mcp` command in Claude Code CLI
  - Check container status with `podman ps` or `docker ps`
  - Check if appropriate ports are open if needed

- **Performance concerns**:
  - Container launches quickly on demand
  - First launch may be slower due to image pulling
  - Consider using a pre-warmed container for faster startup

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
- **Test container Playwright solution**: `node test-final-playwright.js`
- **Run containerized Playwright commands**: `./mcp-playwright-final.sh <command> [arguments]`
- **Set up alternative compatibility layer**: `./setup_mcp_playwright.sh` (fallback option)
- **Test original compatibility layer**: `./mcp-cli.sh node test_mcp_playwright.js` (fallback option)
- **Check container status**: `./mcp-playwright-final.sh status`

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

#### Phase 3.1: Core Resiliency
- Implement API timeout and retry mechanism with exponential backoff
- Add circuit breaker pattern to prevent cascading failures
- Optimize database connection pooling
- Implement rate limiting for all API endpoints
- Create data validation schemas for all endpoints

#### Phase 3.2: Monitoring and Observability
- Implement enhanced structured logging system
- Create comprehensive error tracking system
- Set up performance monitoring and profiling
- Configure alerting system for service disruptions
- Implement frontend error handling

#### Phase 3.3: Deployment and Availability
- Create comprehensive CI/CD pipeline
- Implement blue/green deployment strategy
- Set up automated smoke tests
- Configure scheduled backups and restore procedures
- Implement graceful degradation strategy

#### Phase 3.4: Security and Compliance
- Conduct comprehensive security audit
- Add scanning for vulnerable dependencies
- Enhance data encryption for sensitive information
- Implement API versioning strategy
- Set up secure environment configuration management

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

## Production Readiness Implementation Roadmap

This section outlines our comprehensive production readiness implementation plan with a focus on building a resilient, observable, and maintainable system.

### Core Principles

- **Defense in Depth**: Multiple layers of protection and fallbacks
- **Observability First**: Comprehensive monitoring and logging
- **Graceful Degradation**: Maintain core functionality even during partial failures
- **Progressive Enhancement**: Implement features incrementally
- **Automated Operations**: Minimize manual intervention for routine tasks

### Implementation Priorities

1. **Resiliency Mechanisms**
   - Implement timeout/retry pattern for all external service calls
   - Add circuit breaker pattern to prevent cascading failures
   - Create fallback mechanisms for critical services
   - Ensure proper connection management with pooling

2. **Monitoring and Alerting**
   - Enhanced structured logging with correlation IDs
   - Comprehensive error tracking with deduplication
   - Performance monitoring with automatic profiling for slow operations
   - Multi-channel alerting with appropriate severity levels

3. **Deployment and Testing**
   - Blue-green deployment for zero downtime updates
   - Automated smoke tests to verify deployments
   - Database backup and restore procedures
   - Comprehensive CI/CD pipeline with security scans

4. **Data Integrity and Security**
   - Schema validation for all data operations
   - Rate limiting to prevent abuse
   - API versioning strategy
   - Secure configuration management

### Component Dependencies

```
┌───────────────────┐     ┌───────────────────┐
│ Feature Flags     │────>│ API Versioning    │
└───────────────────┘     └───────────────────┘
        │                         ▲
        ▼                         │
┌───────────────────┐     ┌───────────────────┐
│ Error Handling    │────>│ Graceful          │
│ & Retry           │     │ Degradation       │
└───────────────────┘     └───────────────────┘
        │                         ▲
        ▼                         │
┌───────────────────┐     ┌───────────────────┐
│ Circuit Breaking  │────>│ Monitoring &      │
│                   │     │ Alerting          │
└───────────────────┘     └───────────────────┘
        │                         ▲
        ▼                         │
┌───────────────────┐     ┌───────────────────┐
│ Rate Limiting     │────>│ Deployment        │
│                   │     │ Pipeline          │
└───────────────────┘     └───────────────────┘
```

### Testing Strategy

For each component:
1. **Unit Tests**: Verify component behavior in isolation
2. **Integration Tests**: Verify component interactions
3. **Chaos Tests**: Deliberately cause failures to test resilience
4. **Load Tests**: Verify performance under load
5. **Smoke Tests**: Quick verification of critical paths