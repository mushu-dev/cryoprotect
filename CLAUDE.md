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
- **Best Practices**: 
  - Check Context7 when encountering React 18 and Next.js 14 migration issues
  - Use for Convex Client API documentation validation
  - Leverage for Mol* library implementation guidance
  - Consult for shadcn/ui component integration patterns

### Sequential-thinking MCP Server ✅
- **Purpose**: Provides structured, step-by-step reasoning for complex problems
- **Status**: Functioning correctly
- **Usage**: Use for breaking down complex problems, planning, and tool recommendation
- **Best Practices**:
  - Use for complex implementation planning (especially resiliency patterns)
  - Leverage for React 18 migration planning
  - Apply when implementing Mol* library integration
  - Use when implementing circuit breaker patterns in frontend

### Playwright MCP Server ✅
- **Purpose**: Browser automation for testing and web scraping
- **Status**: Fully functional with containerized solution
- **Usage**: Use for UI testing, particularly for experimental data visualization features
- **Linux Compatibility**: Use our containerized solution for maximum compatibility
- **Commands**:
  - `browser_navigate`: Navigate to protocol and experiment pages for testing
  - `browser_take_screenshot`: Capture screenshots for visual regression tests
  - `browser_click`: Test UI interactions with Protocol Step Editor
  - `browser_type`: Test form inputs in Protocol and Experiment forms
  - `browser_snapshot`: Test accessibility of new UI components
- **Best Practices**:
  - Use the containerized solution (`mcp-playwright-final.sh`) for all tests
  - Create automated tests for all critical protocol workflows
  - Test Mol* visualization across different browsers
  - Verify accessibility standards compliance with snapshots

## Technology Stack Overview

### Frontend
- **Framework**: Next.js (migrating to version 14.x)
- **UI Library**: React (migrating to version 18.x)
- **Component Library**: shadcn/ui (based on Radix UI primitives)
- **Styling**: Tailwind CSS
- **State Management**: Combination of React Query and Zustand
- **Visualization**: Mol* library for 3D molecule visualization
- **Database Client**: Convex Client (for real-time data)
- **Testing**: Jest for unit tests, Playwright for E2E tests

### Backend
- **API**: Flask API with resiliency patterns
- **Database**: Convex
- **Deployment**: Netlify for frontend, Heroku for API
- **Specialized Services**: RDKit service deployed on Fly.io
- **Container Technology**: Docker/Podman for service isolation
- **API Endpoints**:
  - Main API (Heroku): `https://cryoprotect-8030e4025428.herokuapp.com/api`
  - RDKit Service (Fly.io): `https://cryoprotect-rdkit.fly.dev`
  - Convex API: `https://upbeat-parrot-866.convex.cloud`
- **Working Endpoints**:
  - Get all molecules: `GET /api/molecules`
  - Get molecule by ID: `GET /api/molecules/{id}`

## Build/Run/Test Commands

### Development
- **Run frontend development**: `cd frontend && npm run dev` 
- **Run frontend with Convex**: `cd frontend && npm run dev:with-convex`
- **Build frontend**: `cd frontend && npm run build`
- **Run API locally**: `python app.py`
- **Run RDKit service**: `./run_rdkit_service.sh`

### Testing
- **Run all tests**: `cd frontend && npm run test:all`
- **Run end-to-end tests**: `cd frontend && npm run test:e2e`
- **Run experimental data tests**: `cd frontend && npm run test:experimental`
- **Run protein visualizer tests**: `cd frontend && npm run test:protein-visualizer`
- **Run API mock tests**: `cd frontend && npm run test:api:verify`
- **Run API live tests**: `cd frontend && npm run test:api:live`
- **Run API tests (Python)**: `python -m pytest tests/api/`
- **Run resiliency tests**: `python examples/resiliency_demo.py`
- **Run protocol template tests**: `cd frontend && npm run test:protocol-templates`

### Deployment
- **Deploy to Netlify**: 
  - Automated: `git push origin netlify-autodeploy` (pushes to the netlify-autodeploy branch)
  - Manual: `cd minimal-frontend && netlify deploy --prod` (using Netlify CLI)
- **Deploy to Heroku**: `git push heroku main`
- **Deploy RDKit service**: `flyctl deploy`

### Playwright Containers
- **Run containerized Playwright tests**: `cd frontend && npx playwright test`
- **Test container Playwright solution**: `node test-final-playwright.js`
- **Run containerized Playwright commands**: `./mcp-playwright-final.sh <command> [arguments]`
- **Check container status**: `./mcp-playwright-final.sh status`

## Code Style Guidelines

- **Components**: Functional components with hooks (no class components)
- **Docstrings**: JSDoc comments for all modules, components, and functions
- **Imports**: Group imports (React, third-party, local) with empty line between groups
- **Formatting**: 2-space indentation for JS/TS, 4-space for Python, 80-character line length
- **Typing**: TypeScript types for all components, functions, and state
- **Naming**: camelCase for variables/functions, PascalCase for components/classes, UPPER_CASE for constants
- **Error handling**: React error boundaries for components, try/catch with specific error types
- **Logging**: Structured logging with appropriate log levels (using our new logging system)
- **Documentation**: Include detailed JSDoc with parameters and return values
- **Testing**: Component tests, integration tests, and E2E tests with Playwright

## Project Completion Plan

### Current Phase: Production-Ready Convex Integration 🚀

**PRIORITY: Full production deployment with real scientific data in 4 days**

#### Phase 1: Database Migration & Population (Days 1-2) 🔄
- ✅ Analyze existing ChEMBL/PubChem population scripts
- 🔄 Create enhanced ConvexAdapter for production workloads
- 🔄 Convert Supabase-based scripts to Convex-compatible versions
- 🔄 Execute full ChEMBL dataset population (5000+ molecules)
- 🔄 Populate PubChem naming and molecular properties
- 🔄 Implement data validation and integrity checks

#### Phase 2: Integration & Calculations (Days 2-3) 🔄
- 🔄 Enhance RDKit service integration with Convex backend
- 🔄 Update Flask API to use ConvexAdapter exclusively
- 🔄 Implement calculation result caching in Convex
- 🔄 Add comprehensive error handling and retry logic
- 🔄 Test end-to-end calculation workflow

#### Phase 3: Frontend & Real-time Features (Days 3-4) 🔄
- 🔄 Connect frontend to Convex for real-time data
- 🔄 Implement live molecule search and property display
- 🔄 Add real-time calculation updates
- 🔄 Enable collaborative experiment editing
- 🔄 Optimize performance for production loads

#### Phase 4: Production Deployment & Testing (Day 4) 🔄
- 🔄 Coordinate deployments across Netlify, Heroku, Fly.io
- 🔄 Execute comprehensive Playwright testing
- 🔄 Implement monitoring and error tracking
- 🔄 Verify system reliability and performance
- ✅ **DELIVERABLE**: Production-ready scientific research tool

### Priority Tasks (Production Ready in 4 Days)

1. **🚨 CRITICAL: Enhanced ConvexAdapter Development**
   - Create production-ready ConvexAdapter with batch operations
   - Implement connection pooling and retry logic  
   - Add transaction support and error handling
   - Enable bulk inserts for population scripts

2. **🚨 CRITICAL: Database Population Scripts**
   - Convert ChEMBL_CryoProtectants_Supabase.py → ChEMBL_CryoProtectants_Convex.py
   - Convert PubChem_CryoProtectants_Supabase.py → PubChem_CryoProtectants_Convex.py
   - Execute full scientific dataset population (5000+ molecules)
   - Implement progress tracking and resumability

3. **🚨 CRITICAL: RDKit-Convex Integration**
   - Create RDKit calculation service bridge to Convex
   - Implement result caching and error handling
   - Test molecular property calculations end-to-end
   - Optimize calculation performance

4. **🔥 HIGH: Flask API Convex Migration**
   - Replace all Supabase calls with ConvexAdapter
   - Add real-time subscriptions for data changes
   - Implement proper authentication flow
   - Add comprehensive API error handling

5. **🔥 HIGH: Frontend Real-time Features**
   - Enable live molecule search with Convex subscriptions
   - Implement real-time calculation updates
   - Add collaborative experiment editing
   - Optimize frontend performance for production

6. **⚡ MEDIUM: Production Deployment**
   - Coordinate deployments across all services
   - Execute comprehensive Playwright testing
   - Implement monitoring and error tracking
   - Verify end-to-end system functionality

## Testing Strategy

### Comprehensive Testing Approach
We've implemented a multi-layered testing strategy:

1. **Unit Tests**:
   - Test individual components in isolation
   - Mock dependencies and external services
   - Target 90%+ code coverage for critical paths

2. **Integration Tests**:
   - Test component interactions
   - Verify service integration points
   - Test data flows across components

3. **End-to-End Tests**:
   - Test complete user workflows with Playwright
   - Verify critical paths in production-like environment
   - Include visual regression tests

4. **Resilience Tests**:
   - Test system behavior under failure conditions
   - Validate retry mechanisms and circuit breakers
   - Verify graceful degradation

5. **Accessibility Tests**:
   - Verify WCAG compliance across components
   - Test keyboard navigation
   - Verify screen reader compatibility

## Frontend Architecture

### Key Components

```
┌────────────────────────────────────────────────────┐
│ Next.js Application                                │
│                                                    │
│  ┌────────────────────┐    ┌────────────────────┐  │
│  │ UI Components      │    │ Data Layer         │  │
│  │                    │    │                    │  │
│  │  • shadcn/ui       │    │  • Convex Client   │  │
│  │  • Tailwind CSS    │    │  • React Query     │  │
│  │  • Radix UI        │    │  • Zustand         │  │
│  └────────────────────┘    └────────────────────┘  │
│                                                    │
│  ┌────────────────────┐    ┌────────────────────┐  │
│  │ Domain Components  │    │ Visualization      │  │
│  │                    │    │                    │  │
│  │  • Experiments     │    │  • MolstarViewer   │  │
│  │  • Protocols       │    │  • Data Charts     │  │
│  │  • Molecules       │    │  • Interactive UI  │  │
│  └────────────────────┘    └────────────────────┘  │
│                                                    │
│  ┌────────────────────┐    ┌────────────────────┐  │
│  │ Resiliency Layer   │    │ Monitoring         │  │
│  │                    │    │                    │  │
│  │  • Retry Mechanism │    │  • Error Tracking  │  │
│  │  • Circuit Breaker │    │  • Performance     │  │
│  │  • Timeout Handling│    │  • Logging         │  │
│  └────────────────────┘    └────────────────────┘  │
│                                                    │
└────────────────────────────────────────────────────┘
```

### Database Integration (Convex)

Convex serves as our primary database with these key capabilities:

1. **Real-time Data**: Automatic subscriptions to data changes with reactive UI updates
2. **Type Safety**: Full TypeScript integration with generated types from schema
3. **Authentication**: Integrated Auth with role-based access and JWT support
4. **Functions**: Backend functions with secure database access and validation
5. **Performance**: Optimized for React applications with efficient caching
6. **Resiliency**: Connection pooling, retry mechanisms, and error handling
7. **Collaborative Features**: Real-time presence, conflict resolution, and shared editing
8. **Offline Support**: Works offline with automatic synchronization when reconnected
9. **Primary Database**: Now the exclusive database for this project

#### Convex Connection Architecture

```
┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
│ Frontend            │    │ Middleware          │    │ Backend             │
│                     │    │                     │    │                     │
│ ┌─────────────────┐ │    │ ┌─────────────────┐ │    │ ┌─────────────────┐ │
│ │ React Components│ │    │ │ AuthBridge      │ │    │ │ ConvexAdapter   │ │
│ └────────┬────────┘ │    │ └────────┬────────┘ │    │ └────────┬────────┘ │
│          │          │    │          │          │    │          │          │
│ ┌────────┴────────┐ │    │ ┌────────┴────────┐ │    │ ┌────────┴────────┐ │
│ │ Convex Hooks    │ │    │ │ JWT Validator   │ │    │ │ DB Factory      │ │
│ └────────┬────────┘ │    │ └────────┬────────┘ │    │ └────────┬────────┘ │
│          │          │    │          │          │    │          │          │
│ ┌────────┴────────┐ │    │ ┌────────┴────────┐ │    │ ┌────────┴────────┐ │
│ │ Convex Client   │◄├────┼─┤ Sync Manager    │◄├────┼─┤ Flask API       │ │
│ └─────────────────┘ │    │ └─────────────────┘ │    │ └─────────────────┘ │
└─────────────────────┘    └─────────────────────┘    └─────────────────────┘
```

## Experimental Data Enhancement

The experimental data enhancement system is a core feature being developed in the current phase. Key components include:

1. **Protocol Designer**:
   - Comprehensive step editor with advanced features
   - Equipment and parameter management
   - Alert system for critical steps
   - Version management and comparison tools

2. **Experiment Tracking**:
   - Data collection and visualization
   - Results analysis tools
   - Cross-reference with molecular structures
   - Integration with external data sources

3. **Integration Points**:
   - Molecular visualization with MolstarViewer
   - Protocol execution tracking
   - Data export and sharing capabilities
   - Laboratory equipment integration

4. **Protocol Template Management**:
   - Template creation and editing
   - Version control and comparison
   - Template to protocol conversion
   - Access control for sharing templates
   - Component routes:
     - `/protocols/templates`: Template management
     - `/protocols/templates/create`: Create template
     - `/protocols/templates/[id]`: Edit template
     - `/protocols/templates/versions/[id]`: Version comparison

5. **Current Implementation Status**:
   - ✅ Enhanced Protocol Step Editor is completed
   - ✅ Lab Verification System is implemented
   - ✅ Protocol Template Management System is implemented
   - ✅ Time Series Visualization is completed
   - 🔄 Advanced Uncertainty Quantification is pending
   - 🔄 Direct connection to Convex is implemented

## Workflow Optimization Guidelines

### API Cost Management
- Prefer delta context mode: Only analyze changed files when possible
- Use chunking: Focus on one directory/component at a time
- Minimize document generation: Don't create unnecessary documentation files
- Use precise file:line references for targeted changes
- Leverage MCP servers for documentation and database queries to reduce token usage

### Task Management Approach
- Handle simple tasks directly through Claude Code CLI
- Delegate complex tasks to Roo agent team through Master Orchestrator
- Provide detailed instructions with file paths and line numbers
- Track completion status of each delegated task
- Use Context7 for documentation-related tasks instead of regenerating documentation

### Production Deployment Strategy

#### Netlify Deployment Configuration
- **Automated Deployments**: 
  - Push to `netlify-autodeploy` branch to trigger automatic deployments
  - Netlify is configured to deploy the `minimal-frontend` directory
  - Deployment settings are defined in the repo's `netlify.toml` and `minimal-frontend/netlify.toml`
  - Git submodules are disabled with `submodules = false` in the config
  - Content Security Policy is configured to allow connections to API endpoints
  - Environment variables set `NEXT_PUBLIC_USE_CONVEX = "false"` for the minimal frontend

- **Troubleshooting Netlify Deployments**:
  - Check Netlify build logs for errors
  - Ensure no references to Convex client in minimal-frontend
  - Verify all Next.js pages can be built statically
  - Use redirects for API endpoints (configured in netlify.toml)
  - For manual deployment, use `netlify deploy --prod` command

#### Other Deployment Strategies
- Deploy backend API to Heroku with proper scaling configuration
- Use Fly.io for specialized services requiring custom containerization
- Implement proper monitoring and alerting for all deployed services
- Use blue/green deployment for zero-downtime updates

## Conclusion

This guidance document aligns with our production environment and workflow goals. It reflects our current technology choices, priorities, and implementation plans. As we progress through the phases of the project, we will continue to update this document to ensure it remains accurate and valuable for development work.