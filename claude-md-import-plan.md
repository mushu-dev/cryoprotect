# Claude.md Import System Implementation Plan

## Overview

This plan outlines the implementation of file import capabilities for CLAUDE.md files, allowing modular organization of context information for Claude Code.

## Implementation Approach

### 1. Syntax Definition

We'll implement the simple `@path/to/file.md` syntax with these features:
- Basic import: `@api/claude.md`
- Relative path import: `@./database/claude.md`
- Absolute path import: `@/home/mushu/Projects/CryoProtect/api/claude.md`

### 2. Import Structure

The import system will:
- Process imports at Claude Code launch time
- Include content from referenced files inline
- Support nested imports (files can import other files)
- Preserve headings hierarchy from imported files
- Maintain proper markdown formatting

### 3. Module-Specific Files

Each major directory will have its own claude.md with targeted information:
- `/api/claude.md`: API architecture, endpoints, patterns
- `/database/claude.md`: Database schemas, connection patterns, RLS policies
- `/chembl/claude.md`: ChEMBL integration, rate limiting, caching
- `/pubchem/claude.md`: PubChem client, data standardization, import process
- `/monitoring/claude.md`: Logging, metrics, alerting
- `/security/claude.md`: Auth systems, token management, RBAC
- `/tests/claude.md`: Testing frameworks, CI/CD, test data
- `/migrations/claude.md`: Migration patterns, versioning, Supabase integration
- `/docs/claude.md`: Documentation standards and organization

### 4. Main CLAUDE.md Structure

The main CLAUDE.md will be restructured to:
- Provide project-level overview
- Import module-specific files
- Focus on cross-cutting concerns
- Reduce redundancy

## Implementation Steps

1. **Update Module-Specific Files**
   - Review and standardize each module's claude.md
   - Ensure they have consistent heading structure
   - Add module-specific guidance and patterns

2. **Modify Main CLAUDE.md**
   - Add import directives for all module-specific files
   - Remove redundant content now covered in module files
   - Add section explaining the import system

3. **Testing**
   - Test with Claude Code to ensure imports are processed correctly
   - Verify heading levels are maintained properly
   - Check that nested imports work as expected

4. **Documentation**
   - Create a dedicated section in CLAUDE.md explaining the import system
   - Document best practices for maintaining modular documentation
   - Provide examples for adding new module-specific files

## Proposed CLAUDE.md Structure

```markdown
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Claude.md Import System

CLAUDE.md files can now import other files using the `@path/to/file.md` syntax. This allows modular organization of context information.

## Project Overview
[Project overview content]

## Module-Specific Guidance

### API Layer
@api/claude.md

### Database Layer
@database/claude.md

### ChEMBL Integration
@chembl/claude.md

### PubChem Integration
@pubchem/claude.md

### Monitoring
@monitoring/claude.md

### Security
@security/claude.md

### Testing
@tests/claude.md

### Migrations
@migrations/claude.md

### Documentation
@docs/claude.md

## Cross-Cutting Concerns
[Content relevant across all modules]
```

## Benefits

- **Modularity**: Changes to module-specific guidance don't require modifying the main CLAUDE.md
- **Maintenance**: SMEs can update their module's guidance independently
- **Clarity**: Clearer separation of concerns in documentation
- **Scalability**: Easier to add new modules and guidance
- **Context Management**: Better token usage by loading only relevant content when needed

## Risks and Mitigations

- **Path Resolution**: Ensure relative paths work consistently
- **Circular Imports**: Implement detection to prevent infinite loops
- **Heading Levels**: May need to adjust heading levels in imported files
- **Discoverability**: Need clear documentation on the import system