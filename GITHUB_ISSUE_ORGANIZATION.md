# GitHub Issue Organization Plan

## Current Analysis

Our GitHub issue backlog has several critical organizational problems:

1. **Missing descriptions** - Many issues have only titles with no detailed descriptions
2. **Duplicate and fragmented issues** - Multiple issues covering the same topic (e.g., ChEMBL imports)
3. **Inconsistent labeling** - Many issues lack proper type/status labels
4. **No milestone tracking** - Issues aren't grouped into logical milestones
5. **Unclear dependencies** - Dependencies between issues aren't clearly marked
6. **Validation issues** - Separate validation issues instead of using PR reviews

## Proposed Organization Structure

### 1. Standardized Issue Templates with Required Fields

Create templates with required fields to ensure complete information:
- **Description** - Detailed explanation of the issue/task
- **Acceptance Criteria** - Clear definition of what "done" means
- **Dependencies** - Explicit listing of issues this depends on
- **Technical Details** - Relevant technical information
- **Implementation Notes** - Guidance for implementation

### 2. Core Issue Categories

Organize issues into these main categories using labels:
- `area:database` - Database-related issues
- `area:api` - API-related issues
- `area:documentation` - Documentation-related issues
- `area:ui` - UI-related issues
- `area:testing` - Testing infrastructure issues
- `area:devops` - Development/deployment related issues

### 3. Milestone Structure

Create strategic milestones based on the project phases:
- Phase 1: Database Architecture & Foundation
- Phase 2: API Implementation & Core Features
- Phase 3: Deployment & Integration
- Phase 4: Documentation & Testing

### 4. Issue Dependency Management

- Use explicit "#XX depends on #YY" syntax in issue descriptions
- Use the "subissue" label for tasks that are part of a larger issue
- Create parent tracking issues for large features with multiple components

### 5. Status Tracking

Consolidate status labels to:
- `status:planning` - In planning/design phase
- `status:ready` - Ready for implementation
- `status:in-progress` - Currently being worked on
- `status:needs-review` - Implementation complete, needs review
- `status:blocked` - Blocked on some dependency
- `status:completed` - Fully completed and validated

### 6. Priority System

Add priority labels:
- `priority:critical` - Must be addressed immediately
- `priority:high` - Important for current milestone
- `priority:medium` - Should be addressed in current milestone
- `priority:low` - Can be deferred to later milestone

## Immediate Action Items

1. **Add Missing Descriptions to Existing Issues**
   - Identify all issues with insufficient descriptions
   - Add detailed descriptions including context, requirements, acceptance criteria
   - Update with current status information

2. **Consolidate ChEMBL/PubChem Import Issues**
   - Keep issue #243 (Unified Import Pipeline) as the master issue
   - Close or consolidate duplicate/related issues after transferring relevant information
   - Ensure #243 has proper labels and is assigned to an appropriate milestone

3. **Database Connection Issues**
   - Review and consolidate database connection related issues
   - Ensure all have proper area and status labels
   - Update with current status and dependencies

4. **Create Missing Milestones**
   - Create milestone structure based on project phases
   - Assign existing issues to appropriate milestones

5. **Update Issue Templates**
   - Create standardized templates for new issues
   - Include sections for: description, acceptance criteria, technical details, dependencies
   - Make critical fields required

6. **Improve Project Board**
   - Configure GitHub project board to use the new label system
   - Create views filtered by area, status, and priority
   - Set up automation rules based on labels and PR status

## Standardized Label Usage

### Area Labels
- `area:database` - Database schema, connections, queries
- `area:api` - API endpoints, authentication, validation
- `area:ui` - User interface components and design
- `area:docs` - Documentation updates
- `area:testing` - Test infrastructure and implementation
- `area:devops` - CI/CD, deployment, infrastructure

### Type Labels
- `type:feature` - New functionality
- `type:bugfix` - Bug fixes
- `type:refactor` - Code improvements without functional changes
- `type:chore` - Maintenance tasks
- `type:documentation` - Documentation updates

### Status Labels
- `status:planning` - In design/planning phase
- `status:ready` - Ready for implementation
- `status:in-progress` - Currently being worked on
- `status:needs-review` - Ready for review
- `status:blocked` - Blocked by dependencies
- `status:completed` - Fully implemented and verified

### Priority Labels
- `priority:critical` - Urgent issues blocking progress
- `priority:high` - Important for current milestone
- `priority:medium` - Should be addressed in current sprint
- `priority:low` - Nice to have

## Issue Management Guidelines

1. **Complete Descriptions** - Every issue must have a detailed description
2. **One Issue, One Purpose** - Each issue should have a single, clear purpose
3. **Explicit Dependencies** - Clearly state what the issue depends on
4. **Regular Updates** - Update status and add comments on progress
5. **Appropriate Assignment** - Assign to the right person or leave unassigned
6. **Proper Closure** - Close with a reference to the PR that resolved it

## Implementation Timeline

1. **Week 1: Setup and Organization**
   - Create label structure
   - Set up milestones
   - Create issue templates with required fields
   - Configure project board

2. **Week 2: Description Completion and Issue Cleanup**
   - Add/complete missing descriptions for existing issues
   - Consolidate duplicate issues
   - Apply correct labels to existing issues
   - Assign issues to milestones

3. **Week 3: Documentation and Communication**
   - Document the new system
   - Train team on proper issue management
   - Set up automation where possible
   - Review and adjust as needed

## Script for Adding Missing Descriptions

Create a script to help identify issues with missing or minimal descriptions:

```bash
#!/bin/bash
# Script to identify GitHub issues with missing or minimal descriptions

# Get all open issues
gh issue list --limit 100 --json number,title,body > issues.json

# Filter issues with missing or minimal descriptions
jq -r '.[] | select(.body == null or .body == "" or .body | test("^#\\s+[Tt]ask")) | "\(.number): \(.title)"' issues.json > missing_descriptions.txt

echo "Issues with missing or minimal descriptions:"
cat missing_descriptions.txt
```