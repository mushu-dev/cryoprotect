#!/bin/bash

# Setup Experiment Enhancement Project
# This script creates GitHub issues for the Experimental Data Enhancement project

# Make sure gh CLI is available
if ! command -v gh &> /dev/null; then
    echo "GitHub CLI (gh) is not installed. Please install it and try again."
    exit 1
fi

# Check if authenticated with GitHub
if ! gh auth status &> /dev/null; then
    echo "Not authenticated with GitHub. Please run 'gh auth login' first."
    exit 1
fi

# Create a new milestone for the project
MILESTONE_RESPONSE=$(gh api \
  --method POST \
  -H "Accept: application/vnd.github+json" \
  -H "X-GitHub-Api-Version: 2022-11-28" \
  /repos/$(gh repo view --json nameWithOwner -q .nameWithOwner)/milestones \
  -f title="Experimental Data Enhancement" \
  -f state="open" \
  -f description="Enhancement of experimental data capabilities for improved scientific accuracy and usability" \
  -f due_on="$(date -d "+4 months" "+%Y-%m-%dT%H:%M:%SZ")")

MILESTONE_NUMBER=$(echo $MILESTONE_RESPONSE | jq .number)

echo "Created milestone #$MILESTONE_NUMBER: Experimental Data Enhancement"

# Create labels if they don't exist
echo "Creating necessary labels..."
gh label create "phase-1:schema" --color "0366d6" --description "Phase 1: Core Schema and Foundation" || true
gh label create "phase-2:protocols" --color "28a745" --description "Phase 2: Protocol System & Scientific Features" || true
gh label create "phase-3:analytics" --color "6f42c1" --description "Phase 3: Advanced Analytics & Validation" || true
gh label create "phase-4:integration" --color "d73a4a" --description "Phase 4: External Integration & Advanced Features" || true
gh label create "experiments" --color "fbca04" --description "Related to experimental data system" || true
gh label create "schema" --color "1d76db" --description "Database schema changes" || true
gh label create "scientific" --color "0075ca" --description "Scientific functionality" || true
gh label create "api" --color "e4e669" --description "API-related changes" || true
gh label create "ui" --color "fef2c0" --description "User interface components" || true

# Create epic issues for each phase
echo "Creating epic issues for each phase..."

# Phase 1 Epic
PHASE1_BODY="# Phase 1: Core Schema and Foundation

This epic covers the implementation of the core schema and foundation for the experimental data enhancement project.

## Objectives
- Consolidate existing experiment schemas into a unified, comprehensive schema
- Extend database connection management for experimental data
- Implement basic API integration for experiments

## Tasks
- [ ] Design unified experimental schema
- [ ] Create database migration scripts
- [ ] Implement database adapters for experiments
- [ ] Create basic API resources for experiments
- [ ] Implement validation and error handling

## Acceptance Criteria
- Unified schema properly references existing molecules/mixtures
- Migration scripts preserve existing data
- Database adapters follow established connection patterns
- API endpoints conform to established standards
- All changes have appropriate tests

## Timeline
Estimated completion: 4-6 weeks"

gh issue create \
  --title "Epic: Phase 1 - Core Schema and Foundation" \
  --body "$PHASE1_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-1:schema" \
  --label "experiments" \
  --label "schema"

PHASE1_EPIC_NUMBER=$(gh issue list --label "phase-1:schema" --json number --jq '.[0].number')
echo "Created Phase 1 Epic: #$PHASE1_EPIC_NUMBER"

# Phase 2 Epic
PHASE2_BODY="# Phase 2: Protocol System & Scientific Features

This epic covers the implementation of the protocol system and scientific features for experimental data.

## Objectives
- Enhance protocol designer with templates and versioning
- Integrate scientific calculation capabilities
- Implement foundation for laboratory integration

## Tasks
- [ ] Extend protocol designer with template system
- [ ] Implement protocol versioning and comparison
- [ ] Add validation against physical models
- [ ] Integrate thermodynamic property calculations
- [ ] Implement uncertainty propagation
- [ ] Create equipment configuration structures
- [ ] Implement protocol export for lab equipment

## Acceptance Criteria
- Protocol designer supports templates and versioning
- Scientific calculations produce accurate results
- Uncertainty is properly propagated through calculations
- Laboratory protocols can be exported in standard formats
- All features have appropriate tests

## Timeline
Estimated completion: 4-6 weeks

## Dependencies
- Phase 1 Epic (#$PHASE1_EPIC_NUMBER)"

gh issue create \
  --title "Epic: Phase 2 - Protocol System & Scientific Features" \
  --body "$PHASE2_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-2:protocols" \
  --label "experiments" \
  --label "scientific"

PHASE2_EPIC_NUMBER=$(gh issue list --label "phase-2:protocols" --json number --jq '.[0].number')
echo "Created Phase 2 Epic: #$PHASE2_EPIC_NUMBER"

# Phase 3 Epic
PHASE3_BODY="# Phase 3: Advanced Analytics & Validation

This epic covers the implementation of advanced analytics and validation for experimental data.

## Objectives
- Create a results analytics engine for time-series and comparative analysis
- Implement a multi-level data validation framework
- Develop visualization components for data exploration

## Tasks
- [ ] Implement time-series analysis for experimental data
- [ ] Create comparative analysis tools
- [ ] Add statistical process control for quality monitoring
- [ ] Implement multi-level validation pipeline
- [ ] Add statistical outlier detection
- [ ] Create data integrity verification system
- [ ] Develop interactive data exploration components
- [ ] Implement result comparison dashboards

## Acceptance Criteria
- Analytics engine provides insightful analysis of experimental results
- Validation framework ensures data integrity
- Visualization components make data exploration intuitive
- All features have appropriate tests

## Timeline
Estimated completion: 4-6 weeks

## Dependencies
- Phase 1 Epic (#$PHASE1_EPIC_NUMBER)
- Phase 2 Epic (#$PHASE2_EPIC_NUMBER)"

gh issue create \
  --title "Epic: Phase 3 - Advanced Analytics & Validation" \
  --body "$PHASE3_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-3:analytics" \
  --label "experiments" \
  --label "scientific"

PHASE3_EPIC_NUMBER=$(gh issue list --label "phase-3:analytics" --json number --jq '.[0].number')
echo "Created Phase 3 Epic: #$PHASE3_EPIC_NUMBER"

# Phase 4 Epic
PHASE4_BODY="# Phase 4: External Integration & Advanced Features

This epic covers the implementation of external integrations and advanced features for experimental data.

## Objectives
- Connect to external scientific databases and literature sources
- Implement machine learning for protocol optimization and anomaly detection
- Add collaboration and sharing capabilities

## Tasks
- [ ] Connect to literature databases
- [ ] Implement integration with public repositories
- [ ] Add direct import from lab equipment
- [ ] Implement protocol optimization algorithms
- [ ] Add predictive models for experimental outcomes
- [ ] Create anomaly detection for results
- [ ] Add experiment sharing capabilities
- [ ] Implement data export in standard formats

## Acceptance Criteria
- System can integrate with external data sources
- Machine learning models provide valuable insights
- Collaboration features make sharing experiments easy
- All features have appropriate tests

## Timeline
Estimated completion: 4-6 weeks

## Dependencies
- Phase 1 Epic (#$PHASE1_EPIC_NUMBER)
- Phase 2 Epic (#$PHASE2_EPIC_NUMBER)
- Phase 3 Epic (#$PHASE3_EPIC_NUMBER)"

gh issue create \
  --title "Epic: Phase 4 - External Integration & Advanced Features" \
  --body "$PHASE4_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-4:integration" \
  --label "experiments" \
  --label "scientific"

PHASE4_EPIC_NUMBER=$(gh issue list --label "phase-4:integration" --json number --jq '.[0].number')
echo "Created Phase 4 Epic: #$PHASE4_EPIC_NUMBER"

# Create Phase 1 task issues
echo "Creating Phase 1 task issues..."

# Schema design task
SCHEMA_DESIGN_BODY="## Description
Create a comprehensive schema design document for the experimental data system that unifies and extends the current schema.

## Technical Details
- Base design on the schema in \`create_experimental_linkage.py\` but enhance with additional capabilities
- Ensure proper relationships with molecules and mixtures tables
- Add support for uncertainty measurements
- Include metadata fields (provenance, versioning, etc.)
- Design schema for protocol templates and versions

## Acceptance Criteria
- [ ] Schema diagram with all tables and relationships
- [ ] Detailed descriptions of all tables and fields
- [ ] Migration strategy that preserves existing data
- [ ] SQL scripts for schema creation
- [ ] Performance considerations documented

## Related Components
- Database Layer
- Data Access Layer

## Dependencies
None - this is a foundation task"

gh issue create \
  --title "[EXP] Design unified experimental schema" \
  --body "$SCHEMA_DESIGN_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-1:schema" \
  --label "experiments" \
  --label "schema" \
  --assignee "@me"

# Migration scripts task
MIGRATION_SCRIPTS_BODY="## Description
Create database migration scripts to implement the unified experimental schema while preserving existing data.

## Technical Details
- Use the established migration pattern (SQL files in migrations directory)
- Ensure scripts handle data transformation from old to new schema
- Include rollback scripts for each migration
- Add proper indexes for performance

## Acceptance Criteria
- [ ] Migration scripts for all schema changes
- [ ] Rollback scripts for reverting changes
- [ ] Successful migration with data preservation verified
- [ ] Performance impact testing completed
- [ ] Documentation of migration process

## Related Components
- Database Layer

## Dependencies
- Schema design task"

gh issue create \
  --title "[EXP] Create database migration scripts for experiment schema" \
  --body "$MIGRATION_SCRIPTS_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-1:schema" \
  --label "experiments" \
  --label "schema"

# Database adapters task
DB_ADAPTERS_BODY="## Description
Implement database adapters for the experimental data system that follow established connection patterns.

## Technical Details
- Create ExperimentAdapter class that follows existing adapter patterns
- Integrate with ConnectionManager for connection pooling
- Implement caching strategy appropriate for experimental data
- Add batch processing support for large datasets

## Acceptance Criteria
- [ ] ExperimentAdapter implemented and tested
- [ ] TissueTypeAdapter implemented and tested
- [ ] ProtocolAdapter implemented and tested
- [ ] ResultAdapter implemented and tested
- [ ] Integration with connection pooling verified
- [ ] Performance testing completed

## Related Components
- Data Access Layer

## Dependencies
- Schema design task
- Migration scripts task"

gh issue create \
  --title "[EXP] Implement database adapters for experimental data" \
  --body "$DB_ADAPTERS_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-1:schema" \
  --label "experiments" \
  --label "api"

# API resources task
API_RESOURCES_BODY="## Description
Create RESTful API resources for the experimental data system that follow established API patterns.

## Technical Details
- Implement ExperimentResource class with standard CRUD operations
- Add proper validation and error handling
- Follow standardized response format
- Implement pagination for list operations
- Document API endpoints

## Acceptance Criteria
- [ ] ExperimentResource implemented and tested
- [ ] TissueTypeResource implemented and tested
- [ ] ProtocolResource implemented and tested
- [ ] ResultResource implemented and tested
- [ ] API documentation updated
- [ ] Performance testing completed

## Related Components
- API Layer

## Dependencies
- Database adapters task"

gh issue create \
  --title "[EXP] Create API resources for experimental data" \
  --body "$API_RESOURCES_BODY" \
  --milestone "$MILESTONE_NUMBER" \
  --label "phase-1:schema" \
  --label "experiments" \
  --label "api"

echo "GitHub project setup complete!"