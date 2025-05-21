# GitHub Issue Management Guide for CryoProtect

This guide provides detailed instructions for managing GitHub issues in the CryoProtect repository. It outlines best practices, workflows, and tools to help maintain an organized and efficient development process.

## Table of Contents

1. [Epic-Based Organization](#epic-based-organization)
2. [Issue Structure](#issue-structure)
3. [Issue Templates](#issue-templates)
4. [Labels and Categorization](#labels-and-categorization)
5. [Project Board](#project-board)
6. [Issue Lifecycle](#issue-lifecycle)
7. [Automation with GitHub Actions](#automation-with-github-actions)
8. [Issue Management Scripts](#issue-management-scripts)
9. [Best Practices](#best-practices)
10. [Issue Consolidation](#issue-consolidation)
11. [Issue Comments and Communication](#issue-comments-and-communication)

## Epic-Based Organization

CryoProtect uses an epic-based organizational structure to manage issues. All work is categorized into 9 epic categories:

1. **Database Implementation and Optimization**: Database schema, migrations, queries, connection pooling, and performance optimization
2. **ChEMBL Integration and Data Pipeline**: Integration with ChEMBL database, data import processes, and ETL pipelines
3. **PubChem Integration and Molecule Management**: Integration with PubChem, molecule properties, data consolidation
4. **API Development and Standardization**: API endpoints, standardization, documentation, and gateway services
5. **Frontend Implementation and User Experience**: UI components, responsive design, user workflows, and visual styling
6. **Authentication and Security Implementation**: JWT auth, RLS policies, RBAC, and security enhancements
7. **RDKit Integration and Chemical Functionality**: RDKit services, chemical calculations, and molecular analysis
8. **Infrastructure and Deployment**: Deployment processes, container orchestration, CI/CD, and hosting configuration
9. **Testing and Quality Assurance**: Test frameworks, CI validation, and quality assurance procedures

Each epic has a master issue that tracks all related work. When creating a new issue, it should be associated with the appropriate epic by referencing the epic issue number in the "Related Epic" section.

## Issue Structure

All issues should include:

1. **Clear Title**: Brief but descriptive, with a prefix indicating type (e.g., "Bugfix:", "Feature:")
2. **Detailed Description**: What needs to be done and why
3. **Related Epic**: Which epic category the issue belongs to
4. **Acceptance Criteria**: Specific, measurable conditions that define "done"
5. **Area**: Which component of the project is affected
6. **Priority**: How important the issue is
7. **Dependencies**: Issues that must be completed first
8. **Technical Details**: Implementation guidance if relevant

### DO:
- Use the appropriate issue template for your issue type
- Include clear, specific titles that identify the issue
- Follow the format outlined in the template
- Link to related issues with `#issue-number`
- Include relevant area labels
- Set appropriate priority labels
- Add acceptance criteria that define when the issue is complete

### DON'T:
- Create duplicate issues (search first)
- Create vague "omnibus" issues covering multiple unrelated problems
- Leave required template fields blank
- Create issues titled "Test" or with obviously AI-generated content

## Issue Templates

CryoProtect uses the following issue templates:

1. **Epic**: For creating new epic issues that group related tasks
2. **Feature Request**: For new functionality or enhancements
3. **Bug Fix**: For reporting and fixing bugs
4. **Validation**: For validation tasks and testing
5. **Documentation**: For documentation requests and updates

When creating a new issue, always select the appropriate template and fill out all required sections.

## Labels and Categorization

### Area Labels
- `area:database` - Database schema, connections, queries
- `area:api` - API endpoints, authentication, validation
- `area:ui` - User interface components and interactions
- `area:auth` - Authentication, authorization, and security features
- `area:testing` - Test infrastructure and implementation
- `area:chembl` - ChEMBL integration and data handling
- `area:docs` - Documentation updates
- `area:devops` - Infrastructure setup and configuration

### Type Labels
- `type:epic` - Epic-level issue that groups related tasks
- `type:feature` - New functionality
- `type:bugfix` - Bug fixes
- `type:refactor` - Code improvements without functional changes
- `type:documentation` - Documentation updates
- `type:validation` - Testing and validation tasks

### Status Labels
- `status:planning` - In design/planning phase
- `status:ready` - Ready for implementation
- `status:in-progress` - Currently being worked on
- `status:needs-review` - Implementation complete, needs review
- `status:blocked` - Blocked on some dependency
- `status:completed` - Fully completed and validated

### Priority Labels
- `priority:high` - Critical for current milestone
- `priority:medium` - Important, should be addressed soon
- `priority:low` - Nice to have, can be deferred

### Special Labels
- `subissue` - Indicates this is a sub-issue of another issue
- `stale` - Has had no activity for 30+ days

## Project Board

The CryoProtect project uses a GitHub project board with the following columns:

1. **Todo**: Initial backlog and planning
2. **Ready**: Ready for implementation
3. **In Progress**: Currently being worked on
4. **Review**: Complete but pending review
5. **Done**: Fully implemented and validated
6. **Blocked**: Blocked on some dependency

Issues automatically move between columns based on status labels through GitHub Actions automation.

### Project Board Automation

Our project board uses automation to:
- Add new issues to the "Todo" column
- Move issues to appropriate columns based on status labels
- Add pull requests to the "In Progress" column
- Move issues to "Done" when closed

## Issue Lifecycle

1. **Creation**
   - Use the appropriate issue template
   - Add relevant labels (area, type, priority)
   - Associate with the relevant epic category
   - Add to the project board (automated)

2. **Planning**
   - Refine description and acceptance criteria
   - Update labels to `status:planning`
   - Identify dependencies

3. **Ready for Work**
   - Update to `status:ready`
   - Ensure all necessary details are provided

4. **In Progress**
   - Assign to the person working on it
   - Update to `status:in-progress`
   - Make regular updates in comments

5. **Review**
   - Update to `status:needs-review`
   - Reference pull request(s)
   - Tag reviewers in comments

6. **Completion**
   - Update to `status:completed`
   - Ensure all acceptance criteria are met
   - Close the issue

## Automation with GitHub Actions

CryoProtect uses GitHub Actions workflows to automate issue management:

1. **Project Board Assignment**: Automatically adds new issues to the project board
2. **Status Updates**: Moves issues between project board columns based on status labels
3. **Epic Suggestion**: Suggests which epic an issue belongs to
4. **Auto-labeling PRs**: Adds area labels to PRs based on changed files
5. **PR Linking**: Links PRs to issues when mentioned in commit messages
6. **Stale Issue Detection**: Identifies and labels stale issues

These automations are configured in the `.github/workflows/automated-issue-workflow.yml` file.

## Issue Management Scripts

The `.github/scripts/` directory contains several tools for issue management:

1. **setup-github-issues.sh**: Sets up all issue templates, labels, and epics
2. **consolidate-issues.sh**: Helps identify and consolidate duplicate issues
3. **create-labels.sh**: Creates all the standard labels for the repository
4. **create-epics-only.sh**: Creates just the epic issues without other setup
5. **mass-consolidation.py**: Aggressively consolidates issues into epics (use with caution)
6. **step-by-step-consolidation.sh**: More controlled issue consolidation
7. **setup-branch-protection.sh**: Sets up branch protection rules
8. **detect-ai-spam.py**: Identifies and labels AI-generated test issues

### Using the Setup Script

To set up the complete GitHub issue management system:

```bash
cd .github/scripts
chmod +x setup-github-issues.sh
./setup-github-issues.sh
```

Choose option 5 to run the complete setup.

## Best Practices

1. **One Issue, One Purpose**: Each issue should have a single, clear purpose
2. **Always Link to Epics**: Every issue should be linked to one of the 9 epics
3. **Update Status Labels**: Keep issue status up to date as work progresses
4. **Detailed Descriptions**: Include all relevant details in the issue description
5. **Acceptance Criteria**: Always define clear acceptance criteria
6. **Regular Updates**: Provide updates on progress in comments
7. **Close with References**: When closing an issue, reference the PR that resolved it
8. **Keep the Board Clean**: Regularly audit and update issues
9. **Meaningful Titles**: Use clear, descriptive titles that include type prefixes
10. **Follow Templates**: Always use and complete the provided templates

## Issue Consolidation

To maintain a clean and organized issue tracker:

1. **Identify Related Issues**: Use the `consolidate-issues.sh` script to help identify:
   - Duplicate issues
   - Stale issues (30+ days without activity)
   - Test issues that can be closed

2. **Master Issue Selection**: When finding duplicates:
   - Keep the issue with the most complete information
   - Choose the most comprehensive issue as the master
   - Enhance the master issue with information from related issues

3. **Consolidation Process**:
   - Add information from related issues to the master
   - Add comments to related issues pointing to the master
   - Close duplicate issues after consolidating information

4. **Periodic Cleanup**:
   - Run the consolidation script monthly
   - Aim to keep the total number of open issues under 50
   - Use the epic structure to maintain clear organization

For more aggressive consolidation, the `mass-consolidation.py` script can be used with caution:

```bash
cd .github/scripts
python mass-consolidation.py --dry-run
```

## Issue Comments and Communication

- Keep comments focused and relevant
- Use code blocks for code or logs: \`\`\`code\`\`\`
- Use screenshots for visual issues
- @mention relevant people when feedback is needed
- Use checkboxes for sub-tasks: `- [ ] Task`
- Mark your own issues as `in-progress` when you start working on them
- Update status labels as the issue progresses
- Provide clear, specific updates when blocked

For questions or suggestions about this issue management system, please open an issue with the title "Meta: Issue Management System Feedback".