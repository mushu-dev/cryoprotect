# GitHub Issue Management Guidelines

This document provides the updated guidelines for managing issues in the CryoProtect repository.

## Epic-Based Organization

All work is now organized into 9 epic categories:

1. **Database Implementation and Optimization**: Database schema, migrations, queries, connection pooling, and performance optimization
2. **ChEMBL Integration and Data Pipeline**: Integration with ChEMBL database, data import processes, and ETL pipelines
3. **PubChem Integration and Molecule Management**: Integration with PubChem, molecule properties, data consolidation
4. **API Development and Standardization**: API endpoints, standardization, documentation, and gateway services
5. **Frontend Implementation and User Experience**: UI components, responsive design, user workflows, and visual styling
6. **Authentication and Security Implementation**: JWT auth, RLS policies, RBAC, and security enhancements
7. **RDKit Integration and Chemical Functionality**: RDKit services, chemical calculations, and molecular analysis
8. **Infrastructure and Deployment**: Deployment processes, container orchestration, CI/CD, and hosting configuration
9. **Testing and Quality Assurance**: Test frameworks, CI validation, and quality assurance procedures

## Issue Structure

All issues should include:

1. **Clear Title** - Brief but descriptive, with a prefix indicating type (e.g., "Bugfix:", "Feature:")
2. **Detailed Description** - What needs to be done and why
3. **Acceptance Criteria** - Specific, measurable conditions that define "done"
4. **Dependencies** - Issues that must be completed first
5. **Technical Details** - Implementation guidance if relevant

## Labels

### Area Labels
- `area:database` - Database schema, connections, queries
- `area:api` - API endpoints, authentication, validation
- `area:frontend` - User interface components and interactions
- `area:deployment` - Deployment processes and configurations
- `area:security` - Authentication, authorization, and security features
- `area:testing` - Test infrastructure and implementation
- `area:documentation` - Documentation updates
- `area:rdkit` - Chemical functionality and RDKit integration
- `area:infrastructure` - Infrastructure setup and configuration

### Type Labels
- `type:feature` - New functionality
- `type:bugfix` - Bug fixes
- `type:refactor` - Code improvements without functional changes
- `type:documentation` - Documentation updates
- `type:validation` - Testing and validation tasks
- `type:chore` - Maintenance tasks
- `type:epic` - Epic-level issue that groups related tasks

### Status Labels
- `status:planning` - In design/planning phase
- `status:ready` - Ready for implementation
- `status:in-progress` - Currently being worked on
- `status:review` - Implementation complete, needs review
- `status:blocked` - Blocked on some dependency
- `status:completed` - Fully completed and validated
- `status:stale` - No activity for 30+ days

### Priority Labels
- `priority:critical` - Urgent issues blocking progress
- `priority:high` - Important for current milestone
- `priority:medium` - Should be addressed in current sprint
- `priority:low` - Nice to have

### Resolution Labels (for closed issues)
- `resolution:completed` - Successfully implemented
- `resolution:wontfix` - Will not be implemented
- `resolution:duplicate` - Duplicate of another issue

## Project Board Integration

All issues should be added to the CryoProtect project board:

1. **To Do**: Initial backlog and planning
2. **In Progress**: Currently being worked on
3. **Review**: Complete but pending review
4. **Done**: Fully implemented and validated

Issues are automatically added to the project board via GitHub Actions.

## Issue Lifecycle

1. **Creation**
   - Use the appropriate issue template
   - Add relevant labels (area, type, priority)
   - Associate with the relevant epic category

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
   - Update to `status:review`
   - Reference pull request(s)
   - Tag reviewers in comments

6. **Completion**
   - Update to `status:completed`
   - Ensure all acceptance criteria are met
   - Add `resolution:completed` label
   - Close the issue

## Tooling and Automation

### Issue Management Scripts

The `scripts/issue-audit` directory contains tools for issue management:

1. **audit_github_issues.sh**: Generate audit report of all issues
2. **update_issue_audit.py**: Update issues to conform to current standards
3. **cleanup_github_issues.sh**: Clean up issues based on audit results
4. **standardize_github_repository.sh**: Set up standard labels and templates
5. **consolidate_github_issues.sh**: Manage duplicate or related issues

Run regular audits to maintain repository health:

```bash
cd scripts/issue-audit
python update_issue_audit.py --dry-run
```

### GitHub Actions Automation

Automated workflows help maintain repository consistency:

1. **Issue Assignment**: Automatically assigns issues based on labels
2. **Project Board**: Adds issues to project board and updates status
3. **PR Linking**: Links PRs to issues when mentioned in commit messages
4. **Status Updates**: Updates issue status based on PR activity
5. **Label Verification**: Ensures issues have required labels

## TaskMaster Integration

For work in the Cursor IDE environment, use TaskMaster to sync with GitHub issues:

1. **Create TaskMaster Task**: For issues being worked on in Cursor IDE
   ```
   /task-master create --title "Issue #123: Feature Name" --github-issue 123
   ```

2. **Link Existing Task**: To link an existing task to a GitHub issue
   ```
   /task-master link --id TASK_ID --github-issue 123
   ```

3. **Sync Status**: Update TaskMaster when issue status changes
   ```
   /sync-tasks commit
   ```

## Best Practices

1. **One Issue, One Purpose** - Each issue should have a single, clear purpose
2. **Regular Updates** - Provide updates on progress
3. **Detailed Descriptions** - Include all relevant details
4. **Clear Epic Association** - Associate with the appropriate epic category
5. **Consistent Labels** - Use the standard labeling system
6. **Proper Closure** - Close with a reference to the PR that resolved it
7. **Keep the Board Clean** - Regularly audit and update issues
8. **Close Stale Issues** - Don't let inactive issues linger

## Transitioning from Phase-Based to Epic-Based Organization

For historical context, the project previously used a phase-based organization:

1. **Phase 1**: Database Architecture & Foundation
2. **Phase 2**: API Implementation & Core Features
3. **Phase 3**: Deployment & Integration
4. **Phase 4**: Documentation & Testing

Issues with old phase labels (`phase:1`, `phase:2.3`, etc.) should be updated to the new epic-based system. Use the `update_issue_audit.py` script to assist with this transition.