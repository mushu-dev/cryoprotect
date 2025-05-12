# Cursor Migration Guide

This document outlines the process for migrating the CryoProtect project from Claude Code CLI / RooRoo development environment to Cursor IDE.

## Background

Throughout development, we've been using Claude Code CLI along with the RooRoo framework for agent-based development. To streamline development and enable better team collaboration, we're transitioning to Cursor IDE.

## Migration Steps

### 1. Issue Management

Before beginning the migration, we've reorganized GitHub issues with the following structure:

- **Milestones**: Aligned with project phases
  - Phase 1: Technical Foundation
  - Phase 2: Feature Completion
  - Phase 3: Production Readiness
  - Phase 4: Documentation
  - ChEMBL Integration

- **Labels**:
  - **Component-Based**: `component:database`, `component:api`, `component:ui`, etc.
  - **Priority Levels**: `priority:high`, `priority:medium`, `priority:low`
  - **Issue Types**: `type:feature`, `type:bugfix`, `type:refactor`, etc.
  - **Status Tracking**: `status:Implemented`, `status:Validated`, etc.

### 2. Project Setup in Cursor

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/cryoprotect.git
   cd cryoprotect
   ```

2. **Open in Cursor**:
   - Launch Cursor
   - Select "Open Folder" and navigate to the cloned repository

3. **Environment Setup**:
   - Cursor will automatically detect the Python project
   - Use the existing conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate cryoprotect
   ```

### 3. Cursor-Specific Features

#### AI-Assisted Development

Cursor provides built-in AI features that can replace some of the Claude CLI functionality:

- **Code Explanation**: Use CMD+L to explain highlighted code
- **Code Generation**: Use CMD+K to generate code based on comments
- **Code Transformation**: Use CMD+SHIFT+L to transform code

#### Configuration

We've added a `.cursorrules` file to the repository to help Cursor understand our:
- Project structure
- Coding standards
- Testing conventions

### 4. Testing Workflow

Run tests in Cursor's integrated terminal:

```bash
# Run all tests
python tests/run_tests.py

# Run specific tests
python tests/run_tests.py -t test_file_name.py
```

### 5. Adapting RooRoo Scripts

Some scripts were designed specifically for the RooRoo framework. We'll need to:

1. Update scripts in `scripts/` directory to work with Cursor
2. Remove RooRoo-specific prompts and task delegation
3. Adapt tools to work directly with Cursor's AI capabilities

### 6. Documentation Updates

- README.md has been updated with Cursor-specific instructions
- Deprecated RooRoo documentation has been archived in the `docs/archived/` directory
- New Cursor-related documentation is being added to the `docs/cursor/` directory

## Timeline

- **Week 1**: Complete issue reorganization and initial Cursor setup
- **Week 2**: Update scripts and tools to work with Cursor
- **Week 3**: Complete full transition and verify all functionality
- **Week 4**: Provide training and support for team members

## References

- [Cursor Documentation](https://cursor.sh/docs)
- [GitHub Issues in Cursor](https://cursor.sh/docs/github-integration)
- [Cursor AI Commands](https://cursor.sh/docs/ai-commands)

## Support

If you encounter issues during the migration process, please:

1. Check the existing GitHub issues for similar problems
2. Create a new issue with the `cursor-migration` label if needed
3. Reach out to the technical lead for assistance