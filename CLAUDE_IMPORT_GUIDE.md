# Claude.md Import System Guide

This guide explains how to use the import functionality in CLAUDE.md files to create a modular documentation system that enhances Claude Code's context management.

## Overview

The Claude.md import system allows you to:
- Break down large context files into modular components
- Import module-specific documentation where needed
- Maintain clear separation of concerns in documentation
- Keep related information together while reducing duplication
- Optimize token usage by modularizing context

## Basic Syntax

To import another markdown file into a CLAUDE.md file, use the following syntax:

```markdown
@path/to/file.md
```

This directive will be processed when Claude Code loads, and the content of the referenced file will be included inline at that position.

## Import Formats

The system supports several path formats:

1. **Relative paths** (from the main CLAUDE.md location):
   ```markdown
   @api/claude.md
   ```

2. **Absolute paths**:
   ```markdown
   @/home/mushu/Projects/CryoProtect/chembl/claude.md
   ```

3. **Dot notation** for explicit relative paths:
   ```markdown
   @./pubchem/claude.md
   ```

## Usage Patterns

### Module Structure

Each major component of the system should have its own claude.md file:

```
/api/claude.md
/database/claude.md
/chembl/claude.md
/pubchem/claude.md
/monitoring/claude.md
/security/claude.md
/tests/claude.md
/migrations/claude.md
/docs/claude.md
```

### Main CLAUDE.md Organization

The main CLAUDE.md should:
1. Provide high-level project overview
2. Import module-specific documentation in relevant sections
3. Handle cross-cutting concerns

Example:

```markdown
# CLAUDE.md

## Project Overview
[High-level description of the project]

## API Layer
@api/claude.md

## Database Layer
@database/claude.md

## Additional Components
[Other project information]
```

## Best Practices

### File Organization

1. **Consistent structure**: Each module's claude.md should follow a similar structure
   - Start with a clear heading and overview
   - Use consistent heading levels
   - Include key components, patterns, and reference files

2. **Single responsibility**: Each module's file should focus on its specific domain
   - Don't duplicate information across files
   - Cross-reference rather than duplicate where possible

3. **Heading hierarchy**: Maintain proper heading levels
   - Start each module file with a top-level heading (# Heading)
   - Main CLAUDE.md will preserve the heading structure

### Content Guidelines

1. **Self-contained**: Each module file should make sense on its own
   - Avoid references to "above" or "previous section"
   - Include necessary context within each file

2. **Concise**: Keep content focused and relevant
   - Aim for specific, actionable guidance
   - Include code examples where helpful
   - Use bullet points for clarity

3. **Link to source code**: Reference specific files and line numbers
   - Use format: `file_path:line_number`
   - Link examples to real code in the codebase

## Example Module File

Here's a template for a module-specific claude.md file:

```markdown
# Module Name Reference

## Overview

Brief description of the module and its purpose in the system.

## Key Components

### Component 1
- Description of the component
- Key interfaces or patterns
- Example usage

### Component 2
- Description of the component
- Key interfaces or patterns
- Example usage

## Key Files
- `path/to/file1.py`: Description of the file
- `path/to/file2.py`: Description of the file

## Best Practices
1. Specific guidance for working with this module
2. Common patterns to follow
3. Important considerations

## Common Pitfalls
1. Mistake to avoid
2. Another mistake to avoid
```

## Maintenance Guidelines

1. **Review regularly**: Audit all claude.md files when significant changes occur
2. **Update modules first**: When changing functionality, update the module-specific file first
3. **Consistency check**: Ensure imports in main CLAUDE.md match actual files
4. **Version control**: Commit claude.md updates with related code changes

## Implementation Notes

- The import system processes imports at Claude Code launch time
- Nested imports are supported (a file can import other files)
- The system maintains heading hierarchy from imported files
- Imports are processed in the order they appear in the file

## Troubleshooting

If imports aren't working properly:

1. Check file paths for accuracy
2. Ensure the @ symbol is at the beginning of a line
3. Verify that no spaces exist between @ and the file path
4. Confirm files exist at the specified locations
5. Look for any circular imports that might cause issues