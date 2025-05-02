# Phase 3.1: Deployment Infrastructure Optimized Implementation

You're implementing the Deployment Infrastructure phase (Phase 3.1) of CryoProtect v2. I've provided comprehensive resources including exact line references to make your work highly efficient and cost-effective.

## Key Resources - Review These First

1. **PHASE_3.1_RESOURCE_GUIDE.md** - Maps all files and resources
2. **PHASE_3.1_LINE_REFERENCES.md** - Specific line numbers for targeted modifications
3. **PHASE_3.1_OPTIMIZATION_STRATEGY.md** - Strategy for efficient implementation
4. **PHASE_3.1_IMPLEMENTATION_DIRECTIVE.md** - Specific tasks with priorities

## Implementation Approach

For maximum efficiency, follow this approach:

1. **Targeted File Access**
   - Request only the specific line ranges mentioned in PHASE_3.1_LINE_REFERENCES.md
   - Example: "View file: Dockerfile, lines 10-25"

2. **Batched Modifications**
   - Group similar tasks across multiple files
   - Process all CI/CD workflow files together
   - Handle all Docker configurations in one batch

3. **Focus Areas**
   - Enhance CI/CD pipeline configuration (GitHub Actions workflows)
   - Optimize Docker configuration (Dockerfile, docker-compose.yml)
   - Standardize environment configuration (config.py and variants)
   - Implement blue/green deployment scripts

## Example of Optimized Workflow

```
1. Check LINE_REFERENCES.md for specific lines to modify
2. Request only those line ranges: "View file X, lines 100-150"
3. Make targeted changes using the provided code patterns
4. Test the specific component modified
5. Document the changes and move to the next area
```

## Token Usage Optimization

- Never search for files already mapped in the resources
- Only load the specific line ranges you need to modify
- Use the provided code patterns rather than creating new ones
- Update the line references document as you progress

## CI/CD Configuration Examples

GitHub Actions workflow patterns are provided in the RESOURCE_GUIDE.md file:
- Full CI workflow example for testing and validation
- Deployment workflow with environment-specific logic
- Security scanning integration pattern

## Docker Configuration Examples

Docker best practices with examples:
- Multi-stage build pattern for efficient images
- Non-root user setup for security
- Health check implementation
- Resource constraints configuration

## Environment Configuration Pattern

Configuration class hierarchy:
- Base `Config` class with default settings
- Environment-specific subclasses
- Environment variable management

## Deployment Script Example

Blue/green deployment script structure:
- Environment determination
- Container deployment
- Health checking
- Traffic switching
- Rollback procedure

By using these precise line references and focused approach, you'll minimize token usage while maximizing productivity. All the information you need has been mapped and organized - you just need to focus on implementation.

When you've completed the implementation, provide a summary of what was accomplished and which files were updated.