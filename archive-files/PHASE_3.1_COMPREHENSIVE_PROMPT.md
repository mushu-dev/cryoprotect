# Phase 3.1: Deployment Infrastructure Implementation - CryoProtect v2

You are implementing Phase 3.1 (Deployment Infrastructure) of the CryoProtect v2 project. Review all provided resources before beginning implementation.

## Resource Overview - READ FIRST

I've provided comprehensive resources to maximize your efficiency:

1. **PHASE_3.1_RESOURCE_GUIDE.md**: Maps all files and resources
2. **PHASE_3.1_LINE_REFERENCES.md**: Specific line numbers to focus on
3. **PHASE_3.1_IMPLEMENTATION_DIRECTIVE.md**: Detailed tasks with examples
4. **PHASE_3.1_OPTIMIZATION_STRATEGY.md**: Token-efficient implementation approach
5. **PHASE_3.1_DATABASE_GUIDE.md**: Database population strategy

## Critical Problem to Solve

**The Supabase database currently has no scientific data.** In addition to the standard deployment infrastructure, you must implement the database population process to ensure the application works with real cryoprotectant data.

## Key Files for Database Population

- `populate_database_supabase.py` - Main script for database population (exists but needs integration)
- `service_role_helper.py` - Helper for Supabase service role authentication (needs implementation)
- `.env.template` - Template for environment variables (needs creation)
- `.github/workflows/deploy.yml` - Deployment workflow (needs database population integration)

## Implementation Tasks

### 1. Complete CI/CD Pipeline with Database Handling

- Update GitHub Actions workflows to include database schema and data population
- Ensure service role authentication works securely
- Implement environment-specific deployments (dev/staging/prod)

### 2. Optimize Docker Configuration

- Create multi-stage builds for smaller, more secure images
- Configure proper volume mounts for data persistence
- Set up appropriate networking and security

### 3. Standardize Environment Configuration

- Create a unified approach across environments (dev/staging/prod)
- Implement secure handling of Supabase credentials
- Create documentation for environment setup

### 4. Implement Blue/Green Deployment with Database Awareness

- Create deployment scripts that include database readiness checks
- Add database population to deployment process
- Ensure rollback preserves database integrity

## Implementation Approach

For maximum efficiency:

1. Request only specific line ranges from PHASE_3.1_LINE_REFERENCES.md
2. Group similar modifications using BatchTool
3. Use provided code patterns rather than creating new ones
4. Start with the database population integration as it's critical for functionality

## Sample Code Patterns

### Service Role Authentication:

```python
def get_supabase_client():
    """Connect to Supabase using service role key."""
    url = os.getenv("SUPABASE_URL")
    key = os.getenv("SUPABASE_SERVICE_KEY")
    if not url or not key:
        raise ValueError("Missing Supabase configuration")
    return create_client(url, key)
```

### GitHub Actions Database Population Step:

```yaml
- name: Apply database migrations
  run: python apply_migrations.py
  env:
    SUPABASE_URL: ${{ secrets.SUPABASE_URL }}
    SUPABASE_SERVICE_KEY: ${{ secrets.SUPABASE_SERVICE_KEY }}

- name: Populate database
  run: python populate_database_supabase.py
  env:
    SUPABASE_URL: ${{ secrets.SUPABASE_URL }}
    SUPABASE_SERVICE_KEY: ${{ secrets.SUPABASE_SERVICE_KEY }}
    DEFAULT_USER_ID: ${{ secrets.DEFAULT_USER_ID }}
```

### Docker Multi-Stage Build:

```dockerfile
# Build stage
FROM python:3.9-slim AS builder
WORKDIR /app
COPY requirements.txt .
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /app/wheels -r requirements.txt

# Final stage
FROM python:3.9-slim
WORKDIR /app
COPY --from=builder /app/wheels /wheels
RUN pip install --no-cache /wheels/*
COPY . .
```

## Expected Deliverables

1. Complete CI/CD configuration with database integration
2. Service role authentication implementation
3. Docker configuration optimized for production
4. Environment configuration templates
5. Blue/green deployment scripts with database awareness
6. Comprehensive documentation for operations

Begin by examining the database population requirements, as this is critical for a functional deployment. The scientific data must be loaded during deployment for the application to work properly.

When you've completed the implementation, provide a summary of what was accomplished and how the database population is integrated into the deployment process.