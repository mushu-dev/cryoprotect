# Branch Protection Rules for CryoProtect

This document outlines the recommended branch protection rules for the CryoProtect repository to ensure code quality, security, and stability throughout the development lifecycle.

## Overview

Branch protection rules help maintain code quality by enforcing certain conditions before changes can be merged into important branches. These rules are especially important for a scientific application like CryoProtect where data integrity and reproducibility are critical.

## Current Branch Structure

The current repository uses `master` as the primary branch. This document updates the branch protection recommendations to align with the current repository structure.

## Recommended Branch Protection Rules

### Master Branch (`master`)

The `master` branch represents the production-ready state of the codebase and should be protected with the following rules:

- **Require pull request reviews before merging**
  - Require at least 1 approval
  - Dismiss stale pull request approvals when new commits are pushed
  - Require review from Code Owners (once CODEOWNERS file is created)
  - Restrict who can dismiss pull request reviews (limit to repository administrators)

- **Require status checks to pass before merging**
  - Require branches to be up to date before merging
  - Required status checks:
    - CI/CD Pipeline / Run Tests (python)
    - CI/CD Pipeline / Run Tests (node)
    - CI/CD Pipeline / Security Scan
    - Vercel Deployment / deploy (for frontend)
    - Heroku Deployment / build (for backend)

- **Require linear history**
  - Prevents merge commits and enforces a clean, linear commit history

- **Include administrators**
  - Ensure these rules apply to repository administrators as well

- **Restrict who can push to matching branches**
  - Limit direct pushes to repository administrators only

- **Allow force pushes**
  - Disallow force pushes to protect commit history

- **Allow deletions**
  - Disallow branch deletion to preserve history

### Development Branches

For feature branches and other development branches, we recommend the following naming conventions and rules:

- **Naming Conventions**
  - `feature/*` - For new features
  - `bugfix/*` - For bug fixes
  - `hotfix/*` - For urgent fixes to production
  - `release/*` - For release preparation
  - `docs/*` - For documentation updates
  - `refactor/*` - For code refactoring

- **Branch Creation**
  - All feature branches should be created from the latest `master` branch

- **Automated Checks**
  - While not enforced through branch protection, all branches should pass automated tests before being considered for merge

### GitHub Actions Configuration

The repository already has several GitHub Actions workflows:

1. **CI/CD Pipeline** (.github/workflows/ci-cd.yml)
2. **CryoProtect v2 Deployment Pipeline** (.github/workflows/deploy.yml)
3. **Vercel Deployment** (.github/workflows/vercel-deploy.yml) - Added in current PR
4. **Heroku Deployment** (.github/workflows/heroku-deploy.yml) - Added in current PR

These workflows ensure code quality and automated deployment. We recommend enforcing the passing of these checks before merging.

## Protected Tags

In addition to branch protection, we recommend protecting tags that follow the semantic versioning pattern:

- **Protected Tag Pattern**: `v*.*.*`
- **Allow Tag Creation**: Limited to repository administrators
- **Allow Tag Deletion**: Disallow tag deletion to preserve release history

## Implementation

These branch protection rules should be implemented through the GitHub repository settings. Navigate to:

1. Repository Settings
2. Branches
3. Branch protection rules
4. Add rule

For the master branch (`master`), configure the appropriate rules as described above.

## GitHub Project Board Integration

The repository is linked to the "CryoProtect Development" project board. For proper tracking:

1. Ensure all new issues are added to the project board
2. Use appropriate labels for categorization
3. Update issue status in the project board as work progresses

## Secret Management for GitHub Actions

For the newly added deployment workflows, the following secrets need to be configured in the repository:

1. **Vercel Deployment**:
   - `VERCEL_TOKEN`
   - `VERCEL_ORG_ID`
   - `VERCEL_PROJECT_ID`

2. **Heroku Deployment**:
   - `HEROKU_API_KEY`
   - `HEROKU_APP_NAME`
   - `HEROKU_EMAIL`
   - `HEROKU_PIPELINE_ID` (for review apps)

## Conclusion

Following these branch protection rules and GitHub configuration recommendations will help maintain a high-quality, secure codebase for CryoProtect. These rules should be reviewed periodically and updated as the project evolves.