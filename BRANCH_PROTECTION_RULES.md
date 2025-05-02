# Branch Protection Rules for CryoProtect v2

This document outlines the recommended branch protection rules for the CryoProtect v2 repository to ensure code quality, security, and stability throughout the development lifecycle.

## Overview

Branch protection rules help maintain code quality by enforcing certain conditions before changes can be merged into important branches. These rules are especially important for a scientific application like CryoProtect v2 where data integrity and reproducibility are critical.

## Recommended Branch Protection Rules

### Main Branch (`main`)

The `main` branch represents the production-ready state of the codebase and should be protected with the following rules:

- **Require pull request reviews before merging**
  - Require at least 2 approvals
  - Dismiss stale pull request approvals when new commits are pushed
  - Require review from Code Owners (defined in CODEOWNERS file)
  - Restrict who can dismiss pull request reviews (limit to repository administrators)

- **Require status checks to pass before merging**
  - Require branches to be up to date before merging
  - Required status checks:
    - CI/CD Pipeline / Run Tests
    - CI/CD Pipeline / Build and Test Docker
    - Security vulnerability scans
    - Code coverage thresholds (minimum 80%)

- **Require signed commits**
  - All commits must be signed with GPG to verify author identity

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

### Staging Branch (`staging`)

The `staging` branch represents the pre-production environment and should have similar but slightly less restrictive rules:

- **Require pull request reviews before merging**
  - Require at least 1 approval
  - Dismiss stale pull request approvals when new commits are pushed

- **Require status checks to pass before merging**
  - Require branches to be up to date before merging
  - Required status checks:
    - CI/CD Pipeline / Run Tests
    - CI/CD Pipeline / Build and Test Docker

- **Require linear history**
  - Prevents merge commits and enforces a clean commit history

- **Include administrators**
  - Ensure these rules apply to repository administrators as well

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
  - All feature branches should be created from the latest `staging` branch
  - Hotfix branches should be created from `main`

- **Automated Checks**
  - While not enforced through branch protection, all branches should pass automated tests before being considered for merge

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

For each branch pattern (`main`, `staging`), configure the appropriate rules as described above.

## Enforcement and Exceptions

While these rules are recommended for most situations, there may be cases where exceptions are needed:

- Emergency hotfixes may require expedited review processes
- Initial repository setup may require temporary relaxation of certain rules

Any exceptions should be documented and temporary, with normal protections restored as soon as possible.

## Conclusion

Following these branch protection rules will help maintain a high-quality, secure codebase for CryoProtect v2. These rules should be reviewed periodically and updated as the project evolves.