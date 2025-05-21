# Contributing to CryoProtect

Thank you for considering contributing to CryoProtect! This document outlines the process for contributing to the project.

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Environment](#development-environment)
4. [Branching Strategy](#branching-strategy)
5. [Issue Tracking](#issue-tracking)
6. [Pull Requests](#pull-requests)
7. [Testing](#testing)
8. [Documentation](#documentation)
9. [Versioning](#versioning)

## Code of Conduct

Our project adheres to the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/yourusername/cryoprotect.git`
3. Setup your development environment as described in [Development Environment](#development-environment)
4. Create a feature branch for your work
5. Make your changes
6. Submit a pull request

## Development Environment

Please follow the instructions in the [README.md](README.md) file to set up your development environment. We recommend using the Cursor IDE for development.

### Prerequisites

- Python 3.9+
- Conda
- Docker (optional)

### Setup

1. Create and activate the Conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate cryoprotect
   ```

2. Configure Supabase:
   - Create a `.env` file with your Supabase credentials
   - Apply the database migration

3. Run the application:
   ```bash
   ./run_app.sh  # Linux/Mac
   run_app.bat   # Windows
   ```

## Branching Strategy

We use a feature branch workflow:

- `master` - Main branch that always contains production-ready code
- `feature/*` - Feature branches for new features or enhancements
- `bugfix/*` - Bug fix branches
- `release/*` - Release preparation branches
- `hotfix/*` - Emergency fixes for production issues

Always branch from `master` for new work and use pull requests to merge changes back.

## Issue Tracking

We use GitHub Issues for tracking work. Before starting any work:

1. Check if an issue already exists for the change you want to make
2. If not, create a new issue describing the problem or feature
3. Discuss the approach in the issue before implementing

### Issue Labels

We use a comprehensive labeling system:

- Type labels (`type:feature`, `type:bugfix`, etc.)
- Area labels (`area:database`, `area:api`, etc.)
- Status labels (`status:planning`, `status:in-progress`, etc.)
- Priority labels (`priority:high`, `priority:medium`, etc.)

## Pull Requests

When submitting a pull request:

1. Reference the issue number in the PR description
2. Include a clear summary of the changes
3. Update relevant documentation
4. Ensure all tests pass
5. Request reviews from appropriate team members

Use the pull request template provided in the repository.

## Testing

All code changes should include appropriate tests:

- Unit tests for individual functions
- Integration tests for API endpoints
- Database migration tests when applicable

Run tests locally before submitting a PR:
```bash
./run_tests.sh
```

## Documentation

Documentation is as important as code. Please update:

- Code comments and docstrings
- README and other markdown files
- API documentation when changing endpoints
- User documentation for feature changes

## Versioning

We use [Semantic Versioning](https://semver.org/) for releases:

- MAJOR.MINOR.PATCH
- MAJOR version for incompatible API changes
- MINOR version for new backwards-compatible functionality
- PATCH version for backwards-compatible bug fixes

---

Thank you for contributing to CryoProtect! Your efforts help make this project better for everyone.