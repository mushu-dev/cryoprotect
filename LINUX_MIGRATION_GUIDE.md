# CryoProtect Linux Migration Guide

This guide provides instructions for migrating the CryoProtect project from Windows to Linux.

## Migration Overview

CryoProtect is designed to be cross-platform compatible, but several Windows-specific components need to be updated for Linux:

1. Shell script equivalents for batch files
2. Path separators in configuration
3. Environment activation commands
4. Permission settings for executable files
5. Database configuration and connections

## Prerequisites

Before starting the migration, ensure you have the following installed on your Linux system:

- Python 3.9+
- Conda (Miniconda or Anaconda)
- PostgreSQL
- Git
- Node.js (for database migrations)

## Migration Steps

### 1. Clone and Prepare Repository

```bash
# Clone the repository
git clone https://github.com/yourusername/cryoprotect-analyzer.git
cd cryoprotect-analyzer

# Make shell scripts executable
chmod +x *.sh
chmod +x batch_scripts/*.sh
```

### 2. Run Linux Environment Setup

We've created a helper script to prepare your Linux environment:

```bash
./create_linux_env.sh
```

This script will:
- Create a `.env` file from template (if needed)
- Update Docker secrets path for Linux
- Make all shell scripts executable
- Create a local PostgreSQL database (if needed)

### 3. Set Up Conda Environment

```bash
./setup_environment.sh
```

This will:
- Create a conda environment named `cryoprotect`
- Install Python 3.9 and dependencies
- Install RDKit for molecular analysis

### 4. Configure Database Connection

Edit the `.env` file with your Supabase credentials:

```
SUPABASE_URL=your-project-url
SUPABASE_KEY=your-service-role-key
```

### 5. Apply Database Migrations

```bash
# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Apply migrations
node migrations/apply_migration.js
```

### 6. Test the Connection

```bash
python check_supabase_connection.py
```

### 7. Run the Application

```bash
./run_app_with_fix.sh
```

## Linux Shell Script Equivalents

We've created Linux equivalents for the following Windows batch files:

| Windows Batch File | Linux Shell Script |
|--------------------|-------------------|
| `run_app_with_fix.bat` | `run_app_with_fix.sh` |
| `run_tests_conda.bat` | `run_tests_conda.sh` |
| `batch_scripts/start_server.bat` | `batch_scripts/start_server.sh` |

## Platform-Specific Code

The following files contain platform-specific code that's been updated to work on Linux:

1. `ip_resolver.py` - Modified to use proper DNS resolution methods on Linux
2. `config.py` - Updated Docker secrets path for Linux
3. All scripts - Updated to use Linux-style environment activation

## Common Commands

| Task | Command |
|------|---------|
| Run application | `./run_app_with_fix.sh` |
| Run tests | `./run_tests_conda.sh` or `./run_tests.sh` |
| Start server | `./batch_scripts/start_server.sh` |
| Apply migrations | `node migrations/apply_migration.js` |
| Create database backup | `python create_database_backup.py` |
| Check database connection | `python check_supabase_connection.py` |

## Troubleshooting

### Common Issues

1. **Permission Denied for Shell Scripts**
   ```bash
   chmod +x *.sh
   chmod +x batch_scripts/*.sh
   ```

2. **Conda Environment Not Found**
   ```bash
   conda env list
   # If not listed, run:
   ./setup_environment.sh
   ```

3. **PostgreSQL Connection Issues**
   ```bash
   sudo systemctl status postgresql
   # Make sure PostgreSQL is running:
   sudo systemctl start postgresql
   ```

4. **RDKit Import Errors**
   ```bash
   conda activate cryoprotect
   conda install -c conda-forge rdkit
   ```

### Testing Your Migration

Run the full test suite to verify everything works correctly:

```bash
./run_tests.sh
```

## Migration Verification Checklist

- [ ] Environment setup script runs without errors
- [ ] Conda environment created with all dependencies
- [ ] Database connection succeeds
- [ ] Application runs without errors
- [ ] All tests pass
- [ ] API endpoints function correctly
- [ ] Database migrations apply successfully
- [ ] Molecular visualization works