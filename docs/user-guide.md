# User Guide

This User Guide provides instructions for using the CryoProtect system, including how to run the fix scripts, verify changes, and rollback if needed.

## Overview

The CryoProtect system has undergone significant improvements to address various issues with the database schema, security implementation, relationship design, and API integration. This guide will help you understand how to:

1. Run the fix scripts to apply these improvements
2. Verify that the changes have been applied correctly
3. Rollback changes if needed

## Table of Contents

1. [Running Fix Scripts](./user/running-fix-scripts.md)
   - Overview of the fix scripts
   - Running the master integration script
   - Running individual fix scripts
   - Command-line options

2. [Verifying Changes](./user/verifying-changes.md)
   - Verifying database schema changes
   - Verifying security implementation
   - Verifying relationship fixes
   - Verifying API integration fixes
   - Automated verification tools

3. [Rollback Procedures](./user/rollback-procedures.md)
   - When to rollback
   - Using the built-in rollback mechanisms
   - Manual rollback procedures
   - Restoring from backups

## Quick Start

For a quick start, you can run the master integration script to apply all fixes:

```bash
# Windows
run_cryoprotect_fixes.bat

# Unix/Linux/macOS
./run_cryoprotect_fixes.sh
```

This will run all the fix scripts in the correct order and verify the changes.

## System Requirements

Before running the fix scripts, ensure that your system meets the following requirements:

- Python 3.9+
- Supabase account with appropriate permissions
- Environment variables set in a `.env` file:
  ```
  SUPABASE_URL=your-supabase-url
  SUPABASE_KEY=your-supabase-key
  SUPABASE_USER=your-supabase-user
  SUPABASE_PASSWORD=your-supabase-password
  ```

## Important Notes

- Always create a backup before running the fix scripts
- Run the scripts in a test environment before applying to production
- Verify the changes after running each script
- If you encounter any issues, check the log files for error messages

For detailed instructions on each aspect of using the CryoProtect system, please refer to the specific sections linked in the Table of Contents.