# Fix Scripts Consolidation Plan

## Current Situation

There are 13 separate fix scripts with overlapping functionality:

```
fix_api_integration.py
fix_auth_service_role.py
fix_auth_simple.py
fix_database_modular.py
fix_database_tables.py
fix_endpoint_registration.py
fix_foreign_key_relationships.py
fix_missing_tables.py
fix_relationships.py
fix_remaining_api_issues.py
fix_rls_implementation.py
fix_supabase_auth.py
fix_table_names_and_repopulate.py
```

This causes confusion about which scripts to run and in what order.

## Recommended Structure

Create a unified maintenance utility with modules for different categories:

```
maintenance/
├── __init__.py
├── api_fixes.py  (API-related fixes)
├── auth_fixes.py (Authentication fixes)
├── database_fixes.py (Database schema and data fixes)
├── rls_fixes.py (Row Level Security fixes)
├── cli.py (Command-line interface)
└── utils.py (Common utilities for fixes)
```

## Immediate Actions

1. Create the directory structure
2. Group fix functions by category
3. Create a unified CLI entry point
4. Update documentation

## Implementation Plan

1. First, create the directory structure:

```bash
mkdir -p maintenance
touch maintenance/__init__.py
touch maintenance/cli.py
touch maintenance/utils.py
touch maintenance/api_fixes.py
touch maintenance/auth_fixes.py
touch maintenance/database_fixes.py
touch maintenance/rls_fixes.py
```

2. Then, move functionality from existing scripts to the appropriate module:

```python
# maintenance/cli.py
import argparse
import sys
from . import api_fixes, auth_fixes, database_fixes, rls_fixes

def main():
    """Main entry point for the maintenance utility."""
    parser = argparse.ArgumentParser(description="CryoProtect maintenance utility")
    subparsers = parser.add_subparsers(dest="command", help="Maintenance command to run")
    
    # API fixes commands
    api_parser = subparsers.add_parser("api", help="API-related fixes")
    api_subparsers = api_parser.add_subparsers(dest="api_command")
    
    integration_parser = api_subparsers.add_parser("integration", help="Fix API integration issues")
    endpoint_parser = api_subparsers.add_parser("endpoints", help="Fix endpoint registration issues")
    # Add more API fix commands...
    
    # Auth fixes commands
    auth_parser = subparsers.add_parser("auth", help="Authentication-related fixes")
    auth_subparsers = auth_parser.add_subparsers(dest="auth_command")
    
    service_role_parser = auth_subparsers.add_parser("service-role", help="Fix service role authentication")
    simple_auth_parser = auth_subparsers.add_parser("simple", help="Fix simple authentication")
    # Add more auth fix commands...
    
    # Database fixes commands
    db_parser = subparsers.add_parser("db", help="Database-related fixes")
    db_subparsers = db_parser.add_subparsers(dest="db_command")
    
    tables_parser = db_subparsers.add_parser("tables", help="Fix database tables")
    relationships_parser = db_subparsers.add_parser("relationships", help="Fix foreign key relationships")
    # Add more database fix commands...
    
    # RLS fixes commands
    rls_parser = subparsers.add_parser("rls", help="Row Level Security fixes")
    # Add RLS fix commands...
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run the appropriate command
    if args.command == "api":
        if args.api_command == "integration":
            api_fixes.fix_integration()
        elif args.api_command == "endpoints":
            api_fixes.fix_endpoints()
        # Handle other API fix commands...
    
    elif args.command == "auth":
        if args.auth_command == "service-role":
            auth_fixes.fix_service_role()
        elif args.auth_command == "simple":
            auth_fixes.fix_simple_auth()
        # Handle other auth fix commands...
    
    elif args.command == "db":
        if args.db_command == "tables":
            database_fixes.fix_tables()
        elif args.db_command == "relationships":
            database_fixes.fix_relationships()
        # Handle other database fix commands...
    
    elif args.command == "rls":
        rls_fixes.fix_implementation()
    
    else:
        parser.print_help()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

3. Create a main script in the root directory:

```python
#!/usr/bin/env python3
# maintain.py
from maintenance.cli import main

if __name__ == "__main__":
    main()
```

4. Make the script executable:

```bash
chmod +x maintain.py
```

## Commands for Setup

```bash
# Create directory structure
mkdir -p maintenance
touch maintenance/__init__.py
touch maintenance/cli.py
touch maintenance/utils.py
touch maintenance/api_fixes.py
touch maintenance/auth_fixes.py
touch maintenance/database_fixes.py
touch maintenance/rls_fixes.py

# Create the main script
echo '#!/usr/bin/env python3
# maintain.py
from maintenance.cli import main

if __name__ == "__main__":
    main()' > maintain.py

# Make it executable
chmod +x maintain.py

# After implementing the modules, commit the changes
git add maintenance/ maintain.py
git commit -m "Create unified maintenance utility"
```