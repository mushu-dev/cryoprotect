"""
Verification script for migration management module imports.

This script tests importing the migration management functions from both
database.migrations and directly from database.
"""

import sys
import traceback

def test_imports():
    """Test importing migration management functions."""
    results = {
        'status': 'SUCCESS',
        'details': []
    }
    
    # Test 1: Import from database.migrations
    try:
        print("Test 1: Importing from database.migrations...")
        from database.migrations import apply_migrations, rollback_migrations, get_migration_status, initialize_migration_tracking
        results['details'].append({
            'test': 'Import from database.migrations',
            'status': 'SUCCESS',
            'message': 'Successfully imported all functions from database.migrations'
        })
        print("✓ Success: All functions imported from database.migrations")
    except Exception as e:
        results['status'] = 'ERROR'
        error_traceback = traceback.format_exc()
        results['details'].append({
            'test': 'Import from database.migrations',
            'status': 'ERROR',
            'message': f'Error importing from database.migrations: {str(e)}',
            'traceback': error_traceback
        })
        print(f"✗ Error: Failed to import from database.migrations")
        print(f"  Error message: {str(e)}")
        print(f"  Traceback: {error_traceback}")
    
    # Test 2: Import from database
    try:
        print("\nTest 2: Importing from database...")
        from database import apply_migrations, rollback_migrations, get_migration_status, initialize_migration_tracking
        results['details'].append({
            'test': 'Import from database',
            'status': 'SUCCESS',
            'message': 'Successfully imported all functions from database'
        })
        print("✓ Success: All functions imported from database")
    except Exception as e:
        results['status'] = 'ERROR'
        error_traceback = traceback.format_exc()
        results['details'].append({
            'test': 'Import from database',
            'status': 'ERROR',
            'message': f'Error importing from database: {str(e)}',
            'traceback': error_traceback
        })
        print(f"✗ Error: Failed to import from database")
        print(f"  Error message: {str(e)}")
        print(f"  Traceback: {error_traceback}")
    
    return results

if __name__ == "__main__":
    print("Verifying migration management module imports...\n")
    results = test_imports()
    
    print("\nSummary:")
    print(f"Overall status: {results['status']}")
    for detail in results['details']:
        status_symbol = "✓" if detail['status'] == 'SUCCESS' else "✗"
        print(f"{status_symbol} {detail['test']}: {detail['status']}")
    
    # Exit with appropriate code
    sys.exit(0 if results['status'] == 'SUCCESS' else 1)