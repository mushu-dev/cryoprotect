"""
Main entry point for the database migrations package.
Allows running the package as a module: python -m database.migrations
"""

from database.migrations.runner import main

if __name__ == '__main__':
    main()