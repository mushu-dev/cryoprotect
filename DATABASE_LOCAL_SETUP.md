# Local Database Setup Guide

This guide explains how to set up and use the local PostgreSQL database environment for CryoProtect v2 development.

## Prerequisites

Before setting up the local database, ensure you have the following installed:

1. **PostgreSQL** (version 12 or higher)
   - [Download PostgreSQL](https://www.postgresql.org/download/)
   - During installation, note the password you set for the `postgres` user

2. **Python** (version 3.8 or higher)
   - Required packages: `psycopg2`, `python-dotenv` (install via `pip install psycopg2-binary python-dotenv`)

## Configuration

1. **Environment Variables**

   Copy the `.env.template` file to `.env` in the project root directory:

   ```bash
   cp .env.template .env
   ```

2. **Update Local Database Settings**

   Edit the `.env` file and update the following variables:

   ```
   # Local PostgreSQL database configuration
   LOCAL_DB_HOST=localhost
   LOCAL_DB_PORT=5432
   LOCAL_DB_NAME=cryoprotect
   LOCAL_DB_USER=postgres
   LOCAL_DB_PASSWORD=your_postgres_password
   LOCAL_DB_MIN_CONNECTIONS=1
   LOCAL_DB_MAX_CONNECTIONS=10
   ```

   Replace `your_postgres_password` with the password you set during PostgreSQL installation.

## Database Initialization

The `init_local_db.py` script automates the setup of your local database. It performs the following tasks:

1. Creates the database if it doesn't exist
2. Applies all SQL migrations from the `migrations` directory
3. Creates test data for development

### Running the Initialization Script

```bash
python database/init_local_db.py
```

### Command-line Options

- `--reset`: Reset the database if it already exists (drops and recreates it)
- `--skip-migrations`: Skip applying migrations
- `--skip-test-data`: Skip creating test data

Examples:

```bash
# Reset the database and start fresh
python database/init_local_db.py --reset

# Only create the database structure without test data
python database/init_local_db.py --skip-test-data

# Reset the database but skip migrations
python database/init_local_db.py --reset --skip-migrations
```

## Database Schema

The initialization script creates the following schema:

1. **Core Tables**
   - `molecules`: Stores basic molecule information
   - `property_types`: Defines available property types
   - `molecular_properties`: Stores properties for molecules
   - `mixtures`: Defines mixtures of molecules
   - `mixture_components`: Links molecules to mixtures with concentrations

2. **User Tables**
   - `auth.users`: Authentication users
   - `user_profile`: User profile information

3. **Metadata Tables**
   - `migrations`: Tracks applied migrations

## Test Data

The initialization script creates the following test data:

1. **Test User**
   - Email: `test@example.com`
   - Display Name: `Test User`
   - Affiliation: `Test Organization`

2. **Test Molecules**
   - Dimethyl sulfoxide (DMSO)
   - Glycerol
   - Trehalose

3. **Test Mixtures**
   - DMSO-Glycerol Mix
   - Trehalose Solution

## Connecting to the Local Database

### Using the Database Adapter

The recommended way to connect to the database is through the database adapter system:

```python
from database.connection_manager import ConnectionManager

# Get connection manager instance
db = ConnectionManager.get_instance()

# Execute a query
results = db.execute_query("SELECT * FROM molecules LIMIT 10")
```

The connection manager will automatically use the local database configuration from your `.env` file when running in development mode.

### Direct Connection

For debugging or administrative tasks, you can connect directly to the database:

```python
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import os

# Load environment variables
load_dotenv()

# Connect to database
conn = psycopg2.connect(
    host=os.getenv('LOCAL_DB_HOST', 'localhost'),
    port=os.getenv('LOCAL_DB_PORT', '5432'),
    dbname=os.getenv('LOCAL_DB_NAME', 'cryoprotect'),
    user=os.getenv('LOCAL_DB_USER', 'postgres'),
    password=os.getenv('LOCAL_DB_PASSWORD', '')
)

# Create cursor
cursor = conn.cursor(cursor_factory=RealDictCursor)

# Execute query
cursor.execute("SELECT * FROM molecules LIMIT 10")
results = cursor.fetchall()

# Close connection
cursor.close()
conn.close()
```

### Using pgAdmin or Other Tools

You can also connect to the local database using pgAdmin or other PostgreSQL client tools:

- **Host**: localhost (or the value of `LOCAL_DB_HOST`)
- **Port**: 5432 (or the value of `LOCAL_DB_PORT`)
- **Database**: cryoprotect (or the value of `LOCAL_DB_NAME`)
- **Username**: postgres (or the value of `LOCAL_DB_USER`)
- **Password**: Your PostgreSQL password (from `LOCAL_DB_PASSWORD`)

## Troubleshooting

### Connection Issues

1. **Authentication Failed**

   If you see an error like:
   ```
   FATAL: password authentication failed for user "postgres"
   ```

   Ensure the `LOCAL_DB_PASSWORD` in your `.env` file matches the password you set during PostgreSQL installation.

2. **Database Does Not Exist**

   If you see an error like:
   ```
   FATAL: database "cryoprotect" does not exist
   ```

   Run the initialization script to create the database:
   ```bash
   python database/init_local_db.py
   ```

3. **Port Already in Use**

   If PostgreSQL is running on a different port, update the `LOCAL_DB_PORT` in your `.env` file.

### Migration Issues

1. **Missing Migration Files**

   If you see a warning about missing migration files, ensure the `migrations` directory exists and contains SQL migration files.

2. **Migration Errors**

   If a migration fails, check the error message for details. You may need to manually fix the database schema or modify the migration file.

## Best Practices

1. **Regular Backups**

   Before making significant changes, back up your local database:
   ```bash
   pg_dump -U postgres -d cryoprotect > backup.sql
   ```

2. **Reset When Needed**

   If your local database gets into an inconsistent state, reset it:
   ```bash
   python database/init_local_db.py --reset
   ```

3. **Sync Migrations**

   Keep your migrations directory up to date with the latest schema changes.

4. **Test with Real Data**

   For more realistic testing, you can import production data into your local database (if available).

## Additional Resources

- [PostgreSQL Documentation](https://www.postgresql.org/docs/)
- [psycopg2 Documentation](https://www.psycopg.org/docs/)
- [python-dotenv Documentation](https://github.com/theskumar/python-dotenv)