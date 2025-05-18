# CryoProtect Analyzer

A tool for analyzing cryoprotectant molecules using RDKit and Convex.

## Overview

CryoProtect Analyzer is a Flask-based web application that allows users to analyze cryoprotectant molecules, calculate their properties, and store the results in a Convex database. The application uses RDKit for molecular property calculations and visualization.

## Features

- Molecular property calculation using RDKit
- Visualization of molecules
- Substructure search
- Similarity calculation
- Database storage using Convex
- Web interface using Next.js and Netlify
- Real-time data updates

## Development Environment

### Cursor IDE (Recommended)

CryoProtect now supports development using [Cursor IDE](https://cursor.sh/), which provides AI-assisted development features.

1. **Installation**:
   - Download and install Cursor from [cursor.sh](https://cursor.sh/)
   - Open the CryoProtect project folder in Cursor

2. **Features**:
   - AI-assisted code exploration and generation
   - Integrated issue tracking
   - Streamlined development workflow
   
3. **Migration**:
   - Run `./migrate_to_cursor.sh` to prepare the project for Cursor
   - See [CURSOR_MIGRATION_GUIDE.md](CURSOR_MIGRATION_GUIDE.md) for detailed instructions

### Traditional Setup

#### Prerequisites

- Python 3.9+
- Conda
- Docker (optional)

#### Setup

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/cryoprotect-analyzer.git
   cd cryoprotect-analyzer
   ```

2. Set up the environment:
   ```
   ./setup_environment.sh  # Linux/Mac
   setup_environment.bat   # Windows
   ```

3. Configure Convex:
   - Create a `.env` file with your Convex credentials
   - Apply the database configuration using `npx convex deploy`

4. Run the application:
   ```
   ./run_app.sh  # Linux/Mac
   run_app.bat   # Windows
   ```

### Docker

You can also run the application using Docker:

```
docker-compose up
```

## Usage

1. Open your browser and navigate to `http://localhost:5000`
2. Use the web interface to analyze molecules, create mixtures, and view predictions

## Documentation

- [API Documentation](README_API.md)
- [RDKit Integration](README_RDKit.md)
- [Convex Integration Guide](CONVEX_INTEGRATION_GUIDE.md)
- [Convex Integration Complete Report](CONVEX_INTEGRATION_COMPLETE.md)
- [Web Interface](README_Web.md)
- [RDKit Troubleshooting](README_RDKit_Troubleshooting.md)
- [Database Remediation Final Report](reports/DATABASE_REMEDIATION_FINAL_REPORT.md)
- [Database Maintenance Guide](reports/DATABASE_MAINTENANCE_GUIDE.md)
- [Cursor Migration Guide](CURSOR_MIGRATION_GUIDE.md)
- [Netlify Analytics Setup](NETLIFY_ANALYTICS_SETUP.md)
- [Netlify Dynamic Routes](NETLIFY_DYNAMIC_ROUTES.md)
- [Backend Frontend Integration](BACKEND_FRONTEND_INTEGRATION.md)
- [Convex Implementation Plan](CONVEX_IMPLEMENTATION_PLAN.md)

## Database Remediation Project

CryoProtect v2 has undergone a comprehensive database remediation project to address critical issues in the database structure, data integrity, security implementation, and performance. The project successfully remediated multiple issues across several key areas:

1. **Database Integrity**: Fixed empty tables and orphaned records that were causing data inconsistencies
2. **Row Level Security (RLS)**: Corrected incomplete RLS implementation to ensure proper data access controls
3. **Foreign Key Relationships**: Replaced fan traps with proper junction tables and implemented proper 3NF normalization
4. **Performance Optimization**: Added appropriate indexes and optimized query performance
5. **API Integration**: Fixed issues with API endpoints to ensure proper database interaction

The remediation has significantly improved the stability, security, and performance of the CryoProtect v2 database, providing a solid foundation for future development and research activities.

For detailed information about the remediation project, please refer to:
- [Database Remediation Final Report](reports/DATABASE_REMEDIATION_FINAL_REPORT.md) - Comprehensive documentation of all changes made
- [Database Maintenance Guide](reports/DATABASE_MAINTENANCE_GUIDE.md) - Guidelines for maintaining the new database structure

## Scheduled Job Monitoring & Failure Notification

CryoProtect v2 supports automated monitoring and notification for scheduled jobs (e.g., PubChem data updates). If a scheduled job fails (nonzero exit code), maintainers can be notified via email, webhook, or both.

### How It Works

- Use `run_with_notification.py` to wrap any scheduled job (e.g., `update_pubchem_data.py`).
- If the job fails, a notification is sent using the settings below.
- Notification logic is reusable for any automated process.

### Configuration

Set environment variables (in your shell, `.env`, or Task Scheduler) to enable notifications:

#### Email (SMTP) Notification

- `NOTIFY_SMTP_SERVER` — SMTP server address (e.g., `smtp.gmail.com`)
- `NOTIFY_SMTP_PORT` — SMTP port (default: `587`)
- `NOTIFY_SMTP_USER` — SMTP username (optional if server allows anonymous)
- `NOTIFY_SMTP_PASSWORD` — SMTP password
- `NOTIFY_FROM_ADDR` — Sender email address
- `NOTIFY_TO_ADDRS` — Comma-separated recipient emails (e.g., `admin@example.com,maintainer@example.com`)
- `NOTIFY_SMTP_USE_TLS` — `true` (default) or `false`

#### Webhook Notification

- `NOTIFY_WEBHOOK_URL` — Webhook endpoint (e.g., Slack, Discord, custom)
- `NOTIFY_WEBHOOK_HEADERS` — (Optional) JSON string of HTTP headers (e.g., `'{"Authorization": "Bearer ..."}'`)

#### Example `.env` for Email

```
NOTIFY_SMTP_SERVER=smtp.gmail.com
NOTIFY_SMTP_PORT=587
NOTIFY_SMTP_USER=youruser@gmail.com
NOTIFY_SMTP_PASSWORD=yourpassword
NOTIFY_FROM_ADDR=cryoprotect@yourdomain.com
NOTIFY_TO_ADDRS=admin@yourdomain.com,maintainer@yourdomain.com
NOTIFY_SMTP_USE_TLS=true
```

#### Example `.env` for Webhook

```
NOTIFY_WEBHOOK_URL=https://hooks.slack.com/services/XXX/YYY/ZZZ
NOTIFY_WEBHOOK_HEADERS={"Content-Type": "application/json"}
```

### Usage

Wrap your scheduled job with `run_with_notification.py`:

**Cron Example (Linux/macOS):**
```
30 2 * * * /usr/bin/env python3 /path/to/run_with_notification.py "PubChem Update" python /path/to/update_pubchem_data.py --log-file /path/to/logfile.log
```

**Task Scheduler Example (Windows):**
- Program/script: `python`
- Add arguments: `C:\path\to\run_with_notification.py "PubChem Update" python C:\path\to\update_pubchem_data.py --log-file C:\path\to\logfile.log`

### Testing Notifications

- To test email/webhook, run:
  ```
  python run_with_notification.py "Test Job" python -c "import sys; sys.exit(1)"
  ```
- You should receive a notification for the simulated failure.

### Extending the Notification System

- The notification logic is in `notify.py` (email and webhook).
- You can import and use `notify_failure()` in other scripts if needed.
- Supports both environment variables and direct config dicts.

### Script Locations

- `notify.py` — Notification logic (email/webhook)
- `run_with_notification.py` — Wrapper for scheduled jobs

For more details, see comments in each script.

## Database Utilities

- `create_database_backup.py` - Script to create backups of the database
- `create_production_backup.py` - Script to create production database backups
- `run_database_backup.bat` - Windows batch file to run database backup

## License

MIT