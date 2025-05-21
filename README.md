# CryoProtect Analyzer 🚀

**Production-ready scientific research platform for cryoprotectant molecule analysis**

> 🚀 **PRODUCTION READY:** Full Convex integration with 5000+ molecules from ChEMBL/PubChem, real-time calculations via RDKit service, and collaborative experiment management. See [CONVEX_PRODUCTION_INTEGRATION_PLAN.md](CONVEX_PRODUCTION_INTEGRATION_PLAN.md) for complete details.

> 🔬 **SCIENTIFIC DATA:** Comprehensive database of cryoprotectant molecules with ChEMBL identifiers, PubChem naming, molecular properties, and experimental protocols validated by domain experts.

> ⚡ **REAL-TIME:** Live molecular property calculations, collaborative experiment editing, and instant search across the complete scientific dataset.

## Overview

CryoProtect Analyzer is a production-ready scientific research platform that enables researchers to analyze cryoprotectant molecules, calculate properties in real-time, and collaborate on experiments. Built with Convex for real-time data, RDKit for molecular calculations, and a modern React frontend.

## Features

- Molecular property calculation using RDKit
- Visualization of molecules
- Substructure search
- Similarity calculation
- Database storage using Convex
- Web interface using Next.js and Netlify
- Real-time data updates
- **NEW:** Advanced resiliency patterns for robust operation
- **NEW:** Enhanced Protocol Designer with comprehensive step management

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

3. Configure the Database (Supabase or Convex):
   - For Supabase (default):
     - Create a `.env` file with your Supabase credentials
   - For Convex (recommended):
     - Copy `.env.convex.sample` to `.env` and fill in your Convex credentials
     - Configure the main environment variables:
       ```
       CONVEX_DB_ENABLED=true
       CONVEX_URL=https://your-deployment-id.convex.cloud
       CONVEX_DEPLOYMENT_KEY=your-deployment-key
       JWT_SECRET=your-jwt-secret-key
       ```
     - For additional settings, see [README_CONVEX_BACKEND.md](README_CONVEX_BACKEND.md)
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

### Core Documentation
- [API Documentation](README_API.md)
- [Web Interface](README_Web.md)
- [RDKit Integration](README_RDKit.md)
- [RDKit Troubleshooting](README_RDKit_Troubleshooting.md)

### Database Documentation
- [Database Remediation Final Report](reports/DATABASE_REMEDIATION_FINAL_REPORT.md)
- [Database Maintenance Guide](reports/DATABASE_MAINTENANCE_GUIDE.md)

### Convex Integration (Latest)
- [**✅ Enhanced Convex Backend Integration**](README_CONVEX_BACKEND.md) - **NEW**
- [**✅ Convex Backend Integration Documentation**](CONVEX_BACKEND_INTEGRATION.md) - **NEW**
- [**✅ Convex Frontend Integration**](CONVEX_FRONTEND_INTEGRATION.md) - **NEW**
- [Backend Frontend Integration](BACKEND_FRONTEND_INTEGRATION.md)
- [Convex Implementation Plan](CONVEX_IMPLEMENTATION_PLAN.md)

### Production Readiness (Latest)
- [**✅ Production Readiness Plan**](PRODUCTION_READINESS_PLAN.md) - **NEW**
- [Resiliency Patterns](api/resiliency/README.md) - **NEW**
- Demo: `examples/resiliency_demo.py` - **NEW**

### RDKit Service
- [**✅ RDKit Service Deployment**](https://cryoprotect-rdkit.fly.dev/health) - Deployed to fly.io
- Core endpoints:
  - Health check: `https://cryoprotect-rdkit.fly.dev/health`
  - CORS test: `https://cryoprotect-rdkit.fly.dev/test-cors`
  - Property calculation: `https://cryoprotect-rdkit.fly.dev/calculate`

### Deployment
- [**🚀 Netlify Automatic Deployment**](NETLIFY_AUTOMATIC_DEPLOYMENT.md) - **NEW!** Complete guide for GitHub Actions deployment
- [**⚠️ Netlify Deployment Fix**](NETLIFY_DEPLOYMENT_FIX.md) - Fix for static fallback page
- [Netlify Deployment Solution](NETLIFY_DEPLOYMENT_SOLUTION.md) - Comprehensive solution for UI rendering issues
- [Netlify Analytics Setup](NETLIFY_ANALYTICS_SETUP.md)
- [Netlify Dynamic Routes](NETLIFY_DYNAMIC_ROUTES.md)

### Experimental Data
- [**✨ Protocol Designer Documentation**](frontend/src/features/protocols/README.md) - **NEW**
- [Experimental Data Enhancement](EXPERIMENTAL_DATA_ENHANCEMENT_SUMMARY.md)

### Development Environment
- [Cursor Migration Guide](CURSOR_MIGRATION_GUIDE.md)

## Resiliency Patterns

CryoProtect v2 now includes advanced resiliency patterns to ensure robust operation in production environments. These patterns help prevent cascading failures, improve system stability, and provide better user experience during service degradation.

### Core Resiliency Patterns

1. **Retry with Exponential Backoff**
   - Automatically retries failed operations with increasing delays
   - Configurable retry attempts, backoff factors, and jitter
   - Supports specific exception types for targeted retries
   - Built-in logging and observability

2. **Circuit Breaker**
   - Prevents cascading failures by stopping calls to failing services
   - Automatic recovery with half-open state for testing
   - Configurable thresholds and recovery timeouts
   - Integration with monitoring systems

3. **Timeout Management**
   - Enforces time limits on operations to prevent resource blockage
   - Thread-safe implementation for multi-threaded applications
   - Integration with logging and observability

4. **Service Health Tracking**
   - Real-time monitoring of service health based on success rates and response times
   - Automatic detection of degraded or failing services
   - Integration with circuit breaker and retry mechanisms

### Example Usage

```python
from api.resiliency import retry_with_backoff, circuit_breaker, with_timeout, track_service_health

# Basic retry example
@retry_with_backoff(max_retries=3, base_delay=0.5)
def external_api_call():
    # This function will retry up to 3 times with exponential backoff
    pass

# Combined patterns example
@retry_with_backoff(max_retries=3)
@circuit_breaker(name="database")
@with_timeout(seconds=5)
@track_service_health("database")
def database_operation():
    # This function has comprehensive resiliency
    pass
```

For a complete demonstration, run:

```
python examples/resiliency_demo.py
```

See the [Production Readiness Plan](PRODUCTION_READINESS_PLAN.md) for more details.

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