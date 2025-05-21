# Supabase Connection Troubleshooting Guide

This guide provides detailed solutions for common Supabase connection issues encountered in the CryoProtect v2 project.

## Quick Start

We've created several utility scripts to diagnose and fix common Supabase connection issues:

1. **Windows users**: Run `fix_supabase_connection.bat` as Administrator
2. **Linux/WSL users**: Run the following commands:
   ```bash
   python supabase_connection_diagnostic.py
   python fix_supabase_dns.py
   python fix_supabase_mcp.py
   ```

## Common Issues

### DNS Resolution Issues

**Symptoms:**
- "DNS error" or "getaddrinfo EAI_AGAIN" error messages
- "Could not resolve host" errors
- Connection timeouts when attempting to connect to Supabase

**Solutions:**

1. **Run the DNS fix script**:
   ```bash
   python fix_supabase_dns.py
   ```
   This script will:
   - Attempt to resolve Supabase hostnames using alternative DNS servers
   - Add entries to your hosts file (may require administrator privileges)
   - Update your .env file with direct IP addresses if needed

2. **Manual DNS fix**:
   - Find the IP address for your Supabase project:
     ```bash
     nslookup <your-project-id>.supabase.co
     nslookup db.<your-project-id>.supabase.co
     ```
   - Add entries to your hosts file:
     - Windows: `C:\Windows\System32\drivers\etc\hosts`
     - Linux/Mac: `/etc/hosts`
     ```
     XX.XX.XX.XX <your-project-id>.supabase.co
     YY.YY.YY.YY db.<your-project-id>.supabase.co
     ```

3. **Use alternative DNS servers**:
   - Configure your system to use Google DNS (8.8.8.8, 8.8.4.4) or Cloudflare DNS (1.1.1.1, 1.0.0.1)

### MCP (Supabase CLI) Issues

**Symptoms:**
- "Command not found" when trying to use MCP tools
- Authentication errors when running MCP commands
- Errors about missing packages or dependencies

**Solutions:**

1. **Run the MCP fix script**:
   ```bash
   python fix_supabase_mcp.py
   ```
   This script will:
   - Check and install Node.js and npm if needed
   - Install or update the Supabase MCP package
   - Update the access token in supabase_mcp_tools.py
   - Test the MCP connection

2. **Manual MCP fix**:
   - Ensure Node.js and npm are installed: https://nodejs.org/
   - Install the Supabase MCP package:
     ```bash
     npm install -g @supabase/mcp-server-supabase@latest
     ```
   - Generate a new access token at https://supabase.com/dashboard/account/tokens
   - Update the token in supabase_mcp_tools.py

### Authentication Issues

**Symptoms:**
- "Unauthorized" or "Invalid API key" errors
- 401 or 403 status codes when connecting to Supabase
- "Invalid JWT" or "JWT expired" errors

**Solutions:**

1. **Check your API keys**:
   - Visit https://supabase.com/dashboard/project/<your-project-id>/settings/api
   - Copy the `service_role` key (not the `anon` key)
   - Update your .env file:
     ```
     SUPABASE_KEY=your-service-role-key
     ```

2. **Check project status**:
   - Visit https://supabase.com/dashboard/projects
   - Ensure your project is running and not paused
   - If paused, restart it and wait a few minutes

3. **Reset your project password** (if using direct database connections):
   - Visit https://supabase.com/dashboard/project/<your-project-id>/settings/database
   - Click "Reset Database Password"
   - Update your .env file with the new password

## Project-Specific Configuration

For CryoProtect v2, ensure your .env file contains the following Supabase configuration:

```
# Supabase Connection Details
SUPABASE_URL=https://<your-project-id>.supabase.co
SUPABASE_KEY=<your-service-role-key>
SUPABASE_PROJECT_ID=<your-project-id>

# Database credentials for direct connection (if needed)
SUPABASE_DB_HOST=db.<your-project-id>.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=<your-db-password>
```

Replace:
- `<your-project-id>` with your Supabase project ID (e.g., "tsdlmynydfuypiugmkev")
- `<your-service-role-key>` with your Supabase service_role key
- `<your-db-password>` with your database password

## Alternative Connection Methods

If you continue to experience issues, consider these alternative approaches:

### 1. Direct PostgreSQL Connection

Instead of using the Supabase REST API, you can connect directly to the PostgreSQL database:

```python
import psycopg2

conn = psycopg2.connect(
    host="db.<your-project-id>.supabase.co",
    port=5432,
    dbname="postgres",
    user="postgres",
    password="<your-db-password>"
)
```

### 2. Local Development Database

For development, you can run a local PostgreSQL database with Docker:

```bash
docker run --name postgres -e POSTGRES_PASSWORD=postgres -p 5432:5432 -d postgres:14
```

Then update your .env file:
```
SUPABASE_URL=http://localhost:8000
DB_HOST=localhost
DB_PORT=5432
DB_NAME=postgres
DB_USER=postgres
DB_PASSWORD=postgres
```

### 3. Use Supabase's Local Development Environment

Supabase provides a local development environment that you can run with Docker:

```bash
npx supabase start
```

This will give you a local Supabase instance with all features available at `http://localhost:8000`.

## Getting Help

If you're still experiencing issues:

1. Run the diagnostic script and share the results:
   ```bash
   python supabase_connection_diagnostic.py > diagnostic_results.txt
   ```

2. Check the Supabase status page: https://status.supabase.com/

3. Contact Supabase support: https://supabase.com/support

4. Consult the Supabase documentation: https://supabase.com/docs

## Troubleshooting Checklist

- [ ] Verified project ID is correct
- [ ] Confirmed service_role key is valid
- [ ] Checked project is active (not paused)
- [ ] Ensured DNS resolution works for project hostname
- [ ] Verified network allows outbound connections to Supabase
- [ ] Checked Node.js and npm are installed for MCP tools
- [ ] Validated database credentials (if using direct connection)
- [ ] Confirmed RLS policies allow access (if using service_role key)
- [ ] Tested connection with alternative methods