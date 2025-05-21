# Migrating from Simplified App to Full CryoProtect App on Heroku

This guide provides comprehensive instructions for migrating from the simplified test app to the full CryoProtect application on Heroku, including database setup, environment configuration, and troubleshooting.

## Prerequisites

- A working Heroku deployment of the simplified app
- Heroku CLI installed and authenticated
- Git repository access
- PostgreSQL add-on already attached to your Heroku app
- Redis add-on already attached (or you'll need to add it)

## Migration Steps

### Automated Migration

The easiest way to migrate is using the provided migration script:

1. Make sure you're in the root directory of the CryoProtect project
2. Run the migration script:

```bash
chmod +x migrate_to_full_app.sh
./migrate_to_full_app.sh
```

The script will:
- Create a new git branch for the full app deployment
- Update the Procfile to use app.py with a release phase for database setup
- Configure the enhanced database initialization script
- Set up the Heroku-specific configuration
- Commit the changes and push to Heroku
- Monitor the deployment logs

### Manual Migration

If you prefer to perform the migration manually, follow these steps:

1. **Create a new branch for the full app deployment**:

```bash
git checkout -b heroku-full-app
```

2. **Update the Procfile to use the full app with database initialization**:

```bash
echo "release: python setup_database.py" > Procfile
echo "web: gunicorn app:app --log-file -" >> Procfile
```

3. **Ensure the improved database initialization script is ready**:

The project now includes an enhanced `setup_database.py` that:
- Waits for the database to be available
- Sets necessary environment variables
- Runs migrations if available
- Checks and creates required database structure
- Applies RLS policies if needed

4. **Ensure the Heroku configuration is properly set up**:

The `config_heroku.py` file now:
- Automatically configures the app for Heroku environment
- Parses DATABASE_URL to set connection parameters
- Configures connection pooling with reasonable defaults
- Sets up Redis for rate limiting and caching

5. **Update the app.py imports**:

The `app.py` file should already include:
```python
# Import Heroku config
import config_heroku
```

6. **Install additional Python dependencies**:

Ensure all scientific packages that may be needed are properly listed in `requirements.txt`:
```bash
# Update requirements.txt to include numpy and other scientific packages
# if they're not already there
```

7. **Set required environment variables**:

```bash
heroku config:set \
  FLASK_ENV=production \
  FLASK_APP=app.py \
  DEPLOYMENT_COLOR=production \
  RDKIT_ENABLED=false \
  SECRET_KEY=$(openssl rand -hex 32) \
  SUPABASE_URL=https://your-supabase-url.supabase.co \
  SUPABASE_KEY=your-supabase-key \
  SUPABASE_SERVICE_KEY=your-service-role-key \
  LOG_LEVEL=INFO \
  PYTHONUNBUFFERED=1 \
  --app your-app-name
```

8. **Commit and push the changes to Heroku**:

```bash
git add Procfile setup_database.py config_heroku.py app.py requirements.txt
git commit -m "Migrate to full CryoProtect app on Heroku"
git push heroku heroku-full-app:master
```

9. **Monitor the deployment**:

```bash
heroku logs --tail
```

## Database Add-ons and Configuration

### PostgreSQL Configuration

The full app requires PostgreSQL for database storage. Ensure it's properly configured:

```bash
# Check existing PostgreSQL add-on
heroku addons | grep heroku-postgresql

# If not present, add it
heroku addons:create heroku-postgresql:hobby-dev

# Get connection details
heroku pg:info
```

### Redis for Rate Limiting and Caching

The full app uses Redis for rate limiting and caching:

```bash
# Check existing Redis add-on
heroku addons | grep heroku-redis

# If not present, add it
heroku addons:create heroku-redis:hobby-dev

# Get connection details
heroku redis:info
```

## RDKit Configuration

For molecular modeling features, you have two options:

1. **Disable RDKit features** (simplest approach for Heroku):
   ```bash
   heroku config:set RDKIT_ENABLED=false
   ```

2. **Use an external RDKit service**:
   ```bash
   heroku config:set RDKIT_ENABLED=true
   heroku config:set RDKIT_SERVICE_URL=your_rdkit_service_url
   ```

## Verification Steps

After deployment, verify the application is running correctly:

1. **Check the process status**:
   ```bash
   heroku ps
   ```

2. **Verify the web interface**:
   ```bash
   heroku open
   ```

3. **Check the health endpoint**:
   ```bash
   curl $(heroku info -s | grep web_url | cut -d= -f2)health
   ```

4. **Verify database connectivity**:
   ```bash
   curl $(heroku info -s | grep web_url | cut -d= -f2)health/readiness
   ```

5. **Monitor application logs**:
   ```bash
   heroku logs --tail
   ```

## Troubleshooting

### Application Errors

If the application shows errors or won't start:

1. **Check for missing dependencies**:
   - Review logs for import errors
   - Add missing packages to requirements.txt
   - Redeploy with updated requirements

2. **Check for environment variable issues**:
   - Verify all required variables are set
   - Check for typos in variable names or values
   ```bash
   heroku config
   ```

3. **Database connection problems**:
   - Check the database provisioning status
   - Verify connection parameters
   ```bash
   heroku pg:info
   ```

4. **Manual database initialization**:
   If the automatic setup fails, run it manually:
   ```bash
   heroku run python setup_database.py
   ```

### Common Error Messages

1. **ModuleNotFoundError: No module named 'numpy'**:
   - Add numpy and other scientific packages to requirements.txt
   - Redeploy the application

2. **Database Connection Errors**:
   - Check if the PostgreSQL add-on is provisioned
   - Verify the DATABASE_URL environment variable
   - Restart the application: `heroku ps:restart`

3. **Application crashes after deployment**:
   - Check for exceptions in the logs
   - Verify dyno resources are sufficient
   - Restart the application: `heroku ps:restart`

## Scaling and Resource Management

For better performance of the full app:

```bash
# Scale web dynos
heroku ps:scale web=1

# Upgrade to a larger dyno if needed
heroku ps:resize web=standard-1x

# Upgrade database if needed
heroku addons:upgrade heroku-postgresql:hobby-basic
```

## Rollback Procedure

If you need to roll back to the simplified app:

```bash
# Switch back to the simple app branch
git checkout main

# Update the Procfile to use the simple app
echo "web: gunicorn simple_app:app --log-file -" > Procfile

# Commit and push
git add Procfile
git commit -m "Revert to simplified app"
git push heroku main:master
```

## Production Recommendations

For production deployments:

1. **Enable automatic backups**:
   ```bash
   heroku pg:backups:schedule DATABASE_URL --at '04:00 UTC'
   ```

2. **Configure SSL**:
   ```bash
   heroku certs:auto:enable
   ```

3. **Add monitoring**:
   ```bash
   heroku addons:create librato
   heroku addons:create papertrail
   ```

4. **Set up a staging environment**:
   ```bash
   heroku create cryoprotect-staging --remote staging
   ```

5. **Configure health checks**:
   ```bash
   heroku features:enable http-health-checks
   heroku config:set ENABLE_HEALTH_CHECKS=true
   ```