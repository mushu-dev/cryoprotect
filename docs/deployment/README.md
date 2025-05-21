# RDKit Service Deployment Steps

Follow these steps to deploy the RDKit microservice to Fly.io and connect it to your Heroku application.

## 1. Initial Deployment to Fly.io

```bash
# Navigate to the RDKit service directory
cd /home/mushu/Projects/CryoProtect/rdkit-service

# Log in to Fly.io (this will open a browser)
flyctl auth login

# Launch the app on Fly.io
flyctl launch --dockerfile Dockerfile --name cryoprotect-rdkit

# When prompted for organization, select your organization or create a new one
# When asked if you'd like to deploy now, say 'N' (No)

# Deploy the app manually to ensure everything is set up correctly
flyctl deploy
```

## 2. Verify the Deployment

```bash
# Check application status
flyctl status

# View application logs
flyctl logs

# Test the health endpoint
curl https://cryoprotect-rdkit.fly.dev/health
```

## 3. Configure Environment Variables in Heroku

```bash
# Set the RDKit service URL in your Heroku app
heroku config:set RDKIT_SERVICE_URL=https://cryoprotect-rdkit.fly.dev -a your-heroku-app-name

# Optional: Set an API key for security
API_KEY=$(openssl rand -hex 16)
echo "Generated API key: $API_KEY"

# Set the API key in both applications
flyctl secrets set API_KEY=$API_KEY
heroku config:set RDKIT_API_KEY=$API_KEY -a your-heroku-app-name
```

## 4. Integrate the Client in Your Heroku App

1. Copy the `heroku-client.py` file to your Heroku application

2. Update your main app to use the client:

```python
from heroku_client import RDKitServiceClient

# Initialize the client (reads environment variables automatically)
rdkit_client = RDKitServiceClient()

# Use the client in your application code
properties = rdkit_client.calculate_properties("CCO")
```

## 5. Monitor the Service

```bash
# View real-time logs
flyctl logs

# Open the dashboard in your browser
flyctl dashboard
```

## 6. Scale as Needed

```bash
# Scale the number of VMs (machines)
flyctl scale count 2

# Scale memory
flyctl scale memory 512
```

## 7. Troubleshooting

If you encounter issues:

- Check the Fly.io logs: `flyctl logs`
- SSH into the VM for debugging: `flyctl ssh console`
- Test the service directly: `curl https://cryoprotect-rdkit.fly.dev/health`
- If there are issues with the Docker build, try building locally first:
  ```bash
  docker build -t cryoprotect-rdkit .
  docker run -p 8080:8080 cryoprotect-rdkit
  ```

## 8. Production Readiness Checklist

- [ ] Configure auto-scaling for production workloads
- [ ] Add monitoring/alerting for service uptime
- [ ] Set up regular backups if applicable
- [ ] Configure a custom domain if desired:
  ```bash
  flyctl certs create your-custom-domain.com
  ```
- [ ] Consider adding rate limiting for the API