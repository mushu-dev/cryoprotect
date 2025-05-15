# CryoProtect Deployment Steps

This document outlines the steps to complete the repository transfer from `blueprint-house` to `mushu-dev` and fix the deployment issues.

## Repository Transfer

1. **Transfer the repository on GitHub**
   - Go to: https://github.com/blueprint-house/CryoProtect/settings
   - Scroll down to the "Danger Zone" section and click "Transfer"
   - Enter "blueprint-house/CryoProtect" as the repository name
   - Enter "mushu-dev" as the new owner
   - Click "I understand, transfer this repository"

2. **Update Local Repository**
   - After transfer is complete, run the transfer script:
   ```bash
   cd /home/mushu/Projects/CryoProtect
   ./scripts/transfer-repository.sh
   ```
   - This will:
     - Update all GitHub URLs in project files
     - Update GitHub Actions workflow files
     - Update Vercel deployment settings
     - Update Heroku deployment settings

## Deployment Platform Updates

### Vercel Deployment

1. **Update Vercel Project Configuration**
   - Log in to Vercel dashboard: https://vercel.com/dashboard
   - Select the CryoProtect project
   - Go to "Settings" > "Git"
   - Disconnect the current GitHub repository (blueprint-house)
   - Connect to the new repository (mushu-dev/CryoProtect)
   - Update environment variables if necessary
   
2. **Check Frontend-Backend Communication**
   - The frontend is configured to look for the API at `https://api.cryoprotect.app/v1`
   - Make sure this URL is correct and the backend is running at this address
   - Check that CORS headers are correctly set in both frontend and backend

### Heroku Deployment

1. **Update Heroku App Configuration**
   - Log in to Heroku dashboard: https://dashboard.heroku.com/apps
   - Select the CryoProtect app
   - Go to "Deploy" tab
   - Under "Deployment method", connect to the new GitHub repository (mushu-dev/CryoProtect)
   - Enable automatic deploys for the master branch
   
2. **Verify Heroku Environment Variables**
   - Ensure all necessary environment variables are set in Heroku dashboard
   - Key variables include:
     - `DATABASE_URL`: Connection to the Supabase database
     - `SECRET_KEY`: For Flask application security
     - `RDKIT_SERVICE_URL`: URL for the RDKit microservice on Fly.io

### Fly.io Deployment for RDKit Microservice

1. **Verify Fly.io Deployment**
   - Ensure the GitHub Actions workflow for Fly.io deployment is updated with the new repository
   - Verify that the FLY_API_TOKEN secret is set in the GitHub repository settings
   - Run a manual deployment to test the setup:
   ```bash
   cd rdkit-service-files
   flyctl deploy
   ```

2. **Check RDKit Service Connectivity**
   - The RDKit service should be available at the URL specified in Heroku environment variables
   - Test the connectivity from the backend by calling the health endpoint:
   ```bash
   curl https://<rdkit-app-name>.fly.dev/health
   ```

## Testing and Verification

1. **Test Complete Deployment Chain**
   - Frontend (Vercel) → Backend (Heroku) → RDKit Service (Fly.io) → Database (Supabase)
   
2. **Verify Functionality**
   - Test key features that require all components to work together:
     - User authentication
     - Molecule property calculations (requires RDKit service)
     - Data retrieval from database
     
3. **Check Logs and Monitoring**
   - Verify logs in Vercel, Heroku, and Fly.io for any errors
   - Check for any CORS or connectivity issues
   
4. **Troubleshooting Common Issues**
   - **CORS errors**: Ensure CORS headers in frontend and backend match
   - **API connectivity**: Check network requests in browser dev tools
   - **Authentication issues**: Verify JWT token generation and validation
   - **RDKit service errors**: Check logs on Fly.io

## GitHub Repository Configuration

1. **Update GitHub Secrets**
   - Transfer all secrets from the old repository to the new one
   - Key secrets include:
     - `HEROKU_API_KEY`
     - `HEROKU_APP_NAME`
     - `HEROKU_EMAIL`
     - `VERCEL_ORG_ID`
     - `VERCEL_PROJECT_ID`
     - `VERCEL_TOKEN`
     - `FLY_API_TOKEN`
     
2. **Set Up Branch Protection**
   - Configure branch protection rules for the master branch
   - Require pull request reviews before merging
   - Require status checks to pass before merging

## Completion Checklist

- [ ] Repository successfully transferred on GitHub
- [ ] All GitHub URLs updated in project files
- [ ] GitHub Actions workflows updated
- [ ] Vercel deployment reconnected and tested
- [ ] Heroku deployment reconnected and tested
- [ ] Fly.io deployment reconnected and tested
- [ ] Frontend can communicate with backend API
- [ ] Backend can communicate with RDKit service
- [ ] All environment variables properly set
- [ ] GitHub secrets transferred to new repository
- [ ] Branch protection rules reconfigured