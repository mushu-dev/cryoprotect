# GitHub Secrets Setup for Heroku Deployment

To enable automated deployment to Heroku via GitHub Actions, you need to set up the following secrets in your GitHub repository.

## Required Secrets

1. **HEROKU_API_KEY**
   - Value: `HRKU-AAlkch0IahH7s6VdUOKQf3YyujOpCQQc68OT5jDxDrfg_____wYCQ1TUD4eO`
   - Purpose: Authenticates GitHub Actions with Heroku

2. **HEROKU_APP_NAME**
   - Value: `cryoprotect`
   - Purpose: Production application name

3. **HEROKU_STAGING_APP_NAME**
   - Value: `cryoprotect-staging`
   - Purpose: Staging application name

4. **HEROKU_REVIEW_APP_NAME**
   - Value: `cryoprotect-review`
   - Purpose: Review application name

5. **HEROKU_EMAIL**
   - Value: `eluecheelip@gmail.com`
   - Purpose: Heroku account email

6. **HEROKU_PIPELINE_ID** (for PR review apps)
   - Value: `f8ae7d38-92ee-4568-a57b-4e4cc3ab8969`
   - Purpose: Allows creation of review apps for pull requests

## How to Add Secrets

1. Go to your GitHub repository
2. Click on **Settings**
3. Click on **Secrets and variables** → **Actions**
4. Click **New repository secret**
5. Add each secret name and value
6. Click **Add secret**

## Creating a Heroku Pipeline

To set up a pipeline for your apps:

1. Run: `heroku pipelines:create cryoprotect-pipeline --app cryoprotect`
2. Add staging app: `heroku pipelines:add cryoprotect-pipeline --app cryoprotect-staging --stage staging`
3. Add review app: `heroku pipelines:add cryoprotect-pipeline --app cryoprotect-review --stage development`
4. Get pipeline ID: `heroku pipelines:info cryoprotect-pipeline --json`

## Additional Notes

- The API key above is created specifically for GitHub Actions
- Make sure to enable review apps in your pipeline settings
- Consider enabling automatic deployments from GitHub for continuous deployment