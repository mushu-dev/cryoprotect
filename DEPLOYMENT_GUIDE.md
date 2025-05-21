# CryoProtect Deployment Guide

This guide provides instructions for deploying the CryoProtect project to various hosting platforms.

## Backend (Heroku)

The backend is already deployed and available at:
```
https://cryoprotect-8030e4025428.herokuapp.com
```

You can test the backend directly:
- Health check: https://cryoprotect-8030e4025428.herokuapp.com/health
- API connect: https://cryoprotect-8030e4025428.herokuapp.com/api/connect
- List molecules: https://cryoprotect-8030e4025428.herokuapp.com/api/molecules

## Frontend Deployment Options

### Option 1: Deploy to Netlify (Recommended)

Netlify offers free hosting without authentication requirements:

1. Create an account at [Netlify](https://netlify.com/)
2. From the dashboard, click "Add new site" → "Import an existing project"
3. Connect your Git repository or use drag-and-drop:
   - Connect to GitHub OR
   - Drag and drop the `public_connection_test` folder
4. If using Git, ensure the `netlify.toml` file is at the root
5. Click "Deploy site"

**Direct Drop Method**:
1. Drag and drop the `public_connection_test` folder onto the Netlify dashboard
2. Netlify will automatically deploy without requiring configuration

### Option 2: Deploy to GitHub Pages

GitHub Pages is another excellent free option:

1. Create a GitHub repository
2. Push your `public_connection_test` folder contents to the repository
3. Go to repository Settings → Pages
4. Under "Source", select "main" branch and "/ (root)" folder
5. Click "Save"
6. Your site will be available at `https://[username].github.io/[repository-name]/`

### Option 3: Configure Vercel (With Team Account)

Vercel personal accounts often require authentication. To bypass this:

1. Create a team on Vercel (free teams are available)
2. Transfer your project to the team:
   - Go to Project Settings → General
   - Scroll to "Transfer Project"
   - Select your team
3. Deploy using the team scope:
   ```bash
   vercel --scope [team-name] --prod --public
   ```

See `vercel_dashboard_guide.md` for detailed instructions on disabling authentication.

## Frontend Development

When developing your frontend application locally, use these settings:

```javascript
// API configuration
const API_URL = 'https://cryoprotect-8030e4025428.herokuapp.com';

// Example API request
fetch(`${API_URL}/api/molecules`)
  .then(response => response.json())
  .then(data => console.log(data));
```

## Connection Testing

For quick connectivity tests:
1. Open `public_connection_test/index.html` in your browser
2. Click "Run Tests" to verify backend connectivity
3. Check the results for successful API connections

## GitHub Auto-Deployment

To set up auto-deployment with GitHub:

1. Connect your GitHub repository to your chosen platform (Netlify/Vercel/GitHub Pages)
2. Configure the build settings according to platform requirements
3. Each push to the main branch will trigger a new deployment

For Vercel specifically, you may need to set `VERCEL_PROTECTION_BYPASS=1` in your environment variables.

## Need Help?

If you're encountering deployment issues:
- Check the platform's documentation
- Verify CORS settings in the backend
- Test API connectivity directly using curl or Postman
- Try an alternative deployment platform from the options above