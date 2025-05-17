# CryoProtect Netlify Deployment

This branch contains pre-built static files for direct deployment to Netlify.

## Features
- Uses pre-built static files in the /out directory
- Configured with proper redirects for dynamic routes
- Includes security headers
- API routes forwarded to the Heroku backend

## How to deploy

Connect this branch to Netlify and it will automatically deploy the static site.

1. In the Netlify UI, go to "Site settings" > "Build & deploy" > "Continuous deployment"
2. Connect to your GitHub repository and select this branch (netlify-deployment)
3. The deployment should happen automatically without running a build command

## Environment Variables

Make sure to set these environment variables in the Netlify UI:

- NEXT_PUBLIC_NETLIFY=true
- NEXT_PUBLIC_API_URL=https://cryoprotect-8030e4025428.herokuapp.com/api
- NEXTAUTH_URL=https://[your-site-name].netlify.app
- NEXTAUTH_SECRET=[your-auth-secret]
- NEXT_PUBLIC_SUPABASE_URL=[your-supabase-url]
- NEXT_PUBLIC_SUPABASE_ANON_KEY=[your-supabase-anon-key]
