# CryoProtect Minimal Deployment

This is a minimal version of the CryoProtect frontend, created to diagnose deployment issues with Netlify.

## Structure
- Uses basic Next.js setup
- Simple single page with no complex dependencies
- Minimalist Netlify configuration
- No API connections or authentication

## Deployment
Run the deployment script:
```bash
./deploy-minimal.sh
```

This will build and deploy the minimal app to Netlify.

## Verification
Once deployed, verify the site is accessible at the provided Netlify URL.

## Next Steps
After confirming the minimal deployment works:
1. Gradually add back features
2. Test each addition to identify what causes issues
3. Fix the specific component or configuration causing the problem
4. Apply the fix to the main application