#!/bin/bash

echo "Building the Next.js application for static export..."
npm run build

# The export is now done as part of the build with output: 'export' in next.config.js
# So we don't need a separate export step anymore

echo "Deploying to Netlify..."
npx netlify deploy --dir=out --prod

echo "Deployment complete!"