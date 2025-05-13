#!/bin/bash
# Script to deploy the frontend to Vercel with environment variables

# Generate a random NEXTAUTH_SECRET if not provided
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Deploy to Vercel
vercel deploy --prod --yes \
  -e NEXT_PUBLIC_API_URL=https://api.cryoprotect.app/v1 \
  -e NEXTAUTH_URL=https://frontend-cryoprotect.vercel.app \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET"