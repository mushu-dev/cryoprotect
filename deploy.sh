#!/bin/bash
# Script to deploy both frontend and backend to separate Vercel projects

# Set some colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${YELLOW}CryoProtect Deployment Script${NC}"
echo "This script will deploy both frontend and backend to Vercel"
echo

# Ask which component to deploy
echo "Which component do you want to deploy?"
echo "1. Frontend only"
echo "2. Backend only"
echo "3. Both frontend and backend"
read -p "Enter your choice (1-3): " choice

# Check if Vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo -e "${RED}Error: Vercel CLI is not installed${NC}"
    echo "Please install it with: npm install -g vercel"
    exit 1
fi

# Ensure we're logged in to Vercel
vercel whoami &> /dev/null
if [ $? -ne 0 ]; then
    echo -e "${YELLOW}You need to log in to Vercel first${NC}"
    vercel login
    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to log in to Vercel${NC}"
        exit 1
    fi
fi

# Deploy frontend
deploy_frontend() {
    echo -e "${GREEN}Deploying frontend...${NC}"
    cd frontend

    # Ensure NEXTAUTH_SECRET is set
    if [ -z "$NEXTAUTH_SECRET" ]; then
        echo -e "${YELLOW}NEXTAUTH_SECRET is not set, generating a new one...${NC}"
        export NEXTAUTH_SECRET=$(openssl rand -base64 32)
        echo "Using NEXTAUTH_SECRET: $NEXTAUTH_SECRET"
        echo "Store this value securely for future deployments!"
    fi

    # Build the frontend
    echo "Building the frontend..."
    npm run build

    # Deploy to Vercel
    echo "Deploying to Vercel..."
    vercel deploy --prod \
      -e NEXT_PUBLIC_API_URL=https://api.cryoprotect.app/v1 \
      -e NEXTAUTH_URL=https://www.cryoprotect.app \
      -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET"

    cd ..
    echo -e "${GREEN}Frontend deployment completed!${NC}"
}

# Deploy backend
deploy_backend() {
    echo -e "${GREEN}Deploying backend...${NC}"
    
    # Deploy the backend to Vercel
    vercel deploy --prod \
      -e FLASK_ENV=production \
      -e FLASK_APP=app.py \
      -e API_VERSION=v1

    echo -e "${GREEN}Backend deployment completed!${NC}"
}

# Execute the deployment based on user's choice
case $choice in
    1)
        deploy_frontend
        ;;
    2)
        deploy_backend
        ;;
    3)
        deploy_backend
        echo
        deploy_frontend
        ;;
    *)
        echo -e "${RED}Invalid choice. Exiting.${NC}"
        exit 1
        ;;
esac

echo -e "${GREEN}Deployment completed successfully!${NC}"