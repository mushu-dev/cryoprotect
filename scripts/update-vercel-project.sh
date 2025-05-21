#!/bin/bash
# Script to update Vercel project configuration for GitHub organization integration

# Check if the Vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo "Vercel CLI is not installed. Please install it using: npm install -g vercel"
    exit 1
fi

# Ensure user is logged in to Vercel
vercel whoami &> /dev/null
if [ $? -ne 0 ]; then
    echo "Please log in to Vercel first using: vercel login"
    exit 1
fi

# Set variables
ORG_NAME="blueprint-house"
REPO_NAME="CryoProtect"
FRONTEND_DIR="frontend"
PROJECT_NAME="cryoprotect1"  # Use the project name from vercel project ls

echo "Updating Vercel project for GitHub organization repository..."

# Link the project
echo "Linking Vercel project to GitHub repository..."
cd "$FRONTEND_DIR" || (echo "Frontend directory not found"; exit 1)

# Check if project is already linked
if [ -f ".vercel/project.json" ]; then
    echo "Project is already linked. Checking configuration..."
    
    # Get the current project ID
    PROJECT_ID=$(cat .vercel/project.json | grep -o '"projectId": "[^"]*"' | cut -d'"' -f4)
    
    # Print current configuration
    echo "Current project ID: $PROJECT_ID"
    
    # Get project details
    vercel project ls
    
    # Offer to unlink
    read -p "Do you want to unlink and relink the project? (y/N) " UNLINK
    if [[ "$UNLINK" =~ ^[Yy]$ ]]; then
        vercel unlink
        echo "Project unlinked. Proceeding to relink..."
    else
        echo "Keeping current link. Updating settings..."
    fi
fi

# Link project if not already linked or after unlinking
if [ ! -f ".vercel/project.json" ]; then
    echo "Linking project..."
    vercel link -p "$PROJECT_NAME"
fi

# Set up environment variables
echo "Setting up environment variables..."
vercel env add NEXT_PUBLIC_API_URL
vercel env add NEXT_PUBLIC_RDKIT_SERVICE_URL
vercel env add NEXTAUTH_URL
vercel env add NEXTAUTH_SECRET

# Configure project settings for the GitHub organization
echo "Configuring project for GitHub organization..."
vercel git connect "github" "$ORG_NAME/$REPO_NAME"

# Set up build configuration
echo "Setting up build configuration..."
vercel --prod

echo "Vercel project has been updated to use the GitHub organization repository."
echo "Please check the Vercel dashboard to verify the configuration."
echo ""
echo "Next steps:"
echo "1. Go to the Vercel dashboard: https://vercel.com/dashboard"
echo "2. Select the project: $PROJECT_NAME"
echo "3. Go to Settings > Git"
echo "4. Verify the GitHub organization repository is correctly connected"
echo "5. Check the automatic deployments for both production and preview environments"