#!/bin/bash
# CryoProtect Vercel Deployment Script

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
RESET='\033[0m'

# Print section header
section() {
  echo -e "\n${BLUE}========== $1 ==========${RESET}\n"
}

# Print success message
success() {
  echo -e "${GREEN}✓ $1${RESET}"
}

# Print error message
error() {
  echo -e "${RED}✗ $1${RESET}"
}

# Print warning message
warning() {
  echo -e "${YELLOW}! $1${RESET}"
}

# Print info message
info() {
  echo -e "${PURPLE}i $1${RESET}"
}

# Check if Vercel CLI is installed
section "Checking Prerequisites"
if ! command -v vercel &> /dev/null; then
  error "Vercel CLI is not installed. Please install it with 'npm install -g vercel'"
  exit 1
else
  success "Vercel CLI is installed"
fi

# Update vercel.json to fix routing issues
section "Updating Vercel Configuration"
cat > /home/mushu/Projects/CryoProtect/vercel.json << EOL
{
  "version": 2,
  "buildCommand": "cd frontend && npm ci && npm run build",
  "outputDirectory": "frontend/.next",
  "framework": "nextjs",
  "routes": [
    { 
      "src": "/api/v1/(.*)", 
      "dest": "/api/index.py"
    },
    { 
      "src": "/static/(.*)", 
      "dest": "/static/$1"
    },
    {
      "src": "^/app.py$",
      "status": 404,
      "dest": "/404.html"
    },
    {
      "src": "/_next/(.*)",
      "dest": "/frontend/_next/$1"
    },
    {
      "src": "/assets/(.*)",
      "dest": "/frontend/assets/$1"
    },
    {
      "src": "^/?$",
      "dest": "/frontend/index.html"
    },
    {
      "src": "^/(profile|settings|molecules|mixtures|auth).*",
      "dest": "/frontend/$1"
    },
    {
      "src": "/(.*)",
      "dest": "/frontend/$1"
    }
  ]
}
EOL

success "Updated vercel.json with correct routing configuration"

# Check for deployment option
section "Deployment Options"
echo -e "Choose a deployment method:"
echo -e "1. ${GREEN}Standard deployment${RESET} (uses Next.js server for auth/dynamic routes)"
echo -e "2. ${YELLOW}Static export${RESET} (better compatibility but limited auth functionality)"

read -p "Enter option (1 or 2): " deploy_option

# Build and deploy based on the chosen option
if [ "$deploy_option" == "1" ]; then
  section "Preparing Standard Deployment"
  info "This method uses the Next.js server for API routes and server-side rendering."

  cd /home/mushu/Projects/CryoProtect/frontend
  info "Running standard build..."
  npm run build

  cd /home/mushu/Projects/CryoProtect
  success "Build completed successfully"
  
  info "Deploying to Vercel..."
  vercel --prod
  
elif [ "$deploy_option" == "2" ]; then
  section "Preparing Static Export"
  warning "Static exports don't support API routes or server-side rendering."
  warning "Authentication will be client-side only!"

  cd /home/mushu/Projects/CryoProtect/frontend
  info "Running static export build..."
  ./build-static.sh

  # Update vercel.json for static export
  cd /home/mushu/Projects/CryoProtect
  cat > vercel.json << EOL
{
  "version": 2,
  "buildCommand": "cd frontend && npm ci && ./build-static.sh",
  "outputDirectory": "frontend/out",
  "routes": [
    { 
      "src": "/api/v1/(.*)", 
      "dest": "/api/index.py"
    },
    { 
      "src": "/static/(.*)", 
      "dest": "/static/$1"
    },
    {
      "src": "^/app.py$",
      "status": 404,
      "dest": "/404.html"
    },
    {
      "src": "/_next/(.*)",
      "dest": "/frontend/_next/$1"
    },
    {
      "src": "/assets/(.*)",
      "dest": "/frontend/assets/$1"
    },
    {
      "src": "^/?$",
      "dest": "/frontend/index.html"
    },
    {
      "src": "^/(profile|settings|molecules|mixtures|auth).*",
      "dest": "/frontend/$1.html"
    },
    {
      "src": "/(.*)",
      "dest": "/frontend/$1"
    }
  ]
}
EOL

  success "Updated vercel.json for static export"
  info "Deploying to Vercel..."
  vercel --prod
  
else
  error "Invalid option. Please enter 1 or 2."
  exit 1
fi

section "Deployment Complete"
info "Your CryoProtect app should now be deployed to Vercel."
info "Check the Vercel dashboard for your deployment URL and status."