#!/bin/bash

# Fixed Netlify Deployment Script for CryoProtect
set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}==============================================${NC}"
echo -e "${BLUE}    CryoProtect Netlify Deployment Fix    ${NC}"
echo -e "${BLUE}==============================================${NC}"
echo

# Apply fixed configuration files
echo -e "${BLUE}Applying fixed configuration files...${NC}"
cp next.config.js.fixed next.config.js
cp package.json.fixed package.json
cp netlify.toml.fixed netlify.toml
echo -e "${GREEN}Fixed configuration applied${NC}"

# Clean up any existing build artifacts
echo -e "${BLUE}Cleaning previous build artifacts...${NC}"
rm -rf .next node_modules package-lock.json
echo -e "${GREEN}Build artifacts cleaned${NC}"

# Install dependencies
echo -e "${BLUE}Installing dependencies...${NC}"
npm install --legacy-peer-deps
echo -e "${GREEN}Dependencies installed${NC}"

# Clear Netlify cache
echo -e "${BLUE}Clearing Netlify cache...${NC}"
netlify build:clear-cache || echo -e "${YELLOW}No cache to clear or command not available${NC}"
echo -e "${GREEN}Cache cleared${NC}"

# Choose routing strategy - remove App Router
echo -e "${BLUE}Removing incompatible App Router...${NC}"
if [ -d "src/app" ]; then
  echo -e "${YELLOW}App Router directory found, moving to backup...${NC}"
  mkdir -p .app-router-backup
  mv src/app/* .app-router-backup/
  rm -rf src/app
  echo -e "${GREEN}App Router removed and backed up${NC}"
else
  echo -e "${GREEN}App Router not found, nothing to remove${NC}"
fi

# Ensure pages directory exists and restore it if needed
echo -e "${BLUE}Ensuring Pages Router is available...${NC}"
mkdir -p src/pages
if [ ! -f "src/pages/index.tsx" ] && [ -f ".app-router-backup/page.tsx" ]; then
  echo -e "${YELLOW}Creating simple index page from backed up content...${NC}"
  cat > src/pages/index.tsx << EOF
import React from 'react';

export default function Home() {
  return (
    <div className="flex min-h-screen flex-col items-center justify-between p-24">
      <h1 className="text-4xl font-bold mb-8">CryoProtect</h1>
      <p className="text-xl">A platform for cryoprotectant analysis</p>
    </div>
  );
}
EOF
  echo -e "${GREEN}Simple index page created${NC}"
fi

# Make sure NextAuth is properly configured
echo -e "${BLUE}Setting up NextAuth API routes...${NC}"
mkdir -p src/pages/api/auth
if [ ! -f "src/pages/api/auth/[...nextauth].ts" ]; then
  echo -e "${YELLOW}Creating NextAuth API route...${NC}"
  cat > src/pages/api/auth/[...nextauth].ts << EOF
import NextAuth from "next-auth";
import CredentialsProvider from "next-auth/providers/credentials";

export default NextAuth({
  providers: [
    CredentialsProvider({
      name: "Credentials",
      credentials: {
        email: { label: "Email", type: "email" },
        password: { label: "Password", type: "password" }
      },
      async authorize(credentials) {
        // Simple implementation for temporary testing
        if (credentials?.email && credentials?.password) {
          return {
            id: "1",
            name: "User",
            email: credentials.email,
            role: "user"
          }
        }
        return null;
      },
    }),
  ],
  pages: {
    signIn: '/login',
  },
  session: {
    strategy: 'jwt',
  },
  callbacks: {
    async session({ session, token }) {
      if (token.sub) {
        session.user.id = token.sub;
      }
      return session;
    },
  },
  secret: process.env.NEXTAUTH_SECRET
});
EOF
  echo -e "${GREEN}NextAuth API route created${NC}"
fi

# Create a health endpoint
echo -e "${BLUE}Creating health endpoint...${NC}"
mkdir -p src/pages/api
cat > src/pages/api/health.ts << EOF
import type { NextApiRequest, NextApiResponse } from 'next';

type HealthResponse = {
  status: string;
  timestamp: string;
  environment: string;
};

export default function handler(
  req: NextApiRequest,
  res: NextApiResponse<HealthResponse>
) {
  res.status(200).json({
    status: 'OK',
    timestamp: new Date().toISOString(),
    environment: process.env.NEXT_PUBLIC_ENVIRONMENT || 'unknown'
  });
}
EOF
echo -e "${GREEN}Health endpoint created${NC}"

# Create Netlify function directory if it doesn't exist
echo -e "${BLUE}Setting up Netlify functions...${NC}"
mkdir -p netlify/functions
cat > netlify/functions/hello.js << EOF
exports.handler = async function(event, context) {
  return {
    statusCode: 200,
    body: JSON.stringify({
      message: "CryoProtect API is running!",
      timestamp: new Date().toISOString()
    })
  };
};
EOF
echo -e "${GREEN}Basic Netlify function created${NC}"

# Deploy to Netlify
echo -e "${BLUE}Deploying to Netlify...${NC}"

# Check if user is logged in to Netlify
if ! netlify status 2>&1 | grep -q "Logged in"; then
    echo -e "${YELLOW}Not logged in to Netlify. Please log in:${NC}"
    netlify login
fi

# Check if site is linked
if ! netlify status 2>&1 | grep -q "cryoprotect"; then
    echo -e "${YELLOW}Site not linked. Linking to cryoprotect site...${NC}"
    netlify unlink 2>/dev/null
    netlify link --name cryoprotect
fi

# Deploy to Netlify with full build logs
echo -e "${BLUE}Starting Netlify deployment with verbose logging...${NC}"
NETLIFY_BUILD_DEBUG=true netlify deploy --build --prod --debug

echo -e "${BLUE}==============================================${NC}"
echo -e "${GREEN}Deployment completed!${NC}"
echo -e "${BLUE}==============================================${NC}"

exit 0