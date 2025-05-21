# Comprehensive Frontend-Backend Integration Repair Plan

## Overview
This document outlines a complete plan to fix the frontend-backend integration issues for CryoProtect. Our system currently has DNS configured for cryoprotect.app, but the site is showing a 404 error. The goal is to ensure a fully working home page that enables users to view molecules, mixtures, and experiments with proper integration between Heroku, Fly.io, Convex, and Netlify services.

## Current Status
- **DNS**: Properly configured but site returns 404
- **Frontend**: Deployed on Netlify but not properly rendering
- **Backend API**: Running on Heroku (healthy)
- **RDKit Service**: Running on Fly.io (healthy)
- **Convex Database**: Integrated but needs data population
- **Authentication**: NextAuth routes have been deleted and need restoration

## Action Plan

### 1. Fix Netlify Deployment Issues
- [x] Verify DNS configuration with dig/nslookup
- [ ] Run Netlify deployment fix script
- [x] Update Next.js configuration for proper SSR handling
- [x] Configure proper redirects in netlify.toml
- [ ] Deploy and verify site accessibility

### 2. Verify Backend Service Connectivity and Configurations
- [x] Check Heroku backend API health
- [x] Verify RDKit service on fly.io
- [ ] Test Convex database connectivity
- [ ] Update environment variables for all services

### 3. Restore Authentication Functionality
- [ ] Restore NextAuth API routes
- [ ] Configure Convex authentication
- [ ] Test authentication flow end-to-end

### 4. Populate Convex Database with Real Data
- [ ] Create database population script for Convex
- [ ] Migrate data from Supabase to Convex
- [ ] Verify data integrity after migration

### 5. Fix Cross-Origin Resource Sharing (CORS) Configuration
- [ ] Update CSP headers in netlify.toml
- [ ] Configure CORS on backend services
- [ ] Test cross-domain requests

### 6. Implement and Test Dynamic Routes for Molecules and Mixtures
- [ ] Fix generateStaticParams for molecules
- [ ] Fix generateStaticParams for mixtures
- [ ] Test dynamic routes with real data

### 7. End-to-End Testing with Playwright
- [ ] Create comprehensive test suite
- [ ] Test all main user flows
- [ ] Verify integration between all services

## Key Technical Considerations

### Federation of Services
The architecture relies on multiple services:
- Netlify for frontend hosting
- Heroku for API backend
- Fly.io for RDKit service
- Convex for database

Each needs proper configuration and environment variables to work together.

### Configuration Management
- Create a unified config file documenting all service connections
- Use environment variables for all service URLs and API keys
- Ensure all redirects are properly configured in netlify.toml

### Database Factory Pattern
- db_factory.py implements the adapter pattern
- This allows seamless switching between Supabase and Convex
- Use this flexibility to gradually migrate data and validate before full switch

### Authentication Flow
- Ensure the authentication system works with Convex
- NextAuth routes need to be properly restored and configured

## Implementation Timeline
1. **Day 1**: Fix Netlify deployment and verify basic connectivity
2. **Day 2**: Develop and execute database migration scripts
3. **Day 3**: Fix authentication and CORS issues
4. **Day 4**: Implement dynamic routes and comprehensive testing
5. **Day 5**: Final validation and performance optimization