# Convex Integration Guide for CryoProtect

This document provides a comprehensive guide for integrating Convex with the CryoProtect application, including population scripts, frontend setup, and testing.

## Overview

CryoProtect now uses a microservices architecture with these components:

1. **Frontend**: Next.js application deployed on Netlify (https://cryoprotect.netlify.app)
2. **API Backend**: Flask application deployed on Heroku (https://cryoprotect-8030e4025428.herokuapp.com)
3. **RDKit Service**: Specialized microservice (https://cryoprotect-rdkit.fly.dev)
4. **Convex Database**: Database service (https://hallowed-malamute-424.convex.cloud for dev, https://upbeat-parrot-866.convex.cloud for prod)

## What's Been Done

### 1. CORS Configuration for API Backend

We've implemented standardized CORS configuration for the API backend to allow cross-origin requests from the frontend, RDKit service, and Convex. The CORS configuration:

- Uses environment variables to define allowed origins
- Allows all necessary HTTP methods and headers
- Supports credentials for authenticated requests
- Includes debug and test endpoints for verification

The implementation is in `api/cors_config.py` and has been deployed to the Heroku app.

### 2. Convex Database Integration

We've updated the Heroku backend to use Convex as the database instead of Supabase:

- Created a Convex adapter that provides a Supabase-like interface for backward compatibility
- Implemented a database factory that can switch between Supabase and Convex based on configuration
- Set environment variables for Convex configuration (URL, deployment key, etc.)
- Added diagnostic endpoints to verify the database connection

### 3. Convex Data Population

We've created scripts to populate the Convex database with real data:

- Direct import of molecules data using the `convex import` command:
  ```bash
  npx convex import sample_data.json --table molecules
  npx convex import property_types_data.json --table propertyTypes
  ```

- Custom Python script for generating and importing molecular properties data:
  ```bash
  python create_molecular_properties.py
  ```

- Convex dashboard integration for data management and viewing

### 4. Frontend Convex Integration

We've updated the frontend components to use Convex:

- Added a ConvexClientProvider component that wraps the application when enabled
- Created dynamic routing in the molecules page to use either API or Convex based on configuration
- Implemented a molecules page that uses Convex queries directly
- Added a dedicated run script for working with Convex:
  ```bash
  ./run_with_convex.sh
  ```

### 3. Netlify Configuration

We've updated the Netlify configuration in `frontend/netlify.toml` to:

- Point to the correct Heroku app URL (https://cryoprotect-8030e4025428.herokuapp.com)
- Include Convex URL in environment variables
- Set up redirects for API endpoints
- Configure Content Security Policy headers to allow Convex connections

## Implementation Status

### 1. Data Schema ✅

The Convex schema has been defined in `convex/schema/convex_schema.ts` and includes:

- Molecules
- Molecular Properties
- Property Types
- Mixtures
- Mixture Components
- Experiments
- Enhanced Experiments
- Protocols
- Many more tables matching our data model

### 2. Data Population ✅

We have successfully imported data into the following Convex tables:

- **molecules**: 15 records imported
- **propertyTypes**: 6 records imported
- **molecularProperties**: 90 records imported

The data can be viewed in the Convex dashboard or through our test scripts.

### 3. Frontend Integration ✅

The frontend integration has been implemented with these components:

- **ConvexClientProvider**: Wraps the app when `NEXT_PUBLIC_USE_CONVEX=true`
- **Dynamic Page Routing**: Pages like `molecules/index.js` now route to either the API or Convex implementation
- **Convex Molecule Queries**: Implemented in `frontend/src/convex/molecules.ts`

### 4. TypeScript Build Issues ⏳

We're experiencing some TypeScript build issues when running `npx convex dev`:

```
Unexpected Error: SyntaxError: 'from' expected. (59:31)
  57 | import type * as molecules_update from "../molecules/update.js";
  58 | import type * as molecules_validation from "../molecules/validation.js";
> 59 | import type * as node_modules_@commander_js_extra_typings_esm from "../node_modules/@commander-js/extra-typings/esm.js";
```

These need to be resolved for proper development workflow.

### 5. Deployment Script ✅

A script to run the application with Convex has been created:

```bash
./run_with_convex.sh
```

This script:
- Sets the necessary environment variables 
- Starts the Convex development server
- Starts the Next.js development server with Convex enabled

### 6. Testing Utilities ✅

We've created test scripts to verify our Convex integration:

- **test_convex_integration.js**: Tests the Convex HTTP client integration
- **test_convex_http.py**: Tests direct HTTP access to Convex data
- **create_molecular_properties.py**: Tests and populates molecular property data

## Testing and Verification

### 1. Quick Verification

```bash
./verify-integration.sh
```

This script performs a quick check of all services and reports any issues.

### 2. Test Convex Adapter

```bash
python test-convex-adapter.py
```

This script tests the Convex adapter functionality:
- Basic CRUD operations
- Authentication methods
- Error handling

### 3. Full Integration Test

```bash
node test-full-integration.js
```

This script tests all aspects of the integration:
- Basic connectivity between all services
- CORS configuration
- Netlify redirects
- CSP headers
- End-to-end data flow

## Next Steps

### 1. Resolve TypeScript Build Issues

We need to fix the TypeScript compilation errors in Convex to have a smoother development workflow:

- Investigate TypeScript configuration in Convex
- Resolve path import issues for @commander-js/extra-typings
- Consider switching to official Convex Docker development container
- Seek solutions from Convex community forums

### 2. Complete Frontend Integration

After fixing the TypeScript issues, we should:

- Update all relevant pages to use Convex when enabled
- Test the entire application with Convex
- Implement real-time subscriptions for dynamic updates
- Add Convex-specific features like real-time collaboration

### 3. Authentication Integration

We need to configure Convex to work with our authentication system:

- Set up proper authentication with Convex
- Implement Convex identity management
- Ensure backend roles are properly mapped to Convex roles
- Test comprehensive authentication workflows

### 4. Migration Utilities

To fully transition from Supabase to Convex, we need:

- Comprehensive data migration utilities
- Validation and verification of migrated data
- Rollback mechanisms in case of migration failures
- Operational documentation for migration process

## Development Workflow

### Working with Convex

To work with Convex in the development environment:

1. **Start the Development Environment**:
   ```bash
   ./run_with_convex.sh
   ```
   
2. **Access the Convex Dashboard**:
   - Open https://dashboard.convex.dev/d/hallowed-malamute-424
   - View and manage tables, data, and functions
   
3. **Import Data**:
   ```bash
   npx convex import your_data.json --table your_table --append
   ```
   
4. **Run Tests**:
   ```bash
   python create_molecular_properties.py  # Generates and imports properties
   ```

### Troubleshooting

If you encounter issues with Convex:

1. **Build Errors**:
   - Check TypeScript configuration
   - Ensure all dependencies are installed correctly
   - Verify Node.js version (should be >= 18.0.0)

2. **Import Errors**:
   - Use the `--append` flag when adding to existing tables
   - Check JSON file format for syntax errors
   - Verify table names match the schema

3. **Runtime Errors**:
   - Check Convex logs in the dashboard
   - Verify environment variables are correctly set
   - Check that the frontend has the correct Convex URL

## Conclusion

The Convex integration is a significant step forward for the CryoProtect application. It provides:

- A modern document database with real-time capabilities
- Better schema flexibility for our scientific data model
- Improved frontend integration with React and Next.js
- Forward compatibility with our serverless architecture

The integration is mostly complete, with a few remaining tasks to be addressed. Once the TypeScript build issues are resolved, we can fully deploy and leverage Convex's capabilities.

## References

- [Convex Documentation](https://docs.convex.dev/)
- [Convex React Integration](https://docs.convex.dev/react)
- [Data Modeling in Convex](https://docs.convex.dev/database/schemas)
- [Convex Deployment](https://docs.convex.dev/production/deployment)
- [Convex Auth](https://docs.convex.dev/auth)