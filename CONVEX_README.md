# CryoProtect Convex Implementation

This document provides an overview of the CryoProtect implementation using Convex.

## Overview

CryoProtect is being reimplemented using Convex, a cloud-based database and backend platform that provides real-time reactivity. This implementation follows the detailed plan outlined in `CONVEX_IMPLEMENTATION_PLAN.md`.

## Project Structure

The Convex implementation is organized as follows:

```
/convex
  /schema            # Schema definitions
    convex_schema.ts # Main schema file
  /auth              # Authentication functions
    clerk.ts         # Clerk integration
    users.ts         # User management
    client.ts        # Client-side utilities
    index.ts         # Auth exports
  /molecules         # Molecule data functions
    create.ts        # Creation functions
    query.ts         # Query functions
    update.ts        # Update functions
    delete.ts        # Delete functions
    helpers.ts       # Helper utilities
    validation.ts    # Validation functions
    types.ts         # Type definitions
    index.ts         # Molecule exports
  convex.json        # Project configuration
  auth.config.js     # Auth configuration
```

### Configuration Files

**convex.json**
```json
{
  "project": "cryoprotect",
  "team": "cryoprotect-team",
  "functions": "./"
}
```

**auth.config.js** (Development)
```javascript
export default {
  // Allow unauthenticated access during development
  allowAnonymousUsers: true,
}
```

**auth.config.js** (Production with Clerk)
```javascript
export default {
  providers: [
    {
      domain: "https://your-clerk-domain.clerk.accounts.dev/",
      applicationID: "your-clerk-application-id",
    }
  ],
  roles: ["admin", "scientist", "viewer"]
}
```

## Getting Started

### Prerequisites

- Node.js 18+
- Clerk account for authentication
- Convex account for database and functions

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/your-org/cryoprotect.git
   cd cryoprotect
   ```

2. Install dependencies:
   ```bash
   npm install
   ```

3. Configure environment variables:
   ```bash
   cp .env.example .env
   # Edit .env file with your Convex and Clerk credentials
   ```

4. Configure authentication (for production deployment):
   
   Edit the `convex/auth.config.js` file with your Clerk credentials:
   ```javascript
   // convex/auth.config.js
   export default {
     providers: [
       {
         domain: "https://YOUR_CLERK_DOMAIN.clerk.accounts.dev/",
         applicationID: "YOUR_CLERK_APPLICATION_ID",
       }
     ],
     roles: ["admin", "scientist", "viewer"]
   };
   ```

5. Start the development server:
   ```bash
   npm run convex:dev
   ```

### Frontend Integration

To integrate the Convex backend with your frontend:

1. Install required packages in your frontend project:
   ```bash
   cd frontend
   npm install convex @clerk/clerk-react convex-react-clerk
   ```

2. Set up environment variables in your frontend:
   ```
   NEXT_PUBLIC_CONVEX_URL=your_convex_deployment_url
   NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY=your_clerk_publishable_key
   ```

3. Create a Convex client configuration in your frontend:
   ```typescript
   // convex/client.ts
   import { ConvexReactClient } from "convex/react";
   import { ClerkProvider, useAuth } from "@clerk/clerk-react";
   import { ConvexProviderWithClerk } from "convex/react-clerk";

   export function createConvexClient() {
     return new ConvexReactClient(process.env.NEXT_PUBLIC_CONVEX_URL || "");
   }
   ```

4. Wrap your application with the Convex and Clerk providers:
   ```typescript
   // pages/_app.tsx or app/layout.tsx
   import { ClerkProvider } from "@clerk/clerk-react";
   import { ConvexProviderWithClerk } from "convex/react-clerk";
   import { createConvexClient } from "../convex/client";
   import { useAuth } from "@clerk/clerk-react";

   const convex = createConvexClient();

   export default function App({ Component, pageProps }) {
     return (
       <ClerkProvider 
         publishableKey={process.env.NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY}
       >
         <ConvexProviderWithClerk client={convex} useAuth={useAuth}>
           <Component {...pageProps} />
         </ConvexProviderWithClerk>
       </ClerkProvider>
     );
   }
   ```

5. Use Convex functions in your components:
   ```typescript
   // components/MoleculeList.tsx
   import { useQuery, useMutation } from "convex/react";
   import { api } from "../convex/_generated/api";

   export function MoleculeList() {
     const molecules = useQuery(api.molecules.query.searchMolecules) || [];
     const createMolecule = useMutation(api.molecules.create.createMolecule);

     const handleCreateMolecule = async () => {
       await createMolecule({
         input: {
           name: "New Molecule",
           formula: "C6H12O6",
         },
       });
     };

     return (
       <div>
         <button onClick={handleCreateMolecule}>Add Molecule</button>
         <ul>
           {molecules.molecules?.map((molecule) => (
             <li key={molecule._id.toString()}>{molecule.name}</li>
           ))}
         </ul>
       </div>
     );
   }
   ```

## Key Features

### 1. Document-Oriented Database

The Convex implementation uses a document-oriented approach, optimizing for:

- Scientific data integrity
- Performance for critical operations
- Flexible schema evolution
- Real-time reactivity for collaborative features

### 2. Authentication

Authentication is implemented using Clerk, providing:

- Secure user management
- Role-based access control
- JWT-based authentication
- Seamless client integration

### 3. Molecule Management

The core molecule system provides:

- Comprehensive CRUD operations
- Validation for scientific data
- Molecule consolidation capabilities
- Efficient querying and filtering

## Testing

The test suite for the Convex implementation is located in the `/tests/convex` directory. This separation keeps the tests outside of the Convex deployment, preventing any test-related dependencies from interfering with the deployment process.

Run the test suite to verify the implementation:

```bash
npm run test:convex
```

For development with real-time test feedback:

```bash
npm run test:convex:watch
```

### Test Structure

- `/tests/convex/` - Contains all test files for Convex functions
  - `setup.ts` - Test setup and utilities
  - `molecules.test.ts` - Tests for molecule functions
  - `auth.test.ts` - Tests for authentication functions

### Writing New Tests

When adding new tests:

1. Create test files in the `/tests/convex/` directory
2. Import the functions you want to test from the Convex implementation
3. Use the mock context from `setup.ts` to simulate the Convex environment
4. Run tests using the npm scripts

## Implementation Progress

Initial implementation includes:

- [x] Schema design
- [x] Core molecule data model
- [x] Authentication integration
- [x] Basic CRUD operations
- [x] Test suite

Next steps:

- [ ] Property management
- [ ] Mixture composition
- [ ] Scientific calculations
- [ ] Migration tools
- [ ] Frontend integration

## Contributing

See the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on contributing to this project.

## Migration Strategy

The migration from Supabase to Convex will follow these steps:

1. Complete the Convex implementation
2. Implement migration utilities
3. Perform controlled test migration
4. Validate scientific calculations
5. Execute production migration
6. Operate in parallel temporarily
7. Complete cutover

Refer to the `CONVEX_IMPLEMENTATION_PLAN.md` for detailed migration information.