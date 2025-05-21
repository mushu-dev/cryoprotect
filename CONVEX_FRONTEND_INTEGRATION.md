# Convex Frontend Integration Guide

This document provides a comprehensive guide on how to integrate the Convex database with the frontend of CryoProtect. The integration enables real-time collaboration features and improved data synchronization.

## Overview

The CryoProtect frontend now supports two data sources:

1. **Standard API** - The original Flask API that connects to Supabase
2. **Convex API** - Direct connection to Convex for real-time features

The integration is designed to be non-disruptive, allowing for a gradual migration from the original API to Convex. The application can operate in three modes:

- **Standard Mode** - Uses only the original API
- **Convex Mode** - Uses only the Convex API
- **Hybrid Mode** - Uses both APIs, with Convex for collaborative features

## Configuration

### Environment Variables

To enable Convex integration, set the following environment variables in your `.env.local` file:

```bash
# Enable Convex integration
NEXT_PUBLIC_USE_CONVEX=true

# Convex project URL
NEXT_PUBLIC_CONVEX_URL=https://your-instance.convex.cloud

# Convex deployment key (for backend operations)
CONVEX_DEPLOYMENT_KEY=your-deployment-key-here

# Optional: Enable hybrid mode (both APIs)
ENABLE_HYBRID_MODE=true
```

A template file is provided at `.env.convex.template` which you can copy to `.env.local` and customize.

### Running in Convex Mode

To run the application in Convex mode:

1. Create a `.env.local` file with the configuration above
2. Start the frontend development server:

```bash
cd frontend
npm run dev:with-convex
```

This will start the Next.js development server with Convex integration enabled.

## Key Components

### 1. ConvexClientProvider

The `ConvexClientProvider` component provides the Convex client context to the application. It is conditionally included based on the `NEXT_PUBLIC_USE_CONVEX` environment variable.

```jsx
// Located at: frontend/src/convex/ConvexClientProvider.js
import { ConvexProvider } from 'convex/react';
import { convex } from './client';

export function ConvexClientProvider({ children }) {
  return (
    <ConvexProvider client={convex}>
      {children}
    </ConvexProvider>
  );
}
```

### 2. useConvexExperimentData Hook

This custom hook provides real-time access to experiment data from Convex:

```jsx
// Located at: frontend/src/features/experiments/hooks/useConvexExperimentData.js
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';

export default function useConvexExperimentData(initialFilters = {}, options = {}) {
  // Real-time query for experiments with filters
  const experiments = useQuery(
    api.experiments.enhanced_experiments.listEnhancedExperiments,
    { filter: filters, options: queryOptions }
  ) || [];

  // Convex mutations for experiment management
  const createExperiment = useMutation(api.experiments.enhanced_experiments.createEnhancedExperiment);
  const updateExperiment = useMutation(api.experiments.enhanced_experiments.updateEnhancedExperiment);
  
  // ...additional code...
}
```

### 3. CollaborativeExperimentEditor

This component provides real-time collaborative editing for experiment data:

```jsx
// Located at: frontend/src/components/CollaborativeExperimentEditor.js
import { useCollaborativeDocument } from '../hooks/useConvexSubscription';
import { useMutation } from 'convex/react';
import { api } from '../convex/_generated/api';

const CollaborativeExperimentEditor = ({ experimentId, onSave }) => {
  // Use the collaborative document hook
  const {
    document: experiment,
    collaborators,
    updateDocument,
    resolveConflicts
  } = useCollaborativeDocument(experimentId, 'experiments');
  
  // ...additional code...
}
```

### 4. CollaborativeExperimentEditorWrapper

A wrapper component that conditionally renders either the collaborative editor or the standard experiment detail component:

```jsx
// Located at: frontend/src/features/experiments/components/CollaborativeExperimentEditorWrapper.jsx
export default function CollaborativeExperimentEditorWrapper({ 
  experiment, 
  onSave,
  readOnly = false
}) {
  const [isConvexEnabled, setIsConvexEnabled] = useState(false);
  
  // ...additional code...
  
  // Render collaborative editor or standard detail view
  if (isConvexEnabled && isEditing) {
    return <CollaborativeExperimentEditor experimentId={experiment._id} onSave={handleSave} />;
  } else {
    return <ExperimentDetail experiment={experiment} />;
  }
}
```

## Integration Points

### 1. Experiment List Page

The experiment list page now supports both data sources:

```jsx
// frontend/src/features/experiments/components/ExperimentsList.js
export default function ExperimentsList() {
  // Environment configuration for determining data source
  const useConvex = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';
  
  // Standard API data hook
  const standardData = useExperimentData();
  
  // Convex real-time data hook
  const convexData = useConvexExperimentData();

  // Use the appropriate data source based on configuration
  const { 
    experiments, 
    filteredExperiments, 
    // ...other properties
  } = useConvex ? convexData : standardData;
  
  // ...additional code...
}
```

### 2. Experiment Detail Page

The experiment detail page uses the `CollaborativeExperimentEditorWrapper` to handle both viewing and editing:

```jsx
// frontend/src/pages/experiments/[id].js
export default function ExperimentDetailPage() {
  // ...code to fetch experiment...
  
  return (
    <div className="container mx-auto px-4 py-8">
      <CollaborativeExperimentEditorWrapper 
        experiment={experiment}
        onSave={(updatedExperiment) => {
          setExperiment(updatedExperiment);
        }}
        readOnly={false}
      />
    </div>
  );
}
```

### 3. Create Experiment Page

The create experiment page uses the appropriate creation method based on configuration:

```jsx
// frontend/src/pages/experiments/create.js
const handleSubmit = async (e) => {
  // ...validation code...
  
  try {
    // Use appropriate creation method based on configuration
    if (useConvex && createConvexExperiment) {
      // Create experiment in Convex
      const experimentId = await createConvexExperiment(formData);
    } else if (createStandardExperiment) {
      // Create experiment using standard API
      const experimentId = await createStandardExperiment(formData);
    }
    
    // ...additional code...
  } catch (error) {
    // ...error handling...
  }
};
```

## Real-time Collaboration Features

With Convex integration, users can collaborate in real-time on experiment editing. Key features include:

1. **Real-time Updates**: Changes made by one user are instantly visible to others
2. **Conflict Resolution**: If two users edit the same field, the system detects conflicts and provides resolution options
3. **Presence Indicators**: Visual indicators show who else is currently viewing/editing the experiment
4. **Offline Support**: Changes made offline are synchronized when connection is restored
5. **Change History**: A history of changes is maintained for auditing purposes

## Data Flow

### Standard Mode (Original API)

```
User → Next.js Frontend → Flask API → Supabase → PostgreSQL
```

### Convex Mode

```
User → Next.js Frontend → Convex Client → Convex Cloud → Convex Database
```

### Hybrid Mode

```
User → Next.js Frontend → Convex Client → Convex Cloud → Convex Database
      ↓
      Flask API → Supabase → PostgreSQL
```

## Testing the Integration

To test the Convex integration:

1. Set up your environment variables as described above
2. Start the frontend in Convex mode: `npm run dev:with-convex`
3. Open the application in two browser windows
4. Navigate to the same experiment in both windows
5. Make changes in one window and observe real-time updates in the other

## Troubleshooting

### Common Issues

1. **Convex Client Not Initialized**: Ensure `NEXT_PUBLIC_USE_CONVEX` is set to `true` and `NEXT_PUBLIC_CONVEX_URL` is correctly set.

2. **API Errors**: Check the browser console for error messages. Ensure your Convex deployment is running.

3. **Missing Types**: If you see type errors related to Convex, run `npx convex codegen` to generate the latest types.

4. **Authentication Issues**: Ensure your Convex project is properly configured for authentication.

### Debugging

Enable debug logging by adding `DEBUG=convex:*` to your environment:

```bash
DEBUG=convex:* npm run dev:with-convex
```

## Next Steps

1. **Complete Frontend Integration**: Expand Convex integration to all frontend components.
2. **Authentication Integration**: Integrate Convex authentication with the existing auth system.
3. **Data Migration**: Create tools to migrate data from Supabase to Convex.
4. **Performance Optimization**: Optimize Convex queries and mutations for better performance.

## Additional Resources

- [Convex Documentation](https://docs.convex.dev/)
- [Next.js Integration Guide](https://docs.convex.dev/client/react/next)
- [Real-time Collaboration Patterns](https://docs.convex.dev/patterns/realtime-collaboration)