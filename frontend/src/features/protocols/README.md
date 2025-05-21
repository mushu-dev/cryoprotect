# Protocol Management Module

This module provides functionality for creating, editing, and managing scientific protocols within the CryoProtect application.

## Components

### Protocol Builder

The Protocol Builder (`protocol-builder.tsx`) provides a comprehensive interface for creating and editing protocols. It includes:

- Basic protocol information (name, description, version, tags)
- Step management with an intuitive interface
- Preview mode for reviewing the complete protocol
- Validation to ensure protocol integrity

Usage:
```tsx
<ProtocolBuilder 
  protocolId={existingProtocolId} // Optional, for editing existing protocols
  onSave={handleSaveProtocol}
  onCancel={handleCancel}
/>
```

### Protocol Step Editor

The Protocol Step Editor (`protocol-step-editor.tsx`) is a specialized component for managing individual protocol steps. It provides:

- Tabbed interface for organizing step information
- Basic information tab (name, description, duration, temperature)
- Advanced settings tab (equipment, parameters)
- Alerts & warnings tab for defining conditions and their severity
- Validation to prevent invalid step data

Usage:
```tsx
<ProtocolStepEditor
  step={stepToEdit}
  isNew={false} // Set to true for new steps
  onSave={handleSaveStep}
  onCancel={handleCancel}
  onDelete={handleDeleteStep}
/>
```

### Protocol Steps Visualization

The Protocol Steps Visualization (`protocol-steps-visualization.tsx`) displays protocol steps in a structured, user-friendly format:

- Sequential step display with accordion expansion
- Summary information for each step
- Equipment, parameter, and alert indicators
- Interactive editing capabilities when in edit mode
- Total duration and temperature range calculations

Usage:
```tsx
<ProtocolStepsVisualization
  protocol={protocol}
  onStepClick={handleStepClick} // Optional, for editing steps
  editable={true} // Optional, enables editing controls
  activeStepId={currentStepId} // Optional, highlights the active step
/>
```

## Hooks

### useProtocols

Provides functionality for retrieving and filtering multiple protocols.

```tsx
const {
  protocols,
  loading,
  error,
  totalCount,
  page,
  perPage,
  totalPages,
  updateParams,
  changePage,
  changePerPage,
  refreshProtocols
} = useProtocols({
  page: 1,
  per_page: 10,
  sort_by: 'created_at',
  sort_order: 'desc'
});
```

### useProtocol

Manages a single protocol, providing CRUD operations.

```tsx
const {
  protocol,
  loading,
  error,
  updateProtocol,
  validateProtocol,
  createVersion,
  exportProtocol,
  publishToLibrary,
  refreshProtocol
} = useProtocol(protocolId);
```

### useProtocolVersions

Manages protocol versions and comparisons.

```tsx
const {
  versions,
  loading,
  error,
  compareVersions,
  refreshVersions
} = useProtocolVersions(protocolId);
```

### useProtocolTemplates

Provides access to protocol templates.

```tsx
const {
  templates,
  loading,
  error,
  createFromTemplate,
  refreshTemplates
} = useProtocolTemplates();
```

### useProtocolLibrary

Accesses shared protocol library.

```tsx
const {
  libraryProtocols,
  loading,
  error,
  getLibraryProtocol,
  refreshLibraryProtocols
} = useProtocolLibrary();
```

### useProtocolCreation

Simplified interface for creating new protocols.

```tsx
const {
  creating,
  error,
  createProtocol,
  importProtocol
} = useProtocolCreation();
```

## Services

### ProtocolService

The `ProtocolService` interface defines all operations that can be performed on protocols. This interface is implemented by two classes:

- `MockProtocolService`: Provides mock data for development and testing
- `ApiProtocolService`: Connects to the backend API for production use

The service implementation is controlled by the `NEXT_PUBLIC_USE_MOCK_DATA` environment variable.

### ApiProtocolService

The `ApiProtocolService` class implements the `ProtocolService` interface and connects to the backend API. It provides:

- Error handling and request formatting
- Authentication via JWT tokens
- Comprehensive API coverage for all protocol operations
- Additional methods for mixture-specific protocol operations

Usage:
```ts
import { ApiProtocolService } from '../services/protocol-service-api';

const protocolService = new ApiProtocolService('/api/v1');
// Or with authentication:
const authProtocolService = new ApiProtocolService('/api/v1', authToken);

// Use service methods
const protocols = await protocolService.getProtocols({ page: 1, per_page: 10 });
```

## Additional Hooks

### useProtocolDesign

Specialized hook for designing protocols for a specific mixture.

```tsx
const {
  loading,
  error,
  designProtocol,
  saveProtocol,
  getSensitivityProfiles
} = useProtocolDesign(mixtureId);
```

### useProtocolsForMixture

Retrieves protocols associated with a specific mixture.

```tsx
const {
  protocols,
  mixtureName,
  loading,
  error,
  refreshProtocols
} = useProtocolsForMixture(mixtureId);
```

### useProtocolComparison

Compares multiple protocols for analysis.

```tsx
const {
  loading,
  error,
  compareProtocols
} = useProtocolComparison();

// Usage:
const comparison = await compareProtocols([protocolId1, protocolId2]);
```

## Recent Updates

- **API Integration**: Added real API service implementation that connects to the backend
- **Additional Hooks**: Created specialized hooks for protocol design, mixture protocols, and comparison
- **Protocol Step Editor Enhancement**: Added a comprehensive step editor with tabbed interface, improved validation, and support for equipment, parameters, and alerts.
- **Protocol Steps Visualization Improvement**: Enhanced visualization with tooltips, equipment indicators, and improved interaction.
- **Unit Tests**: Added comprehensive test coverage for all protocol components and services.

## Best Practices

1. Always validate protocol data before saving
2. Use the step editor for modifying steps rather than direct object manipulation
3. Keep step names concise and descriptive
4. Provide detailed step descriptions for complex procedures
5. Use alerts for critical conditions that require attention
6. Group related steps together in a logical sequence
7. Add appropriate equipment requirements to each step