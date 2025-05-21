# Protocol Template Management System

This document provides a comprehensive overview of the Protocol Template Management System, a core component of the CryoProtect platform that enables researchers to create, share, and utilize standardized experimental protocols.

## System Overview

The Protocol Template Management System allows researchers to:

1. **Create standardized protocol templates** with detailed steps, parameters, and metadata
2. **Manage template versions** with a comprehensive versioning system
3. **Compare template versions** to understand changes between iterations
4. **Share templates** with colleagues and the wider scientific community
5. **Generate protocols from templates** while customizing parameters for specific experiments

This system enhances research reproducibility and efficiency by providing a structured approach to experimental protocol management.

## Architecture

### Data Model

The Protocol Template Management System is built on the following data model:

#### Protocol Template
- **Name**: Descriptive name of the template
- **Description**: Detailed explanation of the protocol's purpose and usage
- **Steps**: Ordered sequence of protocol steps
- **Parameters**: Global parameters applicable to the entire protocol
- **Version**: Version number of the template
- **Category**: Classification of the protocol (e.g., Vitrification, Slow Freezing)
- **Public**: Flag indicating whether the template is publicly available
- **ParentId**: Reference to the template this version was derived from (for versioning)
- **CreatedBy**: User who created the template
- **CreatedAt**: Timestamp of creation
- **UpdatedAt**: Timestamp of last update

#### Protocol Step
- **Name**: Name of the step
- **Description**: Detailed explanation of what happens in this step
- **Parameters**: Step-specific parameters
- **Duration**: Time required for the step
- **DurationUnit**: Unit of time measurement (seconds, minutes, hours, days)
- **Temperature**: Temperature setting for the step
- **TemperatureUnit**: Unit of temperature measurement (C, F, K)

### Components

The system consists of the following components:

1. **Backend API** (Convex integration):
   - CRUD operations for templates and steps
   - Version management
   - Access control for public/private templates

2. **Frontend Components**:
   - Template Management UI for browsing and managing templates
   - Template Editor for creating and editing templates
   - Template Comparison for analyzing differences between versions
   - Step Editor for defining protocol steps

3. **Integration with Experiment System**:
   - Templates can be used to create actual experiments
   - Experiments maintain a reference to their source template
   - Changes to templates don't affect existing experiments

## Key Features

### 1. Template Creation and Editing

Researchers can create new protocol templates using the Template Editor component, which provides:

- Form-based interface for defining template metadata
- Step-by-step definition of protocol procedures
- Parameter configuration for both global and step-specific settings
- Category assignment and visibility control

### 2. Template Versioning

The system implements a comprehensive versioning approach:

- New versions can be created from existing templates
- Version history is maintained and accessible
- Each version receives a unique identifier and version number
- All versions of a template can be viewed in a chronological list

### 3. Version Comparison

Researchers can compare different versions of a template to understand changes:

- Side-by-side comparison of template metadata
- Identification of added, removed, and modified steps
- Parameter change tracking
- Visual highlighting of differences

### 4. Protocol Generation

Templates can be used to generate actual protocols for experiments:

- Researchers select a template as the starting point
- Template parameters can be customized for the specific experiment
- The resulting protocol is a distinct entity, preserving the original template
- Changes to the template don't affect protocols already generated from it

### 5. Access Control

The system implements access controls for templates:

- Private templates are only visible to their creator
- Public templates are available to all researchers
- Template visibility can be changed at any time
- Original creator maintains edit access

## Implementation Details

### Direct Convex Integration

The system implements direct integration with Convex for real-time data synchronization:

- Custom hooks for template management: `useProtocolTemplates`
- Version tracking with `useTemplateVersions`
- Comparison functionality with `useTemplateComparison`

### User Interface Components

The system provides the following UI components:

- **TemplateManagement**: For listing and filtering available templates
- **TemplateEditor**: For creating and editing templates
- **TemplateComparison**: For comparing different versions of a template
- **Step Editor**: For defining individual protocol steps

### Routing Structure

The system is accessible through the following routes:

- `/protocols/templates`: Displays the template management interface
- `/protocols/templates/create`: Provides the template creation interface
- `/protocols/templates/[id]`: Shows the template editing interface for a specific template
- `/protocols/templates/versions/[id]`: Displays version history and comparison for a template

## Usage Guidelines

### Creating a New Template

1. Navigate to the Templates page and click "Create Template"
2. Fill in the basic information (name, description, category)
3. Add steps by clicking "Add Step" and defining each protocol step
4. Set global parameters if needed
5. Choose whether to make the template public
6. Save the template

### Creating a Protocol from a Template

1. Navigate to the Templates page and find the desired template
2. Click "Use" on the template card
3. Customize the protocol name and description
4. Adjust parameters as needed for your specific experiment
5. Click "Create Protocol" to generate a new protocol based on the template

### Managing Template Versions

1. Edit an existing template and make desired changes
2. The system automatically creates a new version when you save changes
3. Access version history through the template details page
4. Compare different versions to understand changes

## Conclusion

The Protocol Template Management System provides a robust foundation for managing scientific protocols in a standardized, versioned manner. By enabling template sharing and reuse, the system enhances research reproducibility and efficiency in cryopreservation research.

Future enhancements will focus on:

1. Enhanced template categorization and tagging
2. Template recommendation based on experiment parameters
3. Statistical analysis of template usage and success rates
4. Template export to standard formats (JSON, YAML, PDF)
5. Integration with external protocol repositories