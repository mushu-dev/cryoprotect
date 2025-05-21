# CryoProtect Frontend

This is the frontend application for CryoProtect, a platform for cryoprotectant analysis and experiment management.

## Quick Start

```bash
# Install dependencies
npm install

# Run development server
npm run dev

# Build for production
npm run build
```

## UI System

CryoProtect uses a modern UI component system based on Tailwind CSS with a simplified implementation of shadcn/ui components.

### Component Library

To see all available UI components and their usage, start the development server and visit `/components`.

### Dashboard Layout

The application uses a dashboard layout for authenticated pages, which includes:

- Sidebar navigation
- Responsive design for mobile and desktop
- Consistent styling throughout the application

## Testing

```bash
# Run UI tests with Playwright
./run-ui-test.sh
```

## Development

For detailed information on the UI system and how to extend it, see the [UI Setup Guide](./UI_SETUP_GUIDE.md).

## Protein Visualizer Component

The CryoProtect application includes a powerful molecular visualization component based on Mol* (MolStar).

### Features

- Visualize proteins from PDB IDs, SMILES strings, or direct data
- Multiple visualization styles (cartoon, ball-and-stick, surface, etc.)
- Interactive controls for rotation, zoom, and style changes
- Ability to highlight specific residues
- Customizable appearance with background color options
- Screen capture functionality

### Usage

```jsx
import { MolstarViewer } from '@/components/protein-visualizer/MolstarViewer';

// Visualize by PDB ID
<MolstarViewer 
  pdbId="1cbs" 
  name="Cellular Retinoic Acid-Binding Protein Type II" 
  style="cartoon" 
/>

// Visualize by SMILES string
<MolstarViewer 
  smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
  name="Aspirin" 
  style="ball-and-stick" 
/>
```

### Demo Page

A comprehensive demo page is available at `/protein-visualizer-demo` which showcases all the features and provides usage examples.

### Testing

```bash
# Run protein visualizer tests
npm run test:protein-visualizer
```