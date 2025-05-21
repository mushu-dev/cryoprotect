# Molecular Visualization Components

This document provides detailed information about the molecular visualization components available in the CryoProtect frontend application.

## Overview

CryoProtect offers two different molecular visualization components:

1. **MoleculeViewer3D** - A lightweight viewer based on 3DMol.js for basic molecule rendering
2. **MoleculeViewerR3F** - An advanced interactive viewer built with React Three Fiber

Both components are available throughout the application and can be used interchangeably.

## MoleculeViewer3D

This component utilizes the 3DMol.js library to render molecules in 3D. It's lightweight and works well for basic visualization needs.

### Props

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| smiles | string | required | SMILES notation of the molecule to render |
| height | number | 300 | Height of the viewer in pixels |
| name | string | undefined | Name of the molecule to display |
| style | 'stick' \| 'sphere' \| 'cartoon' \| 'surface' | 'stick' | Rendering style for the molecule |
| spin | boolean | false | Whether the molecule should auto-rotate |
| backgroundColor | string | 'transparent' | Background color for the viewer |

### Usage

```tsx
import { MoleculeViewer3D } from '@/features/molecules/components/molecule-viewer-3d'

// Basic usage
<MoleculeViewer3D 
  smiles="C(C(CO)O)O" 
  name="Glycerol" 
/>

// Advanced usage
<MoleculeViewer3D 
  smiles="C(C(CO)O)O" 
  name="Glycerol" 
  height={400} 
  style="sphere" 
  spin={true} 
  backgroundColor="#f5f5f5" 
/>
```

## MoleculeViewerR3F

This component uses React Three Fiber and Three.js to provide a more interactive and customizable 3D visualization experience. It includes features like interactive controls, multiple rendering styles, and screenshot capabilities.

### Props

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| smiles | string | required | SMILES notation of the molecule to render |
| name | string | undefined | Name of the molecule to display |
| height | number | 300 | Height of the viewer in pixels |
| autoRotate | boolean | true | Whether the molecule should auto-rotate |
| backgroundColor | string | 'transparent' | Background color for the viewer |

### Features

- **Interactive Controls**: Pan, zoom, and rotate the molecule
- **Download**: Save the current view as a PNG image
- **Rotation Toggle**: Turn auto-rotation on/off
- **Reset View**: Return to the default camera position
- **Element Colors**: Proper coloring for different atomic elements

### Usage

```tsx
import { MoleculeViewerR3F } from '@/features/molecules/components/molecule-viewer-r3f'

// Basic usage
<MoleculeViewerR3F 
  smiles="C(C(CO)O)O" 
  name="Glycerol" 
/>

// Advanced usage
<MoleculeViewerR3F 
  smiles="C(C(CO)O)O" 
  name="Glycerol" 
  height={400} 
  autoRotate={true} 
  backgroundColor="#f5f5f5" 
/>
```

## Implementation Details

### SMILES Parsing

Both viewers use the 3DMol.js library to parse SMILES notation and generate 3D coordinates for molecules. The MoleculeViewerR3F component extracts the atom and bond information from 3DMol.js and renders it using React Three Fiber.

### React Three Fiber Component

The R3F viewer is implemented as a client component that renders a Three.js scene containing:

- Atoms as spheres with element-appropriate colors and radii
- Bonds as line segments between atoms, with proper coloring and bond order representation
- Interactive orbit controls
- Lighting for proper 3D appearance
- UI controls for user interaction

### Molecule Card Integration

The `MoleculeCard` component integrates both viewers as tabs, allowing users to switch between visualization styles without navigating away from the card.

## Performance Considerations

- The 3DMol.js viewer is more lightweight and may perform better on lower-end devices
- The React Three Fiber viewer provides better interactivity but requires more computational resources
- For large molecules (>100 atoms), consider defaulting to the 3DMol.js viewer
- Both viewers implement optimization techniques to ensure smooth rendering

## Future Enhancements

Planned enhancements for the molecular visualization components include:

1. Surface rendering in the R3F viewer
2. Electron density visualization
3. Measurement tools for bond lengths and angles
4. Support for PDB format and protein visualization
5. Animation of molecular dynamics simulations
6. VR/AR support for immersive visualization

## Troubleshooting

If you encounter issues with the molecular visualization:

1. Ensure the SMILES notation is valid
2. Check browser compatibility (WebGL support is required)
3. For large molecules, increase the height parameter to provide more visual space
4. Use the 3DMol.js viewer for better performance on mobile devices

## Browser Support

These visualization components require:

- WebGL support
- Modern browser (Chrome, Firefox, Safari, Edge)
- JavaScript enabled

Mobile support is provided but may have performance limitations on older devices.