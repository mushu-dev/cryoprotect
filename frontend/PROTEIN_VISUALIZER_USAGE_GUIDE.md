# Protein Visualizer Component Guide

This guide documents how to use the new Protein Visualizer component in the CryoProtect frontend application.

## Overview

The `MolstarViewer` component provides advanced 3D visualization of molecular structures, proteins, and chemical compounds. It supports multiple data sources, visualization styles, and interactive controls.

## Basic Usage

```tsx
import { MolstarViewer } from '@/components/protein-visualizer/MolstarViewer';

// Basic usage with SMILES string
function MyComponent() {
  return (
    <MolstarViewer 
      smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
      name="Aspirin" 
      height={400} 
    />
  );
}

// Using PDB ID
function ProteinView() {
  return (
    <MolstarViewer 
      pdbId="1cbs" 
      name="Cellular retinoic acid-binding protein type II" 
      style="cartoon" 
    />
  );
}

// Using custom data URL
function CustomMolecule() {
  return (
    <MolstarViewer 
      url="https://example.com/molecules/my-structure.pdb" 
      format="pdb"
      showControls={true}
    />
  );
}
```

## Props Reference

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| `pdbId` | string | undefined | PDB ID to load from RCSB PDB |
| `url` | string | undefined | URL to load molecular data from |
| `smiles` | string | undefined | SMILES string representation of a molecule |
| `data` | string | undefined | Raw molecular data as a string |
| `format` | 'pdb' \| 'sdf' \| 'mol2' \| 'mmcif' \| 'cif' | 'pdb' | Format of the data being loaded |
| `height` | number | 400 | Height of the visualization in pixels |
| `width` | string \| number | '100%' | Width of the visualization |
| `style` | 'cartoon' \| 'ball-and-stick' \| 'surface' \| 'ribbon' | 'cartoon' | Visualization style |
| `highlightResidues` | string[] | [] | Array of residue identifiers to highlight |
| `onStructureLoaded` | () => void | undefined | Callback when structure is loaded |
| `showControls` | boolean | true | Whether to show control buttons |
| `name` | string | undefined | Name to display for the molecule |
| `backgroundColor` | string | 'transparent' | Background color of the visualization |

## Data Sources

The component accepts one of the following data sources (in order of precedence):

1. `pdbId` - A PDB ID that will be loaded from RCSB PDB
2. `url` - A URL pointing to molecular data file
3. `smiles` - A SMILES string representation of a molecule
4. `data` - Raw molecular data as a string

## Visualization Styles

### Cartoon
Ribbon-like representation of protein secondary structure. Good for overall protein structure visualization.

```tsx
<MolstarViewer pdbId="1cbs" style="cartoon" />
```

### Ball-and-stick
Shows atoms as spheres and bonds as sticks. Good for small molecules and detailed views.

```tsx
<MolstarViewer smiles="CC(=O)OC1=CC=CC=C1C(=O)O" style="ball-and-stick" />
```

### Surface
Shows molecular surface. Good for visualizing binding pockets and interfaces.

```tsx
<MolstarViewer pdbId="1cbs" style="surface" />
```

### Ribbon
Simple ribbon representation. Good for visualizing protein backbone.

```tsx
<MolstarViewer pdbId="1cbs" style="ribbon" />
```

## Interactive Controls

The component provides several interactive controls:

- **Style Selector**: Change the visualization style
- **Rotation Toggle**: Start/stop automatic rotation
- **Reset View**: Reset the camera to the default view
- **Zoom**: Zoom in and out
- **Download**: Download the current view as an image
- **Fullscreen**: Toggle fullscreen mode

## Highlighting Specific Residues

You can highlight specific residues by passing an array of residue identifiers:

```tsx
<MolstarViewer 
  pdbId="1cbs" 
  highlightResidues={["A:123", "A:124", "A:125"]} 
/>
```

## Event Callbacks

The component provides a callback for when the structure is loaded:

```tsx
<MolstarViewer 
  pdbId="1cbs" 
  onStructureLoaded={() => {
    console.log("Structure loaded successfully!");
  }} 
/>
```

## Accessibility

The component includes accessibility features:

- ARIA labels for the visualization area
- Keyboard-accessible controls
- Screen-reader compatible status messages

## Performance Considerations

For optimal performance:

- Use appropriate visualization styles based on molecule size
- For large proteins (>10,000 atoms), consider using "cartoon" or "ribbon" styles
- For small molecules, "ball-and-stick" provides better detail
- Set appropriate height for your layout to avoid layout shifts

## Examples

### Basic Protein Visualization

```tsx
<MolstarViewer 
  pdbId="1cbs" 
  name="Cellular retinoic acid-binding protein type II" 
  height={500} 
  style="cartoon" 
/>
```

### Small Molecule with Custom Background

```tsx
<MolstarViewer 
  smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
  name="Aspirin" 
  backgroundColor="#f8f9fa" 
  style="ball-and-stick" 
  height={300} 
/>
```

### Protein Surface with Highlighted Active Site

```tsx
<MolstarViewer 
  pdbId="1hsg" 
  name="HIV-1 Protease" 
  style="surface" 
  highlightResidues={["A:25", "A:26", "A:27", "B:25", "B:26", "B:27"]} 
  height={600} 
/>
```

### Loading Custom Data

```tsx
<MolstarViewer 
  data={molecularData} 
  format="mol2" 
  name="Custom Molecule" 
  style="ball-and-stick" 
/>
```

## Troubleshooting

### Common Issues

1. **Molecule doesn't load**:
   - Check that the data source (PDB ID, URL, or SMILES) is correct
   - Verify the format matches the provided data
   - Check browser console for errors

2. **Performance issues**:
   - Reduce molecule size or complexity
   - Use simpler visualization styles for large molecules
   - Ensure adequate system resources (WebGL support)

3. **Style doesn't change**:
   - Some styles may not be applicable to certain molecule types
   - Try different styles to see which works best

### Error Messages

- "Failed to initialize molecule viewer": Check browser WebGL support
- "Failed to load structure": Check data source validity
- "Invalid format": Ensure format matches the provided data

## Integration with Other Components

The Protein Visualizer works well with other components in the CryoProtect application:

```tsx
import { MolstarViewer } from '@/components/protein-visualizer/MolstarViewer';
import { Card } from '@/components/ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';

function MoleculeDetails({ molecule }) {
  return (
    <Card className="p-4">
      <h2 className="text-xl font-bold mb-4">{molecule.name}</h2>
      
      <Tabs defaultValue="visualization">
        <TabsList>
          <TabsTrigger value="visualization">3D Structure</TabsTrigger>
          <TabsTrigger value="properties">Properties</TabsTrigger>
          <TabsTrigger value="data">Data</TabsTrigger>
        </TabsList>
        
        <TabsContent value="visualization">
          <MolstarViewer 
            smiles={molecule.smiles} 
            name={molecule.name} 
            height={400} 
          />
        </TabsContent>
        
        <TabsContent value="properties">
          {/* Properties content */}
        </TabsContent>
        
        <TabsContent value="data">
          {/* Data content */}
        </TabsContent>
      </Tabs>
    </Card>
  );
}
```

## Future Enhancements

Planned future enhancements for the Protein Visualizer include:

- Multiple structure overlay and comparison
- Distance measurement tools
- Sequence viewer integration
- Animation of molecular dynamics
- Enhanced annotation capabilities
- VR/AR support