# Protein Visualizer Component

A powerful 3D molecular visualization component built on top of Mol* (MolStar) for the CryoProtect platform.

## Features

- Multiple data sources support (PDB ID, SMILES, URL, raw data)
- Rich visualization styles (cartoon, ball-and-stick, surface, ribbon, spacefill)
- Advanced measurement tools for distances, angles, and dihedrals
- Support for multiple simultaneous visualization styles
- Clipping plane for cross-sectional views
- Experimental data visualization with color mapping
- Interactive entity selection for data exploration
- Color gradient legends for experimental data
- Interactive controls for manipulation
- Residue highlighting
- Screen capture functionality 
- Responsive design
- Accessibility support

## Usage

```jsx
import { MolstarViewer } from '@/components/protein-visualizer/MolstarViewer';

// Basic usage with PDB ID
<MolstarViewer 
  pdbId="1cbs"
  name="Cellular Retinoic Acid-Binding Protein Type II"
/>

// Using SMILES notation for small molecules
<MolstarViewer 
  smiles="CC(=O)OC1=CC=CC=C1C(=O)O"
  name="Aspirin"
  style="ball-and-stick"
/>

// Advanced usage with highlighted residues
<MolstarViewer 
  pdbId="1hsg"
  name="HIV-1 Protease"
  style="cartoon"
  highlightResidues={["A:25", "A:26", "A:27"]}
  backgroundColor="#f0f9ff"
  spin={true}
/>

// Using measurement tools
<MolstarViewer 
  pdbId="1cbs"
  name="Protein with Measurement Tools"
  showMeasurementTools={true}
  onMeasurement={(type, value, label) => {
    console.log(`Measured ${type}: ${value} (${label})`);
  }}
/>

// Multiple visualization styles
<MolstarViewer 
  pdbId="1cbs"
  name="Multiple Visualization Styles"
  allowMultipleStyles={true}
  style="cartoon"
/>

// Experimental data visualization
<MolstarViewer 
  pdbId="1hsg"
  name="HIV-1 Protease with Temperature Factors"
  experimentalData={[
    { entityId: "A:25", value: 37.5 },
    { entityId: "A:26", value: 42.1 },
    { entityId: "A:27", value: 28.3 },
    // ... more data points
  ]}
  experimentalDataLabel="Temperature Factors (Å²)"
  experimentalDataColorScheme="redWhiteBlue"
  onEntitySelect={(entityId, value) => {
    console.log(`Selected entity ${entityId} with value: ${value}`);
  }}
/>
```

## Props

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| `pdbId` | string | undefined | PDB ID to load from RCSB PDB |
| `url` | string | undefined | URL to load molecular data from |
| `smiles` | string | undefined | SMILES string representation of a molecule |
| `data` | string | undefined | Raw molecular data as a string |
| `format` | 'pdb' \| 'sdf' \| 'mol2' \| 'mmcif' \| 'cif' | 'pdb' | Format of the data |
| `height` | number | 400 | Height of the visualization in pixels |
| `width` | string \| number | '100%' | Width of the visualization |
| `style` | 'cartoon' \| 'ball-and-stick' \| 'surface' \| 'ribbon' \| 'spacefill' | 'cartoon' | Visualization style |
| `highlightResidues` | string[] | [] | Array of residue identifiers to highlight |
| `onStructureLoaded` | () => void | undefined | Callback when structure is loaded |
| `showControls` | boolean | true | Whether to show control buttons |
| `name` | string | undefined | Name to display for the molecule |
| `backgroundColor` | string | 'transparent' | Background color of the visualization |
| `spin` | boolean | false | Whether to auto-rotate the visualization |
| `showMeasurementTools` | boolean | false | Whether to show measurement tools in controls |
| `onMeasurement` | (type: string, value: number, label: string) => void | undefined | Callback when a measurement is made |
| `allowMultipleStyles` | boolean | false | Allow multiple visualization styles simultaneously |
| `experimentalData` | ExperimentalDataPoint[] | [] | Array of experimental data points to visualize |
| `experimentalDataRange` | [number, number] | auto | Min and max values for color scaling |
| `experimentalDataColorScheme` | 'rainbow' \| 'redBlue' \| 'redWhiteBlue' \| 'blueWhiteRed' \| 'custom' | 'rainbow' | Color scheme for experimental data |
| `experimentalDataCustomColors` | [number, number, number][] | undefined | Custom RGB colors for 'custom' color scheme |
| `experimentalDataLabel` | string | undefined | Label for experimental data legend |
| `onEntitySelect` | (entityId: string, value?: number) => void | undefined | Callback when user selects a residue or atom |

## Implementation Details

The MolstarViewer component is built on top of the Mol* (MolStar) library, which provides a comprehensive solution for molecular visualization. The component handles:

- Initialization of the Mol* plugin
- Loading molecular data from various sources
- Applying visualization styles and customizations (including multiple styles)
- Measurement tools for distances, angles, and dihedrals
- Clipping plane functionality for cross-sectional views
- Experimental data visualization with custom color themes
- Interactive entity selection with data retrieval
- Color gradient legends with appropriate scales
- Managing the lifecycle of the visualization
- Providing an interactive UI for manipulation and analysis

### ExperimentalDataPoint Interface

The component accepts experimental data through the `experimentalData` prop, which should be an array of `ExperimentalDataPoint` objects:

```typescript
interface ExperimentalDataPoint {
  // Identifier for the residue, atom, or group (e.g., "A:42" for residue 42 in chain A)
  entityId: string;
  // The data value associated with this entity
  value: number;
  // Optional label to display
  label?: string;
}
```

### Color Schemes

The component supports various color schemes for experimental data visualization:

- `rainbow`: Red → Yellow → Green → Cyan → Blue
- `redBlue`: Red → Blue
- `redWhiteBlue`: Red → White → Blue
- `blueWhiteRed`: Blue → White → Red
- `custom`: Custom colors provided via `experimentalDataCustomColors`

## Demo

A comprehensive demo page is available at `/protein-visualizer-demo` which showcases the component's features and provides usage examples.

## Testing

The component includes Playwright tests to ensure proper rendering and functionality:

```bash
# Run protein visualizer tests
npm run test:protein-visualizer
```

## Dependencies

- molstar: ^4.15.0
- react: 17.0.2
- lucide-react: For UI icons
- @radix-ui/react-select: For style selection dropdown

## Accessibility

The component includes accessibility features:
- Proper ARIA labels
- Keyboard navigation support
- High contrast controls
- Screen reader support for visualization status

## Note on Performance

For optimal performance:
- Prefer smaller structures when possible
- Use appropriate styles (surface rendering is more resource-intensive)
- Consider disabling auto-rotation for lower-end devices