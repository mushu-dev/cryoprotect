# Protein Visualization Enhancement Implementation Plan

This document outlines the comprehensive implementation plan for enhancing the MolstarViewer component in the CryoProtect application, focusing on adding advanced visualization features and integration with experimental data.

## Current Status

The current MolstarViewer component already has:
- Basic Mol* integration with support for PDB, SMILES, and direct data loading
- Multiple visualization styles (cartoon, ball-and-stick, surface, ribbon, spacefill)
- Basic residue highlighting
- Interactive controls for rotation, zoom, and style changes
- Screenshot capability

## Implementation Goals

### Phase 1: Advanced Mol* Integration

1. **Enhanced Interaction Controls**
   - Add measurement tools for bonds and angles
   - Implement clipping plane controls for cross-sectional views
   - Add support for multiple simultaneous visualization styles
   - Implement distance measurement between atoms/residues

2. **Performance Optimization**
   - Implement level-of-detail rendering for large molecules
   - Add asynchronous loading with progress indicators
   - Optimize memory usage for large structures
   - Implement WebGL feature detection and fallbacks

### Phase 2: Cross-Reference with Experimental Data

1. **Data Overlay System**
   - Create a system to map experimental values to structural elements
   - Implement color gradient visualization based on experimental data
   - Add heat map visualization for property distribution
   - Support for toggling between different experimental datasets

2. **Selection and Analysis Tools**
   - Implement advanced selection tools for specific structural elements
   - Add statistical analysis of selected regions
   - Create exportable reports of structure-property relationships
   - Support for custom annotations on the structure

### Phase 3: Integration with Cryoprotection Data

1. **Cryoprotectant Binding Visualization**
   - Visualize binding sites with cryoprotectants
   - Highlight interaction points between cryoprotectants and biomolecules
   - Implement visual comparison between different cryoprotectants
   - Add visualization of hydrogen bonds and hydrophobic interactions

2. **Property Mapping**
   - Map freezing/thawing performance data to structural elements
   - Visualize physical property changes during freezing process
   - Implement time-series visualization for dynamic processes
   - Add correlation views between structure and cryoprotection efficiency

## Technical Implementation Plan

### 1. Enhanced Measurement Tools

We'll extend the MolstarViewer component to support measurement of distances, angles, and dihedrals:

```tsx
// New measurement functionality
const enableDistanceMeasurement = async () => {
  if (!pluginRef.current) return;
  
  await pluginRef.current.managers.interactivity.setProps({
    enableMeasurements: true,
    measurementTypeOptions: [{ 
      type: 'distance', 
      color: Color(0xFF0000), 
      label: true 
    }]
  });
  
  // Enable distance measurement mode
  await pluginRef.current.managers.measurement.startMeasurement({ type: 'distance' });
};

const enableAngleMeasurement = async () => {
  if (!pluginRef.current) return;
  
  await pluginRef.current.managers.interactivity.setProps({
    enableMeasurements: true,
    measurementTypeOptions: [{ 
      type: 'angle', 
      color: Color(0x00FF00), 
      label: true 
    }]
  });
  
  // Enable angle measurement mode
  await pluginRef.current.managers.measurement.startMeasurement({ type: 'angle' });
};

const disableMeasurements = async () => {
  if (!pluginRef.current) return;
  
  // Disable all measurements
  await pluginRef.current.managers.interactivity.setProps({
    enableMeasurements: false
  });
  
  // Clear existing measurements
  await pluginRef.current.managers.measurement.clearMeasurements();
};
```

### 2. Multiple Visualization Styles

We'll implement the ability to show multiple visualization styles simultaneously:

```tsx
const applyMultipleStyles = async (styles: Record<string, boolean>) => {
  if (!pluginRef.current || !structureObject) return;

  try {
    // Reset any existing representations
    await pluginRef.current.builders.structure.representation.removeAll(structureObject);
    
    // Apply each selected style
    for (const [style, enabled] of Object.entries(styles)) {
      if (enabled && stylePresetMap[style]) {
        await pluginRef.current.builders.structure.representation.addRepresentation(
          structureObject,
          { type: stylePresetMap[style] }
        );
      }
    }
    
    // Update the view
    await pluginRef.current.canvas3d?.commit();
  } catch (error) {
    console.error('Failed to apply multiple styles:', error);
  }
};
```

### 3. Data Overlay System

We'll create a system to map experimental data to structural elements:

```tsx
interface ExperimentalDataPoint {
  residueId: string;
  value: number;
  label?: string;
}

const applyExperimentalDataColors = async (
  data: ExperimentalDataPoint[], 
  minValue: number, 
  maxValue: number,
  colorScheme: 'rainbow' | 'redBlue' | 'custom' = 'rainbow',
  customColors?: [number, number, number][]
) => {
  if (!pluginRef.current || !structureObject || data.length === 0) return;

  try {
    // Create a color scale based on the colorScheme
    const colors = getColorScale(colorScheme, customColors);
    
    // Create a selection for each data point
    for (const point of data) {
      // Normalize the value between 0 and 1
      const normalizedValue = (point.value - minValue) / (maxValue - minValue);
      
      // Get the color for this value
      const color = getColorForValue(normalizedValue, colors);
      
      // Create a selection of the specified residue
      const expression = {
        core: 'residue',
        params: { ids: [point.residueId] }
      };
      
      // Create new representation specifically for this data point
      await pluginRef.current.builders.structure.representation.addRepresentation(
        structureObject,
        {
          type: 'ball-and-stick',
          color: 'uniform',
          colorParams: { value: color },
          size: 'uniform',
          sizeParams: { value: 0.5 }
        },
        { selector: expression }
      );
      
      // Add label if specified
      if (point.label) {
        await pluginRef.current.builders.structure.representation.addRepresentation(
          structureObject,
          {
            type: 'label',
            color: 'uniform',
            colorParams: { value: Color(0x000000) },
            text: point.label
          },
          { selector: expression }
        );
      }
    }

    // Update the view
    await pluginRef.current.canvas3d?.commit();
  } catch (error) {
    console.error('Failed to apply experimental data colors:', error);
  }
};

// Helper to get a color scale
const getColorScale = (
  scheme: 'rainbow' | 'redBlue' | 'custom', 
  customColors?: [number, number, number][]
): Color[] => {
  switch (scheme) {
    case 'rainbow':
      return [
        Color(0xFF0000), // Red
        Color(0xFFFF00), // Yellow
        Color(0x00FF00), // Green
        Color(0x00FFFF), // Cyan
        Color(0x0000FF)  // Blue
      ];
    case 'redBlue':
      return [
        Color(0xFF0000), // Red
        Color(0xFFFFFF), // White
        Color(0x0000FF)  // Blue
      ];
    case 'custom':
      if (customColors && customColors.length > 0) {
        return customColors.map(rgb => Color.fromRgb(rgb[0], rgb[1], rgb[2]));
      }
      return [Color(0x000000), Color(0xFFFFFF)];
  }
};

// Helper to get a color for a normalized value
const getColorForValue = (normalizedValue: number, colors: Color[]): Color => {
  if (colors.length === 1) return colors[0];
  
  const segments = colors.length - 1;
  const segmentIndex = Math.min(Math.floor(normalizedValue * segments), segments - 1);
  const segmentPosition = (normalizedValue * segments) - segmentIndex;
  
  const color1 = colors[segmentIndex];
  const color2 = colors[segmentIndex + 1];
  
  return Color.interpolate(color1, color2, segmentPosition);
};
```

### 4. Cross-Reference with Experimental Data

We'll extend the component to accept experimental data and visualize it:

```tsx
// Extend the props interface
interface MolstarViewerProps {
  // ... existing props
  experimentalData?: ExperimentalDataPoint[];
  experimentalDataRange?: [number, number]; // Min and max values
  experimentalDataColorScheme?: 'rainbow' | 'redBlue' | 'custom';
  experimentalDataCustomColors?: [number, number, number][];
  experimentalDataLabel?: string;
  onResidueSelect?: (residueId: string) => void;
}

// In the useEffect that loads the structure
useEffect(() => {
  // ... existing code
  
  // Apply experimental data if provided
  if (experimentalData && experimentalData.length > 0) {
    const [min, max] = experimentalDataRange || [
      Math.min(...experimentalData.map(d => d.value)),
      Math.max(...experimentalData.map(d => d.value))
    ];
    
    applyExperimentalDataColors(
      experimentalData,
      min,
      max,
      experimentalDataColorScheme,
      experimentalDataCustomColors
    );
    
    // Add a legend if a label is provided
    if (experimentalDataLabel) {
      addColorLegend(experimentalDataLabel, min, max, experimentalDataColorScheme);
    }
  }
}, [/* ... existing dependencies */, experimentalData, experimentalDataRange, 
    experimentalDataColorScheme, experimentalDataCustomColors, experimentalDataLabel]);

// Add handler for residue selection
const setupResidueSelection = () => {
  if (!pluginRef.current || !onResidueSelect) return;
  
  pluginRef.current.behaviors.interaction.click.subscribe(({ current }) => {
    if (current.loci && current.loci.kind === 'element-loci') {
      const residueInfo = pluginRef.current.helpers.structureInfo.getResidueInfo(current.loci);
      if (residueInfo) {
        const residueId = `${residueInfo.chainId}:${residueInfo.residueNumber}`;
        onResidueSelect(residueId);
      }
    }
  });
};
```

### 5. User Interface Enhancements

We'll add new UI controls for the enhanced features:

```tsx
// In the component render function
return (
  <Card className="overflow-hidden">
    <div className="protein-visualizer relative" style={{ height, width }}>
      {/* ... existing loading and error indicators */}
      
      {/* Enhanced Controls */}
      {showControls && (
        <div className="absolute top-2 right-2 z-10 flex gap-1 bg-background/40 backdrop-blur-sm rounded-md p-1">
          {/* ... existing controls */}
          
          {/* New Measurements menu */}
          <Select value={measurementMode} onValueChange={handleMeasurementModeChange}>
            <SelectTrigger className="h-8 w-32">
              <SelectValue placeholder="Measure" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="none">No Measurement</SelectItem>
              <SelectItem value="distance">Distance</SelectItem>
              <SelectItem value="angle">Angle</SelectItem>
              <SelectItem value="dihedral">Dihedral</SelectItem>
            </SelectContent>
          </Select>
          
          {/* Clear measurements button */}
          <Button 
            variant="ghost" 
            size="icon" 
            className="h-8 w-8 p-1 text-muted-foreground"
            onClick={handleClearMeasurements}
            title="Clear measurements"
          >
            <Eraser size={16} />
          </Button>
          
          {/* Clipping plane controls */}
          <Button 
            variant="ghost" 
            size="icon" 
            className={`h-8 w-8 p-1 ${isClipping ? 'text-primary' : 'text-muted-foreground'}`}
            onClick={handleToggleClipping}
            title={isClipping ? "Disable clipping" : "Enable clipping"}
          >
            <Scissors size={16} />
          </Button>
          
          {/* ... other existing controls */}
        </div>
      )}
      
      {/* Experimental Data Legend */}
      {experimentalData && experimentalDataLabel && (
        <div className="absolute bottom-2 right-2 z-10 bg-background/40 backdrop-blur-sm rounded-md p-2">
          <p className="text-xs font-medium mb-1">{experimentalDataLabel}</p>
          <div className="flex items-center gap-1">
            <div className="w-24 h-4 rounded-sm" style={{ 
              background: 'linear-gradient(to right, #ff0000, #ffff00, #00ff00, #00ffff, #0000ff)' 
            }} />
            <div className="flex justify-between w-24 text-[10px]">
              <span>{experimentalDataRange?.[0] ?? 'Min'}</span>
              <span>{experimentalDataRange?.[1] ?? 'Max'}</span>
            </div>
          </div>
        </div>
      )}
      
      {/* Visualization container */}
      <div 
        ref={containerRef}
        className="w-full h-full"
        style={{ backgroundColor }}
        aria-label={`3D visualization of ${name || 'molecule'}`}
      />
    </div>
  </Card>
);
```

## Testing Strategy

We will implement comprehensive tests for the enhanced MolstarViewer:

1. **Unit Tests**
   - Test prop validation and default values
   - Test state management (style changes, rotation, etc.)
   - Test UI rendering for different configurations

2. **Integration Tests**
   - Test loading of different data sources (PDB, SMILES, etc.)
   - Test visualization style changes
   - Test measurement tools functionality
   - Test experimental data integration

3. **Visual Regression Tests**
   - Test rendering of different molecule types
   - Test UI appearance with different configurations
   - Test experimental data visualization

4. **Performance Tests**
   - Test loading and rendering of large molecules
   - Test memory usage during interaction
   - Test WebGL performance with complex visualizations

## Implementation Schedule

### Week 1: Core Enhancements
- Implement measurement tools (distance, angle, dihedral)
- Add clipping plane functionality
- Implement support for multiple visualization styles
- Add enhanced selection capabilities

### Week 2: Experimental Data Integration
- Implement data overlay system
- Create color mapping for experimental values
- Add color legend for data visualization
- Implement residue selection callback

### Week 3: Testing and Refinement
- Create comprehensive test suite
- Fix any issues discovered during testing
- Optimize performance for large structures
- Polish user interface and interactions

## Conclusion

This implementation plan provides a comprehensive roadmap for enhancing the MolstarViewer component with advanced visualization features and experimental data integration. The enhanced component will provide scientists with powerful tools for visualizing and analyzing molecular structures in the context of cryoprotection research.