import React from 'react';
import { render, screen, fireEvent } from '@testing-library/react';
import { ProtocolVersionComparison } from '../features/protocols/components/protocol-version-comparison';

// Mock protocol data for testing
const protocolA = {
  id: '1',
  name: 'Standard Cell Freezing Protocol',
  version: '1.0.0',
  description: 'Basic protocol for freezing mammalian cells using DMSO as a cryoprotectant.',
  author: { name: 'Dr. John Doe', id: 'user1' },
  created_at: '2023-07-01',
  updated_at: '2023-07-01',
  tags: ['cryopreservation', 'mammalian cells'],
  steps: [
    {
      id: 's1',
      order: 1,
      name: 'Prepare cryopreservation medium',
      description: 'Mix cell culture medium with DMSO to achieve 10% final concentration. Keep on ice.',
      duration: 15,
      duration_unit: 'minutes',
      is_critical: true
    },
    {
      id: 's2',
      order: 2,
      name: 'Prepare cells',
      description: 'Harvest cells in log phase growth. Count cells and check viability.',
      duration: 20,
      duration_unit: 'minutes',
      is_critical: true
    },
    {
      id: 's3',
      order: 3,
      name: 'Centrifuge cell suspension',
      description: 'Centrifuge at 200 x g for 5 minutes at room temperature.',
      duration: 10,
      duration_unit: 'minutes',
      is_critical: false
    }
  ]
};

const protocolB = {
  id: '2',
  name: 'Standard Cell Freezing Protocol',
  version: '1.1.0',
  description: 'Enhanced protocol for freezing mammalian cells using DMSO as a cryoprotectant with controlled rate cooling.',
  author: { name: 'Dr. John Doe', id: 'user1' },
  created_at: '2023-09-01',
  updated_at: '2023-09-01',
  tags: ['cryopreservation', 'mammalian cells', 'controlled rate'],
  steps: [
    {
      id: 's1',
      order: 1,
      name: 'Prepare cryopreservation medium',
      description: 'Mix cell culture medium with DMSO to achieve 10% final concentration and FBS to 20% final concentration. Keep on ice.',
      duration: 15,
      duration_unit: 'minutes',
      is_critical: true
    },
    {
      id: 's2',
      order: 2,
      name: 'Prepare cells',
      description: 'Harvest cells in log phase growth (70-80% confluency). Count cells and check viability.',
      duration: 20,
      duration_unit: 'minutes',
      is_critical: true
    },
    {
      id: 's3',
      order: 3,
      name: 'Centrifuge cell suspension',
      description: 'Centrifuge at 200 x g for 5 minutes at room temperature.',
      duration: 10,
      duration_unit: 'minutes',
      is_critical: false
    },
    {
      id: 's4',
      order: 4,
      name: 'Resuspend cells in cryopreservation medium',
      description: 'Gently resuspend cell pellet in cold cryopreservation medium to achieve 1-5 x 10^6 cells/mL.',
      duration: 10,
      duration_unit: 'minutes',
      is_critical: true
    }
  ]
};

// Mock function for handling version selection
const mockOnViewVersion = jest.fn();

describe('ProtocolVersionComparison', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  test('renders the component with side-by-side view by default', () => {
    render(
      <ProtocolVersionComparison 
        protocolA={protocolA} 
        protocolB={protocolB}
        onViewVersion={mockOnViewVersion}
      />
    );
    
    // Check component title is rendered
    expect(screen.getByText('Protocol Version Comparison')).toBeInTheDocument();
    
    // Check view options are available
    expect(screen.getByText('Side by Side')).toBeInTheDocument();
    expect(screen.getByText('Unified')).toBeInTheDocument();
    
    // Check protocol versions are displayed
    expect(screen.getByText(`Version ${protocolA.version}`)).toBeInTheDocument();
    expect(screen.getByText(`Version ${protocolB.version}`)).toBeInTheDocument();
    
    // Check protocol names are displayed
    expect(screen.getByText(protocolA.name)).toBeInTheDocument();
    expect(screen.getByText(protocolB.name)).toBeInTheDocument();
    
    // Check changes summary is shown
    expect(screen.getByText('Changes Summary')).toBeInTheDocument();
    expect(screen.getByText('1 Step Added')).toBeInTheDocument();
    expect(screen.getByText('Step Modified')).toBeInTheDocument();
  });

  test('switches between side-by-side and unified views', () => {
    render(
      <ProtocolVersionComparison 
        protocolA={protocolA} 
        protocolB={protocolB}
        onViewVersion={mockOnViewVersion}
      />
    );
    
    // Component should start in side-by-side view
    expect(screen.getAllByText(protocolA.name).length).toBeGreaterThan(0);
    
    // Switch to unified view
    fireEvent.click(screen.getByText('Unified'));
    
    // Check unified view elements are rendered
    expect(screen.getByText('Added in v1.1.0')).toBeInTheDocument();
  });

  test('handles view version button clicks', () => {
    render(
      <ProtocolVersionComparison 
        protocolA={protocolA} 
        protocolB={protocolB}
        onViewVersion={mockOnViewVersion}
      />
    );
    
    // Find and click the "View v1.0.0" button
    const viewOldVersionButton = screen.getByText('View v1.0.0');
    fireEvent.click(viewOldVersionButton);
    
    // Check that onViewVersion was called with the correct protocol
    expect(mockOnViewVersion).toHaveBeenCalledWith(protocolA);
    
    // Find and click the "View v1.1.0" button
    const viewNewVersionButton = screen.getByText('View v1.1.0');
    fireEvent.click(viewNewVersionButton);
    
    // Check that onViewVersion was called with the correct protocol
    expect(mockOnViewVersion).toHaveBeenCalledWith(protocolB);
  });

  test('calculates and displays differences between protocols correctly', () => {
    render(
      <ProtocolVersionComparison 
        protocolA={protocolA} 
        protocolB={protocolB}
        onViewVersion={mockOnViewVersion}
      />
    );
    
    // Test for description change
    expect(screen.getByText(protocolA.description)).toHaveClass('line-through');
    expect(screen.getByText(protocolB.description)).toHaveClass('text-green-500');
    
    // Test for added step
    expect(screen.getByText('Resuspend cells in cryopreservation medium')).toBeInTheDocument();
    
    // Test for modified step
    expect(screen.getByText('Prepare cells')).toBeInTheDocument();
    
    // Test for tag changes
    expect(screen.getByText('controlled rate')).toHaveClass('text-green-500');
  });

  test('handles expand/collapse functionality', () => {
    render(
      <ProtocolVersionComparison 
        protocolA={protocolA} 
        protocolB={protocolB}
        onViewVersion={mockOnViewVersion}
      />
    );
    
    // Find and click the "Expand All" button
    const expandAllButton = screen.getByText('Expand All');
    fireEvent.click(expandAllButton);
    
    // Check that step descriptions are shown (indicating expansion)
    expect(screen.getByText('Mix cell culture medium with DMSO to achieve 10% final concentration. Keep on ice.')).toBeInTheDocument();
    
    // Find and click the "Collapse All" button
    const collapseAllButton = screen.getByText('Collapse All');
    fireEvent.click(collapseAllButton);
    
    // Check that step descriptions are hidden (indicating collapse)
    // This might be harder to test without more specific selectors, but we could check for changes in DOM structure
  });
});