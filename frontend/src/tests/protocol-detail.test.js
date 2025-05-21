import React from 'react';
import { render, screen, fireEvent, waitFor } from '@testing-library/react';
import { useRouter } from 'next/router';
import ProtocolDetailPage from '../pages/protocols/[id]';
import { useProtocol, useProtocolVersions } from '../features/protocols/hooks/use-protocols';

// Mock Next.js router
jest.mock('next/router', () => ({
  useRouter: jest.fn()
}));

// Mock protocol hooks
jest.mock('../features/protocols/hooks/use-protocols', () => ({
  useProtocol: jest.fn(),
  useProtocolVersions: jest.fn()
}));

// Mock protocol data
const mockProtocol = {
  id: '1',
  name: 'Standard Cell Freezing',
  version: '1.2.0',
  description: 'Basic protocol for freezing mammalian cells using DMSO as a cryoprotectant with controlled rate cooling.',
  author: { name: 'Dr. John Doe', id: 'user1' },
  created_at: '2023-09-01',
  updated_at: '2023-10-15',
  duration: '2 hours',
  is_template: false,
  cell_types: ['Human embryonic stem cells', 'Primary fibroblasts', 'CHO cells'],
  equipment: [
    'Controlled-rate freezer',
    'Cryovials',
    'Liquid nitrogen storage tank',
    'Water bath (37°C)',
    'Cell culture materials'
  ],
  cryoprotectants: [
    { id: '1', name: 'DMSO', concentration: '10%' },
    { id: '2', name: 'Fetal bovine serum', concentration: '20%' }
  ],
  freezing_rate: '-1°C/min',
  storage_temperature: '-196°C (liquid nitrogen)',
  thawing_method: 'Rapid thawing in 37°C water bath',
  steps: [
    {
      id: 's1',
      order: 1,
      title: 'Prepare cryopreservation medium',
      description: 'Mix cell culture medium with DMSO to achieve 10% final concentration and FBS to 20% final concentration. Keep on ice.',
      duration: '15 minutes',
      is_critical: true
    },
    {
      id: 's2',
      order: 2,
      title: 'Prepare cells',
      description: 'Harvest cells in log phase growth (70-80% confluency). Count cells and check viability.',
      duration: '20 minutes',
      is_critical: true
    }
  ],
  notes: 'For optimal results, cells should be in log phase of growth and have viability >90% before freezing.',
  references: [
    { title: 'Optimization of cryopreservation procedures for mammalian cells', authors: 'Smith et al.', journal: 'Cryobiology', year: 2020 }
  ],
  used_in_experiments: [
    { id: '1', title: 'DMSO concentration optimization', date: '2023-11-10' }
  ]
};

// Mock version history data
const mockVersions = [
  {
    id: '1',
    version: '1.2.0',
    created_at: '2023-09-01',
    created_by: 'Dr. John Doe',
    description: 'Current version'
  },
  {
    id: '3',
    version: '1.1.0',
    created_at: '2023-08-15',
    created_by: 'Dr. John Doe',
    description: 'Added improved cooling rates'
  },
  {
    id: '4',
    version: '1.0.0',
    created_at: '2023-07-01',
    created_by: 'Dr. John Doe',
    description: 'Initial version'
  }
];

describe('ProtocolDetailPage', () => {
  beforeEach(() => {
    // Reset mocks
    jest.clearAllMocks();
    
    // Mock router
    useRouter.mockReturnValue({
      query: { id: '1' },
      isFallback: false,
      push: jest.fn()
    });
    
    // Mock protocol hooks
    useProtocol.mockReturnValue({
      protocol: mockProtocol,
      loading: false
    });
    
    useProtocolVersions.mockReturnValue({
      versions: mockVersions,
      loading: false,
      compareVersions: jest.fn()
    });
  });
  
  test('renders protocol details correctly', () => {
    render(<ProtocolDetailPage />);
    
    // Check protocol title and metadata
    expect(screen.getByText(mockProtocol.name)).toBeInTheDocument();
    expect(screen.getByText(`v${mockProtocol.version}`)).toBeInTheDocument();
    expect(screen.getByText(mockProtocol.description)).toBeInTheDocument();
    
    // Check author info
    expect(screen.getByText(/Dr. John Doe/)).toBeInTheDocument();
    
    // Check protocol tabs
    expect(screen.getByText('Steps')).toBeInTheDocument();
    expect(screen.getByText('Materials & Equipment')).toBeInTheDocument();
    expect(screen.getByText('Parameters')).toBeInTheDocument();
    expect(screen.getByText('Notes & References')).toBeInTheDocument();
    
    // Check protocol steps
    expect(screen.getByText('Prepare cryopreservation medium')).toBeInTheDocument();
    expect(screen.getByText('Critical Step')).toBeInTheDocument();
  });
  
  test('shows version history when toggle is clicked', async () => {
    render(<ProtocolDetailPage />);
    
    // Version history should not be visible initially
    expect(screen.queryByText('Version History')).not.toBeInTheDocument();
    
    // Click the show version history button
    fireEvent.click(screen.getByText(/Show Version History/));
    
    // Version history should be visible
    expect(screen.getByText('Version History')).toBeInTheDocument();
    expect(screen.getByText('v1.1.0')).toBeInTheDocument();
    expect(screen.getByText('v1.0.0')).toBeInTheDocument();
    expect(screen.getByText('Added improved cooling rates')).toBeInTheDocument();
    expect(screen.getByText('Initial version')).toBeInTheDocument();
    
    // "Compare with Current" button should be disabled initially
    const compareButton = screen.getByText('Compare with Current');
    expect(compareButton).toBeDisabled();
  });
  
  test('enables version comparison when a version is selected', async () => {
    render(<ProtocolDetailPage />);
    
    // Show version history
    fireEvent.click(screen.getByText(/Show Version History/));
    
    // Select a version by clicking the radio button (we can't directly access by ID, so we'll use index)
    const radioButtons = screen.getAllByRole('radio');
    fireEvent.click(radioButtons[1]); // Second radio button (version 1.1.0)
    
    // "Compare with Current" button should be enabled
    const compareButton = screen.getByText('Compare with Current');
    expect(compareButton).not.toBeDisabled();
    
    // Click the compare button
    fireEvent.click(compareButton);
    
    // Version comparison component should be rendered
    await waitFor(() => {
      expect(screen.getByText('Protocol Version Comparison')).toBeInTheDocument();
    });
  });
  
  test('switches tabs when clicked', () => {
    render(<ProtocolDetailPage />);
    
    // Should start with Steps tab active
    expect(screen.getByText('Protocol Steps')).toBeInTheDocument();
    
    // Click Materials & Equipment tab
    fireEvent.click(screen.getByText('Materials & Equipment'));
    expect(screen.getByText('Materials and Equipment')).toBeInTheDocument();
    expect(screen.getByText('Cryoprotectants')).toBeInTheDocument();
    expect(screen.getByText('Equipment')).toBeInTheDocument();
    
    // Click Parameters tab
    fireEvent.click(screen.getByText('Parameters'));
    expect(screen.getByText('Compatible Cell Types')).toBeInTheDocument();
    expect(screen.getByText('Freezing Rate')).toBeInTheDocument();
    
    // Click Notes & References tab
    fireEvent.click(screen.getByText('Notes & References'));
    expect(screen.getByText('Notes')).toBeInTheDocument();
    expect(screen.getByText('References')).toBeInTheDocument();
  });
  
  test('shows loading state when protocol is loading', () => {
    // Mock loading state
    useProtocol.mockReturnValue({
      protocol: null,
      loading: true
    });
    
    render(<ProtocolDetailPage />);
    
    // Should show loading indicator
    expect(screen.getByText('Loading protocol details...')).toBeInTheDocument();
  });
  
  test('navigates to correct version when "View" is clicked', () => {
    const mockRouterPush = jest.fn();
    useRouter.mockReturnValue({
      query: { id: '1' },
      isFallback: false,
      push: mockRouterPush
    });
    
    render(<ProtocolDetailPage />);
    
    // Show version history
    fireEvent.click(screen.getByText(/Show Version History/));
    
    // Click "View" link for version 1.1.0
    const viewLinks = screen.getAllByText('View');
    fireEvent.click(viewLinks[0]);
    
    // Router.push should be called with the correct URL
    expect(mockRouterPush).toHaveBeenCalledWith('/protocols/3');
  });
});