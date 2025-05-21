import { useState, useEffect, useCallback, useMemo } from 'react';
import { filterExperiments } from '../utils/data-transformation';

// Mock data for experiments
const MOCK_EXPERIMENTS = [
  {
    id: '1',
    title: 'DMSO Concentration Optimization',
    status: 'Completed',
    description: 'Evaluating different DMSO concentrations for optimal cryopreservation of HeLa cells.',
    date: 'May 12, 2025',
    createdBy: 'Dr. Jane Smith',
    protocol: {
      id: '1',
      name: 'Standard Cell Freezing Protocol',
      version: '2.1.0'
    },
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '10%' }
    ],
    cellType: 'Human HeLa cells',
    freezingRate: '-1°C/min controlled rate',
    storageTemperature: '-196°C (liquid nitrogen)',
    thawingMethod: 'Rapid thawing in 37°C water bath',
    results: {
      viability: '87%',
      recovery: '92%',
      functionality: '85%',
      notes: 'Cells exhibited normal morphology and growth characteristics following thawing.'
    }
  },
  {
    id: '2',
    title: 'Trehalose DMSO Combination',
    status: 'Completed',
    description: 'Testing combined effect of trehalose and DMSO on fibroblast preservation.',
    date: 'May 15, 2025',
    createdBy: 'Dr. Michael Chen',
    protocol: {
      id: '1',
      name: 'Standard Cell Freezing Protocol',
      version: '2.1.0'
    },
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '5%' },
      { id: '3', name: 'Trehalose', concentration: '100mM' }
    ],
    cellType: 'Human dermal fibroblasts',
    freezingRate: '-1°C/min controlled rate',
    storageTemperature: '-196°C (liquid nitrogen)',
    thawingMethod: 'Rapid thawing in 37°C water bath',
    results: {
      viability: '92%',
      recovery: '95%',
      functionality: '90%',
      notes: 'Combination approach showed significant improvement over DMSO alone.'
    }
  },
  {
    id: '3',
    title: 'Glycerol Cryopreservation',
    status: 'Completed',
    description: 'Evaluation of glycerol as a primary cryoprotectant for plant cells.',
    date: 'May 18, 2025',
    createdBy: 'Dr. Sarah Johnson',
    protocol: {
      id: '3',
      name: 'Plant Tissue Preservation',
      version: '1.5.0'
    },
    cryoprotectants: [
      { id: '2', name: 'Glycerol', concentration: '15%' }
    ],
    cellType: 'Arabidopsis leaf cells',
    freezingRate: '-0.5°C/min controlled rate',
    storageTemperature: '-80°C (ultra-low freezer)',
    thawingMethod: 'Controlled rate thawing',
    results: {
      viability: '72%',
      recovery: '68%',
      functionality: '75%',
      notes: 'Plant cells showed moderate recovery, but maintained good genetic stability.'
    }
  },
  {
    id: '4',
    title: 'Low Temperature DMSO Protocol',
    status: 'In Progress',
    description: 'Testing the effect of reduced temperature during DMSO exposure.',
    date: 'May 20, 2025',
    createdBy: 'Dr. Jane Smith',
    protocol: {
      id: '1',
      name: 'Standard Cell Freezing Protocol',
      version: '2.1.0'
    },
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '10%' }
    ],
    cellType: 'Human T-cells',
    freezingRate: '-1°C/min controlled rate',
    storageTemperature: '-196°C (liquid nitrogen)',
    thawingMethod: 'Rapid thawing in 37°C water bath',
    results: {
      viability: '85%',
      recovery: '82%',
      functionality: '80%',
      notes: 'Preliminary results show promising viability with reduced toxicity.'
    }
  },
  {
    id: '5',
    title: 'Ethylene Glycol Protocol',
    status: 'Completed',
    description: 'Evaluation of ethylene glycol for embryo cryopreservation.',
    date: 'May 22, 2025',
    createdBy: 'Dr. Michael Chen',
    protocol: {
      id: '2',
      name: 'Vitrification Protocol',
      version: '2.0.1'
    },
    cryoprotectants: [
      { id: '5', name: 'Ethylene Glycol', concentration: '30%' }
    ],
    cellType: 'Mouse embryos',
    freezingRate: 'Vitrification (ultra-rapid)',
    storageTemperature: '-196°C (liquid nitrogen)',
    thawingMethod: 'Rapid thawing in 37°C water bath',
    results: {
      viability: '96%',
      recovery: '93%',
      functionality: '94%',
      notes: 'Exceptional results with minimal cryoinjury observed in post-thaw analysis.'
    }
  },
  {
    id: '6',
    title: 'Triple Cryoprotectant Mixture',
    status: 'In Progress',
    description: 'Novel approach using a three-cryoprotectant cocktail for improved cell survival.',
    date: 'May 25, 2025',
    createdBy: 'Dr. Sarah Johnson',
    protocol: {
      id: '2',
      name: 'Vitrification Protocol',
      version: '2.0.1'
    },
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '5%' },
      { id: '4', name: 'Propylene Glycol', concentration: '5%' },
      { id: '5', name: 'Ethylene Glycol', concentration: '5%' }
    ],
    cellType: 'Human iPSCs',
    freezingRate: 'Vitrification (ultra-rapid)',
    storageTemperature: '-196°C (liquid nitrogen)',
    thawingMethod: 'Step-wise thawing protocol',
    results: {
      viability: '91%',
      recovery: '88%',
      functionality: '93%',
      notes: 'Combined approach shows potential for reducing individual cryoprotectant toxicity while maintaining effectiveness.'
    }
  }
];

/**
 * Custom hook for fetching and managing experiment data
 */
export default function useExperimentData() {
  const [experiments, setExperiments] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [filters, setFilters] = useState({});
  const [selectedExperiments, setSelectedExperiments] = useState([]);

  // Fetch experiments (mock implementation)
  const fetchExperiments = useCallback(async () => {
    try {
      setLoading(true);
      // Simulate API call
      await new Promise(resolve => setTimeout(resolve, 800));
      setExperiments(MOCK_EXPERIMENTS);
      setError(null);
    } catch (err) {
      setError('Failed to fetch experiments');
      console.error('Error fetching experiments:', err);
    } finally {
      setLoading(false);
    }
  }, []);

  // Fetch a single experiment by ID
  const fetchExperiment = useCallback(async (id) => {
    try {
      setLoading(true);
      // Simulate API call
      await new Promise(resolve => setTimeout(resolve, 400));
      const experiment = MOCK_EXPERIMENTS.find(exp => exp.id === id);
      
      if (!experiment) {
        throw new Error('Experiment not found');
      }
      
      return experiment;
    } catch (err) {
      setError('Failed to fetch experiment');
      console.error('Error fetching experiment:', err);
      return null;
    } finally {
      setLoading(false);
    }
  }, []);

  // Apply filters to experiments
  const filteredExperiments = useMemo(() => {
    return filterExperiments(experiments, filters);
  }, [experiments, filters]);

  // Toggle experiment selection for comparison
  const toggleExperimentSelection = useCallback((experimentId) => {
    setSelectedExperiments(prev => {
      if (prev.includes(experimentId)) {
        return prev.filter(id => id !== experimentId);
      } else {
        // Limit to comparing max 4 experiments
        if (prev.length >= 4) {
          return [...prev.slice(1), experimentId];
        }
        return [...prev, experimentId];
      }
    });
  }, []);

  // Get selected experiment data
  const selectedExperimentData = useMemo(() => {
    return selectedExperiments.map(id => 
      experiments.find(exp => exp.id === id)
    ).filter(Boolean);
  }, [experiments, selectedExperiments]);

  // Load experiments on initial mount
  useEffect(() => {
    fetchExperiments();
  }, [fetchExperiments]);

  return {
    experiments,
    filteredExperiments,
    loading,
    error,
    filters,
    setFilters,
    fetchExperiment,
    selectedExperiments,
    toggleExperimentSelection,
    selectedExperimentData
  };
}