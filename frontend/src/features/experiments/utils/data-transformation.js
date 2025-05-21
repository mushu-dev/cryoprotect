/**
 * Utility functions for transforming and preparing experiment data for visualization
 */

/**
 * Formats experiment results for chart visualization
 * @param {Object} experiment - The experiment data
 * @returns {Object} Formatted data for charts
 */
export function formatExperimentDataForCharts(experiment) {
  if (!experiment || !experiment.results) {
    return {
      viabilityData: [],
      recoveryData: [],
      functionalityData: []
    };
  }

  // Remove % signs and convert to numbers
  const viability = parseFloat(experiment.results.viability.replace('%', ''));
  const recovery = parseFloat(experiment.results.recovery.replace('%', ''));
  const functionality = parseFloat(experiment.results.functionality.replace('%', ''));

  return {
    viabilityData: [{ name: 'Viability', value: viability }],
    recoveryData: [{ name: 'Recovery', value: recovery }],
    functionalityData: [{ name: 'Functionality', value: functionality }],
    combinedData: [
      { name: 'Viability', value: viability },
      { name: 'Recovery', value: recovery },
      { name: 'Functionality', value: functionality }
    ]
  };
}

/**
 * Prepares data for experiment comparison
 * @param {Array} experiments - Array of experiment data
 * @returns {Object} Formatted data for comparison charts
 */
export function prepareExperimentComparisonData(experiments) {
  if (!experiments || !Array.isArray(experiments) || experiments.length === 0) {
    return {
      comparisonData: [],
      labels: []
    };
  }

  const labels = ['Viability', 'Recovery', 'Functionality'];
  const comparisonData = experiments.map(exp => {
    // Remove % signs and convert to numbers
    const viability = parseFloat(exp.results.viability.replace('%', ''));
    const recovery = parseFloat(exp.results.recovery.replace('%', ''));
    const functionality = parseFloat(exp.results.functionality.replace('%', ''));

    return {
      name: exp.title,
      data: [viability, recovery, functionality]
    };
  });

  return {
    comparisonData,
    labels
  };
}

/**
 * Categorizes experiments by a specific property
 * @param {Array} experiments - Array of experiment data
 * @param {string} property - Property to categorize by (e.g., 'cellType', 'freezingRate')
 * @returns {Object} Categorized experiments
 */
export function categorizeExperiments(experiments, property) {
  if (!experiments || !Array.isArray(experiments) || experiments.length === 0) {
    return {};
  }

  return experiments.reduce((acc, experiment) => {
    const key = experiment[property] || 'Unknown';
    
    if (!acc[key]) {
      acc[key] = [];
    }
    
    acc[key].push(experiment);
    return acc;
  }, {});
}

/**
 * Analyzes cryoprotectant effectiveness across experiments
 * @param {Array} experiments - Array of experiment data
 * @returns {Array} Effectiveness data by cryoprotectant
 */
export function analyzeCryoprotectantEffectiveness(experiments) {
  if (!experiments || !Array.isArray(experiments) || experiments.length === 0) {
    return [];
  }

  const cryoprotectantMap = {};

  // Group experiments by cryoprotectant
  experiments.forEach(exp => {
    if (!exp.cryoprotectants || !Array.isArray(exp.cryoprotectants)) return;

    exp.cryoprotectants.forEach(cp => {
      if (!cryoprotectantMap[cp.name]) {
        cryoprotectantMap[cp.name] = {
          name: cp.name,
          experiments: 0,
          totalViability: 0,
          totalRecovery: 0,
          totalFunctionality: 0
        };
      }

      const viability = parseFloat(exp.results.viability.replace('%', ''));
      const recovery = parseFloat(exp.results.recovery.replace('%', ''));
      const functionality = parseFloat(exp.results.functionality.replace('%', ''));

      cryoprotectantMap[cp.name].experiments += 1;
      cryoprotectantMap[cp.name].totalViability += viability;
      cryoprotectantMap[cp.name].totalRecovery += recovery;
      cryoprotectantMap[cp.name].totalFunctionality += functionality;
    });
  });

  // Calculate averages
  return Object.values(cryoprotectantMap).map(cp => ({
    name: cp.name,
    averageViability: (cp.totalViability / cp.experiments).toFixed(1),
    averageRecovery: (cp.totalRecovery / cp.experiments).toFixed(1),
    averageFunctionality: (cp.totalFunctionality / cp.experiments).toFixed(1),
    experimentCount: cp.experiments
  }));
}

/**
 * Filters experiments based on multiple criteria
 * @param {Array} experiments - Array of experiment data
 * @param {Object} filters - Filter criteria
 * @returns {Array} Filtered experiments
 */
export function filterExperiments(experiments, filters) {
  if (!experiments || !Array.isArray(experiments) || experiments.length === 0) {
    return [];
  }

  return experiments.filter(exp => {
    // Filter by status
    if (filters.status && exp.status !== filters.status) {
      return false;
    }

    // Filter by cell type
    if (filters.cellType && !exp.cellType.toLowerCase().includes(filters.cellType.toLowerCase())) {
      return false;
    }

    // Filter by date range
    if (filters.dateRange) {
      const expDate = new Date(exp.date);
      if (
        (filters.dateRange.start && expDate < new Date(filters.dateRange.start)) ||
        (filters.dateRange.end && expDate > new Date(filters.dateRange.end))
      ) {
        return false;
      }
    }

    // Filter by cryoprotectant
    if (filters.cryoprotectant && 
        (!exp.cryoprotectants || 
         !exp.cryoprotectants.some(cp => cp.name.toLowerCase().includes(filters.cryoprotectant.toLowerCase())))
    ) {
      return false;
    }

    // Filter by performance thresholds
    if (filters.minViability) {
      const viability = parseFloat(exp.results.viability.replace('%', ''));
      if (viability < filters.minViability) {
        return false;
      }
    }

    if (filters.minRecovery) {
      const recovery = parseFloat(exp.results.recovery.replace('%', ''));
      if (recovery < filters.minRecovery) {
        return false;
      }
    }

    // All filters passed
    return true;
  });
}