/**
 * THIS IS A TEMPORARY PLACEHOLDER FILE
 * 
 * This file is a temporary placeholder for the Convex API types.
 * Normally, these types would be generated using the Convex CLI:
 * 
 * npx convex codegen
 * 
 * For now, we're providing a minimal implementation to prevent runtime errors.
 */

// Define a minimal API structure for experiments
export const api = {
  experiments: {
    enhanced_experiments: {
      listEnhancedExperiments: 'experiments/enhanced_experiments:listEnhancedExperiments',
      getEnhancedExperiment: 'experiments/enhanced_experiments:getEnhancedExperiment',
      createEnhancedExperiment: 'experiments/enhanced_experiments:createEnhancedExperiment',
      updateEnhancedExperiment: 'experiments/enhanced_experiments:updateEnhancedExperiment',
      deleteEnhancedExperiment: 'experiments/enhanced_experiments:deleteEnhancedExperiment',
      updateEnhancedExperimentStatus: 'experiments/enhanced_experiments:updateEnhancedExperimentStatus'
    },
    update: 'experiments:update'
  },
  molecules: {
    getAllMolecules: 'molecules:getAllMolecules',
    getMolecule: 'molecules:getMolecule',
    createMolecule: 'molecules:createMolecule',
    updateMolecule: 'molecules:updateMolecule',
    deleteMolecule: 'molecules:deleteMolecule'
  }
};