/**
 * Helper functions for experiments and experiment results
 */

import { Id } from "../_generated/dataModel";
import { DataModel } from "../_generated/dataModel";
import { QueryCtx } from "../_generated/server";
import { 
  Experiment, 
  ExperimentResult, 
  ExperimentWithResults,
  ExperimentWithMixture,
  ExperimentComplete
} from "./types";

/**
 * Expand an experiment with its results
 */
export async function expandExperimentWithResults(
  ctx: QueryCtx,
  experiment: Experiment
): Promise<ExperimentWithResults> {
  const results = await ctx.db
    .query("experimentResults")
    .withIndex("by_experiment", q => q.eq("experimentId", experiment._id))
    .collect();
  
  return {
    ...experiment,
    results
  };
}

/**
 * Expand an experiment with its mixture details
 */
export async function expandExperimentWithMixture(
  ctx: QueryCtx,
  experiment: Experiment
): Promise<ExperimentWithMixture> {
  if (!experiment.mixtureId) {
    return {
      ...experiment,
      mixture: undefined
    };
  }
  
  const mixture = await ctx.db.get(experiment.mixtureId);
  
  if (!mixture) {
    return {
      ...experiment,
      mixture: undefined
    };
  }
  
  // Count the number of components in the mixture
  const componentCount = await ctx.db
    .query("mixtureComponents")
    .withIndex("by_mixture", q => q.eq("mixtureId", mixture._id))
    .count();
  
  return {
    ...experiment,
    mixture: {
      _id: mixture._id,
      name: mixture.name,
      componentCount
    }
  };
}

/**
 * Expand an experiment with complete details (results and mixture with components)
 */
export async function expandExperimentComplete(
  ctx: QueryCtx,
  experiment: Experiment
): Promise<ExperimentComplete> {
  // Get experiment results
  const results = await ctx.db
    .query("experimentResults")
    .withIndex("by_experiment", q => q.eq("experimentId", experiment._id))
    .collect();
  
  // Prepare complete experiment data
  const complete: ExperimentComplete = {
    ...experiment,
    results,
    mixture: undefined
  };
  
  // Add mixture details if available
  if (experiment.mixtureId) {
    const mixture = await ctx.db.get(experiment.mixtureId);
    
    if (mixture) {
      // Get mixture components with molecule names
      const components = await ctx.db
        .query("mixtureComponents")
        .withIndex("by_mixture", q => q.eq("mixtureId", mixture._id))
        .collect();
      
      const componentsWithNames = [];
      
      for (const component of components) {
        const molecule = await ctx.db.get(component.moleculeId);
        
        if (molecule) {
          componentsWithNames.push({
            moleculeId: molecule._id,
            moleculeName: molecule.name,
            concentration: component.concentration,
            units: component.units
          });
        }
      }
      
      complete.mixture = {
        _id: mixture._id,
        name: mixture.name,
        components: componentsWithNames
      };
    }
  }
  
  return complete;
}

/**
 * Format experiment status for display
 */
export function formatExperimentStatus(status: string): string {
  switch (status) {
    case "planned":
      return "Planned";
    case "in-progress":
      return "In Progress";
    case "completed":
      return "Completed";
    case "failed":
      return "Failed";
    default:
      return status;
  }
}

/**
 * Get statistics for experiment results
 * Returns min, max, avg, and count for numeric results
 */
export function getResultStatistics(
  results: ExperimentResult[], 
  parameterName?: string
): { min: number, max: number, avg: number, count: number } | null {
  // Filter results by parameter name if provided
  const filteredResults = parameterName
    ? results.filter(r => r.parameterName === parameterName)
    : results;
  
  // Filter for numeric results only
  const numericResults = filteredResults
    .filter(r => r.numericValue !== undefined)
    .map(r => r.numericValue as number);
  
  if (numericResults.length === 0) {
    return null;
  }
  
  // Calculate statistics
  const min = Math.min(...numericResults);
  const max = Math.max(...numericResults);
  const sum = numericResults.reduce((a, b) => a + b, 0);
  const avg = sum / numericResults.length;
  
  return {
    min,
    max,
    avg,
    count: numericResults.length
  };
}

/**
 * Check if an experiment has results
 */
export async function experimentHasResults(
  ctx: QueryCtx,
  experimentId: Id<"experiments">
): Promise<boolean> {
  const count = await ctx.db
    .query("experimentResults")
    .withIndex("by_experiment", q => q.eq("experimentId", experimentId))
    .count();
  
  return count > 0;
}

/**
 * Group experiment results by parameter name
 */
export function groupResultsByParameter(
  results: ExperimentResult[]
): Record<string, ExperimentResult[]> {
  const grouped: Record<string, ExperimentResult[]> = {};
  
  for (const result of results) {
    if (!grouped[result.parameterName]) {
      grouped[result.parameterName] = [];
    }
    
    grouped[result.parameterName].push(result);
  }
  
  return grouped;
}