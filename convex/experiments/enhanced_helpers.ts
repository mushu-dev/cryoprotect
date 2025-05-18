/**
 * Helper functions for enhanced experiments
 */

import { DatabaseReader } from "../_generated/server";
import { 
  EnhancedExperiment,
  EnhancedExperimentResult,
  Protocol,
  TissueType,
  Equipment,
  ExperimentEquipment,
  TimeSeries,
  TimeSeriesDataPoint,
  EnhancedExperimentWithDetails
} from "./enhanced_types";
import { Id } from "../_generated/dataModel";

/**
 * Expands an enhanced experiment with its results
 */
export async function expandEnhancedExperimentWithResults(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperiment & { results: EnhancedExperimentResult[] }> {
  const results = await ctx.db
    .query("enhancedExperimentResults")
    .withIndex("by_experiment", q => q.eq("experimentId", experiment._id))
    .collect();
  
  return {
    ...experiment,
    results
  };
}

/**
 * Expands an enhanced experiment with its protocol
 */
export async function expandEnhancedExperimentWithProtocol(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperiment & { protocol?: Protocol }> {
  if (!experiment.protocolId) {
    return {
      ...experiment,
      protocol: undefined
    };
  }
  
  const protocol = await ctx.db.get(experiment.protocolId);
  
  return {
    ...experiment,
    protocol: protocol || undefined
  };
}

/**
 * Expands an enhanced experiment with its mixture
 */
export async function expandEnhancedExperimentWithMixture(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperiment & { 
  mixture?: { 
    _id: Id<"mixtures">, 
    name: string, 
    components?: Array<{
      moleculeId: Id<"molecules">;
      moleculeName: string;
      concentration: number;
      units: string;
    }> 
  } 
}> {
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
  
  // Get mixture components
  const components = await ctx.db
    .query("mixtureComponents")
    .withIndex("by_mixture", q => q.eq("mixtureId", experiment.mixtureId!))
    .collect();
  
  // Get molecule names
  const componentDetails = [];
  for (const component of components) {
    const molecule = await ctx.db.get(component.moleculeId);
    if (molecule) {
      componentDetails.push({
        moleculeId: component.moleculeId,
        moleculeName: molecule.name,
        concentration: component.concentration,
        units: component.units
      });
    }
  }
  
  return {
    ...experiment,
    mixture: {
      _id: mixture._id,
      name: mixture.name,
      components: componentDetails
    }
  };
}

/**
 * Expands an enhanced experiment with its tissue types
 */
export async function expandEnhancedExperimentWithTissueTypes(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperiment & { tissueTypes?: TissueType[] }> {
  // Get experiment results to find tissue types
  const results = await ctx.db
    .query("enhancedExperimentResults")
    .withIndex("by_experiment", q => q.eq("experimentId", experiment._id))
    .collect();
  
  if (results.length === 0) {
    return {
      ...experiment,
      tissueTypes: []
    };
  }
  
  // Get unique tissue type IDs
  const tissueTypeIds = [...new Set(results.map(result => result.tissueTypeId))];
  
  // Get tissue type details
  const tissueTypes: TissueType[] = [];
  for (const id of tissueTypeIds) {
    const tissueType = await ctx.db.get(id);
    if (tissueType) {
      tissueTypes.push(tissueType);
    }
  }
  
  return {
    ...experiment,
    tissueTypes
  };
}

/**
 * Expands an enhanced experiment with its equipment
 */
export async function expandEnhancedExperimentWithEquipment(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperiment & { 
  equipment?: Array<{
    _id: Id<"equipment">;
    name: string;
    type: string;
    role?: string;
  }> 
}> {
  // Get experiment equipment links
  const equipmentLinks = await ctx.db
    .query("experimentEquipment")
    .withIndex("by_experiment", q => q.eq("experimentId", experiment._id))
    .collect();
  
  if (equipmentLinks.length === 0) {
    return {
      ...experiment,
      equipment: []
    };
  }
  
  // Get equipment details
  const equipmentDetails = [];
  for (const link of equipmentLinks) {
    const equipment = await ctx.db.get(link.equipmentId);
    if (equipment) {
      const equipmentType = await ctx.db.get(equipment.equipmentTypeId);
      equipmentDetails.push({
        _id: equipment._id,
        name: equipment.name,
        type: equipmentType?.name || "Unknown",
        role: link.role
      });
    }
  }
  
  return {
    ...experiment,
    equipment: equipmentDetails
  };
}

/**
 * Expands an enhanced experiment with its time series
 */
export async function expandEnhancedExperimentWithTimeSeries(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperiment & { timeSeries?: TimeSeries[] }> {
  // Get time series
  const timeSeries = await ctx.db
    .query("timeSeries")
    .withIndex("by_experiment", q => q.eq("experimentId", experiment._id))
    .collect();
  
  return {
    ...experiment,
    timeSeries
  };
}

/**
 * Expands an enhanced experiment with all its details
 */
export async function expandEnhancedExperimentWithDetails(
  ctx: { db: DatabaseReader },
  experiment: EnhancedExperiment
): Promise<EnhancedExperimentWithDetails> {
  // Get results
  const withResults = await expandEnhancedExperimentWithResults(ctx, experiment);
  
  // Get protocol
  const withProtocol = await expandEnhancedExperimentWithProtocol(ctx, withResults);
  
  // Get mixture
  const withMixture = await expandEnhancedExperimentWithMixture(ctx, withProtocol);
  
  // Get tissue types
  const withTissueTypes = await expandEnhancedExperimentWithTissueTypes(ctx, withMixture);
  
  // Get equipment
  const withEquipment = await expandEnhancedExperimentWithEquipment(ctx, withTissueTypes);
  
  // Get time series
  const withTimeSeries = await expandEnhancedExperimentWithTimeSeries(ctx, withEquipment);
  
  return withTimeSeries as EnhancedExperimentWithDetails;
}

/**
 * Expands a time series with its data points
 */
export async function expandTimeSeriesWithData(
  ctx: { db: DatabaseReader },
  timeSeries: TimeSeries
): Promise<TimeSeries & { dataPoints: TimeSeriesDataPoint[] }> {
  const dataPoints = await ctx.db
    .query("timeSeriesData")
    .withIndex("by_series", q => q.eq("timeSeriesId", timeSeries._id))
    .collect();
  
  return {
    ...timeSeries,
    dataPoints
  };
}

/**
 * Checks if an enhanced experiment has results
 */
export async function enhancedExperimentHasResults(
  ctx: { db: DatabaseReader },
  experimentId: Id<"enhancedExperiments">
): Promise<boolean> {
  const result = await ctx.db
    .query("enhancedExperimentResults")
    .withIndex("by_experiment", q => q.eq("experimentId", experimentId))
    .first();
  
  return result !== null;
}

/**
 * Checks if an enhanced experiment has time series data
 */
export async function enhancedExperimentHasTimeSeries(
  ctx: { db: DatabaseReader },
  experimentId: Id<"enhancedExperiments">
): Promise<boolean> {
  const result = await ctx.db
    .query("timeSeries")
    .withIndex("by_experiment", q => q.eq("experimentId", experimentId))
    .first();
  
  return result !== null;
}

/**
 * Checks if an enhanced experiment has equipment
 */
export async function enhancedExperimentHasEquipment(
  ctx: { db: DatabaseReader },
  experimentId: Id<"enhancedExperiments">
): Promise<boolean> {
  const result = await ctx.db
    .query("experimentEquipment")
    .withIndex("by_experiment", q => q.eq("experimentId", experimentId))
    .first();
  
  return result !== null;
}