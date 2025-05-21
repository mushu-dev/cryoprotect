/**
 * Helper functions for predictions and scientific models
 */

import { Id } from "../_generated/dataModel";
import { QueryCtx } from "../_generated/server";
import { 
  Prediction, 
  ScientificModel,
  PredictionWithModel,
  PredictionWithTarget
} from "./types";

/**
 * Expand a prediction with its model data
 */
export async function expandPredictionWithModel(
  ctx: QueryCtx,
  prediction: Prediction
): Promise<PredictionWithModel> {
  const model = await ctx.db.get(prediction.modelId);
  
  if (!model) {
    throw new Error(`Scientific model with ID ${prediction.modelId} not found`);
  }
  
  return {
    ...prediction,
    model: {
      _id: model._id,
      name: model.name,
      version: model.version,
      type: model.type
    }
  };
}

/**
 * Expand a prediction with its target (molecule or mixture) data
 */
export async function expandPredictionWithTarget(
  ctx: QueryCtx,
  prediction: Prediction
): Promise<PredictionWithTarget> {
  // Start with model data
  const expanded = await expandPredictionWithModel(ctx, prediction) as PredictionWithTarget;
  
  // Add molecule data if applicable
  if (prediction.moleculeId) {
    const molecule = await ctx.db.get(prediction.moleculeId);
    
    if (molecule) {
      expanded.molecule = {
        _id: molecule._id,
        name: molecule.name,
        formula: molecule.formula
      };
    }
  }
  
  // Add mixture data if applicable
  if (prediction.mixtureId) {
    const mixture = await ctx.db.get(prediction.mixtureId);
    
    if (mixture) {
      expanded.mixture = {
        _id: mixture._id,
        name: mixture.name
      };
    }
  }
  
  return expanded;
}

/**
 * Format model type for display
 */
export function formatModelType(type: string): string {
  // Format camelCase or snake_case to Title Case
  return type
    .replace(/([A-Z])/g, ' $1') // camelCase to space separated
    .replace(/_/g, ' ')         // snake_case to space separated
    .replace(/^\w/, c => c.toUpperCase()) // capitalize first letter
    .replace(/\s+/g, ' ')       // normalize spaces
    .trim();
}

/**
 * Calculate prediction accuracy when comparing to experimental results
 */
export function calculatePredictionAccuracy(
  predictedValue: number,
  actualValue: number
): number {
  // Return 0 if either value is not usable
  if (typeof predictedValue !== 'number' || typeof actualValue !== 'number') {
    return 0;
  }
  
  // Prevent division by zero
  if (actualValue === 0) {
    return predictedValue === 0 ? 1 : 0;
  }
  
  // Calculate accuracy as 1 - (relative error), bound between 0 and 1
  const relativeError = Math.abs((predictedValue - actualValue) / actualValue);
  const accuracy = Math.max(0, Math.min(1, 1 - relativeError));
  
  return accuracy;
}

/**
 * Group predictions by parameter
 */
export function groupPredictionsByParameter(
  predictions: Prediction[]
): Record<string, Prediction[]> {
  const grouped: Record<string, Prediction[]> = {};
  
  for (const prediction of predictions) {
    if (!grouped[prediction.parameterName]) {
      grouped[prediction.parameterName] = [];
    }
    
    grouped[prediction.parameterName].push(prediction);
  }
  
  return grouped;
}

/**
 * Check if a scientific model has predictions
 */
export async function modelHasPredictions(
  ctx: QueryCtx,
  modelId: Id<"scientificModels">
): Promise<boolean> {
  const count = await ctx.db
    .query("predictions")
    .withIndex("by_model", q => q.eq("modelId", modelId))
    .count();
  
  return count > 0;
}

/**
 * Parse model parameters from string to object
 * Parameters are stored as a JSON string in the database
 */
export function parseModelParameters(parameters?: string): any {
  if (!parameters) {
    return {};
  }
  
  try {
    return JSON.parse(parameters);
  } catch (error) {
    console.error("Failed to parse model parameters:", error);
    return {};
  }
}

/**
 * Stringify model parameters from object to string
 */
export function stringifyModelParameters(parameters?: any): string {
  if (!parameters) {
    return "";
  }
  
  try {
    return JSON.stringify(parameters);
  } catch (error) {
    console.error("Failed to stringify model parameters:", error);
    return "";
  }
}