/**
 * Uncertainty Quantification Functions for Convex
 * 
 * This file contains functions for performing advanced uncertainty analysis
 * on experimental time series data.
 */

import { v } from 'convex/values';
import { query, mutation } from './_generated/server';
import { Id } from './_generated/dataModel';

/**
 * Retrieve uncertainty analysis for a time series
 */
export const getUncertaintyAnalysis = query({
  args: {
    timeSeriesId: v.id('timeSeries'),
    method: v.union(
      v.literal('confidence_interval'),
      v.literal('monte_carlo'),
      v.literal('sensitivity'),
      v.literal('error_propagation')
    ),
    confidenceLevel: v.optional(v.number()),
    monteCarloSamples: v.optional(v.number()),
    sensitivityParameters: v.optional(v.array(v.string())),
  },
  handler: async (ctx, args) => {
    // Get the time series data
    const timeSeriesData = await ctx.db
      .query('timeSeriesData')
      .withIndex('by_timeSeriesId', q => 
        q.eq('timeSeriesId', args.timeSeriesId)
      )
      .order('asc')
      .collect();
    
    // Get the time series metadata
    const timeSeries = await ctx.db.get(args.timeSeriesId);
    if (!timeSeries) {
      throw new Error('Time series not found');
    }
    
    // Get the experiment data if available
    let experiment = null;
    if (timeSeries.experimentId) {
      experiment = await ctx.db.get(timeSeries.experimentId as Id<'enhancedExperiments'>);
    }
    
    // Process based on the selected method
    switch (args.method) {
      case 'confidence_interval':
        return calculateConfidenceIntervals(
          timeSeriesData, 
          args.confidenceLevel || 0.95
        );
      
      case 'monte_carlo':
        return generateMonteCarloSimulation(
          timeSeriesData, 
          args.monteCarloSamples || 100
        );
      
      case 'sensitivity':
        return performSensitivityAnalysis(
          timeSeriesData, 
          experiment, 
          args.sensitivityParameters || []
        );
      
      case 'error_propagation':
        return performErrorPropagation(timeSeriesData);
      
      default:
        // Default to basic confidence intervals
        return calculateConfidenceIntervals(timeSeriesData, 0.95);
    }
  }
});

/**
 * Store uncertainty analysis results for a time series
 */
export const saveUncertaintyAnalysis = mutation({
  args: {
    timeSeriesId: v.id('timeSeries'),
    analysisType: v.string(),
    results: v.any(),
    parameters: v.optional(v.object({
      confidenceLevel: v.optional(v.number()),
      monteCarloSamples: v.optional(v.number()),
      sensitivityParameters: v.optional(v.array(v.string())),
    }))
  },
  handler: async (ctx, args) => {
    // Check if an analysis already exists
    const existingAnalysis = await ctx.db
      .query('uncertaintyAnalysis')
      .withIndex('by_timeSeriesAndType', q => 
        q.eq('timeSeriesId', args.timeSeriesId)
         .eq('analysisType', args.analysisType)
      )
      .first();
    
    // If it exists, update it
    if (existingAnalysis) {
      return await ctx.db.patch(existingAnalysis._id, {
        results: args.results,
        parameters: args.parameters,
        updatedAt: Date.now()
      });
    }
    
    // Otherwise, create a new one
    return await ctx.db.insert('uncertaintyAnalysis', {
      timeSeriesId: args.timeSeriesId,
      analysisType: args.analysisType,
      results: args.results,
      parameters: args.parameters,
      createdAt: Date.now(),
      updatedAt: Date.now()
    });
  }
});

/**
 * Calculate confidence intervals for time series data
 */
function calculateConfidenceIntervals(timeSeriesData: any[], confidenceLevel: number) {
  // Z-score lookup for common confidence levels
  const zScoreMap: Record<number, number> = {
    0.50: 0.674,
    0.80: 1.282,
    0.90: 1.645,
    0.95: 1.960,
    0.99: 2.576
  };
  
  // Get the appropriate z-score or default to 95% confidence
  const zScore = zScoreMap[confidenceLevel] || 1.960;
  
  // Calculate confidence intervals for each data point
  return timeSeriesData.map(point => {
    const result = {
      timestamp: point.timestamp,
      value: point.value,
      standardDeviation: point.uncertainty,
      metadata: point.metadata || {}
    };
    
    // If uncertainty is available, calculate confidence interval
    if (point.uncertainty) {
      result.metadata.confidenceInterval = [
        point.value - zScore * point.uncertainty,
        point.value + zScore * point.uncertainty
      ];
    }
    
    return result;
  });
}

/**
 * Generate Monte Carlo simulation data
 */
function generateMonteCarloSimulation(timeSeriesData: any[], samples: number) {
  return timeSeriesData.map(point => {
    const result = {
      timestamp: point.timestamp,
      value: point.value,
      standardDeviation: point.uncertainty,
      metadata: point.metadata || {}
    };
    
    // If uncertainty is available, generate Monte Carlo samples
    if (point.uncertainty) {
      result.metadata.monteCarlo = Array.from({ length: samples }, () => 
        generateNormalRandom(point.value, point.uncertainty));
    }
    
    return result;
  });
}

/**
 * Perform sensitivity analysis
 */
function performSensitivityAnalysis(
  timeSeriesData: any[], 
  experiment: any | null, 
  parameters: string[]
) {
  // Default parameters to analyze if none specified
  const defaultParameters = [
    'temperature', 
    'concentration', 
    'pressure', 
    'ph', 
    'mixing_rate'
  ];
  
  // Use specified parameters or defaults
  const paramsToAnalyze = parameters.length > 0 
    ? parameters 
    : defaultParameters;
  
  // Get experiment parameters if available
  const experimentParams: Record<string, number> = {};
  if (experiment) {
    if (experiment.temperature) experimentParams.temperature = experiment.temperature;
    if (experiment.concentration) experimentParams.concentration = experiment.concentration;
    if (experiment.pressure) experimentParams.pressure = experiment.pressure;
    if (experiment.ph) experimentParams.ph = experiment.ph;
    if (experiment.mixingRate) experimentParams.mixing_rate = experiment.mixingRate;
    
    // Add any additional parameters from metadata
    if (experiment.parameters) {
      Object.entries(experiment.parameters).forEach(([key, value]) => {
        if (typeof value === 'number') {
          experimentParams[key] = value;
        }
      });
    }
  }
  
  return timeSeriesData.map(point => {
    const result = {
      timestamp: point.timestamp,
      value: point.value,
      standardDeviation: point.uncertainty,
      metadata: { ...(point.metadata || {}) }
    };
    
    // Calculate synthetic sensitivity scores
    // In a real implementation, this would use actual parameter correlations
    const sensitivityScores: Record<string, number> = {};
    
    paramsToAnalyze.forEach(param => {
      // Generate a synthetic sensitivity score
      // This is just for demonstration - in a real system, this would be calculated
      // based on actual data and parameter correlations
      sensitivityScores[param] = generateSyntheticSensitivity(
        param, 
        point.timestamp, 
        experimentParams[param]
      );
    });
    
    result.metadata.sensitivityScores = sensitivityScores;
    return result;
  });
}

/**
 * Perform error propagation analysis
 */
function performErrorPropagation(timeSeriesData: any[]) {
  return timeSeriesData.map(point => {
    const result = {
      timestamp: point.timestamp,
      value: point.value,
      standardDeviation: point.uncertainty,
      metadata: { ...(point.metadata || {}) }
    };
    
    // In a real implementation, this would calculate error propagation
    // based on the mathematical model and input uncertainties
    // For demonstration, we just use the standard deviation
    
    return result;
  });
}

/**
 * Generate a random number from normal distribution
 */
function generateNormalRandom(mean: number, stdDev: number): number {
  // Box-Muller transform for normal distribution
  const u1 = Math.random();
  const u2 = Math.random();
  const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
  return mean + stdDev * z0;
}

/**
 * Generate a synthetic sensitivity score for demonstration
 */
function generateSyntheticSensitivity(
  parameter: string, 
  timestamp: number, 
  paramValue?: number
): number {
  // Generate parameter-specific sensitivity patterns
  // This is purely synthetic for demonstration purposes
  const timeVariation = Math.sin(timestamp / 10000000) * 0.2 + 0.5;
  
  switch (parameter) {
    case 'temperature':
      return 0.8 * timeVariation + 0.1;
    case 'concentration':
      return 0.75 * (1 - timeVariation) + 0.2;
    case 'pressure':
      return 0.5 * timeVariation + 0.3;
    case 'ph':
      return 0.4 * Math.cos(timestamp / 8000000) + 0.5;
    case 'mixing_rate':
      return 0.3 * Math.sin(timestamp / 5000000) + 0.4;
    default:
      // For other parameters, generate a random sensitivity
      return Math.random() * 0.5 + 0.2;
  }
}