/**
 * Supabase to Convex Migration Utility
 * 
 * This module provides utilities for migrating data from Supabase
 * to our enhanced Convex schema. It includes data transformation,
 * enrichment, and validation.
 */

import { createClient } from '@supabase/supabase-js';
import { ConvexClient } from 'convex/browser';
import { Id } from '../_generated/dataModel';
import { api } from '../_generated/api';

// Define interface for migration config
interface MigrationConfig {
  supabaseUrl: string;
  supabaseKey: string;
  convexClient: ConvexClient;
  batchSize?: number;
  logLevel?: 'debug' | 'info' | 'warn' | 'error';
}

// Define interface for migration stats
interface MigrationStats {
  table: string;
  totalRecords: number;
  successfulRecords: number;
  failedRecords: number;
  errors: Array<{
    record: any;
    error: string;
  }>;
  startTime: number;
  endTime?: number;
  durationMs?: number;
}

// Define interface for migration context
interface MigrationContext {
  config: MigrationConfig;
  stats: Record<string, MigrationStats>;
  mappings: Record<string, Record<string, Id<any>>>;
  logMessage: (level: 'debug' | 'info' | 'warn' | 'error', message: string) => void;
}

/**
 * Initialize a migration context
 */
export function initMigrationContext(config: MigrationConfig): MigrationContext {
  const logLevel = config.logLevel || 'info';
  const logLevels = { debug: 0, info: 1, warn: 2, error: 3 };
  
  return {
    config: {
      ...config,
      batchSize: config.batchSize || 100
    },
    stats: {},
    mappings: {},
    logMessage: (level, message) => {
      if (logLevels[level] >= logLevels[logLevel]) {
        const timestamp = new Date().toISOString();
        console.log(`[${timestamp}] [${level.toUpperCase()}] ${message}`);
      }
    }
  };
}

/**
 * Create a Supabase client for data extraction
 */
export function createSupabaseClient(url: string, key: string) {
  return createClient(url, key);
}

/**
 * Initialize migration stats for a table
 */
function initTableStats(table: string): MigrationStats {
  return {
    table,
    totalRecords: 0,
    successfulRecords: 0,
    failedRecords: 0,
    errors: [],
    startTime: Date.now()
  };
}

/**
 * Complete migration stats for a table
 */
function completeTableStats(stats: MigrationStats): MigrationStats {
  return {
    ...stats,
    endTime: Date.now(),
    durationMs: Date.now() - stats.startTime
  };
}

/**
 * Migrate users from Supabase to Convex
 */
export async function migrateUsers(ctx: MigrationContext): Promise<Record<string, Id<"users">>> {
  const { config, logMessage } = ctx;
  const supabase = createSupabaseClient(config.supabaseUrl, config.supabaseKey);
  const convexClient = config.convexClient;
  
  // Initialize stats
  const tableName = 'users';
  ctx.stats[tableName] = initTableStats(tableName);
  const stats = ctx.stats[tableName];
  
  logMessage('info', `Starting migration of ${tableName}`);
  
  // Fetch users from Supabase
  const { data: supabaseUsers, error } = await supabase
    .from('user_profile')
    .select('*');
  
  if (error) {
    logMessage('error', `Error fetching users: ${error.message}`);
    stats.errors.push({ record: null, error: error.message });
    return {};
  }
  
  stats.totalRecords = supabaseUsers.length;
  logMessage('info', `Found ${supabaseUsers.length} users to migrate`);
  
  // Create mapping for user IDs
  const userIdMapping: Record<string, Id<"users">> = {};
  
  // Process users in batches
  const batchSize = config.batchSize!;
  for (let i = 0; i < supabaseUsers.length; i += batchSize) {
    const batch = supabaseUsers.slice(i, i + batchSize);
    logMessage('debug', `Processing batch ${i/batchSize + 1}, size: ${batch.length}`);
    
    // Create Convex users (serially to avoid rate limits)
    for (const user of batch) {
      try {
        // Transform user data to Convex schema
        const convexUser = {
          email: user.email,
          name: user.display_name || user.email.split('@')[0],
          role: user.role || 'viewer',
          lastLogin: user.last_sign_in_at ? new Date(user.last_sign_in_at).getTime() : undefined,
          createdAt: user.created_at ? new Date(user.created_at).getTime() : Date.now(),
          updatedAt: Date.now()
        };
        
        // Insert into Convex
        const convexUserId = await convexClient.mutation(api.users.create, convexUser);
        
        // Store ID mapping
        userIdMapping[user.id] = convexUserId;
        stats.successfulRecords++;
      } catch (error: any) {
        logMessage('error', `Error migrating user ${user.id}: ${error.message}`);
        stats.errors.push({ record: user, error: error.message });
        stats.failedRecords++;
      }
    }
    
    logMessage('info', `Processed ${Math.min((i + batchSize), supabaseUsers.length)} / ${supabaseUsers.length} users`);
  }
  
  // Complete stats
  ctx.stats[tableName] = completeTableStats(stats);
  logMessage('info', `Completed migration of ${stats.successfulRecords} / ${stats.totalRecords} users in ${ctx.stats[tableName].durationMs}ms`);
  
  // Store mapping in context
  ctx.mappings['users'] = userIdMapping;
  
  return userIdMapping;
}

/**
 * Migrate tissue types from Supabase to enhanced Convex schema
 */
export async function migrateTissueTypes(ctx: MigrationContext): Promise<Record<string, Id<"tissueTypes">>> {
  const { config, logMessage } = ctx;
  const supabase = createSupabaseClient(config.supabaseUrl, config.supabaseKey);
  const convexClient = config.convexClient;
  
  // Initialize stats
  const tableName = 'tissueTypes';
  ctx.stats[tableName] = initTableStats(tableName);
  const stats = ctx.stats[tableName];
  
  logMessage('info', `Starting migration of ${tableName}`);
  
  // Fetch tissue types from Supabase (likely found in experiments or a dedicated table)
  // Note: This query might need to be adjusted based on actual Supabase schema
  const { data: supabaseTissueData, error } = await supabase
    .from('tissue_types')
    .select('*');
  
  if (error) {
    logMessage('error', `Error fetching tissue types: ${error.message}`);
    stats.errors.push({ record: null, error: error.message });
    return {};
  }
  
  stats.totalRecords = supabaseTissueData.length;
  logMessage('info', `Found ${supabaseTissueData.length} tissue types to migrate`);
  
  // Create mapping for tissue type IDs
  const tissueTypeIdMapping: Record<string, Id<"tissueTypes">> = {};
  const userIdMapping = ctx.mappings['users'] || {};
  
  // Process tissue types in batches
  const batchSize = config.batchSize!;
  for (let i = 0; i < supabaseTissueData.length; i += batchSize) {
    const batch = supabaseTissueData.slice(i, i + batchSize);
    logMessage('debug', `Processing batch ${i/batchSize + 1}, size: ${batch.length}`);
    
    // Create Convex tissue types (serially to maintain data integrity)
    for (const tissue of batch) {
      try {
        // Transform tissue data to enhanced Convex schema
        const createdBy = tissue.created_by ? userIdMapping[tissue.created_by] : undefined;
        
        // Extract values from tissue data or provide defaults
        // Parse properties from JSON if available
        const properties = tissue.properties ? 
          (typeof tissue.properties === 'string' ? 
            JSON.parse(tissue.properties) : tissue.properties) : {};
            
        // Extract cell type from name or description if available
        let cellType = undefined;
        if (tissue.name) {
          const commonCellTypes = [
            "Hepatocytes", "Oocytes", "Embryos", "Sperm", 
            "Stem cells", "Fibroblasts", "Erythrocytes", "Neurons"
          ];
          
          for (const type of commonCellTypes) {
            if (tissue.name.includes(type) || (tissue.description && tissue.description.includes(type))) {
              cellType = type;
              break;
            }
          }
        }
        
        // Enhanced tissue type with enriched data
        const convexTissueType = {
          name: tissue.name,
          description: tissue.description,
          species: tissue.species,
          taxonomyId: tissue.taxonomy_id ? Number(tissue.taxonomy_id) : undefined,
          
          // Enhanced biospecimen fields (enriched from available data)
          cellType,
          tissueOrigin: properties.origin,
          cellDensity: properties.cell_density,
          cellDensityUnit: properties.cell_density ? 'million/mL' : undefined,
          passageNumber: properties.passage_number,
          cellDiameter: properties.cell_size,
          cellDiameterUnit: properties.cell_size ? 'μm' : undefined,
          waterContent: properties.water_content,
          lipidContent: properties.lipid_content,
          osmolality: properties.osmolality,
          osmolalityUnit: properties.osmolality ? 'mOsm/kg' : undefined,
          preparationMethod: properties.preparation_method,
          storageConditions: properties.storage_conditions,
          cryopreservationHistory: properties.previously_frozen,
          
          properties,
          createdBy,
          createdAt: tissue.created_at ? new Date(tissue.created_at).getTime() : Date.now(),
          updatedAt: Date.now()
        };
        
        // Insert into Convex
        const convexTissueId = await convexClient.mutation(api.tissueTypes.create, convexTissueType);
        
        // Store ID mapping
        tissueTypeIdMapping[tissue.id] = convexTissueId;
        stats.successfulRecords++;
      } catch (error: any) {
        logMessage('error', `Error migrating tissue type ${tissue.id}: ${error.message}`);
        stats.errors.push({ record: tissue, error: error.message });
        stats.failedRecords++;
      }
    }
    
    logMessage('info', `Processed ${Math.min((i + batchSize), supabaseTissueData.length)} / ${supabaseTissueData.length} tissue types`);
  }
  
  // Complete stats
  ctx.stats[tableName] = completeTableStats(stats);
  logMessage('info', `Completed migration of ${stats.successfulRecords} / ${stats.totalRecords} tissue types in ${ctx.stats[tableName].durationMs}ms`);
  
  // Store mapping in context
  ctx.mappings['tissueTypes'] = tissueTypeIdMapping;
  
  return tissueTypeIdMapping;
}

/**
 * Migrate protocols from Supabase to enhanced Convex schema
 */
export async function migrateProtocols(ctx: MigrationContext): Promise<Record<string, Id<"protocols">>> {
  const { config, logMessage } = ctx;
  const supabase = createSupabaseClient(config.supabaseUrl, config.supabaseKey);
  const convexClient = config.convexClient;
  
  // Initialize stats
  const tableName = 'protocols';
  ctx.stats[tableName] = initTableStats(tableName);
  const stats = ctx.stats[tableName];
  
  logMessage('info', `Starting migration of ${tableName}`);
  
  // Fetch protocols from Supabase
  const { data: supabaseProtocols, error } = await supabase
    .from('experimental_protocols')
    .select('*');
  
  if (error) {
    logMessage('error', `Error fetching protocols: ${error.message}`);
    stats.errors.push({ record: null, error: error.message });
    return {};
  }
  
  stats.totalRecords = supabaseProtocols.length;
  logMessage('info', `Found ${supabaseProtocols.length} protocols to migrate`);
  
  // Create mapping for protocol IDs
  const protocolIdMapping: Record<string, Id<"protocols">> = {};
  const userIdMapping = ctx.mappings['users'] || {};
  
  // Process protocols in batches
  const batchSize = config.batchSize!;
  for (let i = 0; i < supabaseProtocols.length; i += batchSize) {
    const batch = supabaseProtocols.slice(i, i + batchSize);
    logMessage('debug', `Processing batch ${i/batchSize + 1}, size: ${batch.length}`);
    
    // Create Convex protocols (serially to maintain data integrity)
    for (const protocol of batch) {
      try {
        // Transform protocol data to enhanced Convex schema
        const createdBy = protocol.created_by ? userIdMapping[protocol.created_by] : undefined;
        
        // Parse steps from protocol content or steps field if available
        let steps: any[] = [];
        if (protocol.steps && typeof protocol.steps === 'string') {
          try {
            steps = JSON.parse(protocol.steps);
          } catch (e) {
            // If steps JSON parsing fails, try to extract steps from content
            steps = extractStepsFromContent(protocol.content || "");
          }
        } else if (Array.isArray(protocol.steps)) {
          steps = protocol.steps;
        } else {
          // Extract steps from content as a fallback
          steps = extractStepsFromContent(protocol.content || "");
        }
        
        // Identify protocol type from name or description
        let protocolType = identifyProtocolType(protocol.name, protocol.description);
        
        // Extract cryopreservation parameters from protocol fields or content
        const cryoParams = extractCryoParameters(protocol, steps);
        
        // Format steps to match Convex schema
        const formattedSteps = steps.map((step, index) => {
          // Ensure step has an id
          const stepId = step.id || `step-${index + 1}`;
          
          // Identify if this is a critical step based on keywords
          const isCritical = isCriticalStep(step);
          
          // Extract substances added in this step
          const substancesAdded = extractSubstancesAdded(step);
          
          return {
            id: stepId,
            name: step.name || `Step ${index + 1}`,
            description: step.description || "",
            parameters: step.parameters || {},
            duration: step.duration,
            durationUnit: step.durationUnit || "minutes",
            temperature: step.temperature,
            temperatureUnit: step.temperatureUnit || "°C",
            
            // Enhanced cryopreservation fields
            rampRate: step.rampRate,
            rampRateUnit: step.rampRateUnit || "°C/min",
            pressureApplied: step.pressure,
            pressureUnit: step.pressureUnit || "atm",
            substancesAdded,
            equipmentRequired: step.equipment || [],
            criticalStep: isCritical,
            qualityControlChecks: isCritical ? 
              ["Visual inspection", "Temperature verification", "Timing verification"] : undefined
          };
        });
        
        // Structured protocol with enhanced fields
        const convexProtocol = {
          name: protocol.name,
          description: protocol.description,
          
          // Cryopreservation-specific fields
          protocolType,
          coolingRate: cryoParams.coolingRate,
          coolingRateUnit: cryoParams.coolingRate ? "°C/min" : undefined,
          warmingRate: cryoParams.warmingRate,
          warmingRateUnit: cryoParams.warmingRate ? "°C/min" : undefined,
          holdTemperature: cryoParams.holdTemperature,
          holdDuration: cryoParams.holdDuration,
          holdDurationUnit: cryoParams.holdDuration ? "minutes" : undefined,
          cryoprotectantAdditionMethod: cryoParams.cryoprotectantAdditionMethod,
          preFreezingTreatment: cryoParams.preFreezingTreatment,
          postThawingTreatment: cryoParams.postThawingTreatment,
          storageTemperature: cryoParams.storageTemperature,
          storageContainerType: cryoParams.storageContainerType,
          cpaEquilibrationTime: cryoParams.cpaEquilibrationTime,
          
          steps: formattedSteps,
          version: protocol.version || 1,
          parentId: protocol.parent_id ? protocolIdMapping[protocol.parent_id] : undefined,
          category: protocol.category,
          isTemplate: protocol.is_template === true,
          parameters: protocol.parameters || {},
          validationStatus: protocol.validated ? "validated" : "unvalidated",
          publications: protocol.references ? 
            extractPublicationReferences(protocol.references) : undefined,
          successRate: protocol.success_rate,
          
          createdBy,
          createdAt: protocol.created_at ? new Date(protocol.created_at).getTime() : Date.now(),
          updatedAt: Date.now(),
          public: protocol.public === true
        };
        
        // Insert into Convex
        const convexProtocolId = await convexClient.mutation(api.protocols.create, convexProtocol);
        
        // Store ID mapping
        protocolIdMapping[protocol.id] = convexProtocolId;
        stats.successfulRecords++;
      } catch (error: any) {
        logMessage('error', `Error migrating protocol ${protocol.id}: ${error.message}`);
        stats.errors.push({ record: protocol, error: error.message });
        stats.failedRecords++;
      }
    }
    
    logMessage('info', `Processed ${Math.min((i + batchSize), supabaseProtocols.length)} / ${supabaseProtocols.length} protocols`);
  }
  
  // Complete stats
  ctx.stats[tableName] = completeTableStats(stats);
  logMessage('info', `Completed migration of ${stats.successfulRecords} / ${stats.totalRecords} protocols in ${ctx.stats[tableName].durationMs}ms`);
  
  // Store mapping in context
  ctx.mappings['protocols'] = protocolIdMapping;
  
  return protocolIdMapping;
}

/**
 * Migrate experiment results to enhanced Convex schema
 */
export async function migrateExperimentResults(
  ctx: MigrationContext, 
  experimentIdMapping: Record<string, Id<"enhancedExperiments">>,
  tissueTypeIdMapping: Record<string, Id<"tissueTypes">>
): Promise<Record<string, Id<"enhancedExperimentResults">>> {
  const { config, logMessage } = ctx;
  const supabase = createSupabaseClient(config.supabaseUrl, config.supabaseKey);
  const convexClient = config.convexClient;
  
  // Initialize stats
  const tableName = 'enhancedExperimentResults';
  ctx.stats[tableName] = initTableStats(tableName);
  const stats = ctx.stats[tableName];
  
  logMessage('info', `Starting migration of ${tableName}`);
  
  // Fetch experiment results from Supabase
  const { data: supabaseResults, error } = await supabase
    .from('experiment_results')
    .select('*');
  
  if (error) {
    logMessage('error', `Error fetching experiment results: ${error.message}`);
    stats.errors.push({ record: null, error: error.message });
    return {};
  }
  
  stats.totalRecords = supabaseResults.length;
  logMessage('info', `Found ${supabaseResults.length} experiment results to migrate`);
  
  // Create mapping for result IDs
  const resultIdMapping: Record<string, Id<"enhancedExperimentResults">> = {};
  const userIdMapping = ctx.mappings['users'] || {};
  
  // Process results in batches
  const batchSize = config.batchSize!;
  for (let i = 0; i < supabaseResults.length; i += batchSize) {
    const batch = supabaseResults.slice(i, i + batchSize);
    logMessage('debug', `Processing batch ${i/batchSize + 1}, size: ${batch.length}`);
    
    // Create Convex experiment results (serially to maintain data integrity)
    for (const result of batch) {
      try {
        const experimentId = result.experiment_id ? 
          experimentIdMapping[result.experiment_id] : undefined;
          
        if (!experimentId) {
          throw new Error(`No experiment mapping found for experiment_id: ${result.experiment_id}`);
        }
        
        // Get tissue type ID from experiment or use a default
        const tissueTypeId = result.tissue_type_id ? 
          tissueTypeIdMapping[result.tissue_type_id] : 
          await getDefaultTissueTypeId(convexClient, result, ctx);
        
        if (!tissueTypeId) {
          throw new Error(`No tissue type mapping found and couldn't create default`);
        }
        
        // Extract method information from result data
        const methodInfo = extractMethodInformation(result);
        
        // Transform to enhanced experiment result
        const convexResult = {
          experimentId,
          tissueTypeId,
          moleculeId: undefined, // Would need molecule mapping
          mixtureId: undefined, // Would need mixture mapping
          concentration: result.concentration,
          concentrationUnit: result.concentration_unit,
          
          // Enhanced viability measurements
          viabilityPercentage: result.viability_percentage,
          viabilityMethod: methodInfo.viabilityMethod,
          recoveryRate: result.recovery_rate,
          recoveryRateUnit: result.recovery_rate ? "%" : undefined,
          functionalityScore: result.functionality_score,
          functionalAssay: methodInfo.functionalAssay,
          integrityMeasure: result.integrity || result.membrane_integrity,
          integrityMethod: methodInfo.integrityMethod,
          
          // Cryopreservation-specific metrics
          iceFormationObserved: getIceFormation(result),
          postThawMorphologyScore: result.morphology_score,
          stressResponseMarkers: getStressMarkers(result),
          timeToRecovery: result.recovery_time,
          
          uncertainty: result.uncertainty || {
            viabilityPercentage: result.viability_uncertainty ? {
              type: "standard_deviation",
              value: result.viability_uncertainty
            } : undefined
          },
          resultDetails: result.details || result.result_details || {},
          notes: result.notes || result.comments,
          createdBy: result.created_by ? userIdMapping[result.created_by] : undefined,
          createdAt: result.created_at ? new Date(result.created_at).getTime() : Date.now(),
          updatedAt: Date.now()
        };
        
        // Insert into Convex
        const convexResultId = await convexClient.mutation(
          api.enhancedExperimentResults.create, 
          convexResult
        );
        
        // Store ID mapping
        resultIdMapping[result.id] = convexResultId;
        stats.successfulRecords++;
      } catch (error: any) {
        logMessage('error', `Error migrating experiment result ${result.id}: ${error.message}`);
        stats.errors.push({ record: result, error: error.message });
        stats.failedRecords++;
      }
    }
    
    logMessage('info', `Processed ${Math.min((i + batchSize), supabaseResults.length)} / ${supabaseResults.length} experiment results`);
  }
  
  // Complete stats
  ctx.stats[tableName] = completeTableStats(stats);
  logMessage('info', `Completed migration of ${stats.successfulRecords} / ${stats.totalRecords} experiment results in ${ctx.stats[tableName].durationMs}ms`);
  
  // Store mapping in context
  ctx.mappings['enhancedExperimentResults'] = resultIdMapping;
  
  return resultIdMapping;
}

// Helper functions for data extraction and transformation

/**
 * Extract steps from protocol content text
 */
function extractStepsFromContent(content: string): any[] {
  if (!content) return [];
  
  // Simple pattern-based extraction (would need refinement based on actual content format)
  const stepPattern = /Step\s+(\d+)[:.-]\s*(.*?)(?=Step\s+\d+[:.-]|$)/gis;
  const steps: any[] = [];
  
  let match;
  while ((match = stepPattern.exec(content)) !== null) {
    const stepNumber = match[1];
    const stepContent = match[2].trim();
    
    // Extract step parameters
    const temperatureMatch = stepContent.match(/(\-?\d+(?:\.\d+)?)\s*°C/);
    const durationMatch = stepContent.match(/(\d+(?:\.\d+)?)\s*(min|hour|sec|h|s)/i);
    
    steps.push({
      id: `step-${stepNumber}`,
      name: `Step ${stepNumber}`,
      description: stepContent,
      parameters: {},
      duration: durationMatch ? parseFloat(durationMatch[1]) : undefined,
      durationUnit: durationMatch ? normalizeTimeUnit(durationMatch[2]) : undefined,
      temperature: temperatureMatch ? parseFloat(temperatureMatch[1]) : undefined,
      temperatureUnit: temperatureMatch ? "°C" : undefined
    });
  }
  
  return steps;
}

/**
 * Normalize time unit to standard format
 */
function normalizeTimeUnit(unit: string): string {
  unit = unit.toLowerCase();
  if (unit === 'min' || unit === 'minute' || unit === 'minutes') return 'minutes';
  if (unit === 'h' || unit === 'hour' || unit === 'hours') return 'hours';
  if (unit === 's' || unit === 'sec' || unit === 'second' || unit === 'seconds') return 'seconds';
  return unit;
}

/**
 * Identify protocol type based on name and description
 */
function identifyProtocolType(name?: string, description?: string): string | undefined {
  const text = `${name || ''} ${description || ''}`.toLowerCase();
  
  if (text.includes('vitrification')) return 'vitrification';
  if (text.includes('slow freezing') || text.includes('slow-freezing')) return 'slow_freezing';
  if (text.includes('controlled rate') || text.includes('controlled-rate')) return 'controlled_rate_freezing';
  if (text.includes('rapid') || text.includes('ultra-rapid')) return 'ultra-rapid_freezing';
  
  // Default to slow freezing as most common
  return text.includes('freez') ? 'slow_freezing' : undefined;
}

/**
 * Extract cryopreservation parameters from protocol data
 */
function extractCryoParameters(protocol: any, steps: any[]): any {
  const params: any = {};
  
  // Try to extract from protocol data first
  if (protocol.parameters) {
    const p = typeof protocol.parameters === 'string' ? 
      JSON.parse(protocol.parameters) : protocol.parameters;
    
    params.coolingRate = p.cooling_rate || p.freezingRate;
    params.warmingRate = p.warming_rate || p.thawingRate;
    params.holdTemperature = p.hold_temperature;
    params.holdDuration = p.hold_duration;
    params.cryoprotectantAdditionMethod = p.cpa_addition_method;
    params.preFreezingTreatment = p.pre_freezing_treatment;
    params.postThawingTreatment = p.post_thawing_treatment;
    params.storageTemperature = p.storage_temperature;
    params.storageContainerType = p.storage_container;
    params.cpaEquilibrationTime = p.cpa_equilibration_time;
  }
  
  // Try to extract from content if fields are still undefined
  const contentText = protocol.content || '';
  
  if (params.coolingRate === undefined) {
    const coolingMatch = contentText.match(/cooling\s+rate\s*(?:of|:)?\s*(\-?\d+(?:\.\d+)?)/i);
    if (coolingMatch) params.coolingRate = parseFloat(coolingMatch[1]);
  }
  
  if (params.warmingRate === undefined) {
    const warmingMatch = contentText.match(/warming\s+rate\s*(?:of|:)?\s*(\d+(?:\.\d+)?)/i);
    if (warmingMatch) params.warmingRate = parseFloat(warmingMatch[1]);
  }
  
  if (params.holdTemperature === undefined) {
    const tempMatch = contentText.match(/hold\s+(?:at|temperature)\s*(?:of|:)?\s*(\-?\d+(?:\.\d+)?)/i);
    if (tempMatch) params.holdTemperature = parseFloat(tempMatch[1]);
  }
  
  // Extract from protocol steps if still undefined
  for (const step of steps) {
    const desc = (step.description || '').toLowerCase();
    
    // Extract storage container
    if (params.storageContainerType === undefined) {
      if (desc.includes('cryovial')) params.storageContainerType = 'Cryovial';
      else if (desc.includes('straw')) params.storageContainerType = 'Straw';
      else if (desc.includes('bag')) params.storageContainerType = 'Cryobag';
    }
    
    // Extract cryoprotectant addition method
    if (params.cryoprotectantAdditionMethod === undefined) {
      if (desc.includes('step') && desc.includes('wise')) params.cryoprotectantAdditionMethod = 'step-wise';
      else if (desc.includes('single step')) params.cryoprotectantAdditionMethod = 'single-step';
      else if (desc.includes('gradient')) params.cryoprotectantAdditionMethod = 'gradient';
    }
    
    // Look for storage temperature in steps
    if (params.storageTemperature === undefined && desc.includes('liquid nitrogen')) {
      params.storageTemperature = -196;
    }
  }
  
  return params;
}

/**
 * Determine if a step is critical to the cryopreservation process
 */
function isCriticalStep(step: any): boolean {
  const text = `${step.name || ''} ${step.description || ''}`.toLowerCase();
  
  // Keywords that indicate critical steps
  const criticalKeywords = [
    'critical', 'important', 'essential', 'key', 'crucial',
    'liquid nitrogen', 'ln2', 'seeding', 'nucleation',
    'cryoprotectant', 'dimethyl sulfoxide', 'dmso', 'glycerol',
    'vitrification', 'plunge'
  ];
  
  return criticalKeywords.some(keyword => text.includes(keyword));
}

/**
 * Extract substances added in a protocol step
 */
function extractSubstancesAdded(step: any): string[] | undefined {
  const text = `${step.name || ''} ${step.description || ''}`.toLowerCase();
  
  // Common cryoprotectants and additives
  const substances = [
    'glycerol', 'dmso', 'dimethyl sulfoxide', 'ethylene glycol', 'propylene glycol',
    'trehalose', 'sucrose', 'mannitol', 'polyvinylpyrrolidone', 'pvp',
    'fetal bovine serum', 'fbs', 'bovine serum albumin', 'bsa'
  ];
  
  const found = substances.filter(substance => text.includes(substance));
  return found.length > 0 ? found : undefined;
}

/**
 * Extract publication references from text
 */
function extractPublicationReferences(text: string): string[] {
  if (!text) return [];
  
  // Extract DOIs
  const doiPattern = /10\.\d{4,}\/[-._;()/:a-zA-Z0-9]+/g;
  const dois = text.match(doiPattern) || [];
  
  // Extract URLs
  const urlPattern = /https?:\/\/[^\s]+/g;
  const urls = text.match(urlPattern) || [];
  
  // Remove duplicates and format
  return [...new Set([
    ...dois.map(doi => doi.startsWith('http') ? doi : `https://doi.org/${doi}`),
    ...urls
  ])];
}

/**
 * Get a default tissue type ID for an experiment result
 */
async function getDefaultTissueTypeId(
  convexClient: ConvexClient,
  result: any,
  ctx: MigrationContext
): Promise<Id<"tissueTypes"> | undefined> {
  try {
    // Check if we can extract tissue info from result
    let name = 'Unknown Tissue';
    let species = 'Unknown';
    
    // Try to extract from notes or details
    const notes = result.notes || result.comments || '';
    const details = typeof result.details === 'string' ? 
      JSON.parse(result.details) : (result.details || {});
    
    // Look for species
    const commonSpecies = ['human', 'mouse', 'rat', 'bovine', 'porcine'];
    for (const s of commonSpecies) {
      if (notes.toLowerCase().includes(s)) {
        species = s.charAt(0).toUpperCase() + s.slice(1);
        break;
      }
    }
    
    // Look for cell type
    const commonCellTypes = [
      'hepatocytes', 'oocytes', 'embryos', 'sperm', 
      'stem cells', 'fibroblasts', 'erythrocytes', 'neurons'
    ];
    for (const c of commonCellTypes) {
      if (notes.toLowerCase().includes(c)) {
        const cellType = c.charAt(0).toUpperCase() + c.slice(1);
        name = `${species} ${cellType}`;
        break;
      }
    }
    
    // Create a default tissue type
    const tissueType = {
      name,
      description: `Default tissue type created for experiment result ${result.id}`,
      species,
      createdAt: Date.now(),
      updatedAt: Date.now()
    };
    
    // Insert into Convex
    return await convexClient.mutation(api.tissueTypes.create, tissueType);
  } catch (error) {
    ctx.logMessage('error', `Failed to create default tissue type: ${error}`);
    return undefined;
  }
}

/**
 * Extract method information from result data
 */
function extractMethodInformation(result: any): {
  viabilityMethod?: string;
  functionalAssay?: string;
  integrityMethod?: string;
} {
  const methods: any = {};
  
  // Extract from result data if available
  if (result.viability_method) {
    methods.viabilityMethod = result.viability_method;
  } else {
    // Try to infer from notes or details
    const notesText = (result.notes || result.comments || '').toLowerCase();
    
    if (notesText.includes('trypan blue')) {
      methods.viabilityMethod = 'trypan_blue';
    } else if (notesText.includes('flow cytometry') || notesText.includes('facs')) {
      methods.viabilityMethod = 'flow_cytometry';
    } else if (notesText.includes('fluorescent') || notesText.includes('fluorescence')) {
      methods.viabilityMethod = 'fluorescent_dye';
    } else if (notesText.includes('mtt') || notesText.includes('metabolic')) {
      methods.viabilityMethod = 'metabolic_assay';
    }
  }
  
  // Extract functional assay information
  if (result.functional_assay) {
    methods.functionalAssay = result.functional_assay;
  } else if (result.functionality_method) {
    methods.functionalAssay = result.functionality_method;
  }
  
  // Extract integrity method information
  if (result.integrity_method) {
    methods.integrityMethod = result.integrity_method;
  }
  
  return methods;
}

/**
 * Determine if ice formation was observed from result data
 */
function getIceFormation(result: any): boolean | undefined {
  if (result.ice_formation !== undefined) return !!result.ice_formation;
  
  // Try to infer from notes
  const notes = (result.notes || result.comments || '').toLowerCase();
  if (notes.includes('ice formation observed')) return true;
  if (notes.includes('no ice formation')) return false;
  
  return undefined;
}

/**
 * Extract stress response markers from result data
 */
function getStressMarkers(result: any): string[] | undefined {
  if (result.stress_markers) {
    return Array.isArray(result.stress_markers) ? 
      result.stress_markers : 
      result.stress_markers.split(',').map((m: string) => m.trim());
  }
  
  // Try to infer from notes
  const notes = (result.notes || result.comments || '').toLowerCase();
  const markers = [];
  
  if (notes.includes('hsp') || notes.includes('heat shock')) markers.push('HSP70');
  if (notes.includes('ros') || notes.includes('reactive oxygen')) markers.push('ROS');
  if (notes.includes('apoptosis') || notes.includes('apoptotic')) markers.push('Apoptosis');
  if (notes.includes('dna damage')) markers.push('DNA damage');
  
  return markers.length > 0 ? markers : undefined;
}

/**
 * Run a complete migration from Supabase to Convex
 */
export async function runFullMigration(
  convexClient: ConvexClient,
  supabaseUrl: string,
  supabaseKey: string,
  config: Partial<MigrationConfig> = {}
): Promise<MigrationContext> {
  // Initialize context
  const ctx = initMigrationContext({
    supabaseUrl,
    supabaseKey,
    convexClient,
    ...config
  });
  
  ctx.logMessage('info', 'Starting full migration from Supabase to Convex');
  
  // Migrate core data
  const userIdMapping = await migrateUsers(ctx);
  const tissueTypeIdMapping = await migrateTissueTypes(ctx);
  const protocolIdMapping = await migrateProtocols(ctx);
  
  // Add more migration steps here for:
  // - Experiments
  // - Lab verifications
  // - Biospecimens
  // - Cryoprotectant effectiveness
  
  // Generate migration summary
  const summary = generateMigrationSummary(ctx);
  ctx.logMessage('info', 'Migration complete');
  ctx.logMessage('info', summary);
  
  return ctx;
}

/**
 * Generate a migration summary
 */
function generateMigrationSummary(ctx: MigrationContext): string {
  const stats = ctx.stats;
  
  let summary = '=== Migration Summary ===\n\n';
  let totalRecords = 0;
  let totalSuccess = 0;
  let totalFailures = 0;
  let totalDuration = 0;
  
  for (const tableName in stats) {
    const tableStats = stats[tableName];
    summary += `Table: ${tableName}\n`;
    summary += `  Records: ${tableStats.totalRecords}\n`;
    summary += `  Success: ${tableStats.successfulRecords}\n`;
    summary += `  Failures: ${tableStats.failedRecords}\n`;
    summary += `  Duration: ${tableStats.durationMs || 0}ms\n`;
    
    if (tableStats.errors.length > 0) {
      summary += `  Error samples: ${tableStats.errors.slice(0, 3).map(e => e.error).join(', ')}${tableStats.errors.length > 3 ? '...' : ''}\n`;
    }
    
    summary += '\n';
    
    totalRecords += tableStats.totalRecords;
    totalSuccess += tableStats.successfulRecords;
    totalFailures += tableStats.failedRecords;
    totalDuration += tableStats.durationMs || 0;
  }
  
  summary += 'Overall:\n';
  summary += `  Total records: ${totalRecords}\n`;
  summary += `  Total success: ${totalSuccess}\n`;
  summary += `  Total failures: ${totalFailures}\n`;
  summary += `  Success rate: ${totalRecords > 0 ? ((totalSuccess / totalRecords) * 100).toFixed(2) : 0}%\n`;
  summary += `  Total duration: ${totalDuration}ms (${(totalDuration / 1000 / 60).toFixed(2)} minutes)\n`;
  
  return summary;
}