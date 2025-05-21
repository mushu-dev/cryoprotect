#!/usr/bin/env node

/**
 * Data Migration Utility: Supabase to Convex
 * 
 * This script migrates data from Supabase to Convex database.
 * It's designed to run in any environment with Node.js installed.
 * 
 * Features:
 * - Table mapping between Supabase and Convex
 * - Data transformation to match schema differences
 * - Resumable migrations with checkpoints
 * - Detailed logging and error handling
 * - Validation of migrated data
 */

const { createClient } = require('@supabase/supabase-js');
const { ConvexClient } = require('convex/browser');
const fs = require('fs');
const path = require('path');
const readline = require('readline');

// Configuration
const CONFIG = {
  // Supabase configuration
  supabase: {
    url: process.env.SUPABASE_URL || '',
    key: process.env.SUPABASE_KEY || '',
    // Tables to migrate and their mapping to Convex tables
    tables: {
      molecules: {
        convexTable: 'molecules',
        keyField: 'id',
        batchSize: 100,
        transform: transformMolecule
      },
      molecular_properties: {
        convexTable: 'moleculeProperties',
        keyField: 'id',
        batchSize: 100,
        transform: transformMoleculeProperty,
        dependencies: ['molecules']
      },
      mixtures: {
        convexTable: 'mixtures',
        keyField: 'id',
        batchSize: 50,
        transform: transformMixture
      },
      mixture_components: {
        convexTable: 'mixtureComponents',
        keyField: 'id',
        batchSize: 100,
        transform: transformMixtureComponent,
        dependencies: ['mixtures', 'molecules']
      },
      tissue_types: {
        convexTable: 'tissueTypes',
        keyField: 'id',
        batchSize: 50,
        transform: transformTissueType
      },
      protocols: {
        convexTable: 'protocols',
        keyField: 'id',
        batchSize: 50,
        transform: transformProtocol
      },
      protocol_steps: {
        convexTable: 'protocolSteps',
        keyField: 'id',
        batchSize: 100,
        transform: transformProtocolStep,
        dependencies: ['protocols']
      },
      experiments: {
        convexTable: 'experiments',
        keyField: 'id',
        batchSize: 50,
        transform: transformExperiment,
        dependencies: ['protocols', 'tissue_types']
      },
      experiment_results: {
        convexTable: 'experimentResults',
        keyField: 'id',
        batchSize: 100,
        transform: transformExperimentResult,
        dependencies: ['experiments', 'molecules', 'mixtures', 'tissue_types', 'protocol_steps']
      },
      time_series: {
        convexTable: 'timeSeries',
        keyField: 'id',
        batchSize: 50,
        transform: transformTimeSeries,
        dependencies: ['experiments', 'experiment_results']
      },
      time_series_data: {
        convexTable: 'timeSeriesDataPoints',
        keyField: 'id',
        batchSize: 200,
        transform: transformTimeSeriesDataPoint,
        dependencies: ['time_series']
      }
    }
  },
  
  // Convex configuration
  convex: {
    url: process.env.CONVEX_URL || ''
  },
  
  // Migration configuration
  migration: {
    checkpointDir: process.env.CHECKPOINT_DIR || './migration_checkpoints',
    logFile: process.env.LOG_FILE || './migration.log',
    validateAfterMigration: true
  }
};

// Initialize clients
let supabase;
let convex;
const idMappings = {};

// Logging
const logger = {
  info: (message) => {
    const logMessage = `[INFO] ${new Date().toISOString()} - ${message}`;
    console.log(logMessage);
    fs.appendFileSync(CONFIG.migration.logFile, logMessage + '\n');
  },
  error: (message, error) => {
    const logMessage = `[ERROR] ${new Date().toISOString()} - ${message}: ${error.message}\n${error.stack}`;
    console.error(logMessage);
    fs.appendFileSync(CONFIG.migration.logFile, logMessage + '\n');
  },
  warn: (message) => {
    const logMessage = `[WARN] ${new Date().toISOString()} - ${message}`;
    console.warn(logMessage);
    fs.appendFileSync(CONFIG.migration.logFile, logMessage + '\n');
  },
  success: (message) => {
    const logMessage = `[SUCCESS] ${new Date().toISOString()} - ${message}`;
    console.log(logMessage);
    fs.appendFileSync(CONFIG.migration.logFile, logMessage + '\n');
  }
};

// Initialization
async function init() {
  // Create checkpoint directory if it doesn't exist
  if (!fs.existsSync(CONFIG.migration.checkpointDir)) {
    fs.mkdirSync(CONFIG.migration.checkpointDir, { recursive: true });
  }
  
  // Initialize log file
  if (!fs.existsSync(path.dirname(CONFIG.migration.logFile))) {
    fs.mkdirSync(path.dirname(CONFIG.migration.logFile), { recursive: true });
  }
  fs.writeFileSync(CONFIG.migration.logFile, `[START] ${new Date().toISOString()} - Migration started\n`);
  
  // Initialize Supabase client
  supabase = createClient(CONFIG.supabase.url, CONFIG.supabase.key);
  logger.info(`Initialized Supabase client for ${CONFIG.supabase.url}`);
  
  // Initialize Convex client
  convex = new ConvexClient(CONFIG.convex.url);
  await convex.waitForInitialization();
  logger.info(`Initialized Convex client for ${CONFIG.convex.url}`);
  
  // Load existing ID mappings
  loadIdMappings();
}

// Load existing ID mappings from checkpoint files
function loadIdMappings() {
  for (const tableName of Object.keys(CONFIG.supabase.tables)) {
    const checkpointFile = path.join(CONFIG.migration.checkpointDir, `${tableName}_mapping.json`);
    if (fs.existsSync(checkpointFile)) {
      try {
        idMappings[tableName] = JSON.parse(fs.readFileSync(checkpointFile, 'utf8'));
        logger.info(`Loaded ID mappings for ${tableName}: ${Object.keys(idMappings[tableName]).length} records`);
      } catch (error) {
        logger.error(`Failed to load ID mappings for ${tableName}`, error);
        idMappings[tableName] = {};
      }
    } else {
      idMappings[tableName] = {};
    }
  }
}

// Save ID mappings to checkpoint files
function saveIdMapping(tableName, supabaseId, convexId) {
  if (!idMappings[tableName]) {
    idMappings[tableName] = {};
  }
  
  idMappings[tableName][supabaseId] = convexId;
  
  // Save to checkpoint file periodically (every 100 records)
  if (Object.keys(idMappings[tableName]).length % 100 === 0) {
    const checkpointFile = path.join(CONFIG.migration.checkpointDir, `${tableName}_mapping.json`);
    fs.writeFileSync(checkpointFile, JSON.stringify(idMappings[tableName], null, 2));
    logger.info(`Checkpoint saved for ${tableName}: ${Object.keys(idMappings[tableName]).length} records`);
  }
}

// Save all ID mappings
function saveAllIdMappings() {
  for (const tableName of Object.keys(idMappings)) {
    const checkpointFile = path.join(CONFIG.migration.checkpointDir, `${tableName}_mapping.json`);
    fs.writeFileSync(checkpointFile, JSON.stringify(idMappings[tableName], null, 2));
    logger.info(`Final checkpoint saved for ${tableName}: ${Object.keys(idMappings[tableName]).length} records`);
  }
}

// Get Convex ID for a Supabase record
function getConvexId(tableName, supabaseId) {
  if (!idMappings[tableName] || !idMappings[tableName][supabaseId]) {
    return null;
  }
  return idMappings[tableName][supabaseId];
}

// Migration functions

// Migrate a single table
async function migrateTable(tableName) {
  const tableConfig = CONFIG.supabase.tables[tableName];
  
  // Check if dependencies are migrated
  if (tableConfig.dependencies) {
    for (const dependency of tableConfig.dependencies) {
      if (!idMappings[dependency] || Object.keys(idMappings[dependency]).length === 0) {
        logger.warn(`Dependency ${dependency} not migrated yet for ${tableName}, skipping...`);
        return false;
      }
    }
  }
  
  logger.info(`Starting migration for table ${tableName} -> ${tableConfig.convexTable}`);
  
  let offset = 0;
  let count = 0;
  let total = 0;
  
  // Get total count
  const { count: totalCount, error: countError } = await supabase
    .from(tableName)
    .select('id', { count: 'exact', head: true });
  
  if (countError) {
    logger.error(`Failed to get count for ${tableName}`, countError);
    return false;
  }
  
  total = totalCount;
  logger.info(`Found ${total} records in ${tableName}`);
  
  // Get already migrated records
  const migratedIds = idMappings[tableName] ? Object.keys(idMappings[tableName]) : [];
  logger.info(`${migratedIds.length} records already migrated for ${tableName}`);
  
  // Calculate remaining records
  const remaining = total - migratedIds.length;
  if (remaining <= 0) {
    logger.info(`All records already migrated for ${tableName}, skipping...`);
    return true;
  }
  
  // Process in batches
  while (true) {
    // Fetch data from Supabase with pagination
    let query = supabase
      .from(tableName)
      .select('*')
      .order(tableConfig.keyField, { ascending: true })
      .range(offset, offset + tableConfig.batchSize - 1);
    
    // Exclude already migrated records
    if (migratedIds.length > 0) {
      query = query.not(tableConfig.keyField, 'in', migratedIds);
    }
    
    const { data, error } = await query;
    
    if (error) {
      logger.error(`Failed to fetch data from ${tableName}`, error);
      return false;
    }
    
    if (!data || data.length === 0) {
      logger.info(`No more data to fetch from ${tableName}`);
      break;
    }
    
    // Transform and insert data into Convex
    for (const record of data) {
      try {
        // Transform data
        const transformedData = tableConfig.transform(record);
        
        // Insert into Convex
        const convexId = await convex.mutation(`${tableConfig.convexTable}:create`, transformedData);
        
        // Save ID mapping
        saveIdMapping(tableName, record[tableConfig.keyField], convexId);
        
        count++;
        
        // Log progress
        if (count % 100 === 0 || count === remaining) {
          logger.info(`Migrated ${count}/${remaining} records from ${tableName}`);
        }
      } catch (error) {
        logger.error(`Failed to migrate record ${record[tableConfig.keyField]} from ${tableName}`, error);
      }
    }
    
    offset += data.length;
    
    // If we got fewer records than the batch size, we're done
    if (data.length < tableConfig.batchSize) {
      break;
    }
  }
  
  logger.success(`Migration completed for ${tableName}: ${count} records migrated`);
  return true;
}

// Transform functions for each table
function transformMolecule(record) {
  return {
    name: record.name || 'Unknown',
    inchikey: record.inchikey || '',
    smiles: record.smiles || '',
    molecularFormula: record.molecular_formula || '',
    molecularWeight: record.molecular_weight || 0,
    isConsolidated: record.is_consolidated || false,
    moleculeStatus: record.molecule_status || 'original',
    primaryMoleculeId: record.primary_molecule_id ? getConvexId('molecules', record.primary_molecule_id) : null,
    primaryMoleculeName: record.primary_molecule_name || null,
    duplicateCount: record.duplicate_count || 0,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformMoleculeProperty(record) {
  return {
    moleculeId: getConvexId('molecules', record.molecule_id),
    propertyTypeId: record.property_type_id,
    propertyName: record.property_name,
    propertyType: record.property_type,
    numericValue: record.numeric_value,
    textValue: record.text_value,
    unit: record.unit,
    source: record.source || 'unknown',
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformMixture(record) {
  return {
    name: record.name,
    description: record.description || '',
    componentCount: record.component_count || 0,
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformMixtureComponent(record) {
  return {
    mixtureId: getConvexId('mixtures', record.mixture_id),
    moleculeId: getConvexId('molecules', record.molecule_id),
    moleculeName: record.molecule_name || 'Unknown',
    concentration: record.concentration || 0,
    concentrationUnit: record.concentration_unit || '%',
    role: record.role,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformTissueType(record) {
  return {
    name: record.name,
    description: record.description,
    species: record.species,
    category: record.category,
    properties: record.properties || {},
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformProtocol(record) {
  return {
    name: record.name,
    description: record.description,
    version: record.version || '1.0',
    parentVersionId: record.parent_version_id ? getConvexId('protocols', record.parent_version_id) : null,
    parameters: record.parameters || {},
    tags: record.tags || [],
    isTemplate: record.is_template || false,
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformProtocolStep(record) {
  return {
    protocolId: getConvexId('protocols', record.protocol_id),
    name: record.name,
    description: record.description,
    order: record.order || 0,
    duration: record.duration,
    durationUnit: record.duration_unit,
    temperature: record.temperature,
    temperatureUnit: record.temperature_unit,
    parameters: record.parameters || {},
    equipment: record.equipment || [],
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformExperiment(record) {
  return {
    name: record.name,
    description: record.description,
    protocolId: getConvexId('protocols', record.protocol_id),
    protocolVersion: record.protocol_version,
    tissueTypeId: getConvexId('tissue_types', record.tissue_type_id),
    experimentType: record.experiment_type,
    startDate: record.start_date,
    endDate: record.end_date,
    status: record.status,
    researcher: record.researcher,
    labId: record.lab_id,
    equipment: record.equipment || [],
    environmentalConditions: record.environmental_conditions || {},
    notes: record.notes,
    tags: record.tags || [],
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformExperimentResult(record) {
  return {
    experimentId: getConvexId('experiments', record.experiment_id),
    tissueTypeId: getConvexId('tissue_types', record.tissue_type_id),
    moleculeId: record.molecule_id ? getConvexId('molecules', record.molecule_id) : null,
    mixtureId: record.mixture_id ? getConvexId('mixtures', record.mixture_id) : null,
    concentration: record.concentration,
    concentrationUnit: record.concentration_unit,
    viabilityPercentage: record.viability_percentage,
    recoveryRate: record.recovery_rate,
    functionalityScore: record.functionality_score,
    resultDetails: record.result_details || {},
    notes: record.notes,
    protocolStepId: record.protocol_step_id ? getConvexId('protocol_steps', record.protocol_step_id) : null,
    timestamp: record.timestamp,
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformTimeSeries(record) {
  return {
    experimentId: getConvexId('experiments', record.experiment_id),
    resultId: record.result_id ? getConvexId('experiment_results', record.result_id) : null,
    parameter: record.parameter,
    unit: record.unit,
    startTime: record.start_time,
    endTime: record.end_time,
    notes: record.notes,
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

function transformTimeSeriesDataPoint(record) {
  return {
    timeSeriesId: getConvexId('time_series', record.time_series_id),
    time: record.time,
    value: record.value,
    uncertainty: record.uncertainty,
    metadata: record.metadata || {},
    createdBy: record.created_by ? getConvexId('users', record.created_by) : null,
    createdAt: new Date(record.created_at).getTime(),
    updatedAt: new Date(record.updated_at).getTime()
  };
}

// Validate migration
async function validateMigration() {
  logger.info('Starting validation of migrated data...');
  
  for (const tableName of Object.keys(CONFIG.supabase.tables)) {
    const tableConfig = CONFIG.supabase.tables[tableName];
    
    // Get count from Supabase
    const { count: supabaseCount, error: supabaseError } = await supabase
      .from(tableName)
      .select('id', { count: 'exact', head: true });
    
    if (supabaseError) {
      logger.error(`Failed to get count from Supabase for ${tableName}`, supabaseError);
      continue;
    }
    
    // Get count from Convex
    const convexCount = await convex.query(`${tableConfig.convexTable}:count`);
    
    // Get count of migrated records
    const migratedCount = idMappings[tableName] ? Object.keys(idMappings[tableName]).length : 0;
    
    logger.info(`Table ${tableName} -> ${tableConfig.convexTable}:`);
    logger.info(`  Supabase count: ${supabaseCount}`);
    logger.info(`  Convex count: ${convexCount}`);
    logger.info(`  Mapped records: ${migratedCount}`);
    
    if (migratedCount < supabaseCount) {
      logger.warn(`  Not all records migrated: ${migratedCount}/${supabaseCount} (${Math.round(migratedCount / supabaseCount * 100)}%)`);
    } else {
      logger.success(`  All records migrated successfully!`);
    }
  }
}

// Main migration function
async function migrate() {
  try {
    await init();
    
    logger.info('Starting migration process...');
    
    // Define migration order based on dependencies
    const migrationOrder = [
      'tissue_types',
      'protocols',
      'molecules',
      'molecular_properties',
      'mixtures',
      'mixture_components',
      'protocol_steps',
      'experiments',
      'experiment_results',
      'time_series',
      'time_series_data'
    ];
    
    // Migrate tables in order
    for (const tableName of migrationOrder) {
      await migrateTable(tableName);
    }
    
    // Save all ID mappings
    saveAllIdMappings();
    
    // Validate migration
    if (CONFIG.migration.validateAfterMigration) {
      await validateMigration();
    }
    
    logger.success('Migration completed successfully!');
  } catch (error) {
    logger.error('Migration failed', error);
  }
}

// Interactive CLI
async function cli() {
  const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
  });
  
  console.log('Supabase to Convex Migration Utility');
  console.log('===================================\n');
  
  if (!CONFIG.supabase.url) {
    CONFIG.supabase.url = await question(rl, 'Enter Supabase URL: ');
  }
  
  if (!CONFIG.supabase.key) {
    CONFIG.supabase.key = await question(rl, 'Enter Supabase service_role key: ');
  }
  
  if (!CONFIG.convex.url) {
    CONFIG.convex.url = await question(rl, 'Enter Convex deployment URL: ');
  }
  
  console.log('\nMigration Configuration:');
  console.log(`- Supabase URL: ${CONFIG.supabase.url}`);
  console.log(`- Convex URL: ${CONFIG.convex.url}`);
  console.log(`- Checkpoint directory: ${CONFIG.migration.checkpointDir}`);
  console.log(`- Log file: ${CONFIG.migration.logFile}`);
  console.log(`- Tables to migrate: ${Object.keys(CONFIG.supabase.tables).join(', ')}`);
  
  const confirm = await question(rl, '\nProceed with migration? (yes/no): ');
  
  if (confirm.toLowerCase() === 'yes' || confirm.toLowerCase() === 'y') {
    await migrate();
  } else {
    console.log('Migration aborted.');
  }
  
  rl.close();
}

// Helper function for CLI
function question(rl, query) {
  return new Promise(resolve => {
    rl.question(query, answer => {
      resolve(answer);
    });
  });
}

// Run CLI if executed directly
if (require.main === module) {
  cli();
}

// Export functions for use in other scripts
module.exports = {
  migrate,
  validateMigration,
  migrateTable
};