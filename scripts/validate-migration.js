#!/usr/bin/env node

/**
 * Migration Validation Utility
 * 
 * This script validates data consistency between Supabase and Convex.
 * It compares record counts, samples random records, and verifies relationships.
 */

const { createClient } = require('@supabase/supabase-js');
const { ConvexClient } = require('convex/browser');
const fs = require('fs');
const path = require('path');

// Configuration
const CONFIG = {
  // Supabase configuration
  supabase: {
    url: process.env.SUPABASE_URL || '',
    key: process.env.SUPABASE_KEY || '',
    // Tables to validate and their mapping to Convex tables
    tables: {
      molecules: { convexTable: 'molecules', keyField: 'id' },
      molecular_properties: { convexTable: 'moleculeProperties', keyField: 'id' },
      mixtures: { convexTable: 'mixtures', keyField: 'id' },
      mixture_components: { convexTable: 'mixtureComponents', keyField: 'id' },
      tissue_types: { convexTable: 'tissueTypes', keyField: 'id' },
      protocols: { convexTable: 'protocols', keyField: 'id' },
      protocol_steps: { convexTable: 'protocolSteps', keyField: 'id' },
      experiments: { convexTable: 'experiments', keyField: 'id' },
      experiment_results: { convexTable: 'experimentResults', keyField: 'id' },
      time_series: { convexTable: 'timeSeries', keyField: 'id' },
      time_series_data: { convexTable: 'timeSeriesDataPoints', keyField: 'id' }
    }
  },
  
  // Convex configuration
  convex: {
    url: process.env.CONVEX_URL || ''
  },
  
  // Validation configuration
  validation: {
    sampleSize: 5, // Number of random records to compare
    outputFile: './validation_report.json',
    idMappingsDir: process.env.CHECKPOINT_DIR || './migration_checkpoints'
  }
};

// Initialize clients
let supabase;
let convex;
const idMappings = {};

// Main validation function
async function validateMigration() {
  // Initialize clients
  supabase = createClient(CONFIG.supabase.url, CONFIG.supabase.key);
  console.log(`Initialized Supabase client for ${CONFIG.supabase.url}`);
  
  convex = new ConvexClient(CONFIG.convex.url);
  await convex.waitForInitialization();
  console.log(`Initialized Convex client for ${CONFIG.convex.url}`);
  
  // Load ID mappings
  loadIdMappings();
  
  const report = {
    timestamp: new Date().toISOString(),
    tables: {},
    relationships: {
      valid: true,
      details: {}
    },
    summary: {
      totalTables: 0,
      tablesWithCorrectCounts: 0,
      tablesWithDataDiscrepancies: 0,
      relationshipsValid: true
    }
  };
  
  // Validate each table
  for (const tableName of Object.keys(CONFIG.supabase.tables)) {
    const tableConfig = CONFIG.supabase.tables[tableName];
    console.log(`Validating table ${tableName} -> ${tableConfig.convexTable}...`);
    
    report.tables[tableName] = await validateTable(tableName, tableConfig);
    
    if (report.tables[tableName].countMatch) {
      report.summary.tablesWithCorrectCounts++;
    }
    
    if (report.tables[tableName].dataSample.discrepancies > 0) {
      report.summary.tablesWithDataDiscrepancies++;
    }
    
    report.summary.totalTables++;
  }
  
  // Validate relationships
  report.relationships = await validateRelationships();
  report.summary.relationshipsValid = report.relationships.valid;
  
  // Write report to file
  fs.writeFileSync(CONFIG.validation.outputFile, JSON.stringify(report, null, 2));
  console.log(`Validation report written to ${CONFIG.validation.outputFile}`);
  
  // Summary
  console.log('\nValidation Summary:');
  console.log(`Total tables: ${report.summary.totalTables}`);
  console.log(`Tables with correct counts: ${report.summary.tablesWithCorrectCounts}`);
  console.log(`Tables with data discrepancies: ${report.summary.tablesWithDataDiscrepancies}`);
  console.log(`Relationships valid: ${report.summary.relationshipsValid}`);
  
  return report;
}

// Load ID mappings from files
function loadIdMappings() {
  for (const tableName of Object.keys(CONFIG.supabase.tables)) {
    const mappingFile = path.join(CONFIG.validation.idMappingsDir, `${tableName}_mapping.json`);
    if (fs.existsSync(mappingFile)) {
      try {
        idMappings[tableName] = JSON.parse(fs.readFileSync(mappingFile, 'utf8'));
        console.log(`Loaded ID mappings for ${tableName}: ${Object.keys(idMappings[tableName]).length} records`);
      } catch (error) {
        console.error(`Failed to load ID mappings for ${tableName}:`, error.message);
        idMappings[tableName] = {};
      }
    } else {
      console.warn(`No ID mapping file found for ${tableName}`);
      idMappings[tableName] = {};
    }
  }
}

// Validate a single table
async function validateTable(tableName, tableConfig) {
  const result = {
    supabaseCount: 0,
    convexCount: 0,
    countMatch: false,
    mappedRecords: 0,
    coverage: 0,
    dataSample: {
      sampleSize: CONFIG.validation.sampleSize,
      sampleRecords: [],
      discrepancies: 0
    }
  };
  
  // Get count from Supabase
  const { count: supabaseCount, error: supabaseError } = await supabase
    .from(tableName)
    .select('id', { count: 'exact', head: true });
  
  if (supabaseError) {
    console.error(`Failed to get count from Supabase for ${tableName}:`, supabaseError.message);
    result.error = supabaseError.message;
    return result;
  }
  
  result.supabaseCount = supabaseCount;
  
  // Get count from Convex
  try {
    const convexCount = await convex.query(`${tableConfig.convexTable}:count`);
    result.convexCount = convexCount;
  } catch (error) {
    console.error(`Failed to get count from Convex for ${tableConfig.convexTable}:`, error.message);
    result.error = error.message;
    return result;
  }
  
  // Get count of mapped records
  result.mappedRecords = idMappings[tableName] ? Object.keys(idMappings[tableName]).length : 0;
  result.coverage = result.supabaseCount > 0 ? (result.mappedRecords / result.supabaseCount) * 100 : 0;
  result.countMatch = result.mappedRecords === result.supabaseCount;
  
  // Sample random records for data validation
  if (result.mappedRecords > 0) {
    // Get random sample of mapped IDs
    const mappedIds = Object.keys(idMappings[tableName]);
    const sampleSize = Math.min(CONFIG.validation.sampleSize, mappedIds.length);
    const randomSample = [];
    
    // Select random IDs without replacement
    const tempIds = [...mappedIds];
    for (let i = 0; i < sampleSize; i++) {
      const randomIndex = Math.floor(Math.random() * tempIds.length);
      randomSample.push(tempIds[randomIndex]);
      tempIds.splice(randomIndex, 1);
    }
    
    // Validate each sample record
    for (const supabaseId of randomSample) {
      const convexId = idMappings[tableName][supabaseId];
      
      // Get record from Supabase
      const { data: supabaseRecord, error: supabaseError } = await supabase
        .from(tableName)
        .select('*')
        .eq(tableConfig.keyField, supabaseId)
        .single();
      
      if (supabaseError) {
        console.error(`Failed to get record from Supabase for ${tableName}:`, supabaseError.message);
        continue;
      }
      
      // Get record from Convex
      try {
        const convexRecord = await convex.query(`${tableConfig.convexTable}:getById`, { id: convexId });
        
        // Compare records
        const sampleResult = {
          supabaseId,
          convexId,
          matches: true,
          fields: {}
        };
        
        // Compare important fields based on table type
        if (tableName === 'molecules') {
          compareFields(supabaseRecord, convexRecord, sampleResult, {
            name: 'name',
            inchikey: 'inchikey',
            smiles: 'smiles',
            molecular_formula: 'molecularFormula',
            molecular_weight: 'molecularWeight'
          });
        } else if (tableName === 'molecular_properties') {
          compareFields(supabaseRecord, convexRecord, sampleResult, {
            property_name: 'propertyName',
            property_type: 'propertyType',
            numeric_value: 'numericValue',
            text_value: 'textValue'
          });
        } else if (tableName === 'experiments') {
          compareFields(supabaseRecord, convexRecord, sampleResult, {
            name: 'name',
            status: 'status',
            start_date: 'startDate',
            experiment_type: 'experimentType'
          });
        } else {
          // Generic comparison for other tables
          compareFields(supabaseRecord, convexRecord, sampleResult, {
            name: 'name',
            description: 'description'
          });
        }
        
        result.dataSample.sampleRecords.push(sampleResult);
        
        if (!sampleResult.matches) {
          result.dataSample.discrepancies++;
        }
      } catch (error) {
        console.error(`Failed to get record from Convex for ${tableConfig.convexTable}:`, error.message);
        result.dataSample.sampleRecords.push({
          supabaseId,
          convexId,
          error: error.message
        });
        result.dataSample.discrepancies++;
      }
    }
  }
  
  return result;
}

// Compare fields between Supabase and Convex records
function compareFields(supabaseRecord, convexRecord, sampleResult, fieldMap) {
  for (const [supabaseField, convexField] of Object.entries(fieldMap)) {
    if (supabaseRecord[supabaseField] !== undefined && convexRecord[convexField] !== undefined) {
      const matches = compareValues(supabaseRecord[supabaseField], convexRecord[convexField]);
      sampleResult.fields[supabaseField] = {
        matches,
        supabaseValue: supabaseRecord[supabaseField],
        convexValue: convexRecord[convexField]
      };
      
      if (!matches) {
        sampleResult.matches = false;
      }
    }
  }
}

// Compare values with type checking
function compareValues(supabaseValue, convexValue) {
  // Handle null/undefined
  if (supabaseValue === null || supabaseValue === undefined) {
    return convexValue === null || convexValue === undefined;
  }
  
  // Handle dates
  if (typeof supabaseValue === 'string' && supabaseValue.match(/^\d{4}-\d{2}-\d{2}T/)) {
    // Compare date strings (ignore millisecond precision differences)
    return new Date(supabaseValue).getTime() === new Date(convexValue).getTime();
  }
  
  // Handle numbers
  if (typeof supabaseValue === 'number') {
    return Math.abs(supabaseValue - convexValue) < 0.000001; // Allow for floating point imprecision
  }
  
  // Handle strings, booleans, etc.
  return supabaseValue === convexValue;
}

// Validate relationships between tables
async function validateRelationships() {
  const result = {
    valid: true,
    details: {}
  };
  
  // Check relationships
  // 1. Check molecule_properties -> molecules
  result.details.molecule_properties_to_molecules = await validateRelationship(
    'molecular_properties', 'molecule_id',
    'molecules', 'molecules', 'moleculeId'
  );
  
  // 2. Check mixture_components -> mixtures
  result.details.mixture_components_to_mixtures = await validateRelationship(
    'mixture_components', 'mixture_id',
    'mixtures', 'mixtureComponents', 'mixtureId'
  );
  
  // 3. Check protocol_steps -> protocols
  result.details.protocol_steps_to_protocols = await validateRelationship(
    'protocol_steps', 'protocol_id',
    'protocols', 'protocolSteps', 'protocolId'
  );
  
  // 4. Check experiments -> protocols
  result.details.experiments_to_protocols = await validateRelationship(
    'experiments', 'protocol_id',
    'protocols', 'experiments', 'protocolId'
  );
  
  // 5. Check experiment_results -> experiments
  result.details.experiment_results_to_experiments = await validateRelationship(
    'experiment_results', 'experiment_id',
    'experiments', 'experimentResults', 'experimentId'
  );
  
  // Check if any relationship is invalid
  for (const key of Object.keys(result.details)) {
    if (!result.details[key].valid) {
      result.valid = false;
    }
  }
  
  return result;
}

// Validate a specific relationship
async function validateRelationship(childTable, childForeignKey, parentTable, convexChildTable, convexForeignKey) {
  const result = {
    valid: true,
    brokenRelationships: 0,
    details: []
  };
  
  // Get sample of child records
  const { data: childRecords, error: childError } = await supabase
    .from(childTable)
    .select(`id, ${childForeignKey}`)
    .limit(CONFIG.validation.sampleSize);
  
  if (childError) {
    console.error(`Failed to get records from ${childTable}:`, childError.message);
    result.error = childError.message;
    result.valid = false;
    return result;
  }
  
  // Check each child record
  for (const childRecord of childRecords) {
    if (!childRecord[childForeignKey]) continue;
    
    const parentId = childRecord[childForeignKey];
    const childId = childRecord.id;
    
    // Get convex IDs
    const convexParentId = idMappings[parentTable]?.[parentId];
    const convexChildId = idMappings[childTable]?.[childId];
    
    if (!convexParentId || !convexChildId) {
      result.details.push({
        childId,
        parentId,
        valid: false,
        reason: 'Missing mapping for one or both IDs'
      });
      result.brokenRelationships++;
      result.valid = false;
      continue;
    }
    
    // Get child record from Convex
    try {
      const convexChild = await convex.query(`${convexChildTable}:getById`, { id: convexChildId });
      
      // Check if the relationship is correct
      if (convexChild[convexForeignKey] !== convexParentId) {
        result.details.push({
          childId,
          parentId,
          convexChildId,
          convexParentId,
          valid: false,
          reason: 'Foreign key mismatch in Convex',
          expected: convexParentId,
          actual: convexChild[convexForeignKey]
        });
        result.brokenRelationships++;
        result.valid = false;
      } else {
        result.details.push({
          childId,
          parentId,
          convexChildId,
          convexParentId,
          valid: true
        });
      }
    } catch (error) {
      result.details.push({
        childId,
        parentId,
        valid: false,
        reason: `Failed to get child record from Convex: ${error.message}`
      });
      result.brokenRelationships++;
      result.valid = false;
    }
  }
  
  return result;
}

// Run validation if executed directly
if (require.main === module) {
  validateMigration().catch(error => {
    console.error('Validation failed:', error);
    process.exit(1);
  });
}

// Export for use in other scripts
module.exports = {
  validateMigration
};