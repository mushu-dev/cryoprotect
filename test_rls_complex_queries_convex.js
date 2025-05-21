/**
 * Test RLS Complex Queries in Convex
 * 
 * This script is used to test the performance of the optimized RLS implementation
 * by comparing query performance before and after optimization.
 */

const { ConvexClient } = require('convex/browser');
const fs = require('fs');
const path = require('path');

// Configuration
const CONVEX_URL = process.env.CONVEX_URL || 'http://localhost:3000';
const AUTH_TOKEN = process.env.AUTH_TOKEN; // JWT token for authentication
const ITERATIONS = parseInt(process.env.ITERATIONS || '3'); // Default to 3 iterations
const OUTPUT_FILE = process.env.OUTPUT_FILE || `convex_rls_test_${Date.now()}.json`;

// Test queries
const TEST_QUERIES = {
  // Original queries
  standard: {
    // Standard molecule query - retrieves a molecule by ID
    singleMoleculeAccess: {
      args: {
        id: '__MOLECULE_ID__' // Replace with actual molecule ID
      },
      path: 'molecules:getMolecule',
      description: 'Standard query to access a single molecule'
    },
    
    // Standard list query - retrieves a list of molecules
    moleculesList: {
      args: {
        filter: {},
        options: { limit: 10 }
      },
      path: 'molecules:searchMolecules',
      description: 'Standard query to list molecules with basic filtering'
    },
    
    // Complex query - molecules with property filtering
    moleculesWithPropertyFilter: {
      args: {
        propertyTypeId: '__PROPERTY_TYPE_ID__', // Replace with actual property type ID
        minValue: 100,
        maxValue: 500,
        limit: 10
      },
      path: 'molecules:searchMoleculesByPropertyRange',
      description: 'Complex query to find molecules by property range'
    },
    
    // Complex query - mixtures with component filtering
    mixturesWithComponentCount: {
      args: {
        minComponents: 2,
        limit: 10
      },
      path: 'mixtures:searchMixturesByComponentCount',
      description: 'Complex query to find mixtures with a minimum number of components'
    }
  },
  
  // Optimized queries with the new access control system
  optimized: {
    // Optimized molecule query with cached access
    singleMoleculeAccessOptimized: {
      args: {
        id: '__MOLECULE_ID__' // Replace with actual molecule ID
      },
      path: 'molecules:getMoleculeOptimized',
      description: 'Optimized query to access a single molecule with cached access control'
    },
    
    // Optimized list query with batch access filtering
    moleculesListOptimized: {
      args: {
        filter: {},
        options: { limit: 10 }
      },
      path: 'molecules:searchMoleculesOptimized',
      description: 'Optimized query to list molecules with batch access control'
    },
    
    // Optimized property filter query
    moleculesWithPropertyFilterOptimized: {
      args: {
        propertyTypeId: '__PROPERTY_TYPE_ID__', // Replace with actual property type ID
        minValue: 100,
        maxValue: 500,
        limit: 10
      },
      path: 'molecules:searchMoleculesByPropertyRangeOptimized',
      description: 'Optimized complex query for property range with batch access control'
    },
    
    // Optimized mixture component query
    mixturesWithComponentCountOptimized: {
      args: {
        minComponents: 2,
        limit: 10
      },
      path: 'mixtures:searchMixturesByComponentCountOptimized',
      description: 'Optimized complex query for mixtures with batch access control'
    }
  }
};

async function getTestIds(client) {
  // Get test IDs (molecule ID, property type ID, etc.)
  console.log('Getting test IDs...');
  
  try {
    // Get a molecule ID
    const molecules = await client.query('molecules:searchMolecules', { 
      filter: {}, 
      options: { limit: 1 } 
    });
    
    if (!molecules || !molecules.molecules || !molecules.molecules.length) {
      throw new Error('No molecules found in the database');
    }
    
    const moleculeId = molecules.molecules[0].id;
    
    // Get a property type ID
    const propertyTypes = await client.query('properties:listPropertyTypes', {
      limit: 1
    });
    
    if (!propertyTypes || !propertyTypes.length) {
      throw new Error('No property types found in the database');
    }
    
    const propertyTypeId = propertyTypes[0].id;
    
    return {
      moleculeId,
      propertyTypeId
    };
  } catch (error) {
    console.error('Error getting test IDs:', error);
    throw error;
  }
}

async function testQuery(client, query, testIds) {
  const result = {
    name: query.path,
    description: query.description,
    iterations: [],
    averageTime: null,
    success: false,
    error: null
  };
  
  try {
    // Replace placeholder IDs with actual IDs
    const args = JSON.parse(JSON.stringify(query.args));
    for (const key in args) {
      if (typeof args[key] === 'string') {
        if (args[key] === '__MOLECULE_ID__') {
          args[key] = testIds.moleculeId;
        } else if (args[key] === '__PROPERTY_TYPE_ID__') {
          args[key] = testIds.propertyTypeId;
        }
      }
    }
    
    // Run the query multiple times
    for (let i = 0; i < ITERATIONS; i++) {
      const startTime = performance.now();
      await client.query(query.path, args);
      const endTime = performance.now();
      result.iterations.push(endTime - startTime);
    }
    
    // Calculate average time
    result.averageTime = result.iterations.reduce((a, b) => a + b, 0) / result.iterations.length;
    result.success = true;
  } catch (error) {
    result.error = error.message;
  }
  
  return result;
}

async function runTests() {
  console.log('Starting RLS query performance tests...');
  console.log(`Running ${ITERATIONS} iterations for each query`);
  
  const results = {
    timestamp: new Date().toISOString(),
    standard: {},
    optimized: {},
    improvements: {},
    summary: {
      averageImprovement: null,
      averageSpeedupFactor: null
    }
  };
  
  try {
    // Initialize Convex client
    const client = new ConvexClient(CONVEX_URL);
    
    // Authenticate if token is provided
    if (AUTH_TOKEN) {
      await client.setAuth(AUTH_TOKEN);
      console.log('Authenticated with provided token');
    } else {
      console.warn('No auth token provided. Some queries may fail if authentication is required.');
    }
    
    // Get test IDs
    const testIds = await getTestIds(client);
    console.log(`Test IDs: moleculeId=${testIds.moleculeId}, propertyTypeId=${testIds.propertyTypeId}`);
    
    // Test standard queries
    console.log('\nTesting standard queries...');
    for (const [name, query] of Object.entries(TEST_QUERIES.standard)) {
      console.log(`Testing ${name}...`);
      results.standard[name] = await testQuery(client, query, testIds);
      if (results.standard[name].success) {
        console.log(`  ✓ ${name}: ${results.standard[name].averageTime.toFixed(2)}ms average`);
      } else {
        console.log(`  ✗ ${name}: Failed - ${results.standard[name].error}`);
      }
    }
    
    // Test optimized queries
    console.log('\nTesting optimized queries...');
    for (const [name, query] of Object.entries(TEST_QUERIES.optimized)) {
      console.log(`Testing ${name}...`);
      results.optimized[name] = await testQuery(client, query, testIds);
      if (results.optimized[name].success) {
        console.log(`  ✓ ${name}: ${results.optimized[name].averageTime.toFixed(2)}ms average`);
      } else {
        console.log(`  ✗ ${name}: Failed - ${results.optimized[name].error}`);
      }
    }
    
    // Calculate improvements
    console.log('\nCalculating improvements...');
    const improvements = [];
    const speedupFactors = [];
    
    for (const [stdName, stdResult] of Object.entries(results.standard)) {
      const optName = `${stdName}Optimized`;
      if (results.optimized[optName] && stdResult.success && results.optimized[optName].success) {
        const standardTime = stdResult.averageTime;
        const optimizedTime = results.optimized[optName].averageTime;
        const improvement = ((standardTime - optimizedTime) / standardTime) * 100;
        const speedupFactor = standardTime / optimizedTime;
        
        results.improvements[stdName] = {
          standardTime,
          optimizedTime,
          improvement,
          speedupFactor
        };
        
        improvements.push(improvement);
        speedupFactors.push(speedupFactor);
        
        console.log(`  ${stdName}: ${improvement.toFixed(2)}% improvement (${speedupFactor.toFixed(2)}x speedup)`);
      }
    }
    
    // Calculate summary
    if (improvements.length > 0) {
      results.summary.averageImprovement = improvements.reduce((a, b) => a + b, 0) / improvements.length;
      results.summary.averageSpeedupFactor = speedupFactors.reduce((a, b) => a + b, 0) / speedupFactors.length;
      
      console.log(`\nAverage improvement: ${results.summary.averageImprovement.toFixed(2)}%`);
      console.log(`Average speedup factor: ${results.summary.averageSpeedupFactor.toFixed(2)}x`);
    } else {
      console.log('\nNo successful comparison pairs found');
    }
    
    // Save results to file
    fs.writeFileSync(path.resolve(OUTPUT_FILE), JSON.stringify(results, null, 2));
    console.log(`\nResults saved to ${OUTPUT_FILE}`);
    
    return results;
  } catch (error) {
    console.error('Error running tests:', error);
    fs.writeFileSync(path.resolve(OUTPUT_FILE), JSON.stringify({ 
      error: error.message,
      stack: error.stack
    }, null, 2));
    throw error;
  }
}

// Run the tests
runTests().catch(error => {
  console.error('Test execution failed:', error);
  process.exit(1);
});