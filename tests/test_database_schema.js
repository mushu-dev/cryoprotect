/**
 * CryoProtect Analyzer - Database Schema Test
 * 
 * This script tests the database schema to ensure it's set up correctly.
 * It performs basic CRUD operations on all tables and verifies relationships.
 * 
 * Prerequisites:
 * - Node.js installed
 * - Supabase project with the CryoProtect schema applied
 * - @supabase/supabase-js package installed (npm install @supabase/supabase-js)
 */

const { createClient } = require('@supabase/supabase-js');

/**
 * RLS/Access Control Test Setup:
 *
 * To test RLS, set the following environment variables:
 * - SUPABASE_URL: Your Supabase project URL
 * - SUPABASE_KEY_AUTHORIZED: Service role key or JWT for a project member (authorized)
 * - SUPABASE_KEY_UNAUTHORIZED: Anon/public key or JWT for a non-member (unauthorized)
 */
const supabaseUrl = process.env.SUPABASE_URL || 'https://your-project-ref.supabase.co';
const supabaseKeyAuthorized = process.env.SUPABASE_KEY_AUTHORIZED || 'your-authorized-key';
const supabaseKeyUnauthorized = process.env.SUPABASE_KEY_UNAUTHORIZED || 'your-unauthorized-key';
const supabaseAuthorized = createClient(supabaseUrl, supabaseKeyAuthorized);
const supabaseUnauthorized = createClient(supabaseUrl, supabaseKeyUnauthorized);

// For legacy tests, keep the default client as authorized
const supabase = supabaseAuthorized;

// Test results tracking
const testResults = {
  passed: 0,
  failed: 0,
  tests: []
};

// Helper function to run a test
async function runTest(name, testFn) {
  console.log(`\n🧪 Running test: ${name}`);
  try {
    await testFn();
    console.log(`✅ Test passed: ${name}`);
    testResults.passed++;
    testResults.tests.push({ name, status: 'passed' });
  } catch (error) {
    console.error(`❌ Test failed: ${name}`);
    console.error(`   Error: ${error.message}`);
    testResults.failed++;
    testResults.tests.push({ name, status: 'failed', error: error.message });
  }
}

// Test 1: Verify tables exist
async function testTablesExist() {
  const tables = [
    'molecules',
    'property_types',
    'molecular_properties',
    'mixtures',
    'mixture_components',
    'calculation_methods',
    'predictions',
    'experiments'
  ];
  
  for (const table of tables) {
    const { data, error } = await supabase
      .from(table)
      .select('*')
      .limit(1);
    
    if (error && !error.message.includes('permission denied')) {
      throw new Error(`Table '${table}' does not exist: ${error.message}`);
    }
  }
}

// Test 2: Verify property types were created
async function testPropertyTypesExist() {
  const { data, error } = await supabase
    .from('property_types')
    .select('name')
    .in('name', ['Molecular Weight', 'LogP', 'TPSA']);
  
  if (error) {
    throw new Error(`Error fetching property types: ${error.message}`);
  }
  
  if (!data || data.length < 3) {
    throw new Error(`Expected property types not found. Found: ${JSON.stringify(data)}`);
  }
}

// Test 3: Verify calculation methods were created
async function testCalculationMethodsExist() {
  const { data, error } = await supabase
    .from('calculation_methods')
    .select('name')
    .in('name', ['PubChem Properties', 'CryoProtect Scoring']);
  
  if (error) {
    throw new Error(`Error fetching calculation methods: ${error.message}`);
  }
  
  if (!data || data.length < 2) {
    throw new Error(`Expected calculation methods not found. Found: ${JSON.stringify(data)}`);
  }
}

// Test 4: Test molecule creation
async function testMoleculeCreation() {
  // Create a test molecule
  const { data: molecule, error: createError } = await supabase
    .from('molecules')
    .insert({
      cid: 999999, // Use a CID that's unlikely to exist
      name: 'Test Molecule',
      molecular_formula: 'C3H8O3',
      smiles: 'C(C(CO)O)O'
    })
    .select()
    .single();
  
  if (createError) {
    throw new Error(`Error creating molecule: ${createError.message}`);
  }
  
  // Verify molecule was created
  const { data: fetchedMolecule, error: fetchError } = await supabase
    .from('molecules')
    .select('*')
    .eq('id', molecule.id)
    .single();
  
  if (fetchError) {
    throw new Error(`Error fetching created molecule: ${fetchError.message}`);
  }
  
  if (!fetchedMolecule || fetchedMolecule.cid !== 999999) {
    throw new Error(`Created molecule not found or incorrect`);
  }
  
  // Clean up - delete the test molecule
  const { error: deleteError } = await supabase
    .from('molecules')
    .delete()
    .eq('id', molecule.id);
  
  if (deleteError) {
    console.warn(`Warning: Could not delete test molecule: ${deleteError.message}`);
  }
  
  return molecule.id; // Return the ID for use in other tests
}

// Test 5: Test molecular properties
async function testMolecularProperties() {
  // First create a test molecule
  const { data: molecule, error: createMoleculeError } = await supabase
    .from('molecules')
    .insert({
      cid: 888888,
      name: 'Test Molecule for Properties',
      molecular_formula: 'C3H8O3',
      smiles: 'C(C(CO)O)O'
    })
    .select()
    .single();
  
  if (createMoleculeError) {
    throw new Error(`Error creating test molecule: ${createMoleculeError.message}`);
  }
  
  // Get property type IDs
  const { data: propertyTypes, error: propertyTypesError } = await supabase
    .from('property_types')
    .select('id, name, data_type')
    .in('name', ['Molecular Weight', 'LogP']);
  
  if (propertyTypesError) {
    throw new Error(`Error fetching property types: ${propertyTypesError.message}`);
  }
  
  if (!propertyTypes || propertyTypes.length < 2) {
    throw new Error(`Required property types not found`);
  }
  
  // Add properties to the molecule
  const mwProperty = propertyTypes.find(pt => pt.name === 'Molecular Weight');
  const logPProperty = propertyTypes.find(pt => pt.name === 'LogP');
  
  const propertyInserts = [
    {
      molecule_id: molecule.id,
      property_type_id: mwProperty.id,
      numeric_value: 92.09
    },
    {
      molecule_id: molecule.id,
      property_type_id: logPProperty.id,
      numeric_value: -1.76
    }
  ];
  
  const { error: insertPropertiesError } = await supabase
    .from('molecular_properties')
    .insert(propertyInserts);
  
  if (insertPropertiesError) {
    throw new Error(`Error inserting molecular properties: ${insertPropertiesError.message}`);
  }
  
  // Verify properties were added
  const { data: fetchedProperties, error: fetchPropertiesError } = await supabase
    .from('molecular_properties')
    .select('property_type_id, numeric_value')
    .eq('molecule_id', molecule.id);
  
  if (fetchPropertiesError) {
    throw new Error(`Error fetching molecular properties: ${fetchPropertiesError.message}`);
  }
  
  if (!fetchedProperties || fetchedProperties.length !== 2) {
    throw new Error(`Expected 2 properties, found ${fetchedProperties?.length || 0}`);
  }
  
  // Clean up - delete the properties and molecule
  await supabase
    .from('molecular_properties')
    .delete()
    .eq('molecule_id', molecule.id);
  
  await supabase
    .from('molecules')
    .delete()
    .eq('id', molecule.id);
}

// Test 6: Test mixture and components
async function testMixtureAndComponents() {
  // Create two test molecules
  const { data: molecules, error: createMoleculesError } = await supabase
    .from('molecules')
    .insert([
      {
        cid: 777771,
        name: 'Test Molecule 1',
        molecular_formula: 'C3H8O3',
        smiles: 'C(C(CO)O)O'
      },
      {
        cid: 777772,
        name: 'Test Molecule 2',
        molecular_formula: 'H2O',
        smiles: 'O'
      }
    ])
    .select();
  
  if (createMoleculesError) {
    throw new Error(`Error creating test molecules: ${createMoleculesError.message}`);
  }
  
  // Create a mixture
  const { data: mixture, error: createMixtureError } = await supabase
    .from('mixtures')
    .insert({
      name: 'Test Mixture',
      description: 'A test mixture for schema validation'
    })
    .select()
    .single();
  
  if (createMixtureError) {
    throw new Error(`Error creating mixture: ${createMixtureError.message}`);
  }
  
  // Add components to the mixture
  const componentInserts = [
    {
      mixture_id: mixture.id,
      molecule_id: molecules[0].id,
      concentration: 30,
      concentration_unit: '%'
    },
    {
      mixture_id: mixture.id,
      molecule_id: molecules[1].id,
      concentration: 70,
      concentration_unit: '%'
    }
  ];
  
  const { error: insertComponentsError } = await supabase
    .from('mixture_components')
    .insert(componentInserts);
  
  if (insertComponentsError) {
    throw new Error(`Error inserting mixture components: ${insertComponentsError.message}`);
  }
  
  // Verify components were added
  const { data: fetchedComponents, error: fetchComponentsError } = await supabase
    .from('mixture_components')
    .select('molecule_id, concentration, concentration_unit')
    .eq('mixture_id', mixture.id);
  
  if (fetchComponentsError) {
    throw new Error(`Error fetching mixture components: ${fetchComponentsError.message}`);
  }
  
  if (!fetchedComponents || fetchedComponents.length !== 2) {
    throw new Error(`Expected 2 components, found ${fetchedComponents?.length || 0}`);
  }
  
  // Test the mixture_with_components view
  const { data: mixtureView, error: mixtureViewError } = await supabase
    .from('mixture_with_components')
    .select('*')
    .eq('id', mixture.id)
    .single();
  
  if (mixtureViewError) {
    throw new Error(`Error fetching from mixture_with_components view: ${mixtureViewError.message}`);
  }
  
  if (!mixtureView || !mixtureView.components || mixtureView.components.length !== 2) {
    throw new Error(`mixture_with_components view not working correctly`);
  }
  
  // Clean up - delete the components, mixture, and molecules
  await supabase
    .from('mixture_components')
    .delete()
    .eq('mixture_id', mixture.id);
  
  await supabase
    .from('mixtures')
    .delete()
    .eq('id', mixture.id);
  
  await supabase
    .from('molecules')
    .delete()
    .in('id', [molecules[0].id, molecules[1].id]);
}

// Test 7: Test predictions and experiments
async function testPredictionsAndExperiments() {
  // Create a test molecule
  const { data: molecule, error: createMoleculeError } = await supabase
    .from('molecules')
    .insert({
      cid: 666666,
      name: 'Test Molecule for Predictions',
      molecular_formula: 'C3H8O3',
      smiles: 'C(C(CO)O)O'
    })
    .select()
    .single();
  
  if (createMoleculeError) {
    throw new Error(`Error creating test molecule: ${createMoleculeError.message}`);
  }
  
  // Create a mixture
  const { data: mixture, error: createMixtureError } = await supabase
    .from('mixtures')
    .insert({
      name: 'Test Mixture for Predictions',
      description: 'A test mixture for predictions and experiments'
    })
    .select()
    .single();
  
  if (createMixtureError) {
    throw new Error(`Error creating mixture: ${createMixtureError.message}`);
  }
  
  // Add molecule to the mixture
  const { error: insertComponentError } = await supabase
    .from('mixture_components')
    .insert({
      mixture_id: mixture.id,
      molecule_id: molecule.id,
      concentration: 100,
      concentration_unit: '%'
    });
  
  if (insertComponentError) {
    throw new Error(`Error inserting mixture component: ${insertComponentError.message}`);
  }
  
  // Get property type and calculation method
  const { data: propertyType, error: propertyTypeError } = await supabase
    .from('property_types')
    .select('id')
    .eq('name', 'Freezing Point')
    .maybeSingle();
  
  // If Freezing Point doesn't exist, create it
  let propertyTypeId;
  if (!propertyType) {
    const { data: newPropertyType, error: createPropertyTypeError } = await supabase
      .from('property_types')
      .insert({
        name: 'Freezing Point',
        description: 'Temperature at which the mixture freezes',
        units: '°C',
        data_type: 'numeric'
      })
      .select()
      .single();
    
    if (createPropertyTypeError) {
      throw new Error(`Error creating property type: ${createPropertyTypeError.message}`);
    }
    
    propertyTypeId = newPropertyType.id;
  } else {
    propertyTypeId = propertyType.id;
  }
  
  const { data: calculationMethod, error: calculationMethodError } = await supabase
    .from('calculation_methods')
    .select('id')
    .eq('name', 'CryoProtect Scoring')
    .single();
  
  if (calculationMethodError) {
    throw new Error(`Error fetching calculation method: ${calculationMethodError.message}`);
  }
  
  // Add a prediction
  const { error: insertPredictionError } = await supabase
    .from('predictions')
    .insert({
      mixture_id: mixture.id,
      property_type_id: propertyTypeId,
      calculation_method_id: calculationMethod.id,
      numeric_value: -15.3,
      confidence: 0.9
    });
  
  if (insertPredictionError) {
    throw new Error(`Error inserting prediction: ${insertPredictionError.message}`);
  }
  
  // Add an experiment
  const { error: insertExperimentError } = await supabase
    .from('experiments')
    .insert({
      mixture_id: mixture.id,
      property_type_id: propertyTypeId,
      numeric_value: -14.8,
      experimental_conditions: 'Standard pressure, cooling rate 1°C/min',
      date_performed: new Date().toISOString().split('T')[0]
    });
  
  if (insertExperimentError) {
    throw new Error(`Error inserting experiment: ${insertExperimentError.message}`);
  }
  
  // Test the compare_prediction_with_experiment function
  const { data: comparison, error: comparisonError } = await supabase.rpc(
    'compare_prediction_with_experiment',
    {
      p_mixture_id: mixture.id,
      p_property_type_id: propertyTypeId
    }
  );
  
  if (comparisonError) {
    throw new Error(`Error comparing prediction with experiment: ${comparisonError.message}`);
  }
  
  if (!comparison || !comparison.prediction || !comparison.experiment) {
    throw new Error(`compare_prediction_with_experiment function not working correctly`);
  }
  
  // Clean up - delete everything
  await supabase
    .from('experiments')
    .delete()
    .eq('mixture_id', mixture.id);
  
  await supabase
    .from('predictions')
    .delete()
    .eq('mixture_id', mixture.id);
  
  await supabase
    .from('mixture_components')
    .delete()
    .eq('mixture_id', mixture.id);
  
  await supabase
    .from('mixtures')
    .delete()
    .eq('id', mixture.id);
  
  await supabase
    .from('molecules')
    .delete()
    .eq('id', molecule.id);
}

// Test 8: Test views
async function testViews() {
  // Create a test molecule
  const { data: molecule, error: createMoleculeError } = await supabase
    .from('molecules')
    .insert({
      cid: 555555,
      name: 'Test Molecule for Views',
      molecular_formula: 'C3H8O3',
      smiles: 'C(C(CO)O)O'
    })
    .select()
    .single();
  
  if (createMoleculeError) {
    throw new Error(`Error creating test molecule: ${createMoleculeError.message}`);
  }
  
  // Get property type IDs
  const { data: propertyTypes, error: propertyTypesError } = await supabase
    .from('property_types')
    .select('id, name, data_type')
    .in('name', ['Molecular Weight', 'LogP']);
  
  if (propertyTypesError) {
    throw new Error(`Error fetching property types: ${propertyTypesError.message}`);
  }
  
  // Add properties to the molecule
  const mwProperty = propertyTypes.find(pt => pt.name === 'Molecular Weight');
  const logPProperty = propertyTypes.find(pt => pt.name === 'LogP');
  
  const propertyInserts = [
    {
      molecule_id: molecule.id,
      property_type_id: mwProperty.id,
      numeric_value: 92.09
    },
    {
      molecule_id: molecule.id,
      property_type_id: logPProperty.id,
      numeric_value: -1.76
    }
  ];
  
  const { error: insertPropertiesError } = await supabase
    .from('molecular_properties')
    .insert(propertyInserts);
  
  if (insertPropertiesError) {
    throw new Error(`Error inserting molecular properties: ${insertPropertiesError.message}`);
  }
  
  // Test the molecule_with_properties view
  const { data: moleculeView, error: moleculeViewError } = await supabase
    .from('molecule_with_properties')
    .select('*')
    .eq('id', molecule.id)
    .single();
  
  if (moleculeViewError) {
    throw new Error(`Error fetching from molecule_with_properties view: ${moleculeViewError.message}`);
  }
  
  if (!moleculeView || !moleculeView.properties) {
    throw new Error(`molecule_with_properties view not working correctly`);
  }
  
  // Verify the properties are in the view
  if (!moleculeView.properties['Molecular Weight'] || !moleculeView.properties['LogP']) {
    throw new Error(`Properties not found in molecule_with_properties view`);
  }
  
  // Clean up
  await supabase
    .from('molecular_properties')
    .delete()
    .eq('molecule_id', molecule.id);
  
  await supabase
    .from('molecules')
    .delete()
    .eq('id', molecule.id);
}

/**
 * RLS and Access Control Tests
 * For each RLS-protected table, attempt SELECT, INSERT, UPDATE, DELETE as both authorized and unauthorized users.
 */
async function testRLSAccessControl() {
  const tables = [
    {
      name: 'molecular_properties',
      sampleRow: { molecule_id: '00000000-0000-0000-0000-000000000001', property_type_id: 1, numeric_value: 1.23 },
      pk: ['molecule_id', 'property_type_id']
    },
    {
      name: 'predictions',
      sampleRow: { mixture_id: '00000000-0000-0000-0000-000000000002', property_type_id: 1, calculation_method_id: 1, numeric_value: 2.34, confidence: 0.9 },
      pk: ['mixture_id', 'property_type_id', 'calculation_method_id']
    },
    {
      name: 'experiment_property',
      sampleRow: { experiment_id: '00000000-0000-0000-0000-000000000003', property_type: 'Test', value: 42, unit: 'unit', method_id: 1, provenance: 'test', created_by: 'test', data_source: 'test', version: 1, modification_history: '[]', created_at: new Date().toISOString(), updated_at: new Date().toISOString() },
      pk: ['experiment_id', 'property_type']
    }
  ];

  for (const { name, sampleRow, pk } of tables) {
    // INSERT
    await runTest(`[RLS] INSERT as authorized into ${name}`, async () => {
      const { error } = await supabaseAuthorized.from(name).insert(sampleRow);
      if (error) throw new Error(`Authorized insert failed: ${error.message}`);
    });
    await runTest(`[RLS] INSERT as unauthorized into ${name}`, async () => {
      const { error } = await supabaseUnauthorized.from(name).insert(sampleRow);
      if (!error || !/permission denied/i.test(error.message)) {
        throw new Error(`Unauthorized insert should be denied, got: ${error ? error.message : 'no error'}`);
      }
    });

    // SELECT
    await runTest(`[RLS] SELECT as authorized from ${name}`, async () => {
      const { error } = await supabaseAuthorized.from(name).select('*').limit(1);
      if (error) throw new Error(`Authorized select failed: ${error.message}`);
    });
    await runTest(`[RLS] SELECT as unauthorized from ${name}`, async () => {
      const { error } = await supabaseUnauthorized.from(name).select('*').limit(1);
      if (!error || !/permission denied/i.test(error.message)) {
        throw new Error(`Unauthorized select should be denied, got: ${error ? error.message : 'no error'}`);
      }
    });

    // UPDATE
    await runTest(`[RLS] UPDATE as authorized on ${name}`, async () => {
      const { error } = await supabaseAuthorized.from(name).update({ test_update: true }).match(sampleRow);
      // If no rows, error will be null, so we don't fail on 0 updated
      if (error && !/no rows/i.test(error.message)) throw new Error(`Authorized update failed: ${error.message}`);
    });
    await runTest(`[RLS] UPDATE as unauthorized on ${name}`, async () => {
      const { error } = await supabaseUnauthorized.from(name).update({ test_update: true }).match(sampleRow);
      if (!error || !/permission denied/i.test(error.message)) {
        throw new Error(`Unauthorized update should be denied, got: ${error ? error.message : 'no error'}`);
      }
    });

    // DELETE
    await runTest(`[RLS] DELETE as authorized from ${name}`, async () => {
      const { error } = await supabaseAuthorized.from(name).delete().match(sampleRow);
      // If no rows, error will be null, so we don't fail on 0 deleted
      if (error && !/no rows/i.test(error.message)) throw new Error(`Authorized delete failed: ${error.message}`);
    });
    await runTest(`[RLS] DELETE as unauthorized from ${name}`, async () => {
      const { error } = await supabaseUnauthorized.from(name).delete().match(sampleRow);
      if (!error || !/permission denied/i.test(error.message)) {
        throw new Error(`Unauthorized delete should be denied, got: ${error ? error.message : 'no error'}`);
      }
    });
  }
}

// Run all tests
async function runAllTests() {
  await testRLSAccessControl();
  console.log('🧪 Starting CryoProtect Analyzer Database Schema Tests');
  
  await runTest('Tables exist', testTablesExist);
  await runTest('Property types exist', testPropertyTypesExist);
  await runTest('Calculation methods exist', testCalculationMethodsExist);
  await runTest('Molecule creation', testMoleculeCreation);
  await runTest('Molecular properties', testMolecularProperties);
  await runTest('Mixture and components', testMixtureAndComponents);
  await runTest('Predictions and experiments', testPredictionsAndExperiments);
  await runTest('Views', testViews);
  
  // Print summary
  console.log('\n📊 Test Summary:');
  console.log(`   ✅ Passed: ${testResults.passed}`);
  console.log(`   ❌ Failed: ${testResults.failed}`);
  console.log(`   🧪 Total:  ${testResults.passed + testResults.failed}`);
  
  if (testResults.failed > 0) {
    console.log('\n❌ Some tests failed. Please check the errors above.');
    process.exit(1);
  } else {
    console.log('\n✅ All tests passed! Your database schema is set up correctly.');
  }
}

// Check if we're authenticated
async function checkAuth() {
  const { data, error } = await supabase.auth.getSession();
  
  if (error) {
    console.error(`Authentication error: ${error.message}`);
    console.log('Please sign in before running the tests:');
    console.log('1. Set SUPABASE_URL and SUPABASE_KEY environment variables');
    console.log('2. Or modify this script to include your Supabase URL and anon key');
    process.exit(1);
  }
  
  if (!data.session) {
    console.log('No active session. Attempting to continue with anon key...');
  }
}

// Main function
async function main() {
  try {
    await checkAuth();
    await runAllTests();
  } catch (error) {
    console.error('An unexpected error occurred:', error);
    process.exit(1);
  }
}

// Run the main function
main();