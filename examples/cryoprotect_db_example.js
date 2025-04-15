/**
 * CryoProtect Analyzer - Database Usage Example
 * 
 * This example demonstrates how to interact with the CryoProtect database
 * using the Supabase JavaScript client.
 * 
 * Prerequisites:
 * - Node.js installed
 * - Supabase project with the CryoProtect schema applied
 * - @supabase/supabase-js package installed (npm install @supabase/supabase-js)
 */

const { createClient } = require('@supabase/supabase-js');

// Replace with your Supabase URL and anon key
const supabaseUrl = 'https://your-project-ref.supabase.co';
const supabaseKey = 'your-anon-key';
const supabase = createClient(supabaseUrl, supabaseKey);

/**
 * Example 1: Import a molecule from PubChem
 * @param {number} cid - PubChem Compound ID
 * @returns {Promise<string>} - UUID of the imported molecule
 */
async function importMoleculeFromPubChem(cid) {
  console.log(`Importing molecule with CID ${cid} from PubChem...`);
  
  // Call the database function to import the molecule
  const { data, error } = await supabase.rpc('import_molecule_from_pubchem', {
    p_cid: cid,
    p_user_id: supabase.auth.user()?.id
  });
  
  if (error) {
    console.error('Error importing molecule:', error);
    throw error;
  }
  
  console.log(`Successfully imported molecule with ID: ${data}`);
  return data;
}

/**
 * Example 2: Add molecular properties
 * @param {string} moleculeId - UUID of the molecule
 * @param {Object} properties - Object with property name-value pairs
 */
async function addMolecularProperties(moleculeId, properties) {
  console.log(`Adding properties for molecule ${moleculeId}...`);
  
  // Get property types
  const { data: propertyTypes, error: propertyTypesError } = await supabase
    .from('property_types')
    .select('id, name, data_type');
  
  if (propertyTypesError) {
    console.error('Error fetching property types:', propertyTypesError);
    throw propertyTypesError;
  }
  
  // Prepare property inserts
  const propertyInserts = [];
  
  for (const [propertyName, value] of Object.entries(properties)) {
    const propertyType = propertyTypes.find(pt => pt.name === propertyName);
    if (!propertyType) {
      console.warn(`Property type "${propertyName}" not found, skipping`);
      continue;
    }
    
    const propertyInsert = {
      molecule_id: moleculeId,
      property_type_id: propertyType.id,
      created_by: supabase.auth.user()?.id
    };
    
    // Set the appropriate value field based on data type
    switch (propertyType.data_type) {
      case 'numeric':
        propertyInsert.numeric_value = typeof value === 'number' ? value : parseFloat(value);
        break;
      case 'text':
        propertyInsert.text_value = String(value);
        break;
      case 'boolean':
        propertyInsert.boolean_value = Boolean(value);
        break;
      default:
        console.warn(`Unknown data type "${propertyType.data_type}" for property "${propertyName}", skipping`);
        continue;
    }
    
    propertyInserts.push(propertyInsert);
  }
  
  // Insert properties
  const { data, error } = await supabase
    .from('molecular_properties')
    .upsert(propertyInserts, { onConflict: 'molecule_id,property_type_id' });
  
  if (error) {
    console.error('Error adding properties:', error);
    throw error;
  }
  
  console.log(`Successfully added ${propertyInserts.length} properties`);
  return data;
}

/**
 * Example 3: Create a mixture
 * @param {string} name - Name of the mixture
 * @param {string} description - Description of the mixture
 * @param {Array} components - Array of { moleculeId, concentration, unit } objects
 * @returns {Promise<string>} - UUID of the created mixture
 */
async function createMixture(name, description, components) {
  console.log(`Creating mixture "${name}"...`);
  
  // Start a transaction
  const { data: mixture, error: mixtureError } = await supabase
    .from('mixtures')
    .insert({
      name,
      description,
      created_by: supabase.auth.user()?.id
    })
    .select()
    .single();
  
  if (mixtureError) {
    console.error('Error creating mixture:', mixtureError);
    throw mixtureError;
  }
  
  const mixtureId = mixture.id;
  
  // Add components
  const componentInserts = components.map(comp => ({
    mixture_id: mixtureId,
    molecule_id: comp.moleculeId,
    concentration: comp.concentration,
    concentration_unit: comp.unit,
    created_by: supabase.auth.user()?.id
  }));
  
  const { error: componentsError } = await supabase
    .from('mixture_components')
    .insert(componentInserts);
  
  if (componentsError) {
    console.error('Error adding mixture components:', componentsError);
    throw componentsError;
  }
  
  console.log(`Successfully created mixture with ID: ${mixtureId}`);
  return mixtureId;
}

/**
 * Example 4: Add a prediction for a mixture
 * @param {string} mixtureId - UUID of the mixture
 * @param {string} propertyName - Name of the property being predicted
 * @param {any} value - Predicted value
 * @param {number} confidence - Confidence level (0-1)
 * @param {string} methodName - Name of the calculation method
 */
async function addPrediction(mixtureId, propertyName, value, confidence, methodName) {
  console.log(`Adding prediction for mixture ${mixtureId}, property "${propertyName}"...`);
  
  // Get property type ID
  const { data: propertyType, error: propertyTypeError } = await supabase
    .from('property_types')
    .select('id, data_type')
    .eq('name', propertyName)
    .single();
  
  if (propertyTypeError) {
    console.error(`Error fetching property type "${propertyName}":`, propertyTypeError);
    throw propertyTypeError;
  }
  
  // Get calculation method ID
  const { data: calculationMethod, error: methodError } = await supabase
    .from('calculation_methods')
    .select('id')
    .eq('name', methodName)
    .single();
  
  if (methodError) {
    console.error(`Error fetching calculation method "${methodName}":`, methodError);
    throw methodError;
  }
  
  // Prepare prediction insert
  const prediction = {
    mixture_id: mixtureId,
    property_type_id: propertyType.id,
    calculation_method_id: calculationMethod.id,
    confidence,
    created_by: supabase.auth.user()?.id
  };
  
  // Set the appropriate value field based on data type
  switch (propertyType.data_type) {
    case 'numeric':
      prediction.numeric_value = typeof value === 'number' ? value : parseFloat(value);
      break;
    case 'text':
      prediction.text_value = String(value);
      break;
    case 'boolean':
      prediction.boolean_value = Boolean(value);
      break;
    default:
      console.error(`Unknown data type "${propertyType.data_type}"`);
      throw new Error(`Unknown data type "${propertyType.data_type}"`);
  }
  
  // Insert prediction
  const { data, error } = await supabase
    .from('predictions')
    .upsert(prediction, { 
      onConflict: 'mixture_id,property_type_id,calculation_method_id' 
    });
  
  if (error) {
    console.error('Error adding prediction:', error);
    throw error;
  }
  
  console.log(`Successfully added prediction`);
  return data;
}

/**
 * Example 5: Record experimental results
 * @param {string} mixtureId - UUID of the mixture
 * @param {string} propertyName - Name of the property being measured
 * @param {any} value - Measured value
 * @param {string} conditions - Experimental conditions
 * @param {string} date - Date of experiment (YYYY-MM-DD)
 */
async function recordExperiment(mixtureId, propertyName, value, conditions, date) {
  console.log(`Recording experiment for mixture ${mixtureId}, property "${propertyName}"...`);
  
  // Get property type ID
  const { data: propertyType, error: propertyTypeError } = await supabase
    .from('property_types')
    .select('id, data_type')
    .eq('name', propertyName)
    .single();
  
  if (propertyTypeError) {
    console.error(`Error fetching property type "${propertyName}":`, propertyTypeError);
    throw propertyTypeError;
  }
  
  // Prepare experiment insert
  const experiment = {
    mixture_id: mixtureId,
    property_type_id: propertyType.id,
    experimental_conditions: conditions,
    date_performed: date,
    created_by: supabase.auth.user()?.id
  };
  
  // Set the appropriate value field based on data type
  switch (propertyType.data_type) {
    case 'numeric':
      experiment.numeric_value = typeof value === 'number' ? value : parseFloat(value);
      break;
    case 'text':
      experiment.text_value = String(value);
      break;
    case 'boolean':
      experiment.boolean_value = Boolean(value);
      break;
    default:
      console.error(`Unknown data type "${propertyType.data_type}"`);
      throw new Error(`Unknown data type "${propertyType.data_type}"`);
  }
  
  // Insert experiment
  const { data, error } = await supabase
    .from('experiments')
    .insert(experiment);
  
  if (error) {
    console.error('Error recording experiment:', error);
    throw error;
  }
  
  console.log(`Successfully recorded experiment`);
  return data;
}

/**
 * Example 6: Compare prediction with experiment
 * @param {string} mixtureId - UUID of the mixture
 * @param {string} propertyName - Name of the property
 */
async function comparePredictionWithExperiment(mixtureId, propertyName) {
  console.log(`Comparing prediction with experiment for mixture ${mixtureId}, property "${propertyName}"...`);
  
  // Get property type ID
  const { data: propertyType, error: propertyTypeError } = await supabase
    .from('property_types')
    .select('id')
    .eq('name', propertyName)
    .single();
  
  if (propertyTypeError) {
    console.error(`Error fetching property type "${propertyName}":`, propertyTypeError);
    throw propertyTypeError;
  }
  
  // Call the database function to compare
  const { data, error } = await supabase.rpc('compare_prediction_with_experiment', {
    p_mixture_id: mixtureId,
    p_property_type_id: propertyType.id
  });
  
  if (error) {
    console.error('Error comparing prediction with experiment:', error);
    throw error;
  }
  
  console.log('Comparison result:', data);
  return data;
}

/**
 * Example 7: Get molecule with all its properties
 * @param {string} moleculeId - UUID of the molecule
 */
async function getMoleculeWithProperties(moleculeId) {
  console.log(`Getting molecule ${moleculeId} with properties...`);
  
  const { data, error } = await supabase
    .from('molecule_with_properties')
    .select('*')
    .eq('id', moleculeId)
    .single();
  
  if (error) {
    console.error('Error getting molecule with properties:', error);
    throw error;
  }
  
  console.log('Molecule with properties:', data);
  return data;
}

/**
 * Example 8: Get mixture with all its components
 * @param {string} mixtureId - UUID of the mixture
 */
async function getMixtureWithComponents(mixtureId) {
  console.log(`Getting mixture ${mixtureId} with components...`);
  
  const { data, error } = await supabase
    .from('mixture_with_components')
    .select('*')
    .eq('id', mixtureId)
    .single();
  
  if (error) {
    console.error('Error getting mixture with components:', error);
    throw error;
  }
  
  console.log('Mixture with components:', data);
  return data;
}

/**
 * Example usage
 */
async function runExample() {
  try {
    // Sign in (replace with your auth method)
    const { user, error } = await supabase.auth.signIn({
      email: 'user@example.com',
      password: 'password'
    });
    
    if (error) {
      console.error('Authentication error:', error);
      return;
    }
    
    // Example 1: Import a molecule from PubChem (Glycerol, CID: 753)
    const moleculeId = await importMoleculeFromPubChem(753);
    
    // Example 2: Add molecular properties
    await addMolecularProperties(moleculeId, {
      'Molecular Weight': 92.09,
      'LogP': -1.76,
      'TPSA': 60.69,
      'H-Bond Donors': 3,
      'H-Bond Acceptors': 3,
      'Toxicity': 'Low toxicity',
      'Stability': 'Stable under normal conditions',
      'Environmental Safety': 'Biodegradable',
      'Total Score': 180
    });
    
    // Example 3: Create a mixture
    const mixtureId = await createMixture(
      'Glycerol-Water Solution',
      'A 30% glycerol solution in water',
      [
        { moleculeId, concentration: 30, unit: '%' },
        // Water would be another component, but we'd need to import it first
      ]
    );
    
    // Example 4: Add a prediction
    await addPrediction(
      mixtureId,
      'Freezing Point',
      -15.3,
      0.9,
      'CryoProtect Scoring'
    );
    
    // Example 5: Record experimental results
    await recordExperiment(
      mixtureId,
      'Freezing Point',
      -14.8,
      'Standard pressure, cooling rate 1°C/min',
      '2025-04-14'
    );
    
    // Example 6: Compare prediction with experiment
    const comparison = await comparePredictionWithExperiment(
      mixtureId,
      'Freezing Point'
    );
    
    // Example 7: Get molecule with properties
    const moleculeWithProps = await getMoleculeWithProperties(moleculeId);
    
    // Example 8: Get mixture with components
    const mixtureWithComps = await getMixtureWithComponents(mixtureId);
    
    console.log('Example completed successfully!');
  } catch (error) {
    console.error('Example failed:', error);
  }
}

// Run the example
runExample();