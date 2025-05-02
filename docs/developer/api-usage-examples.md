# API Usage Examples

This document provides examples of how to use the CryoProtect API for common operations.

## Overview

The CryoProtect API provides RESTful endpoints for interacting with the system. This guide includes examples for:

1. Authentication
2. Working with molecules
3. Working with mixtures
4. Running predictions
5. Managing experiments
6. Error handling

## Authentication

Before using the API, you need to authenticate. CryoProtect uses Supabase for authentication.

### Signing Up

```javascript
// JavaScript example
const signUp = async (email, password) => {
  try {
    const { data, error } = await supabase.auth.signUp({
      email: email,
      password: password,
    });
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error signing up:', error.message);
    throw error;
  }
};
```

```python
# Python example
def sign_up(email, password):
    try:
        response = supabase.auth.sign_up({
            "email": email,
            "password": password
        })
        
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error signing up: {str(e)}")
        raise
```

### Signing In

```javascript
// JavaScript example
const signIn = async (email, password) => {
  try {
    const { data, error } = await supabase.auth.signInWithPassword({
      email: email,
      password: password,
    });
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error signing in:', error.message);
    throw error;
  }
};
```

```python
# Python example
def sign_in(email, password):
    try:
        response = supabase.auth.sign_in_with_password({
            "email": email,
            "password": password
        })
        
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error signing in: {str(e)}")
        raise
```

### Signing Out

```javascript
// JavaScript example
const signOut = async () => {
  try {
    const { error } = await supabase.auth.signOut();
    if (error) throw error;
  } catch (error) {
    console.error('Error signing out:', error.message);
    throw error;
  }
};
```

```python
# Python example
def sign_out():
    try:
        response = supabase.auth.sign_out()
        if response.error:
            raise Exception(response.error.message)
    except Exception as e:
        print(f"Error signing out: {str(e)}")
        raise
```

## Working with Molecules

### Fetching Molecules

```javascript
// JavaScript example
const getMolecules = async () => {
  try {
    const { data, error } = await supabase
      .from('molecules')
      .select('*');
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error fetching molecules:', error.message);
    throw error;
  }
};
```

```python
# Python example
def get_molecules():
    try:
        response = supabase.table('molecules').select('*').execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error fetching molecules: {str(e)}")
        raise
```

### Creating a Molecule

```javascript
// JavaScript example
const createMolecule = async (moleculeData) => {
  try {
    const { data, error } = await supabase
      .from('molecules')
      .insert([moleculeData])
      .select();
    
    if (error) throw error;
    return data[0];
  } catch (error) {
    console.error('Error creating molecule:', error.message);
    throw error;
  }
};
```

```python
# Python example
def create_molecule(molecule_data):
    try:
        response = supabase.table('molecules').insert(molecule_data).execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data[0]
    except Exception as e:
        print(f"Error creating molecule: {str(e)}")
        raise
```

### Updating a Molecule

```javascript
// JavaScript example
const updateMolecule = async (moleculeId, updates) => {
  try {
    const { data, error } = await supabase
      .from('molecules')
      .update(updates)
      .eq('id', moleculeId)
      .select();
    
    if (error) throw error;
    return data[0];
  } catch (error) {
    console.error('Error updating molecule:', error.message);
    throw error;
  }
};
```

```python
# Python example
def update_molecule(molecule_id, updates):
    try:
        response = supabase.table('molecules').update(updates).eq('id', molecule_id).execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data[0]
    except Exception as e:
        print(f"Error updating molecule: {str(e)}")
        raise
```

### Deleting a Molecule

```javascript
// JavaScript example
const deleteMolecule = async (moleculeId) => {
  try {
    const { error } = await supabase
      .from('molecules')
      .delete()
      .eq('id', moleculeId);
    
    if (error) throw error;
  } catch (error) {
    console.error('Error deleting molecule:', error.message);
    throw error;
  }
};
```

```python
# Python example
def delete_molecule(molecule_id):
    try:
        response = supabase.table('molecules').delete().eq('id', molecule_id).execute()
        if response.error:
            raise Exception(response.error.message)
    except Exception as e:
        print(f"Error deleting molecule: {str(e)}")
        raise
```

## Working with Mixtures

### Fetching Mixtures

```javascript
// JavaScript example
const getMixtures = async () => {
  try {
    const { data, error } = await supabase
      .from('mixtures')
      .select(`
        *,
        mixture_components (
          *,
          molecule:molecules (*)
        )
      `);
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error fetching mixtures:', error.message);
    throw error;
  }
};
```

```python
# Python example
def get_mixtures():
    try:
        response = supabase.table('mixtures').select("""
            *,
            mixture_components (
                *,
                molecule:molecules (*)
            )
        """).execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error fetching mixtures: {str(e)}")
        raise
```

### Creating a Mixture

```javascript
// JavaScript example
const createMixture = async (mixtureData, components) => {
  try {
    // Start a transaction
    const { data: mixture, error: mixtureError } = await supabase
      .from('mixtures')
      .insert([mixtureData])
      .select();
    
    if (mixtureError) throw mixtureError;
    
    // Add components
    const mixtureId = mixture[0].id;
    const componentsWithMixtureId = components.map(component => ({
      ...component,
      mixture_id: mixtureId
    }));
    
    const { error: componentsError } = await supabase
      .from('mixture_components')
      .insert(componentsWithMixtureId);
    
    if (componentsError) throw componentsError;
    
    // Fetch the complete mixture with components
    const { data: completeMixture, error: fetchError } = await supabase
      .from('mixtures')
      .select(`
        *,
        mixture_components (
          *,
          molecule:molecules (*)
        )
      `)
      .eq('id', mixtureId)
      .single();
    
    if (fetchError) throw fetchError;
    return completeMixture;
  } catch (error) {
    console.error('Error creating mixture:', error.message);
    throw error;
  }
};
```

```python
# Python example
def create_mixture(mixture_data, components):
    try:
        # Start a transaction
        response = supabase.table('mixtures').insert(mixture_data).execute()
        if response.error:
            raise Exception(response.error.message)
        
        mixture_id = response.data[0]['id']
        
        # Add components
        components_with_mixture_id = [
            {**component, 'mixture_id': mixture_id}
            for component in components
        ]
        
        component_response = supabase.table('mixture_components').insert(components_with_mixture_id).execute()
        if component_response.error:
            raise Exception(component_response.error.message)
        
        # Fetch the complete mixture with components
        complete_response = supabase.table('mixtures').select("""
            *,
            mixture_components (
                *,
                molecule:molecules (*)
            )
        """).eq('id', mixture_id).single().execute()
        
        if complete_response.error:
            raise Exception(complete_response.error.message)
        
        return complete_response.data
    except Exception as e:
        print(f"Error creating mixture: {str(e)}")
        raise
```

## Running Predictions

### Creating a Prediction

```javascript
// JavaScript example
const createPrediction = async (predictionData) => {
  try {
    const { data, error } = await supabase
      .from('predictions')
      .insert([predictionData])
      .select();
    
    if (error) throw error;
    return data[0];
  } catch (error) {
    console.error('Error creating prediction:', error.message);
    throw error;
  }
};
```

```python
# Python example
def create_prediction(prediction_data):
    try:
        response = supabase.table('predictions').insert(prediction_data).execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data[0]
    except Exception as e:
        print(f"Error creating prediction: {str(e)}")
        raise
```

### Fetching Predictions

```javascript
// JavaScript example
const getPredictions = async (mixtureId = null, moleculeId = null) => {
  try {
    let query = supabase.from('predictions').select('*');
    
    if (mixtureId) {
      query = query.eq('mixture_id', mixtureId);
    } else if (moleculeId) {
      query = query.eq('molecule_id', moleculeId);
    }
    
    const { data, error } = await query;
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error fetching predictions:', error.message);
    throw error;
  }
};
```

```python
# Python example
def get_predictions(mixture_id=None, molecule_id=None):
    try:
        query = supabase.table('predictions').select('*')
        
        if mixture_id:
            query = query.eq('mixture_id', mixture_id)
        elif molecule_id:
            query = query.eq('molecule_id', molecule_id)
        
        response = query.execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error fetching predictions: {str(e)}")
        raise
```

## Managing Experiments

### Creating an Experiment

```javascript
// JavaScript example
const createExperiment = async (experimentData, properties = []) => {
  try {
    // Create the experiment
    const { data: experiment, error: experimentError } = await supabase
      .from('experiments')
      .insert([experimentData])
      .select();
    
    if (experimentError) throw experimentError;
    
    // Add properties if provided
    if (properties.length > 0) {
      const experimentId = experiment[0].id;
      const propertiesWithExperimentId = properties.map(property => ({
        ...property,
        experiment_id: experimentId
      }));
      
      const { error: propertiesError } = await supabase
        .from('experiment_properties')
        .insert(propertiesWithExperimentId);
      
      if (propertiesError) throw propertiesError;
    }
    
    // Fetch the complete experiment with properties
    const { data: completeExperiment, error: fetchError } = await supabase
      .from('experiments')
      .select(`
        *,
        experiment_properties (*)
      `)
      .eq('id', experiment[0].id)
      .single();
    
    if (fetchError) throw fetchError;
    return completeExperiment;
  } catch (error) {
    console.error('Error creating experiment:', error.message);
    throw error;
  }
};
```

```python
# Python example
def create_experiment(experiment_data, properties=None):
    try:
        # Create the experiment
        response = supabase.table('experiments').insert(experiment_data).execute()
        if response.error:
            raise Exception(response.error.message)
        
        experiment_id = response.data[0]['id']
        
        # Add properties if provided
        if properties:
            properties_with_experiment_id = [
                {**property_data, 'experiment_id': experiment_id}
                for property_data in properties
            ]
            
            property_response = supabase.table('experiment_properties').insert(properties_with_experiment_id).execute()
            if property_response.error:
                raise Exception(property_response.error.message)
        
        # Fetch the complete experiment with properties
        complete_response = supabase.table('experiments').select("""
            *,
            experiment_properties (*)
        """).eq('id', experiment_id).single().execute()
        
        if complete_response.error:
            raise Exception(complete_response.error.message)
        
        return complete_response.data
    except Exception as e:
        print(f"Error creating experiment: {str(e)}")
        raise
```

### Fetching Experiments

```javascript
// JavaScript example
const getExperiments = async () => {
  try {
    const { data, error } = await supabase
      .from('experiments')
      .select(`
        *,
        experiment_properties (*)
      `);
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error fetching experiments:', error.message);
    throw error;
  }
};
```

```python
# Python example
def get_experiments():
    try:
        response = supabase.table('experiments').select("""
            *,
            experiment_properties (*)
        """).execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error fetching experiments: {str(e)}")
        raise
```

## Error Handling

### Common Error Patterns

```javascript
// JavaScript example
const handleApiRequest = async (apiCall) => {
  try {
    return await apiCall();
  } catch (error) {
    // Handle different types of errors
    if (error.status === 401) {
      console.error('Authentication error. Please sign in again.');
      // Redirect to login page or refresh token
    } else if (error.status === 403) {
      console.error('Permission denied. You do not have access to this resource.');
    } else if (error.status === 404) {
      console.error('Resource not found.');
    } else if (error.status === 422) {
      console.error('Validation error:', error.message);
    } else {
      console.error('An unexpected error occurred:', error.message);
    }
    throw error;
  }
};

// Usage
const fetchData = async () => {
  return handleApiRequest(async () => {
    const { data, error } = await supabase.from('molecules').select('*');
    if (error) throw error;
    return data;
  });
};
```

```python
# Python example
def handle_api_request(api_call):
    try:
        return api_call()
    except Exception as e:
        # Handle different types of errors
        error_message = str(e)
        if "401" in error_message:
            print("Authentication error. Please sign in again.")
            # Redirect to login page or refresh token
        elif "403" in error_message:
            print("Permission denied. You do not have access to this resource.")
        elif "404" in error_message:
            print("Resource not found.")
        elif "422" in error_message:
            print(f"Validation error: {error_message}")
        else:
            print(f"An unexpected error occurred: {error_message}")
        raise

# Usage
def fetch_data():
    def api_call():
        response = supabase.table('molecules').select('*').execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data
    
    return handle_api_request(api_call)
```

### Validation Error Handling

```javascript
// JavaScript example
const validateMoleculeData = (data) => {
  const errors = {};
  
  if (!data.name) {
    errors.name = 'Name is required';
  }
  
  if (!data.smiles) {
    errors.smiles = 'SMILES string is required';
  }
  
  if (Object.keys(errors).length > 0) {
    throw { status: 422, message: 'Validation failed', errors };
  }
  
  return data;
};

const createValidatedMolecule = async (moleculeData) => {
  try {
    const validatedData = validateMoleculeData(moleculeData);
    
    const { data, error } = await supabase
      .from('molecules')
      .insert([validatedData])
      .select();
    
    if (error) throw error;
    return data[0];
  } catch (error) {
    console.error('Error creating molecule:', error);
    throw error;
  }
};
```

```python
# Python example
def validate_molecule_data(data):
    errors = {}
    
    if not data.get('name'):
        errors['name'] = 'Name is required'
    
    if not data.get('smiles'):
        errors['smiles'] = 'SMILES string is required'
    
    if errors:
        raise ValueError(f"Validation failed: {errors}")
    
    return data

def create_validated_molecule(molecule_data):
    try:
        validated_data = validate_molecule_data(molecule_data)
        
        response = supabase.table('molecules').insert(validated_data).execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data[0]
    except ValueError as e:
        print(f"Validation error: {str(e)}")
        raise
    except Exception as e:
        print(f"Error creating molecule: {str(e)}")
        raise
```

## Advanced API Usage

### Filtering and Sorting

```javascript
// JavaScript example
const searchMolecules = async (filters, sortBy, ascending = true) => {
  try {
    let query = supabase.from('molecules').select('*');
    
    // Apply filters
    if (filters) {
      Object.entries(filters).forEach(([key, value]) => {
        if (value !== undefined && value !== null) {
          if (typeof value === 'string' && value.includes('%')) {
            query = query.ilike(key, value);
          } else {
            query = query.eq(key, value);
          }
        }
      });
    }
    
    // Apply sorting
    if (sortBy) {
      query = query.order(sortBy, { ascending });
    }
    
    const { data, error } = await query;
    
    if (error) throw error;
    return data;
  } catch (error) {
    console.error('Error searching molecules:', error.message);
    throw error;
  }
};
```

```python
# Python example
def search_molecules(filters=None, sort_by=None, ascending=True):
    try:
        query = supabase.table('molecules').select('*')
        
        # Apply filters
        if filters:
            for key, value in filters.items():
                if value is not None:
                    if isinstance(value, str) and '%' in value:
                        query = query.ilike(key, value)
                    else:
                        query = query.eq(key, value)
        
        # Apply sorting
        if sort_by:
            query = query.order(sort_by, ascending=ascending)
        
        response = query.execute()
        if response.error:
            raise Exception(response.error.message)
        return response.data
    except Exception as e:
        print(f"Error searching molecules: {str(e)}")
        raise
```

### Pagination

```javascript
// JavaScript example
const fetchPaginatedMolecules = async (page = 0, pageSize = 10) => {
  try {
    const from = page * pageSize;
    const to = from + pageSize - 1;
    
    const { data, error, count } = await supabase
      .from('molecules')
      .select('*', { count: 'exact' })
      .range(from, to);
    
    if (error) throw error;
    
    return {
      data,
      page,
      pageSize,
      totalCount: count,
      totalPages: Math.ceil(count / pageSize)
    };
  } catch (error) {
    console.error('Error fetching paginated molecules:', error.message);
    throw error;
  }
};
```

```python
# Python example
def fetch_paginated_molecules(page=0, page_size=10):
    try:
        from_val = page * page_size
        to_val = from_val + page_size - 1
        
        response = supabase.table('molecules').select('*', count='exact').range(from_val, to_val).execute()
        if response.error:
            raise Exception(response.error.message)
        
        count = response.count
        
        return {
            'data': response.data,
            'page': page,
            'page_size': page_size,
            'total_count': count,
            'total_pages': math.ceil(count / page_size)
        }
    except Exception as e:
        print(f"Error fetching paginated molecules: {str(e)}")
        raise
```

### Transactions

```javascript
// JavaScript example
const createMoleculeWithProperties = async (moleculeData, properties) => {
  try {
    // Start a transaction
    const { data: molecule, error: moleculeError } = await supabase
      .from('molecules')
      .insert([moleculeData])
      .select();
    
    if (moleculeError) throw moleculeError;
    
    // Add properties
    const moleculeId = molecule[0].id;
    const propertiesWithMoleculeId = properties.map(property => ({
      ...property,
      molecule_id: moleculeId
    }));
    
    const { error: propertiesError } = await supabase
      .from('molecular_properties')
      .insert(propertiesWithMoleculeId);
    
    if (propertiesError) throw propertiesError;
    
    // Fetch the complete molecule with properties
    const { data: completeMolecule, error: fetchError } = await supabase
      .from('molecules')
      .select(`
        *,
        molecular_properties (*)
      `)
      .eq('id', moleculeId)
      .single();
    
    if (fetchError) throw fetchError;
    return completeMolecule;
  } catch (error) {
    console.error('Error creating molecule with properties:', error.message);
    throw error;
  }
};
```

```python
# Python example
def create_molecule_with_properties(molecule_data, properties):
    try:
        # Start a transaction
        response = supabase.table('molecules').insert(molecule_data).execute()
        if response.error:
            raise Exception(response.error.message)
        
        molecule_id = response.data[0]['id']
        
        # Add properties
        properties_with_molecule_id = [
            {**property_data, 'molecule_id': molecule_id}
            for property_data in properties
        ]
        
        property_response = supabase.table('molecular_properties').insert(properties_with_molecule_id).execute()
        if property_response.error:
            raise Exception(property_response.error.message)
        
        # Fetch the complete molecule with properties
        complete_response = supabase.table('molecules').select("""
            *,
            molecular_properties (*)
        """).eq('id', molecule_id).single().execute()
        
        if complete_response.error:
            raise Exception(complete_response.error.message)
        
        return complete_response.data
    except Exception as e:
        print(f"Error creating molecule with properties: {str(e)}")
        raise
```

## RDKit Integration

The CryoProtect API provides powerful molecular property calculation, visualization, and searching capabilities through RDKit integration.

### Calculating Molecular Properties

```javascript
// JavaScript example
const calculateProperties = async (moleculeData, inputFormat = 'smiles') => {
  try {
    const response = await fetch('/api/v1/rdkit/properties', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        molecule_data: moleculeData,
        input_format: inputFormat
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error calculating properties');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error calculating molecular properties:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def calculate_properties(molecule_data, input_format='smiles'):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/rdkit/properties',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'molecule_data': molecule_data,
                'input_format': input_format
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error calculating molecular properties: {str(e)}")
        raise
```

### Generating Molecular Visualization

```javascript
// JavaScript example
const generateVisualization = async (moleculeData, inputFormat = 'smiles', width = 400, height = 300, highlightAtoms = []) => {
  try {
    const response = await fetch('/api/v1/rdkit/visualization', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        molecule_data: moleculeData,
        input_format: inputFormat,
        width: width,
        height: height,
        highlight_atoms: highlightAtoms
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error generating visualization');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error generating molecular visualization:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def generate_visualization(molecule_data, input_format='smiles', width=400, height=300, highlight_atoms=None):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/rdkit/visualization',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'molecule_data': molecule_data,
                'input_format': input_format,
                'width': width,
                'height': height,
                'highlight_atoms': highlight_atoms or []
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error generating molecular visualization: {str(e)}")
        raise
```

### Performing Substructure Search

```javascript
// JavaScript example
const performSubstructureSearch = async (queryMolData, targetMolData, queryFormat = 'smarts', targetFormat = 'smiles') => {
  try {
    const response = await fetch('/api/v1/rdkit/substructure', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        query_mol_data: queryMolData,
        target_mol_data: targetMolData,
        query_format: queryFormat,
        target_format: targetFormat
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error performing substructure search');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error performing substructure search:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def perform_substructure_search(query_mol_data, target_mol_data, query_format='smarts', target_format='smiles'):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/rdkit/substructure',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'query_mol_data': query_mol_data,
                'target_mol_data': target_mol_data,
                'query_format': query_format,
                'target_format': target_format
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error performing substructure search: {str(e)}")
        raise
```

### Calculating Molecular Similarity

```javascript
// JavaScript example
const calculateSimilarity = async (mol1Data, mol2Data, mol1Format = 'smiles', mol2Format = 'smiles', fingerprintType = 'morgan') => {
  try {
    const response = await fetch('/api/v1/rdkit/similarity', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        mol1_data: mol1Data,
        mol2_data: mol2Data,
        mol1_format: mol1Format,
        mol2_format: mol2Format,
        fingerprint_type: fingerprintType
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error calculating similarity');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error calculating molecular similarity:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def calculate_similarity(mol1_data, mol2_data, mol1_format='smiles', mol2_format='smiles', fingerprint_type='morgan'):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/rdkit/similarity',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'mol1_data': mol1_data,
                'mol2_data': mol2_data,
                'mol1_format': mol1_format,
                'mol2_format': mol2_format,
                'fingerprint_type': fingerprint_type
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error calculating molecular similarity: {str(e)}")
        raise
```

## Scoring

The CryoProtect API provides endpoints for scoring molecules and mixtures based on their cryoprotection effectiveness.

### Scoring a Molecule

```javascript
// JavaScript example
const scoreMolecule = async (moleculeData, inputFormat = 'smiles', storeResult = false) => {
  try {
    const response = await fetch('/api/v1/scoring/molecules', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        molecule_data: moleculeData,
        input_format: inputFormat,
        store_result: storeResult
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error scoring molecule');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error scoring molecule:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def score_molecule(molecule_data, input_format='smiles', store_result=False):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/scoring/molecules',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'molecule_data': molecule_data,
                'input_format': input_format,
                'store_result': store_result
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error scoring molecule: {str(e)}")
        raise
```

### Scoring a Molecule by ID

```javascript
// JavaScript example
const scoreMoleculeById = async (moleculeId, storeResult = true) => {
  try {
    const response = await fetch(`/api/v1/molecules/${moleculeId}/score`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        store_result: storeResult
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error scoring molecule');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error scoring molecule by ID:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def score_molecule_by_id(molecule_id, store_result=True):
    try:
        response = requests.post(
            f'http://localhost:5000/api/v1/molecules/{molecule_id}/score',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'store_result': store_result
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error scoring molecule by ID: {str(e)}")
        raise
```

### Scoring a Mixture

```javascript
// JavaScript example
const scoreMixture = async (mixtureId, storeResult = true) => {
  try {
    const response = await fetch(`/api/v1/mixtures/${mixtureId}/score`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        store_result: storeResult
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error scoring mixture');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error scoring mixture:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def score_mixture(mixture_id, store_result=True):
    try:
        response = requests.post(
            f'http://localhost:5000/api/v1/mixtures/{mixture_id}/score',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'store_result': store_result
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error scoring mixture: {str(e)}")
        raise
```

### Batch Scoring

```javascript
// JavaScript example
const batchScore = async (entityType, ids, storeResults = true) => {
  try {
    const response = await fetch('/api/v1/scoring/batch', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        entity_type: entityType,
        ids: ids,
        store_results: storeResults
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error performing batch scoring');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error performing batch scoring:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests

def batch_score(entity_type, ids, store_results=True):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/scoring/batch',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'entity_type': entity_type,
                'ids': ids,
                'store_results': store_results
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error performing batch scoring: {str(e)}")
        raise
```

## Export and Sharing

The CryoProtect API provides endpoints for exporting data in various formats and sharing results.

### Exporting Data

```javascript
// JavaScript example
const exportData = async (dataType, format, id = null, includeRelated = false) => {
  try {
    const response = await fetch('/api/v1/export', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        data_type: dataType,
        format: format,
        id: id,
        include_related: includeRelated
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error exporting data');
    }
    
    // Handle file download
    const blob = await response.blob();
    const filename = response.headers.get('Content-Disposition')
      ? response.headers.get('Content-Disposition').split('filename=')[1].replace(/"/g, '')
      : `${dataType}_export.${format}`;
    
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.style.display = 'none';
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
    
    return { success: true, filename };
  } catch (error) {
    console.error('Error exporting data:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests
import os

def export_data(data_type, format, id=None, include_related=False):
    try:
        response = requests.post(
            'http://localhost:5000/api/v1/export',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json={
                'data_type': data_type,
                'format': format,
                'id': id,
                'include_related': include_related
            },
            stream=True  # Stream the response for file downloads
        )
        
        response.raise_for_status()
        
        # Get filename from Content-Disposition header or create a default one
        if 'Content-Disposition' in response.headers:
            filename = response.headers['Content-Disposition'].split('filename=')[1].replace('"', '')
        else:
            filename = f"{data_type}_export.{format}"
        
        # Save the file
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        return {'success': True, 'filename': filename, 'path': os.path.abspath(filename)}
    except requests.exceptions.RequestException as e:
        print(f"Error exporting data: {str(e)}")
        raise
```

### Sharing Results

```javascript
// JavaScript example
const shareResults = async (shareType, dataType, id, options = {}) => {
  try {
    const requestBody = {
      share_type: shareType,
      data_type: dataType,
      id: id,
      ...options
    };
    
    const response = await fetch('/api/v1/share', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify(requestBody)
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error sharing results');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error sharing results:', error.message);
    throw error;
  }
};

// Example usage for different share types
const shareViaLink = async (dataType, id, passwordProtected = false, password = null, expiration = 86400) => {
  return shareResults('link', dataType, id, {
    password_protected: passwordProtected,
    password: password,
    expiration: expiration
  });
};

const shareViaEmail = async (dataType, id, recipients, message = '', passwordProtected = false, password = null) => {
  return shareResults('email', dataType, id, {
    recipients: recipients,
    message: message,
    password_protected: passwordProtected,
    password: password
  });
};

const getEmbedCode = async (dataType, id) => {
  const result = await shareResults('embed', dataType, id);
  return result.embed_code;
};
```

```python
# Python example
import requests

def share_results(share_type, data_type, id, **options):
    try:
        request_body = {
            'share_type': share_type,
            'data_type': data_type,
            'id': id,
            **options
        }
        
        response = requests.post(
            'http://localhost:5000/api/v1/share',
            headers={
                'Content-Type': 'application/json',
                'Authorization': f'Bearer {get_token()}'
            },
            json=request_body
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error sharing results: {str(e)}")
        raise

# Example usage for different share types
def share_via_link(data_type, id, password_protected=False, password=None, expiration=86400):
    return share_results('link', data_type, id,
                        password_protected=password_protected,
                        password=password,
                        expiration=expiration)

def share_via_email(data_type, id, recipients, message='', password_protected=False, password=None):
    return share_results('email', data_type, id,
                        recipients=recipients,
                        message=message,
                        password_protected=password_protected,
                        password=password)

def get_embed_code(data_type, id):
    result = share_results('embed', data_type, id)
    return result.get('embed_code')
```

## Rate Limiting

The CryoProtect API implements rate limiting to ensure fair usage and system stability. Here's how to handle rate limiting in your applications:

### Handling Rate Limit Errors

```javascript
// JavaScript example
const handleRateLimitedRequest = async (url, options) => {
  try {
    const response = await fetch(url, options);
    
    // Check for rate limit headers
    const rateLimitLimit = response.headers.get('X-RateLimit-Limit');
    const rateLimitRemaining = response.headers.get('X-RateLimit-Remaining');
    const rateLimitReset = response.headers.get('X-RateLimit-Reset');
    
    if (rateLimitRemaining && parseInt(rateLimitRemaining) < 5) {
      console.warn(`Rate limit warning: ${rateLimitRemaining}/${rateLimitLimit} requests remaining. Resets at ${new Date(parseInt(rateLimitReset) * 1000).toLocaleTimeString()}`);
    }
    
    if (response.status === 429) {
      // Rate limit exceeded
      const retryAfter = response.headers.get('Retry-After') || 60;
      console.warn(`Rate limit exceeded. Retrying after ${retryAfter} seconds.`);
      
      // Wait for the specified time and retry
      await new Promise(resolve => setTimeout(resolve, parseInt(retryAfter) * 1000));
      return handleRateLimitedRequest(url, options);
    }
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Request failed');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Request error:', error.message);
    throw error;
  }
};
```

```python
# Python example
import requests
import time
from datetime import datetime

def handle_rate_limited_request(url, method='GET', headers=None, json=None, params=None):
    try:
        if method.upper() == 'GET':
            response = requests.get(url, headers=headers, params=params)
        elif method.upper() == 'POST':
            response = requests.post(url, headers=headers, json=json)
        else:
            raise ValueError(f"Unsupported method: {method}")
        
        # Check for rate limit headers
        rate_limit_limit = response.headers.get('X-RateLimit-Limit')
        rate_limit_remaining = response.headers.get('X-RateLimit-Remaining')
        rate_limit_reset = response.headers.get('X-RateLimit-Reset')
        
        if rate_limit_remaining and int(rate_limit_remaining) < 5:
            reset_time = datetime.fromtimestamp(int(rate_limit_reset)).strftime('%H:%M:%S')
            print(f"Rate limit warning: {rate_limit_remaining}/{rate_limit_limit} requests remaining. Resets at {reset_time}")
        
        if response.status_code == 429:
            # Rate limit exceeded
            retry_after = int(response.headers.get('Retry-After', 60))
            print(f"Rate limit exceeded. Retrying after {retry_after} seconds.")
            
            # Wait for the specified time and retry
            time.sleep(retry_after)
            return handle_rate_limited_request(url, method, headers, json, params)
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Request error: {str(e)}")
        raise
```

### Implementing Exponential Backoff

```javascript
// JavaScript example
const fetchWithBackoff = async (url, options, maxRetries = 5) => {
  let retries = 0;
  
  while (retries < maxRetries) {
    try {
      const response = await fetch(url, options);
      
      if (response.status === 429) {
        // Rate limit exceeded
        const retryAfter = parseInt(response.headers.get('Retry-After') || '60');
        const backoffTime = retryAfter * 1000 * Math.pow(2, retries);
        console.warn(`Rate limit exceeded. Retrying after ${backoffTime/1000} seconds (retry ${retries + 1}/${maxRetries}).`);
        
        // Wait with exponential backoff
        await new Promise(resolve => setTimeout(resolve, backoffTime));
        retries++;
        continue;
      }
      
      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.message || 'Request failed');
      }
      
      return await response.json();
    } catch (error) {
      if (retries >= maxRetries - 1) {
        console.error('Max retries reached:', error.message);
        throw error;
      }
      
      // For other errors, also implement backoff
      const backoffTime = 1000 * Math.pow(2, retries);
      console.warn(`Request failed. Retrying after ${backoffTime/1000} seconds (retry ${retries + 1}/${maxRetries}).`);
      await new Promise(resolve => setTimeout(resolve, backoffTime));
      retries++;
    }
  }
};
```

```python
# Python example
import requests
import time
import random

def fetch_with_backoff(url, method='GET', headers=None, json=None, params=None, max_retries=5):
    retries = 0
    
    while retries < max_retries:
        try:
            if method.upper() == 'GET':
                response = requests.get(url, headers=headers, params=params)
            elif method.upper() == 'POST':
                response = requests.post(url, headers=headers, json=json)
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            if response.status_code == 429:
                # Rate limit exceeded
                retry_after = int(response.headers.get('Retry-After', 60))
                backoff_time = retry_after * (2 ** retries)
                # Add jitter to avoid thundering herd problem
                jitter = random.uniform(0, 0.1 * backoff_time)
                backoff_time += jitter
                
                print(f"Rate limit exceeded. Retrying after {backoff_time:.1f} seconds (retry {retries + 1}/{max_retries}).")
                time.sleep(backoff_time)
                retries += 1
                continue
            
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            if retries >= max_retries - 1:
                print(f"Max retries reached: {str(e)}")
                raise
            
            # For other errors, also implement backoff
            backoff_time = (2 ** retries) + random.uniform(0, 1)
            print(f"Request failed. Retrying after {backoff_time:.1f} seconds (retry {retries + 1}/{max_retries}).")
            time.sleep(backoff_time)
            retries += 1
```

## Conclusion

These examples demonstrate how to use the CryoProtect API for common operations. For more detailed information about specific endpoints, refer to the API documentation. Remember to handle errors appropriately and validate input data to ensure robust applications.

When working with the API, keep these best practices in mind:

1. **Authentication**: Always use secure methods to store and manage authentication tokens
2. **Error Handling**: Implement comprehensive error handling for all API requests
3. **Rate Limiting**: Use backoff strategies to handle rate limiting gracefully
4. **Validation**: Validate input data before sending it to the API
5. **Security**: Follow security best practices, especially when handling sensitive data